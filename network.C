/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ CkGroupID mCastGrpId;
extern /*readonly*/ idx_t npdat;
extern /*readonly*/ idx_t npnet;
extern /*readonly*/ idx_t rngseed;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ tick_t tcheck;
extern /*readonly*/ tick_t trecord;
extern /*readonly*/ tick_t tdisplay;
extern /*readonly*/ idx_t equeue;


/**************************************************************************
* Charm++ Reduction
**************************************************************************/

// Reduction for type tick_t
//
CkReduction::reducerType min_tick;
// Initnode
void registerMinTick(void) {
  min_tick = CkReduction::addReducer(minTick);
}
// Reduction function
CkReductionMsg *minTick(int nMsg, CkReductionMsg **msgs) {
  // Initialize to max
  tick_t ret = TICK_T_MAX;
  for (int i = 0; i < nMsg; ++i) {
    // Sanity check
    CkAssert(msgs[i]->getSize() == sizeof(tick_t));
    // Extract data and reduce 
    tick_t m = *(tick_t *)msgs[i]->getData();
    ret = std::min(ret,m);
  }
  // Return minimum tick_t
  return CkReductionMsg::buildNew(sizeof(tick_t),&ret);
}


/**************************************************************************
* Network
**************************************************************************/

// Network constructor
//
Network::Network(mModel *msg) {
  // Bookkeeping
  prtidx = thisIndex;
  // Network part (prtidx) to in/out file (datidx) conversion
  // (with excess parts located modulo at the beginning)
  idx_t datdiv = npnet/npdat;
  idx_t datrem = npnet%npdat;
  datidx = (prtidx-( ((prtidx+1)/(datdiv+1))>datrem ? datrem:((prtidx+1)/(datdiv+1)) ))/datdiv;
  
  // Set up random number generator
  rngine.seed(rngseed+prtidx);
  unifdist = new std::uniform_real_distribution<real_t> (0.0, 1.0);

  // Vertex Models
  for (std::size_t i = 0; i < netmodel.size(); ++i) {
    delete netmodel[i];
  }
  // Set up containers
  netmodel.clear();
  // "none" model
  netmodel.push_back(NetModelFactory::getNetModel()->Create(0));
  netmodel[0]->setPNGActive(false);
  netmodel[0]->setPNGMother(false);
  netmodel[0]->setPNGAnchor(false);
  // User defined models
  for (idx_t i = 1; i < msg->nmodel+1; ++i) {
    // Create model object
    netmodel.push_back(NetModelFactory::getNetModel()->Create(msg->modtype[i-1]));
    CkAssert(netmodel[i]->getNState() == msg->nstate[i-1]);
    CkAssert(netmodel[i]->getNStick() == msg->nstick[i-1]);
    CkAssert(netmodel[i]->getNParam() == msg->xparam[i] - msg->xparam[i-1]);
    netmodel[i]->setParam(msg->param + msg->xparam[i-1]);
    netmodel[i]->setPort(msg->port + msg->xport[i-1]);
    netmodel[i]->setRandom(unifdist, &rngine);
    netmodel[i]->setPNGActive(msg->pngactive[i-1]);
    netmodel[i]->setPNGMother(msg->pngmother[i-1]);
    netmodel[i]->setPNGAnchor(msg->pnganchor[i-1]);

    // Print out model information
    /*
    if (prtidx == 0) {
      std::string params;
      // collect params
      for (idx_t j = msg->xparam[i-1]; j < msg->xparam[i]; ++j) {
        std::ostringstream param;
        param << " " << msg->param[j];
        params.append(param.str());
      }
      CkPrintf("  Network model: %" PRIidx "   ModType: %" PRIidx "   Params:%s\n", i, netmodel[i]->getModType(), params.c_str());
    }
    */
  }
  delete msg;

  // set up auxiliary state
  edgaux.resize(netmodel.size());
  for (std::size_t j = 1; j < netmodel.size(); ++j) {
    if (netmodel[j]->getNAux()) {
      edgaux[j].resize(netmodel.size());
      std::vector<std::string> auxstate = netmodel[j]->getAuxState();
      std::vector<std::string> auxstick = netmodel[j]->getAuxStick();
      for (std::size_t i = 1; i < netmodel.size(); ++i) {
        edgaux[j][i].resize(1);
        edgaux[j][i][0].index = 0;
        edgaux[j][i][0].stateidx.resize(auxstate.size());
        for (std::size_t s = 0; s < auxstate.size(); ++s) {
          edgaux[j][i][0].stateidx[s] = netmodel[i]->getStateIdx(auxstate[s]);
        }
        edgaux[j][i][0].stickidx.resize(auxstick.size());
        for (std::size_t s = 0; s < auxstick.size(); ++s) {
          edgaux[j][i][0].stickidx[s] = netmodel[i]->getStateIdx(auxstick[s]);
        }
      }
    }
  }
  // set up repeating events
  repevt.clear();
  repmodidx.resize(netmodel.size(), false);
  for (std::size_t i = 1; i < netmodel.size(); ++i) {
    netmodel[i]->addRepeat((idx_t) i, repevt);
  }
  if (repevt.size()) {
    std::sort(repevt.begin(), repevt.end());
    trep = repevt[0].diffuse;
    for (std::size_t e = 0; e < repevt.size(); ++e) {
      repmodidx[repevt[e].index] = true;
    }
  }
  else {
    trep = TICK_T_MAX;
  }

  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// Network migration
//
Network::Network(CkMigrateMessage *msg) {
  delete msg;
}

// Network destructor
//
Network::~Network() {
  for (std::size_t i = 0; i < netmodel.size(); ++i) {
    delete netmodel[i];
  }
}


/**************************************************************************
* Network Load Data
**************************************************************************/

// Receive data from Netdata chare array
//
void Network::LoadNetwork(mPart *msg) {
  /* Bookkeeping */
  idx_t nadjcy;  // edges
  idx_t jmodidx; // modidx
  idx_t jstate;  // state
  idx_t jstick;  // time state
  idx_t nevent;  // events

  // Copy over part data
  CkAssert(msg->prtidx == prtidx);

  // Setup data vectors
  vtxdist.resize(npnet+1);
  vtxidx.resize(msg->nvtx);
  vtxmap.clear();
  vtxmodidx.resize(msg->nvtx);
  xyz.resize(msg->nvtx*3);
  pngs.resize(msg->nvtx);
  pngcan.clear();
  pngaux.clear();
  adjcy.resize(msg->nvtx);
  adjmap.clear();
  edgmodidx.resize(msg->nvtx);
  state.resize(msg->nvtx);
  stick.resize(msg->nvtx);
  vtxaux.resize(msg->nvtx);
  event.resize(msg->nvtx);
  evtaux.resize(msg->nvtx);
  evtlog.clear();
  evtext.clear();
  repidx.clear();
  record.clear();
  recordlist.clear();
  recevt.clear();
  recevtlist.resize(EVENT_TOTAL);
  // Record spikes
  recevtlist[EVENT_SPIKE] = true;

  // Graph distribution information
  for (idx_t i = 0; i < npnet+1; ++i) {
    // vtxdist
    vtxdist[i] = msg->vtxdist[i];
  }

  // Initalize counters
  nadjcy = 0;
  jmodidx = 0;
  jstate = 0;
  jstick = 0;
  nevent = 0;

  // Get adjacency matrix
  for (idx_t i = 0; i < adjcy.size(); ++i) {
    vtxidx[i] = vtxdist[prtidx] + i;
    vtxmap[vtxdist[prtidx] + i] = i;
    vtxmodidx[i] = msg->vtxmodidx[i];
    xyz[i*3+0] = msg->xyz[i*3+0];
    xyz[i*3+1] = msg->xyz[i*3+1];
    xyz[i*3+2] = msg->xyz[i*3+2];
    pngs[i].clear();
    // preallocate sizes
    adjcy[i].resize(msg->xadj[i+1] - msg->xadj[i]);
    edgmodidx[i].resize(msg->xadj[i+1] - msg->xadj[i]);
    state[i].resize(msg->xadj[i+1] - msg->xadj[i] + 1);
    stick[i].resize(msg->xadj[i+1] - msg->xadj[i] + 1);
    event[i].resize(equeue);
    nadjcy += adjcy[i].size();
    // copy over vertex data
    state[i][0] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + netmodel[vtxmodidx[i]]->getNState());
    jstate += netmodel[vtxmodidx[i]]->getNState();
    stick[i][0] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + netmodel[vtxmodidx[i]]->getNStick());
    jstick += netmodel[vtxmodidx[i]]->getNStick();
    if (repmodidx[vtxmodidx[i]]) {
      repidx[vtxmodidx[i]].push_back(std::array<idx_t, 2>{{i, 0}});
    }
    // copy over edge data
    for (idx_t j = 0; j < adjcy[i].size(); ++j) {
      adjcy[i][j] = msg->adjcy[msg->xadj[i] + j];
      // models
      edgmodidx[i][j] = msg->edgmodidx[jmodidx++];
      if (edgmodidx[i][j]) {
        // map of targets (edge model in part)
        adjmap[adjcy[i][j]].push_back(std::array<idx_t, 2>{{i, j+1}});
      }
      state[i][j+1] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + netmodel[edgmodidx[i][j]]->getNState());
      jstate += netmodel[edgmodidx[i][j]]->getNState();
      stick[i][j+1] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + netmodel[edgmodidx[i][j]]->getNStick());
      jstick += netmodel[edgmodidx[i][j]]->getNStick();
      if (repmodidx[edgmodidx[i][j]]) {
        repidx[edgmodidx[i][j]].push_back(std::array<idx_t, 2>{{i, j+1}});
      }
    }
    // set up auxiliary state
    vtxaux[i].clear();
    if (netmodel[vtxmodidx[i]]->getNAux()) {
      std::vector<std::string> auxstate = netmodel[vtxmodidx[i]]->getAuxState();
      std::vector<std::string> auxstick = netmodel[vtxmodidx[i]]->getAuxStick();
      for (idx_t j = 0; j < edgmodidx[i].size(); ++j) {
        if (edgmodidx[i][j]) {
          vtxaux[i].push_back(aux_t());
          vtxaux[i].back().index = j+1;
          vtxaux[i].back().stateidx.resize(auxstate.size());
          for (std::size_t s = 0; s < auxstate.size(); ++s) {
            vtxaux[i].back().stateidx[s] = netmodel[edgmodidx[i][j]]->getStateIdx(auxstate[s]);
          }
          vtxaux[i].back().stickidx.resize(auxstick.size());
          for (std::size_t s = 0; s < auxstick.size(); ++s) {
            vtxaux[i].back().stickidx[s] = netmodel[edgmodidx[i][j]]->getStateIdx(auxstick[s]);
          }
        }
      }
    }
    // copy over event data
    event_t evtpre;
    for (idx_t e = msg->xevent[i]; e < msg->xevent[i+1]; ++e) {
      evtpre.diffuse = msg->diffuse[e];
      evtpre.type = msg->type[e];
      evtpre.source = msg->source[e];
      evtpre.index = msg->index[e];
      evtpre.data = msg->data[e];
      // Add to event queue or spillover
      if (evtpre.diffuse/tstep < equeue) {
        event[i][(evtpre.diffuse/tstep)%equeue].push_back(evtpre);
      }
      else {
        evtaux[i].push_back(evtpre);
      }
      ++nevent;
    }
    // open ports if needed
    if (netmodel[vtxmodidx[i]]->getNPort()) {
      netmodel[vtxmodidx[i]]->OpenPorts();
    }
  }
  CkAssert(msg->nedg == nadjcy);
  CkAssert(msg->nedg == jmodidx);
  CkAssert(msg->nstate == jstate);
  CkAssert(msg->nstick == jstick);
  CkAssert(msg->nevent == nevent);

  // Cleanup;
  delete msg;

  // Create Groups
  CreateGroup();

  // Print some information
  CkPrintf("  Network part %" PRIidx ":   vtx: %d   edg: %d   adjvtx: %d   adjprt: %" PRIidx "\n",
           prtidx, adjcy.size(), nadjcy, adjmap.size(), nadjprt);

  // Set up timing
  tsim = 0;
  tdisp = 0;
  iter = 0;
  // Set up coordination
  cadjprt[0] = 0;
  cadjprt[1] = 0;
  prtiter = 0;
  // Set up computation
  compidx = 0;
  evalidx = 0;
  ccomp = 0;
  ncomp = 0;
  tcomp = 0;

  // Set up checkpointing
  checkiter = (idx_t) (tcheck/tstep);
  // Set up recordning
  reciter = (idx_t) (trecord/tstep);
#ifdef STACS_WITH_YARP
  // Set up synchronization
  synciter = IDX_T_MAX;
#endif

  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}


/**************************************************************************
* Network Save Data
**************************************************************************/

// Send network partition to Netdata chare array
//
void Network::SaveNetwork() {
  // Build network part message for saving
  mPart *mpart = BuildPart();
  netdata(datidx).SaveNetwork(mpart);
    
  // Start a new cycle (checked data sent)
  //thisProxy(prtidx).CycleNetwork();
  cbcycleprt.send();
}

// Send network partition to Netdata chare array (final)
//
void Network::SaveFinalNetwork() {
  // Close ports as needed
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    if (netmodel[vtxmodidx[i]]->getNPort()) {
      netmodel[vtxmodidx[i]]->ClosePorts();
    }
  }

  // Build network part message for saving
  mPart *mpart = BuildPart();
  netdata(datidx).SaveFinalNetwork(mpart);
}

// Clean up Network chare array
//
void Network::FinalizeNetwork() {
  // Close ports as needed
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    if (netmodel[vtxmodidx[i]]->getNPort()) {
      netmodel[vtxmodidx[i]]->ClosePorts();
    }
  }

  // Move to cleanup of netdata
  netdata(datidx).FinalizeNetwork();
}


/**************************************************************************
* Network Helpers
**************************************************************************/

// Create multicast network groups from local connectivity
//
void Network::CreateGroup() {
  /* Bookkeeping */
  CkVec<CkArrayIndex1D> elems;
  std::vector<idx_t> adjprt;

  // Initialize for counting
  adjprt.resize(npnet, 0);

  // bin edges to source parts
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    idx_t k = 1;
    for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
      while (adjcy[i][j] >= vtxdist[k]) { ++k; }
      ++adjprt[k-1];
    }
  }

  // count adjacent and add to group
  for (idx_t k = 0; k < npnet; ++k) {
    if (adjprt[k] > 0) {
      elems.push_back(CkArrayIndex1D(k));
    }
  }

  // set group sizes
  nadjprt = elems.size();
  
  // make a map of source parts
  srcprt.clear();
  for (idx_t k = 0; k < nadjprt; ++k) {
    srcprt[*(elems[k].data())] = k;
  }

  // make list of target parts
  trgprt.resize(nadjprt);
  for (idx_t k = 0; k < nadjprt; ++k) {
    trgprt[k] = *(elems[k].data());
  }

  // Create Charm++ multicast sections
  netgroup = CProxySection_Network::ckNew(thisProxy, elems.getVec(), nadjprt);
  netgroup.ckSectionDelegate(CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch());
}

// Reset Network state
//
void Network::ResetNetwork() {
  for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
    // Clear events
    for (idx_t j = 0; j < equeue; ++j) {
      event[i][j].clear();
    }
    evtaux[i].clear();
    // Reset vertex models
    netmodel[vtxmodidx[i]]->Reset(state[i][0], stick[i][0]);
  }
}


/**************************************************************************
* Build Messages
**************************************************************************/

// Build part message
//
mPart* Network::BuildPart() {
  /* Bookkeeping */
  idx_t nadjcy;
  idx_t nstate;
  idx_t nstick;
  idx_t nevent;
  
  // Get total size of adjcy
  nadjcy = 0;
  nstate = 0;
  nstick = 0;
  nevent = 0;
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    nadjcy += adjcy[i].size();
    for (std::size_t j = 0; j < adjcy[i].size()+1; ++j) {
      nstate += state[i][j].size();
      nstick += stick[i][j].size();
    }
    for (std::size_t j = 0; j < event[i].size(); ++j) {
      nevent += event[i][j].size();
    }
    nevent += evtaux[i].size();
  }

  // Initialize partition data message
  int msgSize[MSG_Part];
  msgSize[0] = npnet+1;       // vtxdist
  msgSize[1] = adjcy.size();  // vtxmodidx
  msgSize[2] = adjcy.size()*3;// xyz
  msgSize[3] = adjcy.size()+1;// xadj
  msgSize[4] = nadjcy;        // adjcy
  msgSize[5] = nadjcy;        // edgmodidx
  msgSize[6] = nstate;        // state
  msgSize[7] = nstick;        // stick
  msgSize[8] = adjcy.size()+1;// xevent
  msgSize[9] = nevent;        // diffuse
  msgSize[10] = nevent;       // type
  msgSize[11] = nevent;       // source
  msgSize[12] = nevent;       // index
  msgSize[13] = nevent;       // data
  mPart *mpart = new(msgSize, 0) mPart;

  // Data sizes
  mpart->nvtx = adjcy.size();
  mpart->nedg = nadjcy;
  mpart->nstate = nstate;
  mpart->nstick = nstick;
  mpart->nevent = nevent;
  mpart->prtidx = prtidx;

  // Graph Information
  for (idx_t i = 0; i < npnet+1; ++i) {
    // vtxdist
    mpart->vtxdist[i] = vtxdist[i];
  }

  // Vertex and Edge Information
  idx_t jstate = 0;
  idx_t jstick = 0;
  idx_t jevent = 0;
  mpart->xadj[0] = 0;
  mpart->xevent[0] = 0;
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    // vtxmodidx
    mpart->vtxmodidx[i] = vtxmodidx[i];
    // xyz
    mpart->xyz[i*3+0] = xyz[i*3+0];
    mpart->xyz[i*3+1] = xyz[i*3+1];
    mpart->xyz[i*3+2] = xyz[i*3+2];
    // vertex state
    for (std::size_t s = 0; s < state[i][0].size(); ++s) {
      mpart->state[jstate++] = state[i][0][s];
    }
    for (std::size_t s = 0; s < stick[i][0].size(); ++s) {
      mpart->stick[jstick++] = stick[i][0][s];
    }

    // xadj
    mpart->xadj[i+1] = mpart->xadj[i] + adjcy[i].size();
    for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
      // adjcy
      mpart->adjcy[mpart->xadj[i] + j] = adjcy[i][j];
      // edgmodidx
      mpart->edgmodidx[mpart->xadj[i] + j] = edgmodidx[i][j];
      // state
      for (std::size_t s = 0; s < state[i][j+1].size(); ++s) {
        mpart->state[jstate++] = state[i][j+1][s];
      }
      for (std::size_t s = 0; s < stick[i][j+1].size(); ++s) {
        mpart->stick[jstick++] = stick[i][j+1][s];
      }
    }

    // events
    for (std::size_t j = 0; j < event[i].size(); ++j) {
      for (std::size_t e = 0; e < event[i][j].size(); ++e) {
        mpart->diffuse[jevent] = event[i][j][e].diffuse - tsim;
        mpart->type[jevent] = event[i][j][e].type;
        mpart->source[jevent] = event[i][j][e].source;
        mpart->index[jevent] = event[i][j][e].index;
        mpart->data[jevent++] = event[i][j][e].data;
      }
    }
    // events spillover
    for (std::size_t e = 0; e < evtaux[i].size(); ++e) {
      mpart->diffuse[jevent] = evtaux[i][e].diffuse - tsim;
      mpart->type[jevent] = evtaux[i][e].type;
      mpart->source[jevent] = evtaux[i][e].source;
      mpart->index[jevent] = evtaux[i][e].index;
      mpart->data[jevent++] = evtaux[i][e].data;
    }
    // xevent
    mpart->xevent[i+1] = jevent;
  }
  CkAssert(jstate == nstate);
  CkAssert(jstick == nstick);
  CkAssert(jevent == nevent);

  return mpart;
}


/**************************************************************************
* Charm++ Definitions
**************************************************************************/
#include "network.def.h"
