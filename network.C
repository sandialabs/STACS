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
extern /*readonly*/ unsigned randseed;
extern /*readonly*/ int netparts;
extern /*readonly*/ int netfiles;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ idx_t nevtday;
extern /*readonly*/ idx_t intdisp;
extern /*readonly*/ idx_t intrec;
extern /*readonly*/ idx_t intbal;
extern /*readonly*/ idx_t intsave;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ idx_t grpvtxmin;


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
  int datdiv = netparts/netfiles;
  int datrem = netparts%netfiles;
  datidx = (prtidx-( ((prtidx+1)/(datdiv+1))>datrem ? datrem:((prtidx+1)/(datdiv+1)) ))/datdiv;
  
  // Set up random number generator
  rngine.seed(randseed+prtidx);
  unifdist = new std::uniform_real_distribution<real_t> (0.0, 1.0);
  
  // Simulation configuration
  plastic = msg->plastic;
  episodic = msg->episodic;
  loadbal = msg->loadbal;

  // Network Models
  for (std::size_t i = 0; i < model.size(); ++i) {
    delete model[i];
  }
  // Set up containers
  model.clear();
  // "none" model
  model.push_back(ModelFactory::newModel()->Create(0));
  // TODO: Roll this into model creation perhaps
  model[0]->setActive(false);
  model[0]->setMother(false);
  model[0]->setAnchor(false);
  model[0]->setPlastic(false);

  idx_t jparamname = 0;
  // User defined models
  for (idx_t i = 1; i < msg->nmodel+1; ++i) {
    // Create model object
    model.push_back(ModelFactory::newModel()->Create(msg->modtype[i-1]));
    // TODO: these asserts will no longer work when models are underspecified
    //CkAssert(model[i]->getNState() == msg->nstate[i-1]);
    //CkAssert(model[i]->getNStick() == msg->nstick[i-1]);
    //CkAssert(model[i]->getNParam() == msg->xparam[i] - msg->xparam[i-1]);
    //model[i]->setParam(msg->param + msg->xparam[i-1]);
    model[i]->setPort(msg->port + msg->xport[i-1]);
    model[i]->setRandom(unifdist, &rngine);
    model[i]->setActive(msg->grpactive[i-1]);
    model[i]->setMother(msg->grpmother[i-1]);
    model[i]->setAnchor(msg->grpanchor[i-1]);
    model[i]->setPlastic(msg->plastic);
    // names (may be in a different order than implemented model)
    // Find the mapping from user-provided param names to the implemented param names
    // Find which paramss were not specified (and will need model-supplied defaults)
    std::vector<idx_t> parammap;
    parammap.resize(msg->nparam[i-1]);
    std::vector<real_t> paramvalues;
    paramvalues.resize(model[i]->getNParam());
    std::vector<bool> paramconfig;
    paramconfig.resize(model[i]->getNParam(),false);
    for (std::size_t j = 0; j < msg->nparam[i-1]; ++j) {
      std::string paramname = std::string(msg->paramname + msg->xparamname[jparamname], msg->paramname + msg->xparamname[jparamname+1]);
      parammap[j] = model[i]->getParamIdx(paramname.c_str());
      // some basic error checking
      if (parammap[j] == -1) {
        CkPrintf("  param name: %s is invalid for model: %" PRIidx "\n", paramname.c_str(), model[i]->getModType());
        CkExit();
      }
      paramvalues[parammap[j]] = msg->param[msg->xparam[i-1] + j];
      paramconfig[parammap[j]] = true;
      ++jparamname;
    }
    // Now go through the parameters that weren't defined by the model config
    // These are the false entries in paramconfig
    for (std::size_t j = 0; j < model[i]->getNParam(); ++j) {
      if (paramconfig[j]) { continue; }
      else {
        paramvalues[j] = model[i]->getDefaultParam(j);
      }
    }
    // Set the parameter values after config and defaults populated
    model[i]->setParam(paramvalues.data());

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
      CkPrintf("  Network model: %" PRIidx "   ModType: %" PRIidx "   Params:%s\n", i, model[i]->getModType(), params.c_str());
    }
    */
  }

  // Recording
  evtlog.clear();
  evtloglist.resize(EVENT_TOTAL);
  record.clear();
  recordlist.clear();
  // TODO: Put this in the yml config file
  // TODO: Actually, fold this into a recording model
  //       that combines into a multi-vertex thing
  // Record spikes
  evtloglist[EVENT_SPIKE] = true;
  evtloglist[EVENT_CHGRATE] = true;

  delete msg;

  // set up auxiliary state
  edgaux.resize(model.size());
  for (std::size_t j = 1; j < model.size(); ++j) {
    if (model[j]->getNAux()) {
      edgaux[j].resize(model.size());
      std::vector<std::string> auxstate = model[j]->getAuxState();
      std::vector<std::string> auxstick = model[j]->getAuxStick();
      for (std::size_t i = 1; i < model.size(); ++i) {
        edgaux[j][i].resize(1);
        edgaux[j][i][0].index = 0;
        edgaux[j][i][0].stateidx.resize(auxstate.size());
        for (std::size_t s = 0; s < auxstate.size(); ++s) {
          edgaux[j][i][0].stateidx[s] = model[i]->getStateIdx(auxstate[s]);
        }
        edgaux[j][i][0].stickidx.resize(auxstick.size());
        for (std::size_t s = 0; s < auxstick.size(); ++s) {
          edgaux[j][i][0].stickidx[s] = model[i]->getStickIdx(auxstick[s]);
        }
      }
    }
  }
  // set up periodic events
  leapevt.clear();
  leaplist.resize(model.size(), false);
  leapidx.resize(model.size());
  events.clear();
  for (std::size_t n = 1; n < model.size(); ++n) {
    model[n]->getLeap(events);
    if (events.size()) {
      leaplist[n] = true;
      for (std::size_t e = 0; e < events.size(); ++e) {
        events[e].source = (idx_t) n;
        leapevt.push_back(events[e]);
      }
      events.clear();
    }
  }
  if (leapevt.size()) {
    std::sort(leapevt.begin(), leapevt.end());
    tleap = leapevt.front().diffuse;
  }
  else {
    tleap = TICK_T_MAX; // never leap
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
  for (std::size_t i = 0; i < model.size(); ++i) {
    delete model[i];
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
  vtxdist.resize(netparts+1);
  vtxidx.resize(msg->nvtx);
  vtxmap.clear();
  vtxmodidx.resize(msg->nvtx);
  xyz.resize(msg->nvtx*3);
  adjcy.resize(msg->nvtx);
  adjmap.clear();
  edgmodidx.resize(msg->nvtx);
  state.resize(msg->nvtx);
  stick.resize(msg->nvtx);
  vtxaux.resize(msg->nvtx);
  evtcal.resize(msg->nvtx);
  evtcol.resize(msg->nvtx);
  events.clear();
  evtext.clear();
  evtrpc.clear();
  // Polychronization
  grpstamps.resize(msg->nvtx);
  grpdur.resize(msg->nvtx);
  grpmap.clear();
  grpwindow.resize(msg->nvtx);
  grplog.clear();
  grpseeds.clear();
  grptraces.resize(msg->nvtx);
  grpleg.clear();
  grproute.clear();
  grproutes.clear();

  // Graph distribution information
  for (int i = 0; i < netparts+1; ++i) {
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
    // preallocate sizes
    adjcy[i].resize(msg->xadj[i+1] - msg->xadj[i]);
    edgmodidx[i].resize(msg->xadj[i+1] - msg->xadj[i]);
    state[i].resize(msg->xadj[i+1] - msg->xadj[i] + 1);
    stick[i].resize(msg->xadj[i+1] - msg->xadj[i] + 1);
    evtcal[i].resize(nevtday);
    nadjcy += adjcy[i].size();
    // copy over vertex data
    state[i][0] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + model[vtxmodidx[i]]->getNState());
    jstate += model[vtxmodidx[i]]->getNState();
    stick[i][0] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + model[vtxmodidx[i]]->getNStick());
    jstick += model[vtxmodidx[i]]->getNStick();
    if (leaplist[vtxmodidx[i]]) {
      leapidx[vtxmodidx[i]].push_back(std::array<idx_t, 2>{{i, 0}});
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
      state[i][j+1] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + model[edgmodidx[i][j]]->getNState());
      jstate += model[edgmodidx[i][j]]->getNState();
      stick[i][j+1] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + model[edgmodidx[i][j]]->getNStick());
      jstick += model[edgmodidx[i][j]]->getNStick();
      if (leaplist[edgmodidx[i][j]]) {
        leapidx[edgmodidx[i][j]].push_back(std::array<idx_t, 2>{{i, j+1}});
      }
    }
    // set up auxiliary state
    vtxaux[i].clear();
    if (model[vtxmodidx[i]]->getNAux()) {
      std::vector<std::string> auxstate = model[vtxmodidx[i]]->getAuxState();
      std::vector<std::string> auxstick = model[vtxmodidx[i]]->getAuxStick();
      for (idx_t j = 0; j < edgmodidx[i].size(); ++j) {
        if (edgmodidx[i][j]) {
          vtxaux[i].push_back(auxidx_t());
          vtxaux[i].back().index = j+1;
          vtxaux[i].back().stateidx.resize(auxstate.size());
          for (std::size_t s = 0; s < auxstate.size(); ++s) {
            vtxaux[i].back().stateidx[s] = model[edgmodidx[i][j]]->getStateIdx(auxstate[s]);
          }
          vtxaux[i].back().stickidx.resize(auxstick.size());
          for (std::size_t s = 0; s < auxstick.size(); ++s) {
            vtxaux[i].back().stickidx[s] = model[edgmodidx[i][j]]->getStateIdx(auxstick[s]);
          }
        }
      }
    }
    // copy over event data
    event_t event;
    for (idx_t e = msg->xevent[i]; e < msg->xevent[i+1]; ++e) {
      event.diffuse = msg->diffuse[e];
      event.type = msg->type[e];
      event.source = msg->source[e];
      event.index = msg->index[e];
      event.data = msg->data[e];
      // Add to event queue or spillover
      if (event.diffuse/tstep < nevtday) {
        evtcal[i][(event.diffuse/tstep)%nevtday].push_back(event);
      }
      else {
        evtcol[i].push_back(event);
      }
      ++nevent;
    }
    // open ports if needed
    if (model[vtxmodidx[i]]->getNPort()) {
      model[vtxmodidx[i]]->OpenPorts();
    }
    // initialize polychronization
    grpstamps[i].clear();
    grpdur[i].clear();
    grpwindow[i].clear();
    ReadGroup(i);
  }
  CkAssert(msg->nedg == nadjcy);
  CkAssert(msg->nedg == jmodidx);
  CkAssert(msg->nstate == jstate);
  CkAssert(msg->nstick == jstick);
  CkAssert(msg->nevent == nevent);

  // Cleanup;
  delete msg;

  // Multicast Communication 
  CreateComm();

  // Print some information
  CkPrintf("  Network part %" PRIidx ":   vtx: %d   edg: %d   adjvtx: %d   adjpart: %" PRIidx "\n",
           prtidx, adjcy.size(), nadjcy, adjmap.size(), nadjpart);

  // Set up timing
  tsim = 0;
  iter = 0;
  commiter = 0;
  dispiter = 0;
  reciter = intrec;
  baliter = loadbal ? intbal : IDX_T_MAX;
  saveiter = plastic ? intsave : IDX_T_MAX;
  teps = episodic ? 0 : TICK_T_MAX;
  epsidx = -1;
  // Set up coordination
  cadjpart[0] = 0;
  cadjpart[1] = 0;
  partiter = 0;

#ifdef STACS_WITH_YARP
  // Set up synchronization
  synciter = IDX_T_MAX;
  syncing = false;
#endif

  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// Receive data from Netdata chare array
//
void Network::ReloadNetwork(mPart *msg) {
  /* Bookkeeping */
  idx_t nadjcy;  // edges
  idx_t jmodidx; // modidx
  idx_t jstate;  // state
  idx_t jstick;  // time state
  idx_t nevent;  // events

  // Copy over part data
  CkAssert(msg->prtidx == prtidx);

  // Setup data vectors
  vtxdist.resize(netparts+1);
  vtxidx.resize(msg->nvtx);
  vtxmap.clear();
  vtxmodidx.resize(msg->nvtx);
  xyz.resize(msg->nvtx*3);
  adjcy.resize(msg->nvtx);
  adjmap.clear();
  edgmodidx.resize(msg->nvtx);
  state.resize(msg->nvtx);
  stick.resize(msg->nvtx);
  vtxaux.resize(msg->nvtx);
  evtcal.clear();
  evtcal.resize(msg->nvtx);
  evtcol.clear();
  evtcol.resize(msg->nvtx);
  events.clear();
  evtext.clear();
  evtrpc.clear();

  // Graph distribution information
  for (int i = 0; i < netparts+1; ++i) {
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
    // preallocate sizes
    adjcy[i].resize(msg->xadj[i+1] - msg->xadj[i]);
    edgmodidx[i].resize(msg->xadj[i+1] - msg->xadj[i]);
    state[i].resize(msg->xadj[i+1] - msg->xadj[i] + 1);
    stick[i].resize(msg->xadj[i+1] - msg->xadj[i] + 1);
    evtcal[i].resize(nevtday);
    nadjcy += adjcy[i].size();
    // copy over vertex data
    state[i][0] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + model[vtxmodidx[i]]->getNState());
    jstate += model[vtxmodidx[i]]->getNState();
    stick[i][0] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + model[vtxmodidx[i]]->getNStick());
    jstick += model[vtxmodidx[i]]->getNStick();
    if (leaplist[vtxmodidx[i]]) {
      leapidx[vtxmodidx[i]].push_back(std::array<idx_t, 2>{{i, 0}});
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
      state[i][j+1] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + model[edgmodidx[i][j]]->getNState());
      jstate += model[edgmodidx[i][j]]->getNState();
      stick[i][j+1] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + model[edgmodidx[i][j]]->getNStick());
      jstick += model[edgmodidx[i][j]]->getNStick();
      if (leaplist[edgmodidx[i][j]]) {
        leapidx[edgmodidx[i][j]].push_back(std::array<idx_t, 2>{{i, j+1}});
      }
    }
    // set up auxiliary state
    vtxaux[i].clear();
    if (model[vtxmodidx[i]]->getNAux()) {
      std::vector<std::string> auxstate = model[vtxmodidx[i]]->getAuxState();
      std::vector<std::string> auxstick = model[vtxmodidx[i]]->getAuxStick();
      for (idx_t j = 0; j < edgmodidx[i].size(); ++j) {
        if (edgmodidx[i][j]) {
          vtxaux[i].push_back(auxidx_t());
          vtxaux[i].back().index = j+1;
          vtxaux[i].back().stateidx.resize(auxstate.size());
          for (std::size_t s = 0; s < auxstate.size(); ++s) {
            vtxaux[i].back().stateidx[s] = model[edgmodidx[i][j]]->getStateIdx(auxstate[s]);
          }
          vtxaux[i].back().stickidx.resize(auxstick.size());
          for (std::size_t s = 0; s < auxstick.size(); ++s) {
            vtxaux[i].back().stickidx[s] = model[edgmodidx[i][j]]->getStateIdx(auxstick[s]);
          }
        }
      }
    }
    // copy over event data
    event_t event;
    for (idx_t e = msg->xevent[i]; e < msg->xevent[i+1]; ++e) {
      event.diffuse = msg->diffuse[e];
      event.type = msg->type[e];
      event.source = msg->source[e];
      event.index = msg->index[e];
      event.data = msg->data[e];
      // Add to event queue or spillover
      if ((event.diffuse/tstep - iter) < nevtday) {
        evtcal[i][(event.diffuse/tstep)%nevtday].push_back(event);
      }
      else {
        evtcol[i].push_back(event);
      }
      ++nevent;
    }
  }
  CkAssert(msg->nedg == nadjcy);
  CkAssert(msg->nedg == jmodidx);
  CkAssert(msg->nstate == jstate);
  CkAssert(msg->nstick == jstick);
  CkAssert(msg->nevent == nevent);

  // Cleanup;
  delete msg;

  // Print some information
  CkPrintf("  Network part %" PRIidx ":   vtx: %d   edg: %d   adjvtx: %d   adjpart: %" PRIidx "\n",
           prtidx, adjcy.size(), nadjcy, adjmap.size(), nadjpart);

  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}


/**************************************************************************
* Network Save Data
**************************************************************************/


// Repartition from network
//
void Network::RebalNetwork() {
  // Build network part message for saving
  mPart *mpart = BuildRepart();
  netdata(datidx).LoadRepart(mpart);
}

// Send network partition to Netdata chare array
//
void Network::SaveNetwork() {
  // Build network part message for saving
  mPart *mpart = BuildPart();
  netdata(datidx).SaveNetwork(mpart);
    
  // Start a new cycle (checked data sent)
  cyclepart.send();
}

// Send network partition to Netdata chare array (final)
//
void Network::SaveCloseNetwork() {
  // Close ports as needed
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    if (model[vtxmodidx[i]]->getNPort()) {
      model[vtxmodidx[i]]->ClosePorts();
    }
  }

  // Build network part message for saving
  mPart *mpart = BuildPart();
  netdata(datidx).SaveCloseNetwork(mpart);
}

// Clean up Network chare array
//
void Network::CloseNetwork() {
  // Close ports as needed
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    if (model[vtxmodidx[i]]->getNPort()) {
      model[vtxmodidx[i]]->ClosePorts();
    }
  }

  // Move to cleanup of netdata
  netdata(datidx).CloseNetwork();
}


/**************************************************************************
* Network Helpers
**************************************************************************/

// Create multicast network groups from local connectivity
//
void Network::CreateComm() {
  /* Bookkeeping */
  CkVec<CkArrayIndex1D> elems;
  std::vector<int> adjpart;

  // Initialize for counting
  adjpart.resize(netparts, 0);

  // bin edges to source parts
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    int k = 1;
    for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
      while (adjcy[i][j] >= vtxdist[k]) { ++k; }
      ++adjpart[k-1];
    }
  }

  // count adjacent and add to group
  for (int k = 0; k < netparts; ++k) {
    if (adjpart[k] > 0) {
      elems.push_back(CkArrayIndex1D(k));
    }
  }

  // set group sizes
  nadjpart = elems.size();
  
  // make a map of source parts
  srcpart.clear();
  for (int k = 0; k < nadjpart; ++k) {
    srcpart[*(elems[k].data())] = k;
  }

  // make list of target parts
  trgpart.resize(nadjpart);
  for (int k = 0; k < nadjpart; ++k) {
    trgpart[k] = *(elems[k].data());
  }

  // Create Charm++ multicast sections
  netcomm = CProxySection_Network::ckNew(thisProxy, elems.getVec(), nadjpart);
  netcomm.ckSectionDelegate(CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch());
}

// Reset Network state
//
void Network::ResetNetwork() {
  // Reset models
  for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
    // Clear events
    for (idx_t j = 0; j < nevtday; ++j) {
      evtcal[i][j].clear();
    }
    evtcol[i].clear();
    // Reset vertices
    model[vtxmodidx[i]]->Reset(state[i][0], stick[i][0]);
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
    for (std::size_t j = 0; j < evtcal[i].size(); ++j) {
      nevent += evtcal[i][j].size();
    }
    nevent += evtcol[i].size();
  }

  // Initialize partition data message
  int msgSize[MSG_Part];
  msgSize[0] = netparts+1;    // vtxdist
  msgSize[1] = 0;             // vtxidx (implicit)
  msgSize[2] = adjcy.size();  // vtxmodidx
  msgSize[3] = adjcy.size()*3;// xyz
  msgSize[4] = adjcy.size()+1;// xadj
  msgSize[5] = nadjcy;        // adjcy
  msgSize[6] = nadjcy;        // edgmodidx
  msgSize[7] = nstate;        // state
  msgSize[8] = nstick;        // stick
  msgSize[9] = adjcy.size()+1;// xevent
  msgSize[10] = nevent;        // diffuse
  msgSize[11] = nevent;       // type
  msgSize[12] = nevent;       // source
  msgSize[13] = nevent;       // index
  msgSize[14] = nevent;       // data
  mPart *mpart = new(msgSize, 0) mPart;

  // Data sizes
  mpart->nvtx = adjcy.size();
  mpart->nedg = nadjcy;
  mpart->nstate = nstate;
  mpart->nstick = nstick;
  mpart->nevent = nevent;
  mpart->prtidx = prtidx;

  // Graph Information
  for (int i = 0; i < netparts+1; ++i) {
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
    for (std::size_t j = 0; j < evtcal[i].size(); ++j) {
      for (std::size_t e = 0; e < evtcal[i][j].size(); ++e) {
        mpart->diffuse[jevent] = evtcal[i][j][e].diffuse - tsim;
        mpart->type[jevent] = evtcal[i][j][e].type;
        mpart->source[jevent] = evtcal[i][j][e].source;
        mpart->index[jevent] = evtcal[i][j][e].index;
        mpart->data[jevent++] = evtcal[i][j][e].data;
      }
    }
    // events spillover
    for (std::size_t e = 0; e < evtcol[i].size(); ++e) {
      mpart->diffuse[jevent] = evtcol[i][e].diffuse - tsim;
      mpart->type[jevent] = evtcol[i][e].type;
      mpart->source[jevent] = evtcol[i][e].source;
      mpart->index[jevent] = evtcol[i][e].index;
      mpart->data[jevent++] = evtcol[i][e].data;
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
