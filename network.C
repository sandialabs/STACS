/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"
#include "pup_stl.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ CkGroupID mCastGrpId;
extern /*readonly*/ idx_t npdat;
extern /*readonly*/ idx_t npnet;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ tick_t tcheck;


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

  // Vertex Models
  for (std::size_t i = 0; i < netmodel.size(); ++i) {
    delete netmodel[i];
  }
  // Set up containers
  netmodel.clear();
  // "none" model
  netmodel.push_back(NetModelFactory::getNetModel()->Create(0));
  // User defined models
  for (idx_t i = 1; i < msg->nmodel+1; ++i) {
    // Create model object
    netmodel.push_back(NetModelFactory::getNetModel()->Create(msg->modtype[i-1]));
    CkAssert(netmodel[i]->getNState() == msg->nstate[i-1]);
    CkAssert(netmodel[i]->getNStick() == msg->nstick[i-1]);
    CkAssert(netmodel[i]->getNParam() == msg->xparam[i] - msg->xparam[i-1]);
    netmodel[i]->setParam(msg->param + msg->xparam[i-1]);

    // Print out model information
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
  }
  delete msg;

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

// Coordination with NetData chare array
//
void Network::LoadNetwork(CProxy_NetData cpdat) {
  // Save proxy
  netdata = cpdat;
  
  // Request network part from input
  CkCallback *cb = new CkCallback(CkIndex_Network::LoadNetwork(NULL), thisIndex, thisProxy);
  netdata(datidx).LoadNetwork(prtidx, *cb);
}

// Receive data from NetData chare array
//
void Network::LoadNetwork(mPart *msg) {
  /* Bookkeeping */
  idx_t nadjcy;  // edges
  idx_t jmodidx; // modidx
  idx_t jstate;  // state
  idx_t jstick;  // time state

  // Copy over part data
  CkAssert(msg->prtidx == prtidx);

  // Setup data vectors
  vtxdist.resize(npnet+1);
  vtxidx.resize(msg->nvtx);
  vtxmodidx.resize(msg->nvtx);
  xyz.resize(msg->nvtx*3);
  adjcy.resize(msg->nvtx);
  adjmap.clear();
  edgmodidx.resize(msg->nvtx);
  state.resize(msg->nvtx);
  stick.resize(msg->nvtx);
  
  // Graph distribution information
  for (idx_t i = 0; i < npnet+1; ++i) {
    // vtxdist
    vtxdist[i] = msg->vtxdist[i];
  }

  // Initalize counters
  nadjcy = 0;
  jmodidx = 0;
  jstate = 0; jstick = 0;
  // Get adjacency matrix
  for (idx_t i = 0; i < adjcy.size(); ++i) {
    vtxidx[i] = vtxdist[prtidx] + i;
    vtxmodidx[i] = msg->vtxmodidx[i];
    xyz[i*3+0] = msg->xyz[i*3+0];
    xyz[i*3+1] = msg->xyz[i*3+1];
    xyz[i*3+2] = msg->xyz[i*3+2];
    // preallocate sizes
    adjcy[i].resize(msg->xadj[i+1] - msg->xadj[i]);
    edgmodidx[i].resize(msg->xadj[i+1] - msg->xadj[i]);
    state[i].resize(msg->xadj[i+1] - msg->xadj[i] + 1);
    stick[i].resize(msg->xadj[i+1] - msg->xadj[i] + 1);
    nadjcy += adjcy[i].size();
    // copy over vertex data
    state[i][0] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + netmodel[vtxmodidx[i]]->getNState());
    jstate += netmodel[vtxmodidx[i]]->getNState();
    stick[i][0] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + netmodel[vtxmodidx[i]]->getNStick());
    jstick += netmodel[vtxmodidx[i]]->getNStick();
    // copy over edge data
    for (idx_t j = 0; j < adjcy[i].size(); ++j) {
      adjcy[i][j] = msg->adjcy[msg->xadj[i] + j];
      adjmap[adjcy[i][j]].push_back(std::array<idx_t, 2>{{i, j}});
      // models
      edgmodidx[i][j] = msg->edgmodidx[jmodidx++];
      state[i][j+1] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + netmodel[edgmodidx[i][j]]->getNState());
      jstate += netmodel[edgmodidx[i][j]]->getNState();
      stick[i][j+1] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + netmodel[edgmodidx[i][j]]->getNStick());
      jstick += netmodel[edgmodidx[i][j]]->getNStick();
    }
  }
  CkAssert(msg->nedg == nadjcy);
  CkAssert(msg->nedg == jmodidx);
  CkAssert(msg->nstate == jstate);
  CkAssert(msg->nstick == jstick);

  // Cleanup;
  delete msg;

  // Create Groups
  CreateGroup();

  // Print some information
  CkPrintf("  Network part %" PRIidx ": vtx: %d   edg: %d   adjprt: %" PRIidx "   adjmap: %d\n",
           prtidx, adjcy.size(), nadjcy, nadjprt, adjmap.size());

  // Set up timing
  tsim = 0;
  iter = 0;
  // Set up coordination
  cadjprt[0] = 0;
  cadjprt[1] = 0;
  prtiter = 0;
  // Set up checkpointing
  cpflag = false;
  checkiter = (idx_t) (tcheck/tstep);
#ifdef STACS_WITH_YARP
  // Set up synchronization
  synciter = IDX_T_MAX;
#endif

  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

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


/**************************************************************************
* Network Save Data
**************************************************************************/

// Send network partition to NetData chare array
//
void Network::SaveNetwork() {
  /* Bookkeeping */
  idx_t nadjcy;
  idx_t nstate;
  idx_t nstick;
  
  // Get total size of adjcy
  nadjcy = 0;
  nstate = 0;
  nstick = 0;
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    nadjcy += adjcy[i].size();
    for (std::size_t j = 0; j < adjcy[i].size()+1; ++j) {
      nstate += state[i][j].size();
      nstick += stick[i][j].size();
    }
  }

  // Initialize partition data message
  int msgSize[MSG_Part];
  msgSize[0] = npnet+1;       // vtxdist
  msgSize[1] = adjcy.size();  // vtxmodidx
  msgSize[2] = adjcy.size()*3;// xyz
  msgSize[3] = adjcy.size();  // xadj
  msgSize[4] = nadjcy;        // adjcy
  msgSize[5] = nadjcy;        // edgmodidx
  msgSize[6] = nstate;        // state
  msgSize[7] = nstick;        // stick
  mPart *part = new(msgSize, 0) mPart;

  // Data sizes
  part->nvtx = adjcy.size();
  part->nedg = nadjcy;
  part->nstate = nstate;
  part->nstick = nstick;
  part->prtidx = prtidx;

  // Graph Information
  for (idx_t i = 0; i < npnet+1; ++i) {
    // vtxdist
    part->vtxdist[i] = vtxdist[i];
  }

  // Vertex and Edge Information
  idx_t jstate = 0;
  idx_t jstick = 0;
  part->xadj[0] = 0;
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    // vtxmodidx
    part->vtxmodidx[i] = vtxmodidx[i];
    // xyz
    part->xyz[i*3+0] = xyz[i*3+0];
    part->xyz[i*3+1] = xyz[i*3+1];
    part->xyz[i*3+2] = xyz[i*3+2];
    // vertex state
    for (std::size_t s = 0; s < state[i][0].size(); ++s) {
      part->state[jstate++] = state[i][0][s];
    }
    for (std::size_t s = 0; s < stick[i][0].size(); ++s) {
      part->stick[jstick++] = stick[i][0][s];
    }

    // xadj
    part->xadj[i+1] = part->xadj[i] + adjcy[i].size();
    for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
      // adjcy
      part->adjcy[part->xadj[i] + j] = adjcy[i][j];
      // edgmodidx
      part->edgmodidx[part->xadj[i] + j] = edgmodidx[i][j];
      // state
      for (std::size_t s = 0; s < state[i][j+1].size(); ++s) {
        part->state[jstate++] = state[i][j+1][s];
      }
      for (std::size_t s = 0; s < stick[i][j+1].size(); ++s) {
        part->stick[jstick++] = stick[i][j+1][s];
      }
    }
  }
  CkAssert(jstate == nstate);
  CkAssert(jstick == nstick);

  // Send network part to output file
  if (cpflag) {
    cpflag = false;
    netdata(datidx).CheckNetwork(part);
    
    // Start a new cycle (checked data sent)
    thisProxy(prtidx).Cycle();
  }
  else {
    netdata(datidx).SaveNetwork(part);
  }
}


/**************************************************************************
* Network Simulation
**************************************************************************/

// Receive go-ahead message
//
void Network::GoAhead(mGo *msg) {
  // Increment coordination
  ++cadjprt[(prtiter + (msg->iter - iter))%2];
  delete msg;

  // Start next cycle
  if (cadjprt[prtiter] == nadjprt) {
    // Bookkeepping
    cadjprt[prtiter] = 0;
    prtiter = (prtiter+1)%2;

    // Increment iteration
    ++iter;

    // Start a new cycle
    thisProxy(prtidx).Cycle();
  }
}

// Network Simulation Cycle (control flow)
//
void Network::Cycle() {
  // Check if simulation time is complete
  if (tsim >= tmax) {
    // return control to main
    contribute(0, NULL, CkReduction::nop);
  }
#ifdef STACS_WITH_YARP
  // Synchronization from RPC
  else if (iter == synciter) {
    // Bookkkeeping
    synciter = IDX_T_MAX;

    // Display synchronization information
    if (prtidx == 0) {
      CkPrintf("  Synchronized at iteration %" PRIidx "\n", iter);
    }

    // move control to sychronization callback
    if (cpflag) {
      thisProxy(prtidx).SaveNetwork();
    }
    contribute(0, NULL, CkReduction::nop);
  }
#endif
  // Checkpointing
  else if (iter == checkiter) {
    // Bookkeeping
    checkiter = checkiter + (idx_t) (tcheck/tstep);

    // Checkpoint
    cpflag = true;
    thisProxy(prtidx).SaveNetwork();
  }
  // Simulate next cycle
  else {
    // Display iteration information
    if (prtidx == 0) {
      CkPrintf("  Simulating iteration %" PRIidx "\n", iter);
    }

    // Perform computation
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
      tick_t tstart = 0;
      tick_t tdrift;
      tdrift = tstep - tstart;
      netmodel[vtxmodidx[i]]->Step(tdrift);
    }

    // Increment simulated time
    tsim += tstep;

    // Send messages to neighbors
    mGo *go = new mGo(iter);
    netgroup.GoAhead(go);

    // Send messages to netdata (file)

#ifdef STACS_WITH_YARP
    // Send messages to streams (yarp)
#endif
  }
}

/**************************************************************************
* Charm++ Definitions
**************************************************************************/
#include "network.def.h"
