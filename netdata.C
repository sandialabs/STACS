/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "stacs.h"
#include "network.h"
#include "pup_stl.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ CProxy_Main mainProxy;
extern /*readonly*/ idx_t npdat;
extern /*readonly*/ idx_t npnet;


/**************************************************************************
* Reduction for network distribution
**************************************************************************/

CkReduction::reducerType net_dist;
/*initnode*/
void registerNetDist(void) {
  net_dist = CkReduction::addReducer(netDist);
}

CkReductionMsg *netDist(int nMsg, CkReductionMsg **msgs) {
  std::vector<dist_t> ret;
  ret.clear();
  for (int i = 0; i < nMsg; i++) {
    for (std::size_t j = 0; j < msgs[i]->getSize()/sizeof(dist_t); ++j) {
      // Extract data and reduce 
      ret.push_back(*((dist_t *)msgs[i]->getData() + j));
    }
  }
  return CkReductionMsg::buildNew(ret.size()*sizeof(dist_t), ret.data());
}


/**************************************************************************
* Network Data
**************************************************************************/

// NetData constructor
//
NetData::NetData(mDist *msg) {
  // Bookkeeping
  datidx = thisIndex;
  idx_t ndiv = npnet/npdat;
  idx_t nrem = npnet%npdat;
  cprt = rprt = 0;
  nprt = ndiv + (datidx < nrem);
  xprt = datidx*ndiv + (datidx < nrem ? datidx : nrem);

  // Persistence
  vtxdist.resize(npnet+1);
  edgdist.resize(npnet+1);
  statedist.resize(npnet+1);
  stickdist.resize(npnet+1);
  eventdist.resize(npnet+1);
  for (idx_t i = 0; i < npnet+1; ++i) {
    vtxdist[i] = msg->vtxdist[i];
    edgdist[i] = msg->edgdist[i];
    statedist[i] = msg->statedist[i];
    stickdist[i] = msg->stickdist[i];
    eventdist[i] = msg->eventdist[i];
  }
  // Models
  for (std::size_t i = 0; i < netmodel.size(); ++i) {
    delete netmodel[i];
  }
  // Set up containers
  netmodel.clear();
  modname.resize(msg->nmodel+1);
  // "none" model
  netmodel.push_back(NetModelFactory::getNetModel()->Create(0));
  modname[0] = std::string("none");
  modmap[modname[0]] = 0;
  // User defined models
  for (idx_t i = 1; i < msg->nmodel+1; ++i) {
    netmodel.push_back(NetModelFactory::getNetModel()->Create(msg->modtype[i-1]));
    modname[i] = std::string(msg->modname + msg->xmodname[i-1], msg->modname + msg->xmodname[i]);
    modmap[modname[i]] = i;
    if (datidx == 0) {
      CkPrintf("  NetData model: %" PRIidx "   NStates: %" PRIidx "   Name: %s\n", i, netmodel[i]->getNState(), modname[i].c_str());
    }
  }
  delete msg;

  // Data
  parts.resize(nprt);
  records.resize(nprt);

  // Read in files
  CkPrintf("  Loading input %" PRIidx "\n", datidx);
  ReadCSR();

#ifdef STACS_WITH_YARP
  // Open yarp
  yarp.init();
#endif

  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// NetData migration
//
NetData::NetData(CkMigrateMessage *msg) {
  delete msg;
}

// NetData destructor
//
NetData::~NetData() {
  for (std::size_t i = 0; i < netmodel.size(); ++i) {
    delete netmodel[i];
  }
}


/**************************************************************************
* Network Persistence
**************************************************************************/

// Send data to network partition
//
void NetData::LoadNetwork(idx_t prtidx, const CkCallback &cb) {
  // Send part to network
  cb.send(parts[prtidx - xprt]);
}

// Receive data from network partition (checkpointing)
//
void NetData::CheckNetwork(mPart *msg) {
  // Stash part
  parts[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++cprt == nprt) {
    cprt = 0;

    // Write data
    CkPrintf("  Checking network data %" PRIidx "\n", datidx);
    WriteCSR(true);

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete parts[i];
    }

    // Return control to main
    CkCallback *cb = new CkCallback(CkIndex_Main::CheckNetwork(NULL), mainProxy);
    contribute(nprt*sizeof(dist_t), netdist.data(), net_dist, *cb);
  }
}

// Receive data from network partition
//
void NetData::SaveNetwork(mPart *msg) {
  // Stash part
  parts[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++cprt == nprt) {
    cprt = 0;

    // Write data
    CkPrintf("  Writing network data %" PRIidx "\n", datidx);
    WriteCSR();

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete parts[i];
    }

#ifdef STACS_WITH_YARP
    // Finalize YARP
    yarp.fini();
#endif

    // Return control to main
    CkCallback *cb = new CkCallback(CkIndex_Main::SaveNetwork(NULL), mainProxy);
    contribute(nprt*sizeof(dist_t), netdist.data(), net_dist, *cb);
  }
}

// Receive data from network partition
//
void NetData::CloseNetwork() {
  // Wait for all parts
  if (++cprt == nprt) {
    cprt = 0;

#ifdef STACS_WITH_YARP
    // Finalize YARP
    yarp.fini();
#endif

    // Return control to main
    contribute(0, NULL, CkReduction::nop);
  }
}
