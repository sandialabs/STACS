/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "stacs.h"
#include "network.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ unsigned randseed;
extern /*readonly*/ CProxy_Main mainProxy;
extern /*readonly*/ int netfiles;
extern /*readonly*/ int netparts;


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

// Netdata constructor
//
Netdata::Netdata(mModel *msg) {
  // Bookkeeping
  fileidx = thisIndex;
  int ndiv = netparts/netfiles;
  int nrem = netparts%netfiles;
  cpart = rpart = 0;
  npart = ndiv + (fileidx < nrem);
  xpart = fileidx*ndiv + (fileidx < nrem ? fileidx : nrem);
  
  // Bookkeeping 2
  datidx = thisIndex;
  nprt = ndiv + (datidx < nrem);
  xprt = datidx*ndiv + (datidx < nrem ? datidx : nrem);

  // Set up random number generator
  rngine.seed(randseed+datidx);
  unifdist = new std::uniform_real_distribution<real_t> (0.0, 1.0);
  normdist = new std::normal_distribution<real_t> (0.0, 1.0);

  // RNG types (for errors)
  rngtype.resize(RNGTYPE_NRNG);
  rngtype[RNGTYPE_CONST] = std::string("constant");
  rngtype[RNGTYPE_UNIF] = std::string("uniform");
  rngtype[RNGTYPE_UNINT] = std::string("uniform interval");
  rngtype[RNGTYPE_NORM] = std::string("normal");
  rngtype[RNGTYPE_BNORM] = std::string("bounded normal");
  rngtype[RNGTYPE_BNORM] = std::string("lower bounded normal");
  rngtype[RNGTYPE_LIN] = std::string("linear");
  rngtype[RNGTYPE_LBLIN] = std::string("lower bounded linear");
  rngtype[RNGTYPE_UBLIN] = std::string("upper bounded linear");
  rngtype[RNGTYPE_BLIN] = std::string("bounded linear");
  rngtype[RNGTYPE_FILE] = std::string("file");

  // Set up counters
  idx_t jstateparam = 0;
  idx_t jstickparam = 0;

  // Set up maps
  modmap.clear();
  modname.resize(msg->nmodel+1);
  modname[0] = std::string("none");
  modmap[modname[0]] = 0;
  // Read in models
  modeldata.resize(msg->nmodel);
  for (std::size_t i = 0; i < modeldata.size(); ++i) {
    // modname
    modeldata[i].modname = std::string(msg->modname + msg->xmodname[i], msg->modname + msg->xmodname[i+1]);
    modname[i+1] = modeldata[i].modname;
    modmap[modeldata[i].modname] = i+1;
    // type
    modeldata[i].graphtype = msg->graphtype[i];
    // prepare containers
    modeldata[i].statetype.resize(msg->xstatetype[i+1] - msg->xstatetype[i]);
    modeldata[i].stateparam.resize(msg->xstatetype[i+1] - msg->xstatetype[i]);
    for (std::size_t j = 0; j < modeldata[i].statetype.size(); ++j) {
      // statetype
      modeldata[i].statetype[j] = msg->statetype[msg->xstatetype[i] + j];
      switch (modeldata[i].statetype[j]) {
        case RNGTYPE_CONST:
          modeldata[i].stateparam[j].resize(RNGPARAM_CONST);
          break;
        case RNGTYPE_UNIF:
          modeldata[i].stateparam[j].resize(RNGPARAM_UNIF);
          break;
        case RNGTYPE_UNINT:
          modeldata[i].stateparam[j].resize(RNGPARAM_UNINT);
          break;
        case RNGTYPE_NORM:
          modeldata[i].stateparam[j].resize(RNGPARAM_NORM);
          break;
        case RNGTYPE_BNORM:
          modeldata[i].stateparam[j].resize(RNGPARAM_BNORM);
          break;
        case RNGTYPE_LBNORM:
          modeldata[i].stateparam[j].resize(RNGPARAM_LBNORM);
          break;
        case RNGTYPE_LIN:
          modeldata[i].stateparam[j].resize(RNGPARAM_LIN);
          break;
        case RNGTYPE_BLIN:
          modeldata[i].stateparam[j].resize(RNGPARAM_BLIN);
          break;
        case RNGTYPE_FILE:
          modeldata[i].stateparam[j].resize(RNGPARAM_FILE);
          break;
        default:
          CkPrintf("Error: unknown statetype\n");
          break;
      }
      for (std::size_t s = 0; s < modeldata[i].stateparam[j].size(); ++s) {
        modeldata[i].stateparam[j][s] = msg->stateparam[jstateparam++];
      }
    }
    // prepare containers
    modeldata[i].sticktype.resize(msg->xsticktype[i+1] - msg->xsticktype[i]);
    modeldata[i].stickparam.resize(msg->xsticktype[i+1] - msg->xsticktype[i]);
    for (std::size_t j = 0; j < modeldata[i].sticktype.size(); ++j) {
      // sticktype
      modeldata[i].sticktype[j] = msg->sticktype[msg->xsticktype[i] + j];
      switch (modeldata[i].sticktype[j]) {
        case RNGTYPE_CONST:
          modeldata[i].stickparam[j].resize(RNGPARAM_CONST);
          break;
        case RNGTYPE_UNIF:
          modeldata[i].stickparam[j].resize(RNGPARAM_UNIF);
          break;
        case RNGTYPE_UNINT:
          modeldata[i].stickparam[j].resize(RNGPARAM_UNINT);
          break;
        case RNGTYPE_NORM:
          modeldata[i].stickparam[j].resize(RNGPARAM_NORM);
          break;
        case RNGTYPE_BNORM:
          modeldata[i].stickparam[j].resize(RNGPARAM_BNORM);
          break;
        case RNGTYPE_LBNORM:
          modeldata[i].stickparam[j].resize(RNGPARAM_LBNORM);
          break;
        case RNGTYPE_LIN:
          modeldata[i].stickparam[j].resize(RNGPARAM_LIN);
          break;
        case RNGTYPE_BLIN:
          modeldata[i].stickparam[j].resize(RNGPARAM_BLIN);
          break;
        case RNGTYPE_FILE:
          modeldata[i].stickparam[j].resize(RNGPARAM_FILE);
          break;
        default:
          CkPrintf("Error: unknown statetype\n");
          break;
      }
      for (std::size_t s = 0; s < modeldata[i].stickparam[j].size(); ++s) {
        modeldata[i].stickparam[j][s] = msg->stickparam[jstickparam++];
      }
    }
  }
  // Sanity check
  CkAssert(jstateparam == msg->nstateparam);
  CkAssert(jstickparam == msg->nstickparam);

  // Read in data files
  datafiles.resize(msg->ndatafiles);
  for (std::size_t i = 0; i < datafiles.size(); ++i) {
    // filename
    datafiles[i].filename = std::string(msg->datafiles + msg->xdatafiles[i], msg->datafiles + msg->xdatafiles[i+1]);
    // read in data (as matrix)
    if (ReadDataCSV(datafiles[i])) {
      CkPrintf("Error reading data file %s...\n", datafiles[i].filename.c_str());
      CkExit();
    }
  }

  // Models
  // TODO: set up models the same way as network
  for (std::size_t i = 0; i < model.size(); ++i) {
    delete model[i];
  }
  // Set up containers
  model.clear();
  modname.resize(msg->nmodel+1);
  // "none" model
  model.push_back(ModelFactory::newModel()->Create(0));
  modname[0] = std::string("none");
  modmap[modname[0]] = 0;
  // User defined models
  for (idx_t i = 1; i < msg->nmodel+1; ++i) {
    model.push_back(ModelFactory::newModel()->Create(msg->modtype[i-1]));
    modname[i] = std::string(msg->modname + msg->xmodname[i-1], msg->modname + msg->xmodname[i]);
    modmap[modname[i]] = i;
    /*
    if (fileidx == 0) {
      CkPrintf("  Netdata model: %" PRIidx "   NStates: %" PRIidx "   Name: %s\n", i, model[i]->getNState(), modname[i].c_str());
    }
    */
  }
  delete msg;

#ifdef STACS_WITH_YARP
  // Open yarp
  yarp.init();
#endif

  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// Netdata migration
//
Netdata::Netdata(CkMigrateMessage *msg) {
  delete msg;
}

// Netdata destructor
//
Netdata::~Netdata() {
  for (std::size_t i = 0; i < model.size(); ++i) {
    delete model[i];
  }
}


/**************************************************************************
* Load Network Data
**************************************************************************/

void Netdata::LoadData(mDist *msg) {
  // Persistence
  vtxdist.resize(netparts+1);
  edgdist.resize(netparts+1);
  statedist.resize(netparts+1);
  stickdist.resize(netparts+1);
  eventdist.resize(netparts+1);
  for (idx_t i = 0; i < netparts+1; ++i) {
    vtxdist[i] = msg->vtxdist[i];
    edgdist[i] = msg->edgdist[i];
    statedist[i] = msg->statedist[i];
    stickdist[i] = msg->stickdist[i];
    eventdist[i] = msg->eventdist[i];
  }

  // Network distribution
  maindist = CkCallback(CkIndex_Main::SaveDist(NULL), mainProxy);
  
  // Data
  parts.resize(npart);
  records.resize(npart);

  // Read in files
  CkPrintf("Reading network data files %" PRIidx "\n", fileidx);
  ReadNetwork();

  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// Send data to network partition
//
void Netdata::LoadNetwork(int partidx, const CkCallback &cbpart) {
  // Send part to network
  cbpart.send(parts[partidx - xpart]);
}


/**************************************************************************
* Save Network Data
**************************************************************************/

// Save data from built network
//
void Netdata::SaveBuild() {
  // Write data
  WriteBuild();

  // Return control to main
  contribute(npart*sizeof(dist_t), netdist.data(), net_dist,
      CkCallback(CkIndex_Main::SaveFinalDist(NULL), mainProxy));
}

// Save data from network partition
//
void Netdata::SaveNetwork(mPart *msg) {
  // Stash part
  parts[msg->partidx - xpart] = msg;
  
  // Wait for all parts
  if (++cpart == npart) {
    cpart = 0;

    // Write data
    WriteNetwork();

    // Cleanup stash
    for (idx_t i = 0; i < npart; ++i) {
      delete parts[i];
    }

    // Return control to main
    contribute(npart*sizeof(dist_t), netdist.data(), net_dist, maindist);
  }
}

// Save data from network partition (final)
//
void Netdata::SaveCloseNetwork(mPart *msg) {
  // Stash part
  parts[msg->partidx - xpart] = msg;
  
  // Wait for all parts
  if (++cpart == npart) {
    cpart = 0;

    // Write data
    WriteNetwork();

    // Cleanup stash
    for (idx_t i = 0; i < npart; ++i) {
      delete parts[i];
    }

#ifdef STACS_WITH_YARP
    // Finalize YARP
    yarp.fini();
#endif

    // Return control to main
    contribute(npart*sizeof(dist_t), netdist.data(), net_dist,
        CkCallback(CkIndex_Main::SaveFinalDist(NULL), mainProxy));
  }
}

// Finalize Netdata chare array
//
void Netdata::CloseNetwork() {
  // Wait for all parts
  if (++cpart == npart) {
    cpart = 0;

#ifdef STACS_WITH_YARP
    // Finalize YARP
    yarp.fini();
#endif

    // Return control to main
    contribute(0, NULL, CkReduction::nop);
  }
}

/**************************************************************************
* Save Network Distribution
**************************************************************************/

// Collect distribution from netdata
//
void Main::SaveDist(CkReductionMsg *msg) {
  //CkPrintf("Checkpointing simulation\n");
  
  // Save network part distribution to local
  netdist.clear();
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(dist_t); ++i) {
    netdist.push_back(*((dist_t *)msg->getData()+i));
  }
  delete msg;

  // Write distribution
  if (WriteDist()) {
    CkPrintf("Error writing distribution...\n");
    CkExit();
  }
}

// Collect distribution from netdata (final)
//
void Main::SaveFinalDist(CkReductionMsg *msg) {
  // Save network part distribution to local
  netdist.clear();
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(dist_t); ++i) {
    netdist.push_back(*((dist_t *)msg->getData()+i));
  }
  delete msg;

  // Write distribution
  if (WriteDist()) {
    CkPrintf("Error writing distribution...\n");
    CkExit();
  }

  // Finished
  thisProxy.Halt();
}


/**************************************************************************
* Estimation Recording
**************************************************************************/

// Collect estimates for writing
//
void Netdata::SaveEstimate(CkReductionMsg *msg) {
  // Add to group log
  grplog.clear();
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(event_t); ++i) {
    grplog.push_back(*((event_t *)msg->getData()+i));
  }
  delete msg;

  // Sorting
  std::sort(grplog.begin(), grplog.end());

  // Write estimate
  WriteEstimate();
}

// Collect estimates for writing (final)
//
void Netdata::SaveFinalEstimate(CkReductionMsg *msg) {
  // Add to group log
  grplog.clear();
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(event_t); ++i) {
    grplog.push_back(*((event_t *)msg->getData()+i));
  }
  delete msg;

  // Sorting
  std::sort(grplog.begin(), grplog.end());
  
  // Write estimate
  WriteEstimate();
    
  // Return control to main (halting)
  mainProxy.Halt();
}
