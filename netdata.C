/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
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
  datidx = thisIndex;
  int ndiv = netparts/netfiles;
  int nrem = netparts%netfiles;
  cprt = rprt = 0;
  nprt = ndiv + (datidx < nrem);
  xprt = datidx*ndiv + (datidx < nrem ? datidx : nrem);

  // Models (implemented)
  // TODO: set up models the same way as network?
  for (std::size_t i = 0; i < model.size(); ++i) {
    delete model[i];
  }
  // Set up containers
  model.clear();
  modname.resize(msg->nmodel+1);
  modmap.clear();
  // "none" model
  model.push_back(ModelFactory::newModel()->Create(0));
  modname[0] = std::string("none");
  modmap[modname[0]] = 0;
  // User defined models
  for (idx_t i = 1; i < msg->nmodel+1; ++i) {
    model.push_back(ModelFactory::newModel()->Create(msg->modtype[i-1]));
    modname[i] = std::string(msg->modname + msg->xmodname[i-1], msg->modname + msg->xmodname[i]);
    modmap[modname[i]] = i;
  }

  // Models (user configured information)

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
  idx_t jstatename = 0;
  idx_t jstickname = 0;

  // Read in models
  modeldata.resize(msg->nmodel);
  for (std::size_t i = 0; i < modeldata.size(); ++i) {
    // modname
    modeldata[i].modname = std::string(msg->modname + msg->xmodname[i], msg->modname + msg->xmodname[i+1]);
    // graph type
    modeldata[i].graphtype = msg->graphtype[i];
    // states and param sizes (user specified)
    modeldata[i].nstate = model[modmap[modeldata[i].modname]]->getNState();
    modeldata[i].nstick = model[modmap[modeldata[i].modname]]->getNStick();
    modeldata[i].nparam = model[modmap[modeldata[i].modname]]->getNParam();
    // names (may be in a different order than implemented model)
    // Find the mapping from user-provided state names to the implemented state names
    // Find which states were not specified (and will need model-supplied defaults)
    std::vector<idx_t> statemap;
    statemap.resize(msg->nstate[i]);
    std::vector<bool> stateconfig;
    stateconfig.resize(modeldata[i].nstate, false);
    modeldata[i].statename.resize(modeldata[i].nstate);
    modeldata[i].statename = model[modmap[modeldata[i].modname]]->getStateList();
    for (std::size_t j = 0; j < msg->nstate[i]; ++j) {
      std::string statename = std::string(msg->statename + msg->xstatename[jstatename], msg->statename + msg->xstatename[jstatename+1]);
      statemap[j] = model[modmap[modeldata[i].modname]]->getStateIdx(statename.c_str());
      // some basic error checking
      if (statemap[j] == -1) {
        CkPrintf("  state name: %s is invalid for model: %s\n", statename.c_str(), modeldata[i].modname.c_str());
        CkExit();
      }
      // TODO: We had set the statename earlier, so change this to error checking instead
      modeldata[i].statename[statemap[j]] = statename;
      stateconfig[statemap[j]] = true;
      ++jstatename;
    }
    std::vector<idx_t> stickmap;
    stickmap.resize(msg->nstick[i]);
    std::vector<bool> stickconfig;
    stickconfig.resize(modeldata[i].nstick, false);
    modeldata[i].stickname.resize(modeldata[i].nstick);
    modeldata[i].stickname = model[modmap[modeldata[i].modname]]->getStickList();
    for (std::size_t j = 0; j < msg->nstick[i]; ++j) {
      std::string stickname = std::string(msg->stickname + msg->xstickname[jstickname], msg->stickname + msg->xstickname[jstickname+1]);
      stickmap[j] = model[modmap[modeldata[i].modname]]->getStickIdx(stickname.c_str());
      // some basic error checking
      if (stickmap[j] == -1) {
        CkPrintf("  state name: %s is invalid for model: %s\n", stickname.c_str(), modeldata[i].modname.c_str());
        CkExit();
      }
      modeldata[i].stickname[stickmap[j]] = stickname;
      stickconfig[stickmap[j]] = true;
      ++jstickname;
    }
    // prepare containers for parameters
    modeldata[i].statetype.resize(modeldata[i].nstate);
    modeldata[i].stateparam.resize(modeldata[i].nstate);
    for (std::size_t j = 0; j < msg->nstate[i]; ++j) {
      // statetype
      modeldata[i].statetype[statemap[j]] = msg->statetype[msg->xstatetype[i] + j];
      switch (modeldata[i].statetype[statemap[j]]) {
        case RNGTYPE_CONST:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_CONST);
          break;
        case RNGTYPE_UNIF:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_UNIF);
          break;
        case RNGTYPE_UNINT:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_UNINT);
          break;
        case RNGTYPE_NORM:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_NORM);
          break;
        case RNGTYPE_BNORM:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_BNORM);
          break;
        case RNGTYPE_LBNORM:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_LBNORM);
          break;
        case RNGTYPE_LBLOGNORM:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_LBLOGNORM);
          break;
        case RNGTYPE_LIN:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_LIN);
          break;
        case RNGTYPE_BLIN:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_BLIN);
          break;
        case RNGTYPE_FILE:
          modeldata[i].stateparam[statemap[j]].resize(RNGPARAM_FILE);
          break;
        default:
          CkPrintf("Error: unknown statetype\n");
          break;
      }
      for (std::size_t s = 0; s < modeldata[i].stateparam[statemap[j]].size(); ++s) {
        modeldata[i].stateparam[statemap[j]][s] = msg->stateparam[jstateparam++];
      }
    }
    // prepare containers
    modeldata[i].sticktype.resize(modeldata[i].nstick);
    modeldata[i].stickparam.resize(modeldata[i].nstick);
    for (std::size_t j = 0; j < msg->nstick[i]; ++j) {
      // sticktype
      modeldata[i].sticktype[stickmap[j]] = msg->sticktype[msg->xsticktype[i] + j];
      switch (modeldata[i].sticktype[stickmap[j]]) {
        case RNGTYPE_CONST:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_CONST);
          break;
        case RNGTYPE_UNIF:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_UNIF);
          break;
        case RNGTYPE_UNINT:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_UNINT);
          break;
        case RNGTYPE_NORM:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_NORM);
          break;
        case RNGTYPE_BNORM:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_BNORM);
          break;
        case RNGTYPE_LBNORM:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_LBNORM);
          break;
        case RNGTYPE_LBLOGNORM:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_LBLOGNORM);
          break;
        case RNGTYPE_LIN:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_LIN);
          break;
        case RNGTYPE_BLIN:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_BLIN);
          break;
        case RNGTYPE_FILE:
          modeldata[i].stickparam[stickmap[j]].resize(RNGPARAM_FILE);
          break;
        default:
          CkPrintf("Error: unknown statetype\n");
          break;
      }
      for (std::size_t s = 0; s < modeldata[i].stickparam[stickmap[j]].size(); ++s) {
        modeldata[i].stickparam[stickmap[j]][s] = msg->stickparam[jstickparam++];
      }
    }
    // Now go through the states that weren't defined by the model config
    // These are the false entries in stateconfig
    for (std::size_t j = 0; j < modeldata[i].nstate; ++j) {
      if (stateconfig[j]) { continue; }
      else {
        modeldata[i].statetype[j] = model[modmap[modeldata[i].modname]]->getDefaultStateType(j);
        modeldata[i].stateparam[j] = model[modmap[modeldata[i].modname]]->getDefaultStateParam(j);
      }
    }
    for (std::size_t j = 0; j < modeldata[i].nstick; ++j) {
      if (stickconfig[j]) { continue; }
      else {
        modeldata[i].sticktype[j] = model[modmap[modeldata[i].modname]]->getDefaultStickType(j);
        modeldata[i].stickparam[j] = model[modmap[modeldata[i].modname]]->getDefaultStickParam(j);
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
    // sparse flag
    datafiles[i].filetype = msg->datatypes[i];
    // TODO: mabye combine the two functions in ReadDataCSV?
    if (datafiles[i].filetype == FT_CSV_SPARSE) {
      // read in data (as csr if sparse flag set)
      if (ReadDataCSVSparse(datafiles[i])) {
        CkPrintf("Error reading data file %s...\n", datafiles[i].filename.c_str());
        CkExit();
      }
    }
    else if (datafiles[i].filetype == FT_CSV_DENSE) {
      // read in data (as matrix)
      if (ReadDataCSV(datafiles[i])) {
        CkPrintf("Error reading data file %s...\n", datafiles[i].filename.c_str());
        CkExit();
      }
    }
  }
  
  delete msg;

#ifdef STACS_WITH_YARP
  // Open yarp
  yarp.init();
#endif
  
  // Preparing network data
  parts.resize(nprt);
  records.resize(nprt);
  cpprt = 0;

  // Network distribution
  maindist = CkCallback(CkIndex_Main::SaveDist(NULL), mainProxy);
  
  // Preparing network build (if needed)
  connvtxreq.clear();

  // Repartitioning
  vtxidxreprt.resize(netparts);
  vtxmodidxreprt.resize(netparts);
  xyzreprt.resize(netparts);
  edgmodidxreprt.resize(netparts);
  adjcyreprt.resize(netparts);
  statereprt.resize(netparts);
  stickreprt.resize(netparts);
  eventreprt.resize(netparts);


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
  // cleanup
  delete msg;

  // Read in files
  CkPrintf("Reading network data files %" PRIidx "\n", datidx);
  ReadNetwork();

  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// Send data to network partition
//
void Netdata::LoadNetwork(int prtidx, const CkCallback &cbpart) {
  // Send part to network
  cbpart.send(parts[prtidx - xprt]);
}


/**************************************************************************
* Save Network Data
**************************************************************************/

// Save data from built network
//
void Netdata::SaveBuild() {
  // Write data
  WriteNetwork(0);

  // Return control to main
  contribute(nprt*sizeof(dist_t), netdist.data(), net_dist,
      CkCallback(CkIndex_Main::SaveInitDist(NULL), mainProxy));
}

// Save data from built network
//
void Netdata::SaveCloseBuild() {
  // Write data
  WriteNetwork(0);

  // Return control to main
  contribute(nprt*sizeof(dist_t), netdist.data(), net_dist,
      CkCallback(CkIndex_Main::SaveFinalDist(NULL), mainProxy));
}


// Save data from network partition
//
void Netdata::SaveNetwork(mPart *msg) {
  // Stash part
  parts[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++cprt == nprt) {
    cprt = 0;

    // Write data
    WriteNetwork();

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete parts[i];
    }

    // Return control to main
    contribute(nprt*sizeof(dist_t), netdist.data(), net_dist, maindist);
  }
}

// Save data from network partition (final)
//
void Netdata::SaveCloseNetwork(mPart *msg) {
  // Stash part
  parts[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++cprt == nprt) {
    cprt = 0;

    // Write data
    WriteNetwork();

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete parts[i];
    }

#ifdef STACS_WITH_YARP
    // Finalize YARP
    yarp.fini();
#endif

    // Return control to main
    contribute(nprt*sizeof(dist_t), netdist.data(), net_dist,
        CkCallback(CkIndex_Main::SaveFinalDist(NULL), mainProxy));
  }
}

// Finalize Netdata chare array
//
void Netdata::CloseNetwork() {
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

// Collect distribution from netdata
//
void Main::SaveInitDist(CkReductionMsg *msg) {
  //CkPrintf("Checkpointing simulation\n");
  
  // Save network part distribution to local
  netdist.clear();
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(dist_t); ++i) {
    netdist.push_back(*((dist_t *)msg->getData()+i));
  }
  delete msg;

  // Write distribution
  if (WriteDist(0)) {
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


/**************************************************************************
* Convert built network state to part messages
**************************************************************************/

void Netdata::BuildParts() {
  /* Bookkeeping */
  idx_t nadjcy;
  idx_t nstate;
  idx_t nstick;
  idx_t nevent;
  idx_t xvtxidx;

  //parts.clear();
  //parts.resize(nprt);

  for (int k = 0; k < nprt; ++k) {
    // Get total size of adjcy
    nadjcy = 0;
    nstate = 0;
    nstick = 0;
    nevent = 0;
    xvtxidx = vtxdist[xprt+k] - vtxdist[xprt];
    for (idx_t i = 0; i < norderprt[k]; ++i) {
      nadjcy += adjcy[xvtxidx+i].size();
      for (std::size_t j = 0; j < adjcy[xvtxidx+i].size()+1; ++j) {
        nstate += state[xvtxidx+i][j].size();
        nstick += stick[xvtxidx+i][j].size();
      }
      nevent += event[xvtxidx+i].size();
    }

    // Initialize partition data message
    int msgSize[MSG_Part];
    msgSize[0] = netparts+1;    // vtxdist
    msgSize[1] = 0;             // vtxidx (implicit)
    msgSize[2] = norderprt[k];  // vtxmodidx
    msgSize[3] = norderprt[k]*3;// xyz
    msgSize[4] = norderprt[k]+1;// xadj
    msgSize[5] = nadjcy;        // adjcy
    msgSize[6] = nadjcy;        // edgmodidx
    msgSize[7] = nstate;        // state
    msgSize[8] = nstick;        // stick
    msgSize[9] = norderprt[k]+1;// xevent
    msgSize[10] = nevent;        // diffuse
    msgSize[11] = nevent;       // type
    msgSize[12] = nevent;       // source
    msgSize[13] = nevent;       // index
    msgSize[14] = nevent;       // data
    parts[k] = new(msgSize, 0) mPart;

    // Data sizes
    parts[k]->nvtx = norderprt[k];
    parts[k]->nedg = nadjcy;
    parts[k]->nstate = nstate;
    parts[k]->nstick = nstick;
    parts[k]->nevent = nevent;
    parts[k]->prtidx = xprt + k;
  
    // Graph Information
    for (int i = 0; i < netparts+1; ++i) {
      // vtxdist
      parts[k]->vtxdist[i] = vtxdist[i];
    }
    // Vertex and Edge Information
    idx_t jstate = 0;
    idx_t jstick = 0;
    idx_t jevent = 0;
    parts[k]->xadj[0] = 0;
    parts[k]->xevent[0] = 0;
    for (idx_t i = 0; i < norderprt[k]; ++i) {
      // vtxmodidx
      parts[k]->vtxmodidx[i] = vtxmodidx[xvtxidx+i];
      // xyz
      parts[k]->xyz[i*3+0] = xyz[xvtxidx+i*3+0];
      parts[k]->xyz[i*3+1] = xyz[xvtxidx+i*3+1];
      parts[k]->xyz[i*3+2] = xyz[xvtxidx+i*3+2];
      // vertex state
      for (std::size_t s = 0; s < state[xvtxidx+i][0].size(); ++s) {
        parts[k]->state[jstate++] = state[xvtxidx+i][0][s];
      }
      for (std::size_t s = 0; s < stick[xvtxidx+i][0].size(); ++s) {
        parts[k]->stick[jstick++] = stick[xvtxidx+i][0][s];
      }

      // xadj
      parts[k]->xadj[i+1] = parts[k]->xadj[i] + adjcy[xvtxidx+i].size();
      for (std::size_t j = 0; j < adjcy[xvtxidx+i].size(); ++j) {
        // adjcy
        parts[k]->adjcy[parts[k]->xadj[i] + j] = adjcy[xvtxidx+i][j];
        // edgmodidx
        parts[k]->edgmodidx[parts[k]->xadj[i] + j] = edgmodidx[xvtxidx+i][j];
        // state
        for (std::size_t s = 0; s < state[xvtxidx+i][j+1].size(); ++s) {
          parts[k]->state[jstate++] = state[xvtxidx+i][j+1][s];
        }
        for (std::size_t s = 0; s < stick[xvtxidx+i][j+1].size(); ++s) {
          parts[k]->stick[jstick++] = stick[xvtxidx+i][j+1][s];
        }
      }
      
      // events
      for (std::size_t s = 0; s < event[xvtxidx+i].size(); ++s) {
        parts[k]->diffuse[jevent] = event[xvtxidx+i][s].diffuse;
        parts[k]->type[jevent] = event[xvtxidx+i][s].type;
        parts[k]->source[jevent] = event[xvtxidx+i][s].source;
        parts[k]->index[jevent] = event[xvtxidx+i][s].index;
        parts[k]->data[jevent++] = event[xvtxidx+i][s].data;
      }
      // xevent
      parts[k]->xevent[i+1] = jevent;
    }
    CkAssert(jstate == nstate);
    CkAssert(jstick == nstick);
    CkAssert(jevent == nevent);
  }
  // Clear out file-based network information
  vtxmodidx.clear();
  vtxordidx.clear();
  xyz.clear();
  adjcy.clear();
  edgmodidx.clear();
  state.clear();
  stick.clear();
  event.clear();
  // Clear out partitioning
  vtxprted.clear();
  xyzprted.clear();
  adjcyprted.clear();
  edgmodidxprted.clear();
  stateprted.clear();
  stickprted.clear();
  eventprted.clear();
  // also allocate for reordering
  adjcyreord.clear();
  edgmodidxreord.clear();
  statereord.clear();
  stickreord.clear();
  eventsourcereord.clear();
  eventindexreord.clear();
}
