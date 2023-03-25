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
  normdist = new std::normal_distribution<real_t> (0.0, 1.0);
  
  // Simulation configuration
  plastic = msg->plastic;
  episodic = msg->episodic;
  loadbal = msg->loadbal;
  selfconn = msg->selfconn;

  // Network Models
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
  // TODO: Roll this into model creation perhaps
  model[0]->setPlastic(false);

  idx_t jparamname = 0;
  // User defined models
  for (idx_t i = 1; i < msg->nmodel+1; ++i) {
    // Create model object
    model.push_back(ModelFactory::newModel()->Create(msg->modtype[i-1]));
    modname[i] = std::string(msg->modname + msg->xmodname[i-1], msg->modname + msg->xmodname[i]);
    modmap[modname[i]] = i;
    model[i]->setPort(msg->port + msg->xport[i-1]);
    model[i]->setRandom(unifdist, &rngine);
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
    for (idx_t j = 0; j < msg->nparam[i-1]; ++j) {
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
    for (idx_t j = 0; j < model[i]->getNParam(); ++j) {
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
  modelconf.resize(msg->nmodel);
  for (std::size_t i = 0; i < modelconf.size(); ++i) {
    // modname
    modelconf[i].modname = std::string(msg->modname + msg->xmodname[i], msg->modname + msg->xmodname[i+1]);
    // graph type
    modelconf[i].graphtype = msg->graphtype[i];
    // states and param sizes (user specified)
    modelconf[i].nstate = model[modmap[modelconf[i].modname]]->getNState();
    modelconf[i].nstick = model[modmap[modelconf[i].modname]]->getNStick();
    modelconf[i].nparam = model[modmap[modelconf[i].modname]]->getNParam();
    // names (may be in a different order than implemented model)
    // Find the mapping from user-provided state names to the implemented state names
    // Find which states were not specified (and will need model-supplied defaults)
    std::vector<idx_t> statemap;
    statemap.resize(msg->nstate[i]);
    std::vector<bool> stateconfig;
    stateconfig.resize(modelconf[i].nstate, false);
    modelconf[i].statename.resize(modelconf[i].nstate);
    modelconf[i].statename = model[modmap[modelconf[i].modname]]->getStateList();
    for (idx_t j = 0; j < msg->nstate[i]; ++j) {
      std::string statename = std::string(msg->statename + msg->xstatename[jstatename], msg->statename + msg->xstatename[jstatename+1]);
      statemap[j] = model[modmap[modelconf[i].modname]]->getStateIdx(statename.c_str());
      // some basic error checking
      if (statemap[j] == -1) {
        CkPrintf("  state name: %s is invalid for model: %s\n", statename.c_str(), modelconf[i].modname.c_str());
        CkExit();
      }
      modelconf[i].statename[statemap[j]] = statename;
      stateconfig[statemap[j]] = true;
      ++jstatename;
    }
    std::vector<idx_t> stickmap;
    stickmap.resize(msg->nstick[i]);
    std::vector<bool> stickconfig;
    stickconfig.resize(modelconf[i].nstick, false);
    modelconf[i].stickname.resize(modelconf[i].nstick);
    modelconf[i].stickname = model[modmap[modelconf[i].modname]]->getStickList();
    for (idx_t j = 0; j < msg->nstick[i]; ++j) {
      std::string stickname = std::string(msg->stickname + msg->xstickname[jstickname], msg->stickname + msg->xstickname[jstickname+1]);
      stickmap[j] = model[modmap[modelconf[i].modname]]->getStickIdx(stickname.c_str());
      // some basic error checking
      if (stickmap[j] == -1) {
        CkPrintf("  state name: %s is invalid for model: %s\n", stickname.c_str(), modelconf[i].modname.c_str());
        CkExit();
      }
      modelconf[i].stickname[stickmap[j]] = stickname;
      stickconfig[stickmap[j]] = true;
      ++jstickname;
    }
    // prepare containers for parameters
    modelconf[i].stateinit.resize(modelconf[i].nstate);
    modelconf[i].stateparam.resize(modelconf[i].nstate);
    for (idx_t j = 0; j < msg->nstate[i]; ++j) {
      // stateinit
      modelconf[i].stateinit[statemap[j]] = msg->stateinit[msg->xstateinit[i] + j];
      switch (modelconf[i].stateinit[statemap[j]]) {
        case RNGTYPE_CONST:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_CONST);
          break;
        case RNGTYPE_UNIF:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_UNIF);
          break;
        case RNGTYPE_UNINT:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_UNINT);
          break;
        case RNGTYPE_NORM:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_NORM);
          break;
        case RNGTYPE_BNORM:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_BNORM);
          break;
        case RNGTYPE_LBNORM:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_LBNORM);
          break;
        case RNGTYPE_LBLOGNORM:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_LBLOGNORM);
          break;
        case RNGTYPE_LIN:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_LIN);
          break;
        case RNGTYPE_BLIN:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_BLIN);
          break;
        case RNGTYPE_FILE:
          modelconf[i].stateparam[statemap[j]].resize(RNGPARAM_FILE);
          break;
        default:
          CkPrintf("Error: unknown stateinit\n");
          break;
      }
      for (std::size_t s = 0; s < modelconf[i].stateparam[statemap[j]].size(); ++s) {
        modelconf[i].stateparam[statemap[j]][s] = msg->stateparam[jstateparam++];
      }
    }
    // prepare containers
    modelconf[i].stickinit.resize(modelconf[i].nstick);
    modelconf[i].stickparam.resize(modelconf[i].nstick);
    for (idx_t j = 0; j < msg->nstick[i]; ++j) {
      // stickinit
      modelconf[i].stickinit[stickmap[j]] = msg->stickinit[msg->xstickinit[i] + j];
      switch (modelconf[i].stickinit[stickmap[j]]) {
        case RNGTYPE_CONST:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_CONST);
          break;
        case RNGTYPE_UNIF:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_UNIF);
          break;
        case RNGTYPE_UNINT:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_UNINT);
          break;
        case RNGTYPE_NORM:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_NORM);
          break;
        case RNGTYPE_BNORM:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_BNORM);
          break;
        case RNGTYPE_LBNORM:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_LBNORM);
          break;
        case RNGTYPE_LBLOGNORM:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_LBLOGNORM);
          break;
        case RNGTYPE_LIN:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_LIN);
          break;
        case RNGTYPE_BLIN:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_BLIN);
          break;
        case RNGTYPE_FILE:
          modelconf[i].stickparam[stickmap[j]].resize(RNGPARAM_FILE);
          break;
        default:
          CkPrintf("Error: unknown stateinit\n");
          break;
      }
      for (std::size_t s = 0; s < modelconf[i].stickparam[stickmap[j]].size(); ++s) {
        modelconf[i].stickparam[stickmap[j]][s] = msg->stickparam[jstickparam++];
      }
    }
    // Now go through the states that weren't defined by the model config
    // These are the false entries in stateconfig
    for (idx_t j = 0; j < modelconf[i].nstate; ++j) {
      if (stateconfig[j]) { continue; }
      else {
        modelconf[i].stateinit[j] = model[modmap[modelconf[i].modname]]->getDefaultStateType(j);
        modelconf[i].stateparam[j] = model[modmap[modelconf[i].modname]]->getDefaultStateParam(j);
      }
    }
    for (idx_t j = 0; j < modelconf[i].nstick; ++j) {
      if (stickconfig[j]) { continue; }
      else {
        modelconf[i].stickinit[j] = model[modmap[modelconf[i].modname]]->getDefaultStickType(j);
        modelconf[i].stickparam[j] = model[modmap[modelconf[i].modname]]->getDefaultStickParam(j);
      }
    }
  }
  // Sanity check
  CkAssert(jstateparam == msg->nstateparam);
  CkAssert(jstickparam == msg->nstickparam);
  
  // Pointers to data files
  datafiles.resize(msg->ndatafiles);
  for (std::size_t i = 0; i < datafiles.size(); ++i) {
    // filename
    datafiles[i].filename = std::string(msg->datafiles + msg->xdatafiles[i], msg->datafiles + msg->xdatafiles[i+1]);
    // sparse flag
    datafiles[i].filetype = msg->datatypes[i];
  }

  // Logging events
  evtlog.clear();
  evtloglist.resize(EVENT_TOTAL);
  // Always log spikes?
  evtloglist[EVENT_SPIKE] = true;
  for (idx_t i = 0; i < msg->nevtlog; ++i) {
    evtloglist[msg->evtloglist[i]] = true;
  }
  // Recording values
  record.clear();
  recordlist.clear();
  recordmodset.resize(model.size());
  for (std::size_t i = 0; i < model.size(); ++i) {
    recordmodset[i].clear();
  }
  for (idx_t r = 0; r < msg->nrecord; ++r) {
    recordlist.push_back(track_t());
    recordlist[r].trec = 0;
    recordlist[r].tfreq = msg->rectfreq[r];
    recordlist[r].recmodidx = msg->recmodidx[r];
    recordmodset[msg->recmodidx[r]].insert(r);
    // Get the state/stick indices
    std::string statename = std::string(msg->recstate + msg->xrecstate[r], msg->recstate + msg->xrecstate[r+1]);
    idx_t sttidx = model[msg->recmodidx[r]]->getStateIdx(statename);
    if (sttidx >= 0) {
      recordlist[r].rectype = RECORD_STATE;
      recordlist[r].recsttidx = sttidx;
    } else {
      sttidx = model[msg->recmodidx[r]]->getStickIdx(statename);
      if (sttidx >= 0) {
        recordlist[r].rectype = RECORD_STICK;
        recordlist[r].recsttidx = sttidx;
      } else {
        CkPrintf("  error: state name %s not valid for model %" PRIidx "\n", statename.c_str(), msg->recmodidx[r]);
        CkExit();
      }
    }
    // determining vertices will have to wait until models are loaded
  }

  delete msg;

  // Preparing network build (if needed)
  connvtxreq.clear();
  reordlist.clear();
  cpprt = 0;
  cphnd = 0;
  tsim = 0;
  iter = 0;

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
  vtxmodidx.resize(msg->nvtx);
  xyz.resize(msg->nvtx*3);
  adjcy.resize(msg->nvtx);
  edgmodidx.resize(msg->nvtx);
  state.resize(msg->nvtx);
  stick.resize(msg->nvtx);
  evtcal.clear();
  evtcal.resize(msg->nvtx);
  evtcol.clear();
  evtcol.resize(msg->nvtx);
  norderprt = msg->nvtx;
  
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
  for (idx_t i = 0; (std::size_t) i < adjcy.size(); ++i) {
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
    state[i][0] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + model[vtxmodidx[i]]->getNState());
    jstate += model[vtxmodidx[i]]->getNState();
    stick[i][0] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + model[vtxmodidx[i]]->getNStick());
    jstick += model[vtxmodidx[i]]->getNStick();
    // copy over edge data
    for (idx_t j = 0; (std::size_t) j < adjcy[i].size(); ++j) {
      adjcy[i][j] = msg->adjcy[msg->xadj[i] + j];
      // models
      edgmodidx[i][j] = msg->edgmodidx[jmodidx++];
      state[i][j+1] = std::vector<real_t>(msg->state + jstate, msg->state + jstate + model[edgmodidx[i][j]]->getNState());
      jstate += model[edgmodidx[i][j]]->getNState();
      stick[i][j+1] = std::vector<tick_t>(msg->stick + jstick, msg->stick + jstick + model[edgmodidx[i][j]]->getNStick());
      jstick += model[edgmodidx[i][j]]->getNStick();
    }
    // copy over event data
    event_t event;
    idx_t arrival;
    // copy over event data
    for (idx_t e = msg->xevent[i]; e < msg->xevent[i+1]; ++e) {
      event.diffuse = msg->diffuse[e];
      event.type = msg->type[e];
      event.source = msg->source[e];
      event.index = msg->index[e];
      event.data = msg->data[e];
      // Add to event queue or spillover
      arrival = (idx_t) (event.diffuse/tstep);
      if (arrival < nevtday) {
        evtcal[i][(arrival)%nevtday].push_back(event);
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
  
  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

/**************************************************************************
* Network Simulation Initialization
**************************************************************************/

// Coordination with NetData chare array
//
void Network::InitProxy(CProxy_Netdata cpdata) {
  // Set proxy
  netdata = cpdata;
  
  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// Initialize network for cycling
//
void Network::InitNetwork(std::string runmode, const CkCallback &cbcheck) {
  // Set compute cycle depending on runmode
  if (runmode == std::string(RUNMODE_SIMULATE)) {
    cyclepart = CkCallback(CkIndex_Network::CycleSim(), thisProxy(prtidx));
  }
  else if (runmode == std::string(RUNMODE_SIMGPU)) {
    cyclepart = CkCallback(CkIndex_Network::CycleSimGPU(), thisProxy(prtidx));
  }
  else if (runmode == std::string(RUNMODE_FINDGROUP)) {
    cyclepart = CkCallback(CkIndex_Network::CycleGroup(), thisProxy(prtidx));
  }
  else if (runmode == std::string(RUNMODE_ESTIMATE)) {
    cyclepart = CkCallback(CkIndex_Network::CycleEst(), thisProxy(prtidx));
  }
  maincheck = cbcheck;

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
  
  // open ports if needed
  for (idx_t i = 0; (std::size_t) i < vtxmodidx.size(); ++i) {
    if (model[vtxmodidx[i]]->getNPort()) {
      model[vtxmodidx[i]]->OpenPorts();
    }
  }
  
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
// Load data from file
void Network::LoadData() {
  // Request network part from netdata
  netdata(datidx).LoadNetwork(prtidx,
      CkCallback(CkIndex_Network::LoadNetwork(NULL), thisProxy(prtidx)));
}

// Initialize Network for simulation
//
void Network::StartNetwork() {
  // Network data should already be loaded/built/repartitioned
  idx_t nadjcy;

  // Setup data vectors
  vtxidx.clear();
  vtxidx.resize(norderprt);
  vtxmap.clear();
  adjmap.clear();
  vtxaux.clear();
  vtxaux.resize(norderprt);
  events.clear();
  evtext.clear();
  evtrpc.clear();
  leapidx.clear();
  leapidx.resize(model.size());

  // Initalize counters
  nadjcy = 0;

  // Get adjacency matrix
  for (idx_t i = 0; (std::size_t) i < adjcy.size(); ++i) {
    vtxidx[i] = vtxdist[prtidx] + i;
    vtxmap[vtxdist[prtidx] + i] = i;
    // preallocate sizes
    evtcal[i].resize(nevtday);
    nadjcy += adjcy[i].size();
    // copy over vertex data
    if (leaplist[vtxmodidx[i]]) {
      leapidx[vtxmodidx[i]].push_back(std::array<idx_t, 2>{{i, 0}});
    }
    // copy over edge data
    for (idx_t j = 0; (std::size_t) j < adjcy[i].size(); ++j) {
      if (edgmodidx[i][j]) {
        // map of targets (edge model in part)
        adjmap[adjcy[i][j]].push_back(std::array<idx_t, 2>{{i, j+1}});
      }
      if (leaplist[edgmodidx[i][j]]) {
        leapidx[edgmodidx[i][j]].push_back(std::array<idx_t, 2>{{i, j+1}});
      }
    }
    // set up auxiliary state
    vtxaux[i].clear();
    if (model[vtxmodidx[i]]->getNAux()) {
      std::vector<std::string> auxstate = model[vtxmodidx[i]]->getAuxState();
      std::vector<std::string> auxstick = model[vtxmodidx[i]]->getAuxStick();
      for (idx_t j = 0; (std::size_t) j < edgmodidx[i].size(); ++j) {
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
  }
  
  // Multicast Communication 
  CreateComm();
  
  // Print some information
  CkPrintf("  Network part %d:   vtx: %zu   edg: %" PRIidx "   adjvtx: %zu   adjpart: %" PRIidx "\n",
           prtidx, adjcy.size(), nadjcy, adjmap.size(), nadjpart);
  
  // Set up recording
  // Header information (vtx/edg corresponding to data)
  for (idx_t r = 0; (std::size_t) r < recordlist.size(); ++r) {
    record.push_back(record_t());
    record.back().recidx = r;
    record.back().trec = 0;
    record.back().state.clear();
    record.back().stick.clear();
    record.back().index.clear();
    recordlist[r].recvtxidx.clear();
    recordlist[r].recedgidx.resize(vtxmodidx.size());
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
      recordlist[r].recedgidx[i].clear();
    }
  }
  // Populate recording indices
  std::set<idx_t>::iterator jrec;
  for (idx_t i = 0; (std::size_t) i < vtxmodidx.size(); ++i) {
    // vertices
    for (jrec = recordmodset[vtxmodidx[i]].begin(); jrec != recordmodset[vtxmodidx[i]].end(); ++jrec) {
      recordlist[*jrec].recvtxidx.push_back(i);
      recordlist[*jrec].recedgidx[i].push_back(0);
      record[*jrec].index.push_back(i);
      record[*jrec].index.push_back(0);
    }
    // edges
    for (idx_t j = 0; (std::size_t) j < edgmodidx[i].size(); ++j) {
      for (jrec = recordmodset[edgmodidx[i][j]].begin(); jrec != recordmodset[edgmodidx[i][j]].end(); ++jrec) {
        if (recordlist[*jrec].recvtxidx.size() == 0 || recordlist[*jrec].recvtxidx.back() != i) {
          recordlist[*jrec].recvtxidx.push_back(i);
        }
        recordlist[*jrec].recedgidx[i].push_back(j+1);
        record[*jrec].index.push_back(i);
        record[*jrec].index.push_back(j+1);
      }
    }
  }

  // Start the cycles
  cyclepart.send();
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
  cyclepart.send();
}

void Network::SaveBuild() {
  // Build network part message for saving
  mPart *mpart = BuildPart();
  netdata(datidx).SaveBuild(mpart);
  
  // Return control to main
  contribute(0, NULL, CkReduction::nop);
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
