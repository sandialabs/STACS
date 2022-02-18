/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "genet.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ unsigned int randseed;
/*readonly*/ std::string netwkdir;
/*readonly*/ idx_t netparts;
/*readonly*/ int netfiles;
/*readonly*/ std::string filebase;
/*readonly*/ std::string filesave;


/**************************************************************************
* Main
**************************************************************************/

// Main entry point
//
Main::Main(CkArgMsg* msg) {
  // Display title
  CkPrintf("Generate Network for STACS (genet)\n");

  // MPI Interoperate bookkeeping kludge
  bool initok = true;

  // Command line arguments
  std::string configfile;
  if (msg->argc < 2) {
    configfile = "configs/config.yml"; // default
    mode = std::string("build");
  }
  else if (msg->argc == 2) {
    mode = msg->argv[1];
    if (mode != "build" && mode != "part" && mode != "order") {
      configfile = msg->argv[1];
      mode = std::string("build");
    }
    else {
      configfile = "config.yml"; // default
    }
  }
  else if (msg->argc == 3) {
    configfile = msg->argv[1];
    mode = msg->argv[2];
    if (mode != "build" && mode != "part" && mode != "order") {
      CkPrintf("Error: mode %s not valid\n"
               "       valid modes: build, order\n", mode.c_str());
      //CkExit();
      initok = false;
    }
  }
  else {
    CkPrintf("Usage: [config file] [mode]\n");
    //CkExit();
    initok = false;
  }
  delete msg;

  if (initok) {
    // Parsing config
    if (ParseConfig(configfile)) {
      CkPrintf("Error loading config...\n");
      //CkExit();
      initok = false;
    }
  }

  if (initok) {
    // Basic error check
    if (netfiles != CkNumPes()) {
      CkPrintf("Error: netfiles (%d) does not match CkNumPes (%d)\n"
               "       Use '+p%d' to set %d PEs in Charm++\n",
               netfiles, CkNumPes(), netfiles, netfiles);
      //CkExit();
      initok = false;
    }
  }

  if (initok) {
    // Display configuration information
    CkPrintf("Loaded config from %s\n"
             "  Data Files (netfiles):     %" PRIidx "\n"
             "  Network Parts (netparts):  %" PRIidx "\n",
             configfile.c_str(), netfiles, netparts);

    CkPrintf("Initializing models\n");
  }

  if (initok) {
    // Read model information
    if (ReadModel()) {
      CkPrintf("Error loading models...\n");
      //CkExit();
      initok = false;
    }
  }

  if (initok) {
    // Print out model information
    for (std::size_t i = 0; i < models.size(); ++i) {
      CkPrintf("  Model: %" PRIidx "   ModName: %s   Type: %s   States: %d   Sticks: %d\n",
               i+1, models[i].modname.c_str(), graphtype[models[i].type].c_str(), models[i].statetype.size(), models[i].sticktype.size());
    }
    // TODO: concatenate these onto one line
    for (std::size_t i = 0; i < datafiles.size(); ++i) {
      CkPrintf("  Datafiles: %" PRIidx "   Filename: %s\n", i, datafiles[i].c_str());
    }

    // Set up control flags
    buildflag = true;
    partsflag = true;
    metisflag = true;
    orderflag = true;
    writeflag = true;
    if (mode == "build") {
      partsflag = false;
      metisflag = false;
      orderflag = false;
    }
    else if (mode == "part") {
      buildflag = false;
      metisflag = false;
      orderflag = false;
    }
    else if (mode == "order") {
      buildflag = false;
      partsflag = false;
    }
    // MPI Glue
    mainProxy = thisProxy;

    // Build model message
    mModel *mmodel = BuildModel();

    // Set Round Robin Mapping
    CkArrayOptions opts(netfiles);
    CProxy_RRMap rrMap = CProxy_RRMap::ckNew();
    opts.setMap(rrMap);

    // Create chare array
    CkCallback *cb = new CkCallback(CkReductionTarget(Main, ReturnControl), thisProxy);
    genet = CProxy_GeNet::ckNew(mmodel, opts);
    genet.ckSetReductionClient(cb);
  }
  else {
    // Initialization not okay
    CkExit();
  }
}

// Main migration
//
Main::Main(CkMigrateMessage* msg) {
  delete msg;
}

// Main control
//
void Main::Control() {
  if (mode == "build") {
    if (buildflag) {
      CkPrintf("Building network\n");
      buildflag = false;

      // Read graph information
      if (ReadGraph()) {
        CkPrintf("Error loading graph...\n");
        CkExit();
      }
      mGraph *mgraph = BuildGraph();

      // Build Network
      CkCallback *cb = new CkCallback(CkReductionTarget(Main, Control), thisProxy);
      genet.Build(mgraph);
      genet.ckSetReductionClient(cb);
    }
    else if (writeflag) {
      CkPrintf("Writing network\n");
      writeflag = false;

      CkCallback *cb = new CkCallback(CkIndex_Main::Halt(NULL), thisProxy);
      genet.Write(*cb);
    }
  }
  else if (mode == "part") {
    if (partsflag) {
      CkPrintf("Partitioning network\n");
      partsflag = false;

      CkCallback *cb = new CkCallback(CkReductionTarget(Main, ReturnControl), thisProxy);
      genet.SetPartition();
      genet.ckSetReductionClient(cb);
    }
    else {
      // Already partitioned return control to MPI
      CkExit();
    }
  }
  else if (mode == "order") {
    if (metisflag) {
      CkPrintf("Reading network\n");
      metisflag = false;

      if (ReadMetis()) {
        CkPrintf("Error loading metis...\n");
        CkExit();
      }
      mMetis *mmetis = BuildMetis();

      CkCallback *cb = new CkCallback(CkReductionTarget(Main, Control), thisProxy);
      genet.Read(mmetis);
      genet.ckSetReductionClient(cb);
    }
    else if (orderflag) {
      CkPrintf("Reordering network\n");
      orderflag = false;

      CkCallback *cb = new CkCallback(CkReductionTarget(Main, Control), thisProxy);
      genet.ScatterPart();
      genet.ckSetReductionClient(cb);
    }
    else if (writeflag) {
      CkPrintf("Writing network\n");
      writeflag = false;

      CkCallback *cb = new CkCallback(CkIndex_Main::Halt(NULL), thisProxy);
      genet.Write(*cb);
    }
  }
}

// Main Return Rontrol to MPI
//
void Main::ReturnControl() {
  CkExit();
}

// Main Stop
//
void Main::Halt(CkReductionMsg *msg) {
  CkPrintf("Finalizing network\n");
 
  // Save network part distribution to local
  netdist.clear();
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(dist_t); ++i) {
    netdist.push_back(*((dist_t *)msg->getData()+i));
  }
  CkAssert(netdist.size() == netparts);
  // cleanup
  delete msg;
  
  // Write distribution
  if (WriteDist()) {
    CkPrintf("Error writing distribution...\n");
    CkExit();
  }

  CkExit();
}


/**************************************************************************
* Generate Network
**************************************************************************/

// GeNet Constructor
//
GeNet::GeNet(mModel *msg) {
  // Bookkeeping
  datidx = thisIndex;
  idx_t ndiv = netparts/netfiles;
  idx_t nrem = netparts%netfiles;
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
  models.resize(msg->nmodel);
  for (std::size_t i = 0; i < models.size(); ++i) {
    // modname
    models[i].modname = std::string(msg->modname + msg->xmodname[i], msg->modname + msg->xmodname[i+1]);
    modname[i+1] = models[i].modname;
    modmap[models[i].modname] = i+1;
    // type
    models[i].type = msg->type[i];
    // prepare containers
    models[i].statetype.resize(msg->xstatetype[i+1] - msg->xstatetype[i]);
    models[i].stateparam.resize(msg->xstatetype[i+1] - msg->xstatetype[i]);
    for (std::size_t j = 0; j < models[i].statetype.size(); ++j) {
      // statetype
      models[i].statetype[j] = msg->statetype[msg->xstatetype[i] + j];
      switch (models[i].statetype[j]) {
        case RNGTYPE_CONST:
          models[i].stateparam[j].resize(RNGPARAM_CONST);
          break;
        case RNGTYPE_UNIF:
          models[i].stateparam[j].resize(RNGPARAM_UNIF);
          break;
        case RNGTYPE_UNINT:
          models[i].stateparam[j].resize(RNGPARAM_UNINT);
          break;
        case RNGTYPE_NORM:
          models[i].stateparam[j].resize(RNGPARAM_NORM);
          break;
        case RNGTYPE_BNORM:
          models[i].stateparam[j].resize(RNGPARAM_BNORM);
          break;
        case RNGTYPE_LBNORM:
          models[i].stateparam[j].resize(RNGPARAM_LBNORM);
          break;
        case RNGTYPE_LIN:
          models[i].stateparam[j].resize(RNGPARAM_LIN);
          break;
        case RNGTYPE_BLIN:
          models[i].stateparam[j].resize(RNGPARAM_BLIN);
          break;
        case RNGTYPE_FILE:
          models[i].stateparam[j].resize(RNGPARAM_FILE);
          break;
        default:
          CkPrintf("Error: unknown statetype\n");
          break;
      }
      for (std::size_t s = 0; s < models[i].stateparam[j].size(); ++s) {
        models[i].stateparam[j][s] = msg->stateparam[jstateparam++];
      }
    }
    // prepare containers
    models[i].sticktype.resize(msg->xsticktype[i+1] - msg->xsticktype[i]);
    models[i].stickparam.resize(msg->xsticktype[i+1] - msg->xsticktype[i]);
    for (std::size_t j = 0; j < models[i].sticktype.size(); ++j) {
      // sticktype
      models[i].sticktype[j] = msg->sticktype[msg->xsticktype[i] + j];
      switch (models[i].sticktype[j]) {
        case RNGTYPE_CONST:
          models[i].stickparam[j].resize(RNGPARAM_CONST);
          break;
        case RNGTYPE_UNIF:
          models[i].stickparam[j].resize(RNGPARAM_UNIF);
          break;
        case RNGTYPE_UNINT:
          models[i].stickparam[j].resize(RNGPARAM_UNINT);
          break;
        case RNGTYPE_NORM:
          models[i].stickparam[j].resize(RNGPARAM_NORM);
          break;
        case RNGTYPE_BNORM:
          models[i].stickparam[j].resize(RNGPARAM_BNORM);
          break;
        case RNGTYPE_LBNORM:
          models[i].stickparam[j].resize(RNGPARAM_LBNORM);
          break;
        case RNGTYPE_LIN:
          models[i].stickparam[j].resize(RNGPARAM_LIN);
          break;
        case RNGTYPE_BLIN:
          models[i].stickparam[j].resize(RNGPARAM_BLIN);
          break;
        case RNGTYPE_FILE:
          models[i].stickparam[j].resize(RNGPARAM_FILE);
          break;
        default:
          CkPrintf("Error: unknown statetype\n");
          break;
      }
      for (std::size_t s = 0; s < models[i].stickparam[j].size(); ++s) {
        models[i].stickparam[j][s] = msg->stickparam[jstickparam++];
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

  // cleanup
  delete msg;

  // Try to do MPI stuff
  // MPI initialized by Charm++ if using MPI for communication
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  CkPrintf("GeNet PE: %d   MPI: %d/%d (%s)\n", datidx, world_rank, world_size, processor_name);
  
  // Initialize coordination lists
  adjcyreq.clear();
  ordering.clear();

  // Everything initialized correctly
  GeNet_UnSetDoneFlag();
  // return control to main
  contribute(0, NULL, CkReduction::nop);
}

// GeNet Migration
GeNet::GeNet(CkMigrateMessage* msg) {
  delete msg;
}

// GeNet Destructor
//
GeNet::~GeNet() {
}


/**************************************************************************
* MPI Glue Code
**************************************************************************/

// Main control loop
//
void GeNet_MainControl() {
  // Main is run on PE 0
  if (CkMyPe() == 0) {
    mainProxy.Control();
  }
  // Turn control over to Charm++
  StartCharmScheduler();
}

// Control Flags
//
void GeNet::SetPartition() {
  GeNet_UnSetDoneFlag();
  GeNet_SetPartFlag();

  // return control to main
  contribute(0, NULL, CkReduction::nop);
}


/**************************************************************************
* Charm++ Definitions
**************************************************************************/
#include "genet.def.h"
