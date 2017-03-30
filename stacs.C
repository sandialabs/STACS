/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "stacs.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CkGroupID mCastGrpId;
/*readonly*/ idx_t npdat;
/*readonly*/ idx_t npnet;
/*readonly*/ std::string netdir;
/*readonly*/ std::string recdir;
/*readonly*/ std::string filebase;
/*readonly*/ std::string fileout;
/*readonly*/ tick_t tmax;
/*readonly*/ tick_t tstep;
/*readonly*/ tick_t tqueue;
/*readonly*/ tick_t tcheck;
/*readonly*/ tick_t trecord;
/*readonly*/ tick_t tdisplay;
/*readonly*/ idx_t equeue;
/*readonly*/ idx_t rngseed;
/*readonly*/ idx_t runmode;


/**************************************************************************
* Main
**************************************************************************/

// Main entry point
//
Main::Main(CkArgMsg *msg) {
  // Display title
  CkPrintf("Simulation Tool for Asynchrnous Cortical Streams (stacs)\n");

  // Command line arguments
  std::string configfile;
  if (msg->argc < 2) {
    configfile = "config.yml"; // default
  }
  else {
    configfile = msg->argv[1];
  }
  delete msg;

  // Parsing config
  if (ParseConfig(configfile)) {
    CkPrintf("Error loading config...\n");
    CkExit();
  }

  // Charm information
  real_t netpe = (real_t)npnet/CkNumPes();
  if (netpe < 1) { netpe = 1; }

  // Display configuration information
  if (runmode == RUNMODE_SIM) {
    CkPrintf("Loaded config from %s\n"
             "  Data Files (npdat):     %" PRIidx "\n"
             "  Network Parts (npnet):  %" PRIidx "\n"
             "  Processing Elements:    %d\n"
             "  Network Parts per PE:   %.2g\n"
             "  Simulation run mode:    sim\n"
             "  Total Simulation Time (tmax): %" PRItick "\n"
             "  Simulation Time Step (tstep): %" PRItick "\n"
             "  Checkpoint Interval (tcheck): %" PRItick "\n",
             configfile.c_str(), npdat, npnet,
             CkNumPes(), netpe, tmax, tstep, tcheck);
  }
  else if (runmode == RUNMODE_PNG) {
    std::string pngmodstr;
    // collect active models
    for (std::size_t i = 0; i < pngmods.size(); ++i) {
      std::ostringstream pngmod;
      pngmod << " " << pngmods[i];
      pngmodstr.append(pngmod.str());
    }
    CkPrintf("Loaded config from %s\n"
             "  Data Files (npdat):     %" PRIidx "\n"
             "  Network Parts (npnet):  %" PRIidx "\n"
             "  Processing Elements:    %d\n"
             "  Network Parts per PE:   %.2g\n"
             "  Simulation run mode:    png\n"
             "  Polychronizing models: %s\n",
             configfile.c_str(), npdat, npnet,
             CkNumPes(), netpe, pngmodstr.c_str());
  }

  // Read vertex distribution
  CkPrintf("Initializing simulation\n");
  if (ReadDist()) {
    CkPrintf("Error loading distribution...\n");
    CkExit();
  }
  // Read model information
  if (ReadModel()) {
    CkPrintf("Error loading models...\n");
    CkExit();
  }
  
  // Setup Charm++ variables
  mainProxy = thisProxy;
  mCastGrpId = CProxy_CkMulticastMgr::ckNew();

  // Initialize coordination
  cinit = 0;
  ninit = 0;

  // Setup chare arrays
  CkCallback *cb = new CkCallback(CkReductionTarget(Main, Init), mainProxy);
  // netdata
  ++ninit;
  mDist *mdist = BuildDist();
  netdata = CProxy_Netdata::ckNew(mdist, npdat);
  netdata.ckSetReductionClient(cb);
  // network
  ++ninit;
  mModel *mmodel = BuildModel();
  network = CProxy_Network::ckNew(mmodel, npnet);
  network.ckSetReductionClient(cb);
#ifdef STACS_WITH_YARP
  // stream
  ++ninit;
  mVtxs *mvtxs = BuildVtxs();
  stream = CProxy_Stream::ckNew(mvtxs);
#endif
}

// Main migration
//
Main::Main(CkMigrateMessage *msg) {
  delete msg;
}


/**************************************************************************
* Simulation Startup
**************************************************************************/

// Coordination for file input, initialized chare arrays
//
void Main::Init() {
  // Wait on initialization
  if (++cinit == ninit) {
    CkPrintf("Setting up network parts\n");

    // Load data from input files to network parts
    CkCallback *cb = new CkCallback(CkReductionTarget(Main, Start), mainProxy);
    network.ckSetReductionClient(cb);

    if (runmode == RUNMODE_SIM) {
      cbcycle = CkCallback(CkIndex_Network::CycleSim(), network);
      network.InitSim(netdata);
    }
    else if (runmode == RUNMODE_PNG) {
      cbcycle = CkCallback(CkIndex_Network::CyclePNG(), network);
      network.InitPNG(netdata);
    }
    
#ifdef STACS_WITH_YARP
    // Open RPC port
    stream.OpenRPC(network, cbcycle, startpaused);
#endif
  }
}

// Coordination for starting simulation, network partition setup
//
void Main::Start() {
  // Start simulation
  CkCallback *cb = new CkCallback(CkReductionTarget(Main, Stop), mainProxy);
  network.ckSetReductionClient(cb);

#ifdef STACS_WITH_YARP
  if (startpaused) {
    // Start paused
    CkPrintf("Starting simulation (paused)\n");
  }
  else {
    CkPrintf("Starting simulation\n");
    cbcycle.send();
  }
#else
  CkPrintf("Starting simulation\n");
  cbcycle.send();
#endif

  // Start timer
  tstart = std::chrono::system_clock::now();
}


/**************************************************************************
* Simulation Shutdown
**************************************************************************/

// Coordination for stopping simulation
//
void Main::Stop() {
  CkPrintf("Stopping simulation\n");
  
  // Stop timer
  tstop = std::chrono::system_clock::now();
  // Print timing
  std::chrono::duration<real_t> tduration = std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
  CkPrintf("Elapsed time (wall clock): %" PRIrealsec " seconds\n", tduration.count());

#ifdef STACS_WITH_YARP
  // Close RPC port
  stream.CloseRPC();
#endif

  // Save data from network parts to output files
  chalt = nhalt = 0;
  if (runmode == RUNMODE_SIM) {
    network.SaveNetwork();
    ++nhalt;
    network.SaveRecord();
    ++nhalt;
  }
  else if (runmode == RUNMODE_PNG) {
    network.CloseNetwork();
    ++nhalt;
  }
  
  // Set callback for halting
  CkCallback *cb = new CkCallback(CkReductionTarget(Main, Halt), mainProxy);
  netdata.ckSetReductionClient(cb);
}

// Finish
//
void Main::Halt() {
  // Everything checks back to here
  if (++chalt == nhalt) {
  
    // Simulation complete
    CkExit();
  }
}

/**************************************************************************
* Charm++ Definitions
**************************************************************************/
#include "stacs.def.h"
