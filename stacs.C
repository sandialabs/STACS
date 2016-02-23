/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "stacs.h"
#include "pup_stl.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CkGroupID mCastGrpId;
/*readonly*/ std::string filebase;
/*readonly*/ std::string filein;
/*readonly*/ std::string fileout;
/*readonly*/ idx_t npdat;
/*readonly*/ idx_t npnet;
/*readonly*/ tick_t tmax;
/*readonly*/ tick_t tstep;
/*readonly*/ tick_t tdisplay;
/*readonly*/ tick_t tcheck;
/*readonly*/ tick_t trecord;
/*readonly*/ idx_t evtcal;
/*readonly*/ idx_t rngseed;
/*readonly*/ std::string rpcport;


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
  CkPrintf("Loaded config from %s\n"
           "  Data Files (npdat):     %" PRIidx "\n"
           "  Network Parts (npnet):  %" PRIidx "\n"
           "  Processing Elements:    %d\n"
           "  Network Parts per PE:   %.2g\n"
           "  Total Simulation Time (tmax): %" PRItick "\n"
           "  Simulation Time Step (tstep): %" PRItick "\n"
           "  Checkpoint Interval (tcheck): %" PRItick "\n",
           configfile.c_str(), npdat, npnet,
           CkNumPes(), netpe, tmax, tstep, tcheck);

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

#ifdef STACS_WITH_YARP
  // Initialize YARP
  yarp.init();
#endif

  // Setup chare arrays
  CkCallback *cb = new CkCallback(CkReductionTarget(Main, InitSim), mainProxy);
  // netdata
  ++ninit;
  mDist *mdist = BuildDist();
  netdata = CProxy_NetData::ckNew(mdist, npdat);
  netdata.ckSetReductionClient(cb);
  // network
  ++ninit;
  mModel *mmodel = BuildModel();
  network = CProxy_Network::ckNew(mmodel, npnet);
  network.ckSetReductionClient(cb);
#ifdef STACS_WITH_YARP
  // streamrpc
  ++ninit;
  mVtxDist *mvtxdist = BuildVtxDist();
  streamrpc = CProxy_StreamRPC::ckNew(mvtxdist);
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
void Main::InitSim() {
  // Wait on initialization
  if (++cinit == ninit) {
    CkPrintf("Setting up network parts\n");

    // Load data from input files to network parts
    CkCallback *cb = new CkCallback(CkReductionTarget(Main, StartSim), mainProxy);
    network.ckSetReductionClient(cb);
    network.LoadNetwork(netdata);
    
#ifdef STACS_WITH_YARP
    // Open RPC port
    streamrpc.Open(network);
#endif
  }
}

// Coordination for starting simulation, network partition setup
//
void Main::StartSim() {
  CkPrintf("Starting simulation\n");

  // Start simulation
  CkCallback *cb = new CkCallback(CkReductionTarget(Main, SaveSim), mainProxy);
  network.ckSetReductionClient(cb);
  network.Cycle();
}


/**************************************************************************
* Simulation Checkpointing
**************************************************************************/

// Coordination for checkpointing simulation
//
void Main::CheckSim(CkReductionMsg *msg) {
  CkPrintf("Checkpointing simulation\n");
  
  // Save network part distribution to local
  netdist.clear();
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(dist_t); ++i) {
    netdist.push_back(*((dist_t *)msg->getData()+i));
  }
  delete msg;

  // Write distribution
  if (WriteDist(true)) {
    CkPrintf("Error writing distribution...\n");
    CkExit();
  }
}

/**************************************************************************
* Simulation Shutdown
**************************************************************************/

// Coordination for stopping simulation
//
void Main::SaveSim() {
  CkPrintf("Saving simulation\n");

  // Save data from network parts to output files
  chalt = nhalt = 0;
  network.SaveNetwork();
  ++nhalt;
  network.SaveRecord();
  ++nhalt;
  
  // Set callback for halting
  CkCallback *cb = new CkCallback(CkReductionTarget(Main, Halt), mainProxy);
  netdata.ckSetReductionClient(cb);
}

// Coordination for file output
//
void Main::FiniSim(CkReductionMsg *msg) {
  CkPrintf("Stopping simulation\n");

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

// Halt
//
void Main::Halt() {
  // Everything checks back to here
  if (++chalt == nhalt) {
#ifdef STACS_WITH_YARP
    // Close RPC port
    streamrpc.Close();

    // Finalize YARP
    yarp.fini();
#endif
  
    // Simulation complete
    CkExit();
  }
}

/**************************************************************************
* Charm++ Definitions
**************************************************************************/
#include "stacs.def.h"
