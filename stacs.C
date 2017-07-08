/**
 * Copyright (C) 2017 Felix Wang
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
/*readonly*/ std::string filebase;
/*readonly*/ std::string fileload;
/*readonly*/ std::string filesave;
/*readonly*/ std::string modeldir;
/*readonly*/ std::string recordir;
/*readonly*/ std::string groupdir;
/*readonly*/ idx_t rngseed;
/*readonly*/ tick_t tmax;
/*readonly*/ tick_t tstep;
/*readonly*/ tick_t tqueue;
/*readonly*/ tick_t tcheck;
/*readonly*/ tick_t trecord;
/*readonly*/ tick_t tdisplay;
/*readonly*/ idx_t nevtday;
/*readonly*/ int pnglength;
/*readonly*/ idx_t comprtmin;
/*readonly*/ idx_t comprtmax;
/*readonly*/ idx_t ntrials;
/*readonly*/ tick_t ttrial;


/**************************************************************************
* Main
**************************************************************************/

// Main entry point
//
Main::Main(CkArgMsg *msg) {
  // Display title
  CkPrintf("\nSimulation Tool for Asynchrnous Cortical Streams (stacs)\n");

  // Command line arguments
  std::string configfile;
  if (msg->argc < 2) {
    configfile = "networks/dummy.yml"; // default
  }
  else {
    configfile = msg->argv[1];
  }
  delete msg;

  // Read configuration
  if (ReadConfig(configfile)) {
    CkPrintf("Error reading configuration...\n");
    CkExit();
  }

  // Charm information
  real_t netpe = (real_t)npnet/CkNumPes();
  if (netpe < 1) { netpe = 1; }

  // Display configuration information
  CkPrintf("  STACS run mode: %s\n"
           "  Data Files (npdat):     %" PRIidx "\n"
           "  Network Parts (npnet):  %" PRIidx "\n"
           "  Processing Elements:    %d\n"
           "  Network Parts per PE:   %.2g\n",
           runmode.c_str(), npdat, npnet, CkNumPes(), netpe);

  // Read graph distribution
  if (ReadDist()) {
    CkPrintf("Error reading graph distribution...\n");
    CkExit();
  }
  // Read model information
  if (ReadModel()) {
    CkPrintf("Error reading model information...\n");
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
* STACS Startup
**************************************************************************/

// Coordination for file input, initialized chare arrays
//
void Main::Init() {
  // Wait on initialization
  if (++cinit == ninit) {
    CkPrintf("Initializing network\n");

    // Load data from input files to network parts
    CkCallback *cb = new CkCallback(CkReductionTarget(Main, Start), mainProxy);
    network.ckSetReductionClient(cb);

    if (runmode == RUNMODE_SIM) {
      CkPrintf("  Random Generator Seed (rngseed): %" PRIidx "\n"
               "  Total Simulation Time    (tmax): %" PRIrealms "ms\n"
               "  Simulation Time Step    (tstep): %" PRIrealms "ms\n"
               "  Event Queue Length     (tqueue): %" PRIrealms "ms\n"
               "  Checkpointing Interval (tcheck): %" PRIrealms "ms\n"
               "  Recording Interval    (trecord): %" PRIrealms "ms\n"
               "  Display Interval     (tdisplay): %" PRIrealms "ms\n"
               "  Network Plasticity (plasticity): %s\n",
               rngseed, ((real_t)(tmax/TICKS_PER_MS)),
               ((real_t)(tstep/TICKS_PER_MS)), ((real_t)(tqueue/TICKS_PER_MS)),
               ((real_t)(tcheck/TICKS_PER_MS)), ((real_t)(trecord/TICKS_PER_MS)),
               ((real_t)(tdisplay/TICKS_PER_MS)), (plasticity ? "on" : "off"));
      // Set compute cycle
      if (plasticity) {
        cbcycle = CkCallback(CkIndex_Network::CycleSimPlastic(), network);
        network.InitSimPlastic(netdata);
      }
      else {
        cbcycle = CkCallback(CkIndex_Network::CycleSimStatic(), network);
        network.InitSimStatic(netdata);
      }
    }
    else if (runmode == RUNMODE_EPS) {
      CkPrintf("  Random Generator Seed (rngseed): %" PRIidx "\n"
               "  Simulation Time Step    (tstep): %" PRIrealms "ms\n"
               "  Event Queue Length     (tqueue): %" PRIrealms "ms\n"
               "  Time per Episode       (ttrial): %" PRIrealms "ms\n"
               "  Number of Episodes    (ntrials): %" PRIidx "\n",
               rngseed, ((real_t)(tstep/TICKS_PER_MS)),
               ((real_t)(tqueue/TICKS_PER_MS)), ((real_t)(ttrial/TICKS_PER_MS)), ntrials);
      // Set compute cycle
      cbcycle = CkCallback(CkIndex_Network::CycleEpsPlastic(), network);
      network.InitEpsPlastic(netdata);
    }
    else if (runmode == RUNMODE_PNG) {
      std::string pnginfo;
      // collect active models
      for (std::size_t i = 0; i < pngactives.size(); ++i) {
        std::ostringstream pngactive;
        pngactive << " " << pngactives[i];
        pnginfo.append(pngactive.str());
      }
      // collect mother models
      for (std::size_t i = 0; i < pngmothers.size(); ++i) {
        std::ostringstream pngmother;
        pngmother << " (" << pngmothers[i] << ")";
        pnginfo.append(pngmother.str());
      }
      // collect anchor models
      for (std::size_t i = 0; i < pnganchors.size(); ++i) {
        std::ostringstream pnganchor;
        pnganchor << " <" << pnganchors[i] << ">";
        pnginfo.append(pnganchor.str());
      }
      CkPrintf("  Random Generator Seed (rngseed): %" PRIidx "\n"
               "  Simulation Time Step    (tstep): %" PRIrealms "ms\n"
               "  Event Queue Length     (tqueue): %" PRIrealms "ms\n"
               "  Computation range (inclusive)  : %" PRIidx " to %" PRIidx "\n"
               "  Polychronization Information   :%s\n",
               rngseed, ((real_t)(tstep/TICKS_PER_MS)),
               ((real_t)(tqueue/TICKS_PER_MS)),
               comprtmin, comprtmax, pnginfo.c_str());
      // Set compute cycle
      cbcycle = CkCallback(CkIndex_Network::CyclePNG(), network);
      network.InitPNG(netdata);
    }
    else if (runmode == RUNMODE_EST) {
      CkPrintf("  Random Generator Seed (rngseed): %" PRIidx "\n"
               "  Simulation Time Step    (tstep): %" PRIrealms "ms\n"
               "  Event Queue Length     (tqueue): %" PRIrealms "ms\n"
               "  Classes estimated     (ntrials): %" PRIidx "\n",
               rngseed, ((real_t)(tstep/TICKS_PER_MS)),
               ((real_t)(tqueue/TICKS_PER_MS)), ntrials);
      // Set compute cycle
      cbcycle = CkCallback(CkIndex_Network::CycleEstStatic(), network);
      network.InitEstStatic(netdata);
    }
    else if (runmode == RUNMODE_MON) {
      CkPrintf("  Random Generator Seed (rngseed): %" PRIidx "\n"
               "  Total Simulation Time    (tmax): %" PRIrealms "ms\n"
               "  Simulation Time Step    (tstep): %" PRIrealms "ms\n"
               "  Event Queue Length     (tqueue): %" PRIrealms "ms\n"
               "  Recording Interval    (trecord): %" PRIrealms "ms\n"
               "  Display Interval     (tdisplay): %" PRIrealms "ms\n",
               rngseed, ((real_t)(tmax/TICKS_PER_MS)),
               ((real_t)(tstep/TICKS_PER_MS)), ((real_t)(tqueue/TICKS_PER_MS)),
               ((real_t)(trecord/TICKS_PER_MS)), ((real_t)(tdisplay/TICKS_PER_MS)));
      // Set compute cycle
      cbcycle = CkCallback(CkIndex_Network::CycleMonStatic(), network);
      network.InitMonStatic(netdata);
    }
    
#ifdef STACS_WITH_YARP
    // Open RPC port
    stream.OpenRPC(network, cbcycle, startpaused);
#endif
  }
}

// Coordination for starting simulation, initialized network partitions
//
void Main::Start() {
  // Start simulation
  CkCallback *cb = new CkCallback(CkReductionTarget(Main, Stop), mainProxy);
  network.ckSetReductionClient(cb);

#ifdef STACS_WITH_YARP
  if (startpaused) {
    // Start paused
    CkPrintf("Starting (paused)\n");
  }
  else {
    CkPrintf("Starting\n");
    cbcycle.send();
  }
#else
  CkPrintf("Starting\n");
  cbcycle.send();
#endif

  // Start timer
  tstart = std::chrono::system_clock::now();
}


/**************************************************************************
* STACS Shutdown
**************************************************************************/

// Coordination for stopping simulation
//
void Main::Stop() {
  CkPrintf("Stopping\n");
  
  // Stop timer
  tstop = std::chrono::system_clock::now();
  // Print timing
  std::chrono::duration<real_t> tduration = std::chrono::duration_cast<std::chrono::milliseconds>(tstop - tstart);
  CkPrintf("  Elapsed time (wall clock): %" PRIrealsec " seconds\n", tduration.count());

  CkPrintf("Finalizing network\n");
#ifdef STACS_WITH_YARP
  // Close RPC port
  stream.CloseRPC();
#endif

  // Save data from network parts to output files
  chalt = nhalt = 0;
  if (runmode == RUNMODE_SIM) {
    if (plasticity) {
      network.SaveFinalNetwork();
      ++nhalt;
    }
    else {
      network.FinalizeNetwork();
      ++nhalt;
    }
    network.SaveFinalRecord();
    ++nhalt;
  }
  else if (runmode == RUNMODE_EPS) {
    network.SaveFinalNetwork();
    ++nhalt;
  }
  else if (runmode == RUNMODE_PNG || runmode == RUNMODE_EST) {
    network.FinalizeNetwork();
    ++nhalt;
  }
  else if (runmode == RUNMODE_MON) {
    network.FinalizeNetwork();
    ++nhalt;
    network.SaveFinalRecord();
    ++nhalt;
  }
  
  // Set callback for halting
  CkCallback *cb = new CkCallback(CkReductionTarget(Main, Halt), mainProxy);
  netdata.ckSetReductionClient(cb);
}

// Exit STACS
//
void Main::Halt() {
  // Everything checks back to here
  if (++chalt == nhalt) {
  
    // Halt
    CkExit();
  }
}


/**************************************************************************
* Charm++ Definitions
**************************************************************************/
#include "stacs.def.h"
