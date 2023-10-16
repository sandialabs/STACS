/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "stacs.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CkGroupID mCastGrpId;
/*readonly*/ unsigned int randseed;
/*readonly*/ std::string netwkdir;
/*readonly*/ int netparts;
/*readonly*/ int netfiles;
/*readonly*/ std::string filebase;
/*readonly*/ std::string fileload;
/*readonly*/ std::string filesave;
/*readonly*/ std::string recordir;
/*readonly*/ std::string groupdir;
/*readonly*/ tick_t tstep;
/*readonly*/ idx_t nevtday;
/*readonly*/ idx_t intdisp;
/*readonly*/ idx_t intrec;
/*readonly*/ idx_t intbal;
/*readonly*/ idx_t intsave;
/*readonly*/ tick_t tmax;
/*readonly*/ tick_t tepisode;
/*readonly*/ idx_t episodes;
/*readonly*/ int grpminlen;
/*readonly*/ tick_t grpmaxdur;
/*readonly*/ idx_t grpvtxmin;
/*readonly*/ idx_t grpvtxmax;


/**************************************************************************
* Main
**************************************************************************/

// Main entry point
//
Main::Main(CkArgMsg *msg) {
  // Display title
  CkPrintf("\nSimulation Tool for Asynchronous Cortical Streams (stacs)\n");

  // Command line arguments
  std::string configfile;
  if (msg->argc < 2) {
    configfile = std::string(CONFIG_DEFAULT);
    runmode = std::string(RUNMODE_EMPTY);
  }
  else if (msg->argc == 2) {
    configfile = std::string(msg->argv[1]);
    runmode = std::string(RUNMODE_EMPTY);
  }
  else if (msg->argc == 3) {
    configfile = msg->argv[1];
    runmode = msg->argv[2];
  }
  else {
    CkPrintf("Usage: [config file] [run mode]\n");
    CkExit();
  }
  delete msg;
    
  // Read configuration
  if (ReadConfig(configfile)) {
    CkPrintf("Error reading configuration...\n");
    CkExit();
  }

  // Read model information
  if (ReadModel()) {
    CkPrintf("Error reading model information...\n");
    CkExit();
  }
      
  // Read graph information
  if (ReadGraph()) {
    CkPrintf("Error reading graph specification...\n");
    CkExit();
  }

  //Read graph distribution files (for built networks)
  if (runmode != std::string(RUNMODE_BUILD) &&
      runmode != std::string(RUNMODE_BUILDSIM)) {
    if (ReadDist()) {
      CkPrintf("Error reading graph distribution...\n");
      CkExit();
    }

    // Convert group percentage to index
    grpvtxmin = (idx_t)(std::floor(grpvtxminreal * netdist[netparts].nvtx));
    grpvtxmax = (idx_t)(std::floor(grpvtxmaxreal * netdist[netparts].nvtx));
  }

  // Charm information
  mainProxy = thisProxy;
  mCastGrpId = CProxy_CkMulticastMgr::ckNew();
  real_t netpe = (real_t)netparts/CkNumPes();
  if (netpe < 1) { netpe = 1; }
  
  // Display configuration information
  CkPrintf("  STACS Run Mode      (runmode): %s\n"
           "  Network Data Files (netfiles): %d\n"
           "  Network Partitions (netparts): %d\n"
           "  Charm++ Processing Elements  : %d\n"
           "  Network Partitions per PE    : %.2g\n",
           runmode.c_str(), netfiles, netparts, CkNumPes(), netpe);
  
  // Netdata and network chare arrays
  ninit = 0;
  cinit = 0;
  // Set Round Robin Mapping for data
  // Ideally set this to one chare per compute node
  CkArrayOptions netdataopts(netfiles);
  CProxy_RRMap rrMap = CProxy_RRMap::ckNew();
  netdataopts.setMap(rrMap);
  // Create chare arrays with model information
  ++ninit;
  //mModname *mmodname = BuildModname();
  mModname *mmodname = BuildModname();
  netdata = CProxy_Netdata::ckNew(mmodname, netdataopts);
  // Initialize network chares as well
  ++ninit;
  mModel *mmodel = BuildModel();
  network = CProxy_Network::ckNew(mmodel, netparts);
  // Set callback to return control to main
  CkCallback cbinit(CkReductionTarget(Main, Init), mainProxy);
  netdata.ckSetReductionClient(&cbinit);
  network.ckSetReductionClient(&cbinit);

  // bookkeeping for program control
  readflag = true;
  loadflag = true;
  buildflag = false;
  moveflag = false;
  writeflag = false;
  lbflag = true;
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
  // Determine program control based on the run mode
  if (++cinit == ninit) {
    ninit = 0;
    cinit = 0;
    
    // Populate network (build or from disk)
    if (readflag) {
      readflag = false;

      // Pass some metadata to network chares
      // now that netdata has been initialized
      ++ninit;
      network.InitProxy(netdata);
      
      // Network needs to be built
      if (runmode == std::string(RUNMODE_BUILD)     ||
          runmode == std::string(RUNMODE_BUILDSIM)) {
        buildflag = true;
        CkPrintf("Reading datafiles\n");
        // Load data files (if any) from disk
        ++ninit;
        netdata.LoadFile();
        ++ninit;
        mGraph *mgraph = BuildGraph();
        network.OrderGraph(mgraph);
      }
      // Network needs to be read from disk
      else {
        CkPrintf("Reading network\n");

        // Load network from disk
        ++ninit;
        mDist *mdist = BuildDist();
        netdata.LoadData(mdist);

        // Optional loading of information based on runmode
        // Partitioning information
        if (runmode == std::string(RUNMODE_MIGRATE) ||
            runmode == std::string(RUNMODE_REPART)) {
          moveflag = true;
          ++ninit;
          netdata.LoadRepart();
        }
      }
    }

    // Load network partitions
    else if (loadflag) {
      loadflag = false;
      
      // Load any data files
      if (buildflag) {
        CkPrintf("Loading datafiles\n");
        ++ninit;
        network.LoadFile();
      }

      // Regular loading of network
      else {
        CkPrintf("Loading network\n");
        ++ninit;
        network.LoadData();

        // Optional loading of information based on runmode
        // Partitioning information
        if (runmode == std::string(RUNMODE_MIGRATE) ||
            runmode == std::string(RUNMODE_REPART)) {
          ++ninit;
          network.LoadRepart();
        }
        // Polychronous group information
        else if (runmode == std::string(RUNMODE_FINDGROUP) ||
                 runmode == std::string(RUNMODE_ESTIMATE)) {
          ++ninit;
          mGroup *mgroup = BuildGroup();
          // TODO: Should this stay as network?
          network.LoadGroup(mgroup);
        }
      }
    }

    // Build network (from graph and files)
    else if (buildflag) {
      buildflag = false;
      writeflag = true;
      CkPrintf("Building network\n");

      // Start timer
      wcstart = std::chrono::system_clock::now();

      // Build Network
      ++ninit;
      network.Build();
    }

    // Reorder an already built network
    else if (moveflag) {
      moveflag = false;
      writeflag = true;
      CkPrintf("Repartitioning network\n");
      
      // Start timer
      wcstart = std::chrono::system_clock::now();
      
      ++ninit;
      network.Repart();
    }

    // Write out any generated network
    else if (writeflag) {
      writeflag = false;
      CkPrintf("Writing network\n");
      
      // Stop timer
      wcstop = std::chrono::system_clock::now();
      // Print timing
      std::chrono::duration<real_t> wctime = std::chrono::duration_cast<std::chrono::milliseconds>(wcstop - wcstart);
      CkPrintf("  Elapsed time (wall clock): %" PRIrealsec " seconds\n", wctime.count());

      // Determine if exiting afterwards or not
      if (runmode == std::string(RUNMODE_BUILD)   ||
          runmode == std::string(RUNMODE_MIGRATE) ||
          runmode == std::string(RUNMODE_REPART)) {
        // Halting coordination
        chalt = 0;
        nhalt = 0;
        // Set callback for halting (actually not needed)
        CkCallback cbhalt(CkReductionTarget(Main, Halt), mainProxy);
        network.ckSetReductionClient(&cbhalt);
        // Write network to disk
        ++nhalt;
        network.SaveCloseNetwork();
      }
      else {
        if (runmode == std::string(RUNMODE_BUILDSIM)) {
          // Write network to disk
          ++ninit;
          network.SaveBuild();
          // Convert runmode to simulate
          runmode = std::string(RUNMODE_SIMULATE);
        }
      }
    }

    // Start running the network
    else {
      CkPrintf("Initializing network\n");
        
      // Initialize coordination
      CkCallback cbstart(CkReductionTarget(Main, Start), mainProxy);
      CkCallback cbcheck(CkReductionTarget(Main, Check), mainProxy);

      // Initialize network
      ++ninit;
      network.InitNetwork(runmode, cbcheck);
      network.ckSetReductionClient(&cbstart);
#ifdef STACS_WITH_YARP
      // stream
      ++ninit;
      mVtxs *mvtxs = BuildVtxs();
      stream = CProxy_Stream::ckNew(mvtxs);
      stream.ckSetReductionClient(&cbstart);
#endif
    }
  }
}

// Coordination for starting simulation, initialized network partitions
//
void Main::Start() {
  if (++cinit == ninit) {
    ninit = 0;
    cinit = 0;
    
    // Start simulation
    CkCallback cbstop(CkReductionTarget(Main, Stop), mainProxy);
    network.ckSetReductionClient(&cbstop);

    if (runmode == RUNMODE_SIMULATE) {
      if (episodic) {
        CkPrintf("  Random Number Seed  (randseed): %u\n"
                 "  Network Plasticity   (plastic): %s\n"
                 "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
                 "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
                 "  Time per Episode    (tepisode): %" PRIrealms "ms\n"
                 "  Number of Episodes  (episodes): %" PRIidx "\n",
                 randseed, (plastic ? "yes" : "no"),
                 (((real_t) tstep)/TICKS_PER_MS), teventq,
                 (((real_t) tepisode)/TICKS_PER_MS), episodes);
      }
      else {
        CkPrintf("  Random Number Seed  (randseed): %u\n"
                 "  Network Plasticity   (plastic): %s\n"
                 "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
                 "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
                 "  Display Interval    (tdisplay): %" PRIrealms "ms\n"
                 "  Recording Interval   (trecord): %" PRIrealms "ms\n"
                 "  Rebalance Interval  (tbalance): %" PRIrealms "ms\n"
                 "  Save State Interval    (tsave): %" PRIrealms "ms\n"
                 "  Max Simulation Time     (tmax): %" PRIrealms "ms\n",
                 randseed, (plastic ? "yes" : "no"),
                 (((real_t) tstep)/TICKS_PER_MS), teventq, tdisplay,
                 trecord, tbalance, tsave, (((real_t) tmax)/TICKS_PER_MS));
      }
      // Set compute cycle
      netcycle = CkCallback(CkIndex_Network::CycleSim(), network);
    }
    else if (runmode == RUNMODE_SIMGPU) {
      if (episodic) {
        CkPrintf("  Random Number Seed  (randseed): %u\n"
                 "  Network Plasticity   (plastic): %s\n"
                 "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
                 "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
                 "  Time per Episode    (tepisode): %" PRIrealms "ms\n"
                 "  Number of Episodes  (episodes): %" PRIidx "\n",
                 randseed, (plastic ? "yes" : "no"),
                 (((real_t) tstep)/TICKS_PER_MS), teventq,
                 (((real_t) tepisode)/TICKS_PER_MS), episodes);
      }
      else {
        CkPrintf("  Random Number Seed  (randseed): %u\n"
                 "  Network Plasticity   (plastic): %s\n"
                 "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
                 "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
                 "  Display Interval    (tdisplay): %" PRIrealms "ms\n"
                 "  Recording Interval   (trecord): %" PRIrealms "ms\n"
                 "  Save State Interval    (tsave): %" PRIrealms "ms\n"
                 "  Max Simulation Time     (tmax): %" PRIrealms "ms\n",
                 randseed, (plastic ? "yes" : "no"),
                 (((real_t) tstep)/TICKS_PER_MS), teventq, tdisplay,
                 trecord, tsave, (((real_t) tmax)/TICKS_PER_MS));
      }
      // Set compute cycle
      netcycle = CkCallback(CkIndex_Network::CycleSimGPU(), network);
    }
    else if (runmode == RUNMODE_FINDGROUP) {
      // collect active models
      std::string grpactivestring;
      for (std::size_t i = 0; i < grpactives.size(); ++i) {
        std::ostringstream grpactive;
        grpactive << " " << grpactives[i];
        grpactivestring.append(grpactive.str());
      }
      // collect mother models
      std::string grpmotherstring;
      for (std::size_t i = 0; i < grpmothers.size(); ++i) {
        std::ostringstream grpmother;
        grpmother << " " << grpmothers[i];
        grpmotherstring.append(grpmother.str());
      }
      // collect anchor models
      std::string grpanchorstring;
      for (std::size_t i = 0; i < grpanchors.size(); ++i) {
        std::ostringstream grpanchor;
        grpanchor << " " << grpanchors[i];
        grpanchorstring.append(grpanchor.str());
      }
      // Compute vertices evaluated
      CkPrintf("  Random Number Seed  (randseed): %u\n"
               "  Group Directory     (groupdir): %s\n"
               "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
               "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
               "  Polychronization   (grpactive):%s\n"
               "                     (grpmother):%s\n"
               "                     (grpanchor):%s\n"
               "                     (grpminlen): %d\n"
               "                     (grpmaxdur): %" PRIrealms "ms\n"
               "  Evaluated Vertices (grpvtxmin): %.6g (%" PRIidx ")\n"
               "                     (grpvtxmax): %.6g (%" PRIidx ")\n",
               randseed, groupdir.c_str(), (((real_t) tstep)/TICKS_PER_MS), teventq,
               grpactivestring.c_str(), grpmotherstring.c_str(), grpanchorstring.c_str(),
               grpminlen, (((real_t) grpmaxdur)/TICKS_PER_MS),
               ((real_t)grpvtxminreal), grpvtxmin, ((real_t)grpvtxmaxreal), grpvtxmax);
      // Set compute cycle
      netcycle = CkCallback(CkIndex_Network::CycleGroup(), network);
    }
    else if (runmode == RUNMODE_ESTIMATE) {
      if (episodic) {
        CkPrintf("  Random Number Seed  (randseed): %u\n"
                 "  Group Directory     (groupdir): %s\n"
                 "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
                 "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
                 "  Time per Episode    (tepisode): %" PRIrealms "ms\n"
                 "  Number of Episodes  (episodes): %" PRIidx "\n",
                 randseed, groupdir.c_str(), (((real_t) tstep)/TICKS_PER_MS), teventq,
                 (((real_t) tepisode)/TICKS_PER_MS), episodes);
      }
      else {
        CkPrintf("  Random Number Seed  (randseed): %u\n"
                 "  Group Directory     (groupdir): %s\n"
                 "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
                 "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
                 "  Display Interval    (tdisplay): %" PRIrealms "ms\n"
                 "  Recording Interval   (trecord): %" PRIrealms "ms\n"
                 "  Max Simulation Time     (tmax): %" PRIrealms "ms\n",
                 randseed, groupdir.c_str(), (((real_t) tstep)/TICKS_PER_MS), teventq,
                 tdisplay, trecord, (((real_t) tmax)/TICKS_PER_MS));
      }
      // Set compute cycle
      netcycle = CkCallback(CkIndex_Network::CycleEst(), network);
    }

#ifdef STACS_WITH_YARP
    // Open RPC port
    stream.OpenRPC(network, netcycle, rpcpause);

    // Start computing
    if (rpcpause) {
      // Start paused (wait for rpc)
      CkPrintf("Starting (paused)\n");
    }
    else {
      CkPrintf("Starting\n");
      network.StartNetwork();
    }
#else
    CkPrintf("Starting\n");
    network.StartNetwork();
#endif

    // Start timer
    wcstart = std::chrono::system_clock::now();
  }
}

// Coordination of continuing simulation after repartitioning
//
void Main::Check() {
  if (lbflag) {
    lbflag = false;
    CkPrintf("Repartitioning\n");
    network.Repart();
  }
  else {
    lbflag = true;
    network.StartNetwork();
  }
}


/**************************************************************************
* STACS Shutdown
**************************************************************************/

// Coordination for stopping simulation
//
void Main::Stop() {
  CkPrintf("Stopping\n");
  
  // Stop timer
  wcstop = std::chrono::system_clock::now();
  // Print timing
  std::chrono::duration<real_t> wctime = std::chrono::duration_cast<std::chrono::milliseconds>(wcstop - wcstart);
  CkPrintf("  Elapsed time (wall clock): %" PRIrealsec " seconds\n", wctime.count());

  CkPrintf("Finalizing network\n");
#ifdef STACS_WITH_YARP
  // Close RPC port
  stream.CloseRPC();
#endif

  // Halting coordination
  chalt = 0;
  nhalt = 0;
  // Set callback for halting
  CkCallback cbhalt(CkReductionTarget(Main, Halt), mainProxy);
  netdata.ckSetReductionClient(&cbhalt);

  // Save data from network parts to output files
  if (runmode == RUNMODE_SIMULATE ||
      runmode == RUNMODE_SIMGPU) {
    network.SaveFinalRecord();
    ++nhalt;
    if (plastic) {
      network.SaveCloseNetwork();
      ++nhalt;
    }
    else {
      network.CloseNetwork();
      ++nhalt;
    }
  }
  else if (runmode == RUNMODE_ESTIMATE) {
    network.SaveFinalEstimate();
    ++nhalt; ++nhalt;
    network.CloseNetwork();
    ++nhalt;
  }
  else {
    network.CloseNetwork();
    ++nhalt;
  }
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
