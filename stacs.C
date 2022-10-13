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
  CkPrintf("\nSimulation Tool for Asynchrnous Cortical Streams (stacs)\n");

  // Command line arguments
  // TODO: distinguish between runmode (stacs, genet) and simmode (sim, fg, est)
  // For right now, the following still depends on the config.yml to get runmode
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

  
  // Netdata chare array (used in all runmodes)
  // Set Round Robin Mapping
  // TODO: Ideally set this to one chare per compute node
  CkArrayOptions netdataopts(netfiles);
  CProxy_RRMap rrMap = CProxy_RRMap::ckNew();
  netdataopts.setMap(rrMap);
  // Create chare array with model information
  mModel *mmodel = BuildModel();
  netdata = CProxy_Netdata::ckNew(mmodel, netdataopts);
  // Set callback to return control to main
  CkCallback cbcontrol(CkReductionTarget(Main, Control), mainProxy);
  netdata.ckSetReductionClient(&cbcontrol);

  // bookkeeping for program control
  buildflag = true;
  orderflag = true;
  readflag = true;
  writeflag = true;
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
void Main::Control() {
  // Determine program control based on the run mode
  // TODO: consolidate build and simulate into a single runmode (i.e. for smaller networks)

  // Network needs to be built
  if (runmode == std::string(RUNMODE_BUILD)) {
    if (buildflag) {
      CkPrintf("Building network\n");
      buildflag = false;

      // Start timer
      wcstart = std::chrono::system_clock::now();

      // Build Network
      CkCallback cbcontrol(CkReductionTarget(Main, Control), mainProxy);
      mGraph *mgraph = BuildGraph();
      netdata.Build(mgraph);
      netdata.ckSetReductionClient(&cbcontrol);
    }
    else if (writeflag) {
      CkPrintf("Writing network\n");
      writeflag = false;
      
      // Stop timer
      wcstop = std::chrono::system_clock::now();
      // Print timing
      std::chrono::duration<real_t> wctime = std::chrono::duration_cast<std::chrono::milliseconds>(wcstop - wcstart);
      CkPrintf("  Elapsed time (wall clock): %" PRIrealsec " seconds\n", wctime.count());

      // Halting coordination
      chalt = 0;
      nhalt = 0;
      // Set callback for halting (actually not needed)
      CkCallback cbhalt(CkReductionTarget(Main, Halt), mainProxy);
      netdata.ckSetReductionClient(&cbhalt);
      // Write network to disk
      netdata.SaveBuild();
      ++nhalt;
    }
  }

  // Reorder an already built network
  else if (runmode == std::string(RUNMODE_ORDER)) {
    if (readflag) {
      CkPrintf("Reading network\n");
      readflag = false;

      // Start timer
      wcstart = std::chrono::system_clock::now();

      // Read graph distribution files
      if (ReadDist()) {
        CkPrintf("Error reading graph distribution...\n");
        CkExit();
      }
      CkCallback cbcontrol(CkReductionTarget(Main, Control), mainProxy);
      mDist *mdist = BuildDist();
      netdata.ReadPart(mdist);
      netdata.ckSetReductionClient(&cbcontrol);
    }
    else if (orderflag) {
      CkPrintf("Reordering network\n");
      orderflag = false;

      // Scatter and gather part information
      CkCallback cbcontrol(CkReductionTarget(Main, Control), mainProxy);
      netdata.ScatterPart();
      netdata.ckSetReductionClient(&cbcontrol);
    }
    else if (writeflag) {
      CkPrintf("Writing network\n");
      writeflag = false;
      
      // Stop timer
      wcstop = std::chrono::system_clock::now();
      // Print timing
      std::chrono::duration<real_t> wctime = std::chrono::duration_cast<std::chrono::milliseconds>(wcstop - wcstart);
      CkPrintf("  Elapsed time (wall clock): %" PRIrealsec " seconds\n", wctime.count());

      // Halting coordination
      chalt = 0;
      nhalt = 0;
      // Set callback for halting (actually not needed)
      CkCallback cbhalt(CkReductionTarget(Main, Halt), mainProxy);
      netdata.ckSetReductionClient(&cbhalt);
      // Write network to disk
      netdata.SaveBuild();
      ++nhalt;
    }
  }

  // Simulate an already built network
  else if (runmode == std::string(RUNMODE_SIMULATE) ||
           runmode == std::string(RUNMODE_SIMGPU) ||
           runmode == std::string(RUNMODE_FINDGROUP) ||
           runmode == std::string(RUNMODE_ESTIMATE)) {
    // Read graph distribution files
    if (ReadDist()) {
      CkPrintf("Error reading graph distribution...\n");
      CkExit();
    }
    // Convert group percentage to index
    grpvtxmin = (idx_t)(std::floor(grpvtxminreal * netdist[netparts].nvtx));
    grpvtxmax = (idx_t)(std::floor(grpvtxmaxreal * netdist[netparts].nvtx));

    // Initialize coordination
    cinit = 0;
    ninit = 0;
    CkCallback cbinit(CkReductionTarget(Main, Init), mainProxy);

    // Loading network data from file
    ++ninit;
    mDist *mdist = BuildDist();
    netdata.LoadData(mdist);
    netdata.ckSetReductionClient(&cbinit);
    // Network chare array is used in simulation
    ++ninit;
    mModel *mmodel = BuildModel();
    network = CProxy_Network::ckNew(mmodel, netparts);
    network.ckSetReductionClient(&cbinit);
#ifdef STACS_WITH_YARP
    // stream
    ++ninit;
    mVtxs *mvtxs = BuildVtxs();
    stream = CProxy_Stream::ckNew(mvtxs);
#endif
  }
}
  
// Coordination for file input, initialized chare arrays
//
void Main::Init() {
  // Wait on initialization
  if (++cinit == ninit) {
    CkPrintf("Initializing network\n");

    // Load data from input files to network parts
    CkCallback cbstart(CkReductionTarget(Main, Start), mainProxy);
    network.ckSetReductionClient(&cbstart);

    if (runmode == RUNMODE_SIMULATE) {
      if (episodic) {
        CkPrintf("  Random Number Seed  (randseed): %u\n"
                 "  Network Plasticity   (plastic): %s\n"
                 "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
                 "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
                 "  Time per Episode    (tepisode): %" PRIrealms "ms\n"
                 "  Number of Episodes  (episodes): %" PRIidx "\n",
                 randseed, (plastic ? "yes" : "no"),
                 ((real_t)(tstep/TICKS_PER_MS)), teventq,
                 ((real_t)(tepisode/TICKS_PER_MS)), episodes);
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
                 ((real_t)(tstep/TICKS_PER_MS)), teventq, tdisplay,
                 trecord, tsave, ((real_t)(tmax/TICKS_PER_MS)));
      }
      // Set compute cycle
      netcycle = CkCallback(CkIndex_Network::CycleSim(), network);
      network.InitSim(netdata);
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
                 ((real_t)(tstep/TICKS_PER_MS)), teventq,
                 ((real_t)(tepisode/TICKS_PER_MS)), episodes);
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
                 ((real_t)(tstep/TICKS_PER_MS)), teventq, tdisplay,
                 trecord, tsave, ((real_t)(tmax/TICKS_PER_MS)));
      }
      // Set compute cycle
      netcycle = CkCallback(CkIndex_Network::CycleSimGPU(), network);
      network.InitSimGPU(netdata);
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
      CkPrintf("  Random Number Seed  (randseed): %" PRIidx "\n"
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
               randseed, groupdir.c_str(), ((real_t)(tstep/TICKS_PER_MS)), teventq,
               grpactivestring.c_str(), grpmotherstring.c_str(), grpanchorstring.c_str(),
               grpminlen, ((real_t)(grpmaxdur/TICKS_PER_MS)),
               ((real_t)grpvtxminreal), grpvtxmin, ((real_t)grpvtxmaxreal), grpvtxmax);
      // Set compute cycle
      netcycle = CkCallback(CkIndex_Network::CycleGroup(), network);
      network.InitGroup(netdata);
    }
    else if (runmode == RUNMODE_ESTIMATE) {
      if (episodic) {
        CkPrintf("  Random Number Seed  (randseed): %" PRIidx "\n"
                 "  Group Directory     (groupdir): %s\n"
                 "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
                 "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
                 "  Time per Episode    (tepisode): %" PRIrealms "ms\n"
                 "  Number of Episodes  (episodes): %" PRIidx "\n",
                 randseed, groupdir.c_str(), ((real_t)(tstep/TICKS_PER_MS)), teventq,
                 ((real_t)(tepisode/TICKS_PER_MS)), episodes);
      }
      else {
        CkPrintf("  Random Number Seed  (randseed): %" PRIidx "\n"
                 "  Group Directory     (groupdir): %s\n"
                 "  Simulation Time Step   (tstep): %" PRIrealms "ms\n"
                 "  Event Queue Length   (teventq): %" PRIrealms "ms\n"
                 "  Display Interval    (tdisplay): %" PRIrealms "ms\n"
                 "  Recording Interval   (trecord): %" PRIrealms "ms\n"
                 "  Max Simulation Time     (tmax): %" PRIrealms "ms\n",
                 randseed, groupdir.c_str(), ((real_t)(tstep/TICKS_PER_MS)), teventq,
                 tdisplay, trecord, ((real_t)(tmax/TICKS_PER_MS)));
      }
      // Set compute cycle
      netcycle = CkCallback(CkIndex_Network::CycleEst(), network);
      network.InitEst(netdata);
    }
    
#ifdef STACS_WITH_YARP
    // Open RPC port
    stream.OpenRPC(network, netcycle, rpcpause);
#endif
  }
}

// Coordination for starting simulation, initialized network partitions
//
void Main::Start() {
  // Start simulation
  CkCallback cbstop(CkReductionTarget(Main, Stop), mainProxy);
  network.ckSetReductionClient(&cbstop);

#ifdef STACS_WITH_YARP
  if (rpcpause) {
    // Start paused
    CkPrintf("Starting (paused)\n");
  }
  else {
    CkPrintf("Starting\n");
    netcycle.send();
  }
#else
  CkPrintf("Starting\n");
  netcycle.send();
#endif

  // Start timer
  wcstart = std::chrono::system_clock::now();
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
