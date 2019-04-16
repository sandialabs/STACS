/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifdef STACS_WITH_YARP

#include "stacs.h"
#include "stream.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ CProxy_Main mainProxy;
extern /*readonly*/ int netparts;
extern /*readonly*/ tick_t tstep;


/**************************************************************************
* Stream (Remote Procedure Call)
**************************************************************************/

// Stream constructor
//
Stream::Stream(mVtxs *msg) {
  // Read in distribution
  vtxdist.resize(netparts+1);
  for (idx_t i = 0; i < netparts+1; ++i) {
    // vtxdist
    vtxdist[i] = msg->vtxdist[i];
  }
  // Get RPC port name
  rpcport = std::string(msg->rpcport, msg->rpcport + msg->xrpcport);

  delete msg;

  // Return control to main
  mainProxy.Init();
}

// Stream migration
//
Stream::Stream(CkMigrateMessage *msg) {
  delete msg;
}

// Stream destructor
Stream::~Stream() {
}

// Open RPC port
//
void Stream::OpenRPC(CProxy_Network cpnet, const CkCallback &cbcycle, bool paused) {
  // Charm++ Coordination
  network = cpnet;
  netcycle = cbcycle;

  // Open RPC and attach callback
  CkPrintf("Opening port %s\n", rpcport.c_str());
  this->open(rpcport.c_str());
  // Test port
  if (this->where().getPort() > 0) {
    // Setup callback
    netpause = CkCallback(CkReductionTarget(Stream, Pause), thisProxy);
    netstop = CkCallback(CkReductionTarget(Main, Stop), mainProxy);
    // Set RPC reader
    rpcreader = new RPCReader(vtxdist, network, netcycle, netstop,
        CkCallback(CkReductionTarget(Stream, Pause), thisProxy),
        CkCallback(CkReductionTarget(Stream, Sync), thisProxy));
    if (rpcreader != NULL) {
      this->setReader(*rpcreader);
      if (paused) {
        rpcreader->SetSyncFlag(RPCSYNC_SYNCED);
      }
      else {
        rpcreader->SetSyncFlag(RPCSYNC_UNSYNCED);
      }
    }
    else {
      CkPrintf("rpcreader failed to initialize...\n");
    }
  }
  else {
    CkPrintf("%s failed to open...\n", rpcport.c_str());
  }
}

// Close RPC port
//
void Stream::CloseRPC() {
  // Close RPC
  if (this->where().getPort() > 0) {
    CkPrintf("Closing port %s\n", rpcport.c_str());
    this->close();
    if (rpcreader != NULL) {
      delete rpcreader;
    }
  }
}


/**************************************************************************
* Stream synchronization callbacks
**************************************************************************/


// Network callback (syncing)
//
void Stream::Sync(idx_t synciter) {
  //idx_t synciter = *(idx_t *)msg;
  //delete msg;
  // Display some information
  CkPrintf("Syncing at %" PRIidx" \n", synciter+1);

  // Set next callback
  network.ckSetReductionClient(&netpause);
  
  // Broadcast command to network
  mRPC *mrpc = BuildRPCSync(synciter+1);
  if (mrpc != NULL) {
    network.CommRPC(mrpc);
  }
}


// Network callback (paused/synced)
//
void Stream::Pause() {
  // Display some information
  CkPrintf("Paused\n");

  // Set Synchronization flag
  rpcreader->SetSyncFlag(RPCSYNC_SYNCED);

  // Reset callback
  network.ckSetReductionClient(&netstop);
}


/**************************************************************************
* RPC Reader
**************************************************************************/

// Read input message
//
bool RPCReader::read(yarp::os::ConnectionReader& connection) {
  /* Bottles */
  yarp::os::Bottle in, out;

  // Test if we read something
  if (!(in.read(connection))) {
    return false;
  }

  // Extract message
  std::string command(in.get(0).asString().c_str());
  int cmdid = RPCCOMMAND_NONE;
  // Prepare output
  out.clear();

  // Process message, with callbacks
  // TODO: return a clean help page
  // TODO: p - pause/play
  if (command == "help") {
    CkPrintf("RPC Command: Help\n");
    out.addVocab(yarp::os::Vocab::encode("many"));
    out.add("received command: help"
            "  commands to control simulation are:\n"
            "   - help: show this help page\n"
            "   - pause: un/pause the simulation\n"
            "   - stop: stop the simulation\n"
            "   - save: save simulation state\n"
            "   - step <t>: step the simulation <t> ms\n"
            "   - stim <t o a d>: apply stimulation\n");
  }

  // Pausing
  // (also syncronizes the simulation to the same iteration)
  //
  else if (command == "pause") {
    // Check for synchronization (toggles behavior)
    if (syncflag == RPCSYNC_SYNCING) {
      CkPrintf("RPC Error: Simulation is Syncing\n");
    }
    else if (syncflag == RPCSYNC_UNSYNCED) {
      syncflag = RPCSYNC_SYNCING;

      // Pause
      cmdid = RPCCOMMAND_PAUSE;
      CkPrintf("RPC Command: Pausing Simulation\n");

      // Modify reduction client
      network.ckSetReductionClient(&netsync);
    }
    else if (syncflag == RPCSYNC_SYNCED) {
      syncflag = RPCSYNC_UNSYNCED;

      // Unpause
      cmdid = RPCCOMMAND_UNPAUSE;
      CkPrintf("RPC Command: Unpausing Simulation\n");
    }
    out.add("received command: pause");
  }

  // Stop Simulation
  //
  else if (command == "stop") {
    // Check for synchronization
    if (syncflag == RPCSYNC_SYNCING) {
      CkPrintf("RPC Error: Simulation is Syncing\n");
    }
    else if (syncflag == RPCSYNC_UNSYNCED) {
      CkPrintf("RPC Error: Simulation not Paused\n");
    }
    else if (syncflag == RPCSYNC_SYNCED) {
      // Stop while paused
      CkPrintf("RPC Command: Stopping Simulation\n");
      netstop.send();
    }
    out.add("received command: stop");
  }

  // Saving network state
  //
  else if (command == "save") {
    // Check for synchronization
    if (syncflag == RPCSYNC_SYNCING) {
      CkPrintf("RPC Error: Simulation is Syncing\n");
    }
    else if (syncflag == RPCSYNC_UNSYNCED) {
      CkPrintf("RPC Error: Simulation not Paused\n");
    }
    else if (syncflag == RPCSYNC_SYNCED) {
      syncflag = RPCSYNC_SYNCING;

      // Checkpoint while paused
      cmdid = RPCCOMMAND_CHECK;
      CkPrintf("RPC Command: Checkpointing Simulation\n");
      
      // Modify reduction client
      network.ckSetReductionClient(&netpause);
    }
    out.add("received command: check");
  }
  
  // Stepping
  //
  else if (command == "step") {
    // Check for synchronization
    if (syncflag == RPCSYNC_SYNCING) {
      CkPrintf("RPC Error: Simulation is Syncing\n");
    }
    else if (syncflag == RPCSYNC_UNSYNCED) {
      CkPrintf("RPC Error: Simulation not Paused\n");
    }
    else if (syncflag == RPCSYNC_SYNCED) {
      syncflag = RPCSYNC_SYNCING;

      // Checkpoint while paused
      cmdid = RPCCOMMAND_STEP;
      
      real_t step = 0.0;
      if (in.size() > 1) {
        step = (real_t) ((((tick_t) (in.get(1).asDouble())*TICKS_PER_MS))/tstep);
      }
      if (step <= 0.0) {
        step = (real_t) (tstep/TICKS_PER_MS);
      }
      // Get step
      CkPrintf("RPC Command: Stepping Simulation for %" PRIrealms " ms\n", step);
      
      // Modify reduction client
      network.ckSetReductionClient(&netpause);
    }
    out.add("received command: step");
  }
  // Stimulation (
  //
  else if (command == "stim") {
    // The variables of <t o a d> are
    // target (vertex), offset (in time from current iteration)
    // amplitude (of stimulation), and duration
    // Check for synchronization
    if (syncflag == RPCSYNC_SYNCING) {
      CkPrintf("RPC Error: Simulation is Syncing\n");
    }
    else if (syncflag == RPCSYNC_UNSYNCED) {
      CkPrintf("RPC Error: Simulation not Paused\n");
    }
    else if (syncflag == RPCSYNC_SYNCED) {
      syncflag = RPCSYNC_SYNCING;

      // Stimulation while paused
      cmdid = RPCCOMMAND_PSTIM;
      CkPrintf("RPC Command: Stimulating Simulation\n");
      
      // Modify reduction client
      network.ckSetReductionClient(&netpause);
    }
    out.add("received command: stim");
  }
  // Open YARP stream
  //
  else if (command == "open") {
    // Check for synchronization
    if (syncflag == RPCSYNC_SYNCING) {
      CkPrintf("RPC Error: Simulation is Syncing\n");
    }
    else if (syncflag == RPCSYNC_UNSYNCED) {
      CkPrintf("RPC Error: Simulation not Paused\n");
    }
    else if (syncflag == RPCSYNC_SYNCED) {
      syncflag = RPCSYNC_SYNCING;

      // Stimulation while paused
      cmdid = RPCCOMMAND_OPEN;
      CkPrintf("RPC Command: Opening YARP Connection\n");
      
      // Modify reduction client
      network.ckSetReductionClient(&netpause);
    }
    out.add("received command: open");
  }
  // Close YARP stream
  //
  else if (command == "close") {
    // Check for synchronization
    if (syncflag == RPCSYNC_SYNCING) {
      CkPrintf("RPC Error: Simulation is Syncing\n");
    }
    else if (syncflag == RPCSYNC_UNSYNCED) {
      CkPrintf("RPC Error: Simulation not Paused\n");
    }
    else if (syncflag == RPCSYNC_SYNCED) {
      syncflag = RPCSYNC_SYNCING;

      // Stimulation while paused
      cmdid = RPCCOMMAND_CLOSE;
      CkPrintf("RPC Command: Closing YARP Connection\n");
      
      // Modify reduction client
      network.ckSetReductionClient(&netpause);
    }
    out.add("received command: close");
  }
  else {
    CkPrintf("RPC Message: %s\n", in.toString().c_str());
    out.add("message recieved");
  }

  // Broadcast command to network
  mRPC *mrpc = BuildRPC(cmdid, in.tail());
  if (mrpc != NULL) {
    network.CommRPC(mrpc);
  }

  // Reply with output
  yarp::os::ConnectionWriter *reply = connection.getWriter();
  if (reply != NULL) {
    out.write(*reply);
  }
  return true;
}

#endif //STACS_WITH_YARP
