/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "stacs.h"
#include "network.h"

#ifdef STACS_WITH_YARP

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ CProxy_Main mainProxy;
extern /*readonly*/ tick_t tstep;


/**************************************************************************
* Network Remote Procedure Call
**************************************************************************/

// Receive and handle RPC messages
//
void Network::RPCMsg(mRPC *msg) {
  //CkPrintf("Received %" PRIidx" on %" PRIidx " in iteration %" PRIidx "\n", msg->command, prtidx, iter);

  // Pausing/Checkpointing/Stepping
  //
  if (msg->command == RPCCOMMAND_PAUSE ||
      msg->command == RPCCOMMAND_STOP) {
    // Coordinate the synchronization iteration
    synciter = iter;
    if (prtiter % 2 == 0 && cadjprt[0] > cadjprt[1]) { ++synciter; }
    if (prtiter % 2 == 1 && cadjprt[1] > cadjprt[0]) { ++synciter; }
    // TODO make this print with debugging flag
    //CkPrintf("Messages on %" PRIidx ": c0: %d, c1: %d, pi: % " PRIidx "\n", prtidx, cadjprt[0], cadjprt[1], prtiter);
    //CkPrintf("Pausing %" PRIidx " in iteration %" PRIidx "\n", prtidx, synciter);
  }
  else if (msg->command == RPCCOMMAND_UNPAUSE) {
    // Resume simulation
    thisProxy(prtidx).Cycle();
  }
  else if (msg->command == RPCCOMMAND_CHECK) {
    // Coordinate the synchronization iteration
    synciter = iter;
    // Perform checkpointing
    cpflag = true;
    thisProxy(prtidx).SaveNetwork();
  }
  else if (msg->command == RPCCOMMAND_STEP) {
    if (msg->nrpcdata == 0) {
      // Step for one iteration
      synciter = iter + 1;
    }
    else {
      // Coordinate the synchronization iteration
      synciter = iter + (idx_t) (((tick_t) (msg->rpcdata[0]*TICKS_PER_MS))/tstep);
    }
    // Resume simulation until synchronization point
    thisProxy(prtidx).Cycle();
  }

  // Stimulation
  //
  else if (msg->command == RPCCOMMAND_STIM) {
    // Coordinate the syncronization iteration
    synciter = iter;
    if (prtiter % 2 == 0 && cadjprt[0] > cadjprt[1]) { ++synciter; }
    if (prtiter % 2 == 1 && cadjprt[1] > cadjprt[0]) { ++synciter; }

    // Apply Stimuli (with offset)

  }
  else if (msg->command == RPCCOMMAND_PSTIM) {
    // Coordinate the syncronization iteration
    synciter = iter;
    // Apply Stimuli
    //msg->nrpcdata
    // Resync
    thisProxy(prtidx).Cycle();
  }

  // YARP Port Streams
  //
  else if (msg->command == RPCCOMMAND_OPEN) {
    // Coordinate the syncronization iteration
    synciter = iter;
    // Resync
    thisProxy(prtidx).Cycle();
  }
  else if (msg->command == RPCCOMMAND_CLOSE) {
    // Coordinate the syncronization iteration
    synciter = iter;
    // Resync
    thisProxy(prtidx).Cycle();
  }

  // cleanup
  delete msg;
}



/**************************************************************************
* RPC Messages
**************************************************************************/

// Build RPC message (single)
//
mRPC* RPCReader::BuildRPCMsg(idx_t command, yarp::os::Bottle message) {
  /* Message data */
  std::vector<real_t> rpcdata;
  
  // Prepare data
  rpcdata.clear();

  // Don't send message if nothing to send
  if (command == RPCCOMMAND_NONE) {
    return NULL;
  }

  // Send (optional) step size
  if (command == RPCCOMMAND_STEP) {
    real_t step = 0.0;
    if (message.size() > 0) {
      step = (real_t) ((((tick_t) (message.get(0).asDouble())*TICKS_PER_MS))/tstep);
    }
    if (step > 0.0) {
      rpcdata.push_back(step);
    }
  }

  // Send stimuli information
  if (command == RPCCOMMAND_STIM ||
      command == RPCCOMMAND_PSTIM) {
  }

  // Initialize rpc message
  int msgSize[MSG_RPC];
  msgSize[0] = rpcdata.size();    //rpcdata
  mRPC *rpcmsg = new(msgSize, 0) mRPC;
  // Sizes
  rpcmsg->nrpcdata = rpcdata.size();
  rpcmsg->command = command;

  for (std::size_t i = 0; i < rpcdata.size(); ++i) {
    rpcmsg->rpcdata[i] = rpcdata[i];
  }

  // Return built message
  return rpcmsg;
}

// Network callback (continue)
//
void StreamRPC::RPCSync() {
  // Display some information
  CkPrintf("Simulation Synced\n");

  // Set Synchronization flag
  rpcreader->SetSyncFlag(RPCSYNC_UNSYNCED);

  // Reset callback
  network.ckSetReductionClient(cbmain);

  // Restart network
  network.Cycle();
}

// Network callback (paused)
//
void StreamRPC::RPCPause() {
  // Display some information
  CkPrintf("Simulation Paused\n");

  // Set Synchronization flag
  rpcreader->SetSyncFlag(RPCSYNC_SYNCED);

  // Reset callback
  network.ckSetReductionClient(cbmain);
}

#endif //STACS_WITH_YARP
