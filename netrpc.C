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
extern /*readonly*/ idx_t evtcal;


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
  else if (msg->command == RPCCOMMAND_STIM || msg->command == RPCCOMMAND_PSTIM) {
    // Coordinate the syncronization iteration
    synciter = iter;
    if (msg->command == RPCCOMMAND_STIM) {
      if (prtiter % 2 == 0 && cadjprt[0] > cadjprt[1]) { ++synciter; }
      if (prtiter % 2 == 1 && cadjprt[1] > cadjprt[0]) { ++synciter; }
    }

    // Apply Stimuli (with offset)
    if (msg->nrpcdata > 0) {
      idx_t stimtype = (idx_t) msg->rpcdata[0];
      // per neuron
      if (stimtype == RPCSTIM_POINT) {
        idx_t numvtx = (idx_t) msg->rpcdata[1];
        idx_t pulses = (idx_t) msg->rpcdata[2];
        // build stim event
        std::vector<event_t> evtpre;
        evtpre.resize(pulses*2);
        for (idx_t i = 0; i < pulses; ++i) {
          evtpre[i*2  ].diffuse = ((tick_t) synciter*tstep) + ((tick_t) msg->rpcdata[3+numvtx+i*3]*TICKS_PER_MS);
          evtpre[i*2+1].diffuse = evtpre[i*2].diffuse + ((tick_t) msg->rpcdata[3+numvtx+i*3+1]*TICKS_PER_MS);
          evtpre[i*2  ].index = 0;
          evtpre[i*2+1].index = 0;
          evtpre[i*2  ].type = EVTYPE_STIM;
          evtpre[i*2+1].type = EVTYPE_STIM;
          evtpre[i*2  ].data = msg->rpcdata[3+numvtx+i*3+2];
          evtpre[i*2+1].data = -msg->rpcdata[3+numvtx+i*3+2];
        }
        // add events to vertices
        for (idx_t i = 3; i < 3 + numvtx; ++i) {
          std::unordered_map<idx_t, idx_t>::iterator target = vtxmap.find((idx_t) msg->rpcdata[i]);
          if (target != vtxmap.end()) {
            for (size_t e = 0; e < evtpre.size(); ++e) {
              if ((evtpre[e].diffuse/tstep - synciter) < evtcal) {
                event[target->second][(evtpre[e].diffuse/tstep)%evtcal].push_back(evtpre[e]);
              }
              else {
                evtaux[target->second].push_back(evtpre[e]);
              }
            }
          }
        }
      }
      // circle
      else if (stimtype == RPCSTIM_CIRCLE) {
      }
      else if (stimtype == RPCSTIM_SPHERE) {
      }
    }
    // Resync if necessary
    if (msg->command == RPCCOMMAND_PSTIM) {
      thisProxy(prtidx).Cycle();
    }
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
    // Decompose stimuli into events
    if (message.size() > 1) {
      idx_t stimtype = message.get(0).asInt();
      // stim file
      if (stimtype == RPCSTIM_FILE) {
        std::string filename(message.get(1).asString().c_str());
        CkPrintf("  Opening stim file: %s\n", filename.c_str());
      }
      // per neuron
      else if (stimtype == RPCSTIM_POINT) {
        if (message.size() > 2) {
          idx_t numvtx = message.get(1).asInt();
          idx_t pulses = message.get(2).asInt();
          if (message.size() == 3 + numvtx + pulses*3) {
            rpcdata.push_back((real_t) stimtype);
            rpcdata.push_back((real_t) numvtx);
            rpcdata.push_back((real_t) pulses);
            // neurons
            for (idx_t i = 3; i < 3 + numvtx; ++i) {
              rpcdata.push_back((real_t) message.get(i).asInt());
            }
            // pulses (offset, duration, amplitude)
            for (idx_t i = 3 + numvtx; i < 3 + numvtx + pulses*3; ++i) {
              rpcdata.push_back((real_t) message.get(i).asDouble());
            }
            CkPrintf("  Stimulating %" PRIidx " neurons with %" PRIidx " pulses\n", numvtx, pulses);
          }
        }
      }
      // circle
      else if (stimtype == RPCSTIM_CIRCLE) {
      }
      // sphere
      else if (stimtype == RPCSTIM_SPHERE) {
      }
    }
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
  CkPrintf("  Simulation Synced\n");

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
  CkPrintf("  Simulation Paused\n");

  // Set Synchronization flag
  rpcreader->SetSyncFlag(RPCSYNC_SYNCED);

  // Reset callback
  network.ckSetReductionClient(cbmain);
}

#endif //STACS_WITH_YARP
