/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"
#include "stream.h"

#ifdef STACS_WITH_YARP

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ idx_t nevtday;


/**************************************************************************
* Network Remote Procedure Call
**************************************************************************/

// Receive and handle RPC messages
//
void Network::CommRPC(mRPC *msg) {
  //CkPrintf("Received %" PRIidx" on %" PRIidx " in iteration %" PRIidx "\n", msg->command, partidx, iter);

  // Pausing/Checkpointing/Stepping
  //
  if (msg->command == RPCCOMMAND_PAUSE ||
      msg->command == RPCCOMMAND_STOP) {
    synciter = iter;
    cyclepart.send();
    // Coordinate the synchronization iteration
    /*
    if (commiter < iter) {
      synciter = iter;
      if (partiter % 2 == 0 && cadjpart[1] > 0) { ++synciter; }
      if (partiter % 2 == 1 && cadjpart[0] > 0) { ++synciter; }
      //if (partiter % 2 == 0 && cadjpart[1] > 0 && cadjpart[0] > 0) { ++synciter; }
      //if (partiter % 2 == 1 && cadjpart[0] > 0 && cadjpart[1] > 0) { ++synciter; }
    }
    else if (commiter == iter) {
      synciter = commiter+1;
      if (partiter % 2 == 0 && cadjpart[1] > cadjpart[0]) { --synciter; }
      if (partiter % 2 == 1 && cadjpart[0] > cadjpart[1]) { --synciter; }
    }
    else {
      synciter = commiter;
    }
    synciter = commiter;
    if (partiter % 2 == 0 && cadjpart[0] > cadjpart[1]) { ++synciter; }
    if (partiter % 2 == 1 && cadjpart[1] > cadjpart[0]) { ++synciter; }
    */
    // TODO make this print with debugging flag
    //CkPrintf("Messages on %" PRIidx ": c0: %d, c1: %d, pi: % " PRIidx "\n", partidx, cadjpart[0], cadjpart[1], partiter);
    //CkPrintf("Pausing %" PRIidx " in iteration %" PRIidx " (comm: %" PRIidx ", sim: %" PRIidx ")\n", partidx, synciter, commiter, iter);
  }
  else if (msg->command == RPCCOMMAND_PAUSED) {
    synciter = msg->nrpcdata;
    cyclepart.send();
  }
  else if (msg->command == RPCCOMMAND_UNPAUSE) {
    // Resume simulation
    syncing = false;
    cyclepart.send();
  }
  else if (msg->command == RPCCOMMAND_CHECK) {
    // Coordinate the synchronization iteration
    synciter = iter;
    // Perform checkpointing
    thisProxy(partidx).SaveNetwork();
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
    cyclepart.send();
  }

  // Stimulation
  //
  else if (msg->command == RPCCOMMAND_STIM || msg->command == RPCCOMMAND_PSTIM) {
    // Coordinate the syncronization iteration
    synciter = iter;
    if (msg->command == RPCCOMMAND_STIM) {
      if (partiter % 2 == 0 && cadjpart[0] > cadjpart[1]) { ++synciter; }
      if (partiter % 2 == 1 && cadjpart[1] > cadjpart[0]) { ++synciter; }
    }

    // Apply Stimuli (with offset)
    if (msg->nrpcdata > 0) {
      idx_t stimtype = (idx_t) msg->rpcdata[0];
      // per neuron
      if (stimtype == RPCSTIM_POINT) {
        idx_t numvtx = (idx_t) msg->rpcdata[1];
        idx_t pulses = (idx_t) msg->rpcdata[2];
        // build stim event
        evtrpc.resize(pulses*2);
        for (idx_t i = 0; i < pulses; ++i) {
          evtrpc[i*2  ].diffuse = ((tick_t) synciter*tstep) + ((tick_t) msg->rpcdata[3+numvtx+i*3]*TICKS_PER_MS);
          evtrpc[i*2+1].diffuse = evtrpc[i*2].diffuse + ((tick_t) msg->rpcdata[3+numvtx+i*3+1]*TICKS_PER_MS);
          evtrpc[i*2  ].type = EVENT_STIM;
          evtrpc[i*2+1].type = EVENT_STIM;
          evtrpc[i*2  ].source = -1;
          evtrpc[i*2+1].source = -1;
          evtrpc[i*2  ].index = 0;
          evtrpc[i*2+1].index = 0;
          evtrpc[i*2  ].data = msg->rpcdata[3+numvtx+i*3+2];
          evtrpc[i*2+1].data = -msg->rpcdata[3+numvtx+i*3+2];
        }
        // add events to vertices
        for (idx_t i = 3; i < 3 + numvtx; ++i) {
          std::unordered_map<idx_t, idx_t>::iterator target = vtxmap.find((idx_t) msg->rpcdata[i]);
          if (target != vtxmap.end()) {
            for (size_t e = 0; e < evtrpc.size(); ++e) {
              evtrpc[e].source -= (idx_t) msg->rpcdata[i];
              if ((evtrpc[e].diffuse/tstep - synciter) < nevtday) {
                evtcal[target->second][(evtrpc[e].diffuse/tstep)%nevtday].push_back(evtrpc[e]);
              }
              else {
                evtcol[target->second].push_back(evtrpc[e]);
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
      //thisProxy(partidx).CycleNetwork();
      cyclepart.send();
    }
  }

  // YARP Port Streams
  //
  else if (msg->command == RPCCOMMAND_OPEN) {
    // Coordinate the syncronization iteration
    synciter = iter;
    // Resync
    //thisProxy(partidx).CycleNetwork();
    cyclepart.send();
  }
  else if (msg->command == RPCCOMMAND_CLOSE) {
    // Coordinate the syncronization iteration
    synciter = iter;
    // Resync
    //thisProxy(partidx).CycleNetwork();
    cyclepart.send();
  }

  // cleanup
  delete msg;
}

#endif //STACS_WITH_YARP
