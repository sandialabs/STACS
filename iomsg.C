/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 *
 */

#include "stacs.h"
#include "network.h"
#include "stream.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ idx_t npnet;
extern /*readonly*/ tick_t tstep;


/**************************************************************************
* Build Messages (Main)
**************************************************************************/

// Build distribution for netdata
//
mDist* Main::BuildDist() {
  /* Bookkeeping */
  idx_t nmodname;

  // Get total name sizes
  nmodname = 0;
  for (std::size_t i = 0; i < models.size(); ++i) {
    nmodname += models[i].modname.size();
  }

  // Initialize distribution message
  int msgSize[MSG_Dist];
  msgSize[0] = npnet+1;         // vtxdist
  msgSize[1] = npnet+1;         // edgdist
  msgSize[2] = npnet+1;         // statedist
  msgSize[3] = npnet+1;         // stickdist
  msgSize[4] = npnet+1;         // eventdist
  msgSize[5] = models.size();   // modtype
  msgSize[6] = models.size()+1; // xmodname
  msgSize[7] = nmodname;        // modname
  mDist *mdist = new(msgSize, 0) mDist;
  mdist->nmodel = models.size();

  // Get distribution info
  for (idx_t i = 0; i < npnet+1; ++i) {
    // vtxdist
    mdist->vtxdist[i] = netdist[i].nvtx;
    // edgdist
    mdist->edgdist[i] = netdist[i].nedg;
    // statedist
    mdist->statedist[i] = netdist[i].nstate;
    // stickdist
    mdist->stickdist[i] = netdist[i].nstick;
    // eventdist
    mdist->eventdist[i] = netdist[i].nevent;
  }

  // Prefixes starts with zero
  mdist->xmodname[0] = 0;

  // Get model info
  for (std::size_t i = 0; i < models.size(); ++i) {
    // modtype
    mdist->modtype[i] = models[i].modtype;
    // xmodname
    mdist->xmodname[i+1] = mdist->xmodname[i] + models[i].modname.size();
    for (std::size_t j = 0; j < models[i].modname.size(); ++j) {
      mdist->modname[mdist->xmodname[i] + j] = models[i].modname[j];
    }
  }

  // Return distribution message
  return mdist;
}

// Build models for network
//
mModel* Main::BuildModel() {
  /* Bookkeeping */
  idx_t nparam;
  idx_t nport;

  // Get total size of param
  nparam = 0;
  nport = 0;
  for (std::size_t i = 0; i < models.size(); ++i) {
    nparam += models[i].param.size();
    for (std::size_t j = 0; j < models[i].port.size(); ++j) {
      nport += models[i].port[j].size() + 1;
    }
  }

  // Initialize model message
  int msgSize[MSG_Model];
  msgSize[0] = models.size();     // modtype
  msgSize[1] = models.size();     // state
  msgSize[2] = models.size();     // stick
  msgSize[3] = models.size()+1;   // xparam
  msgSize[4] = nparam;            // param
  msgSize[5] = models.size()+1;   // xport
  msgSize[6] = nport;             // port
  msgSize[7] = models.size();     // pngactive
  msgSize[8] = models.size();     // pngmother
  msgSize[9] = models.size();     // pnganchor
  mModel *mmodel = new(msgSize, 0) mModel;
  // Sizes
  mmodel->nmodel = models.size();

  // Prefixes starts with zero
  mmodel->xparam[0] = 0;
  mmodel->xport[0] = 0;
  
  // Copy over model information
  for (std::size_t i = 0; i < models.size(); ++i) {
    // modtype
    mmodel->modtype[i] = models[i].modtype;
    // nstate
    mmodel->nstate[i] = models[i].nstate;
    // nstick
    mmodel->nstick[i] = models[i].nstick;
    // xparam
    mmodel->xparam[i+1] = mmodel->xparam[i] + models[i].param.size();
    for (std::size_t j = 0; j < models[i].param.size(); ++j) {
      // vtxparam
      mmodel->param[mmodel->xparam[i]+j] = models[i].param[j];
    }
    // xport
    mmodel->xport[i+1] = mmodel->xport[i];
    for (std::size_t j = 0; j < models[i].port.size(); ++j) {
      // port
      for (std::size_t c = 0; c < models[i].port[j].size(); ++c) {
        mmodel->port[mmodel->xport[i+1]+c] = models[i].port[j][c];
      }
      mmodel->port[mmodel->xport[i+1]+models[i].port[j].size()] = '\0';
      // xport
      mmodel->xport[i+1] += models[i].port[j].size() + 1;
    }
    // pngactive
    mmodel->pngactive[i] = models[i].pngactive;
    // pngmother
    mmodel->pngmother[i] = models[i].pngmother;
    // pnganchor
    mmodel->pnganchor[i] = models[i].pnganchor;
  }

  // Return model
  return mmodel;
}


#ifdef STACS_WITH_YARP
/**************************************************************************
* Build Messages (Stream)
**************************************************************************/

// Build graph adjacency information (just vertices)
//
mVtxs* Main::BuildVtxs() {
  // Initialize distribution message
  int msgSize[MSG_Vtxs];
  msgSize[0] = npnet+1; // vtxdist
  msgSize[1] = rpcport.size(); // rpcport
  mVtxs *mvtxs = new(msgSize, 0) mVtxs;

  // Get distribution info
  for (idx_t i = 0; i < npnet+1; ++i) {
    //vtxdist
    mvtxs->vtxdist[i] = netdist[i].nvtx;
  }
  // RPC port
  mvtxs->xrpcport = rpcport.size();
  for (std::size_t c = 0; c < rpcport.size(); ++c) {
    mvtxs->rpcport[c] = rpcport[c];
  }

  // Return distribution message
  return mvtxs;
}

// Build RPC message (single)
//
mRPC* RPCReader::BuildRPC(idx_t command, yarp::os::Bottle message) {
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
  mRPC *mrpc = new(msgSize, 0) mRPC;
  // Sizes
  mrpc->nrpcdata = rpcdata.size();
  mrpc->command = command;

  for (std::size_t i = 0; i < rpcdata.size(); ++i) {
    mrpc->rpcdata[i] = rpcdata[i];
  }

  // Return built message
  return mrpc;
}

#endif
