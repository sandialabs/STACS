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
extern /*readonly*/ int netparts;
extern /*readonly*/ tick_t tstep;


/**************************************************************************
* Build Messages (Main)
**************************************************************************/

// Build distribution for netdata
//
mDist* Main::BuildDist() {
  // Initialize distribution message
  int msgSize[MSG_Dist];
  msgSize[0] = netparts+1;      // vtxdist
  msgSize[1] = netparts+1;      // edgdist
  msgSize[2] = netparts+1;      // statedist
  msgSize[3] = netparts+1;      // stickdist
  msgSize[4] = netparts+1;      // eventdist
  mDist *mdist = new(msgSize, 0) mDist;

  // Get distribution info
  for (int i = 0; i < netparts+1; ++i) {
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

  // Return distribution message
  return mdist;
}

// Build models for network
//
mModel* Main::BuildModel() {
  /* Bookkeeping */
  // TODO: probably better to make these ints instead of idx_t
  idx_t nmodname;
  idx_t nstatename;
  idx_t nstatenamestrings;
  idx_t nstatetype;
  idx_t nstateparam;
  idx_t jstateparam;
  idx_t nstickname;
  idx_t nsticknamestrings;
  idx_t nsticktype;
  idx_t nstickparam;
  idx_t jstickparam;
  idx_t jstatename;
  idx_t jstickname;
  idx_t jparamname;
  idx_t nparamname;
  idx_t nparamnamestrings;
  idx_t nparam;
  idx_t nport;
  idx_t ndatafile;

  // Get total name sizes
  nmodname = 0;
  nstatename = 0;
  nstatenamestrings = 0;
  nstatetype = 0;
  nstateparam = 0;
  nstickname = 0;
  nsticknamestrings = 0;
  nsticktype = 0;
  nstickparam = 0;
  for (std::size_t i = 0; i < models.size(); ++i) {
    nmodname += models[i].modname.size();
    // TODO: consolidate these counters (many are duplicated)
    nstatename += models[i].statename.size();
    nstickname += models[i].stickname.size();
    nstatetype += models[i].statetype.size();
    nsticktype += models[i].sticktype.size();
    for (std::size_t j = 0; j < models[i].statetype.size(); ++j) {
      // strings have an extra character for the null termination
      nstatenamestrings += models[i].statename[j].size() + 1;
      nstateparam += models[i].stateparam[j].size();
    }
    for (std::size_t j = 0; j < models[i].sticktype.size(); ++j) {
      nsticknamestrings += models[i].stickname[j].size() + 1;
      nstickparam += models[i].stickparam[j].size();
    }
  }
  // Get total size of param
  nparamname = 0;
  nparamnamestrings = 0;
  nparam = 0;
  nport = 0;
  for (std::size_t i = 0; i < models.size(); ++i) {
    nparam += models[i].param.size();
    nparamname += models[i].paramname.size();
    for (std::size_t j = 0; j < models[i].paramname.size(); ++j) {
      nparamnamestrings += models[i].paramname[j].size() + 1;
    }
    for (std::size_t j = 0; j < models[i].port.size(); ++j) {
      nport += models[i].port[j].size() + 1;
    }
  }
  // Get size of datafiles
  ndatafile = 0;
  for (std::size_t i = 0; i < datafiles.size(); ++i) {
    ndatafile += datafiles[i].size();
  }

  // Initialize model message
  int msgSize[MSG_Model];
  msgSize[0] = models.size();     // modtype
  msgSize[1] = models.size();     // graphtype
  msgSize[2] = models.size()+1;   // xmodname
  msgSize[3] = nmodname;          // modname
  msgSize[4] = models.size();     // state
  msgSize[5] = models.size();     // stick
  msgSize[6] = models.size();     // param
  msgSize[7] = nstatename+1;        // xstatename
  msgSize[8] = nstickname+1;        // xstickname
  msgSize[9] = nstatenamestrings;   // statename
  msgSize[10] = nsticknamestrings;  // stickname
  msgSize[11] = models.size()+1;   // xstatetype
  msgSize[12] = models.size()+1;   // xsticktype
  msgSize[13] = nstatetype;        // statetype
  msgSize[14] = nsticktype;        // sticktype
  msgSize[15] = nstateparam;       // stateparam
  msgSize[16] = nstickparam;       // stickparam
  msgSize[17] = nparamname+1;      // xparamname
  msgSize[18] = nparamnamestrings; // paramname
  msgSize[19] = models.size()+1;   // xparam
  msgSize[20] = nparam;            // param
  msgSize[21] = models.size()+1;   // xport
  msgSize[22] = nport;             // port
  msgSize[23] = datafiles.size()+1;  // xdatafiles
  msgSize[24] = ndatafile;           // datafiles
  msgSize[25] = ndatafile;           // datatypes
  msgSize[26] = models.size();     // grpactive
  msgSize[27] = models.size();     // grpmother
  msgSize[28] = models.size();     // grpanchor
  mModel *mmodel = new(msgSize, 0) mModel;
  // Sizes
  mmodel->nmodel = models.size();
  mmodel->nstateparam = nstateparam;
  mmodel->nstickparam = nstickparam;
  mmodel->ndatafiles = datafiles.size();
  // Configuration
  mmodel->plastic = plastic;
  mmodel->episodic = episodic;

  // Prefixes starts with zero
  mmodel->xmodname[0] = 0;
  mmodel->xstatename[0] = 0;
  mmodel->xstickname[0] = 0;
  mmodel->xstatetype[0] = 0;
  mmodel->xsticktype[0] = 0;
  mmodel->xparamname[0] = 0;
  mmodel->xparam[0] = 0;
  mmodel->xport[0] = 0;
  mmodel->xdatafiles[0] = 0;

  // Set up counters
  jstateparam = 0;
  jstickparam = 0;
  jstatename = 0;
  jstickname = 0;
  jparamname = 0;

  // Copy over model information
  for (std::size_t i = 0; i < models.size(); ++i) {
    // modtype
    mmodel->modtype[i] = models[i].modtype;
    // graphtype
    mmodel->graphtype[i] = models[i].graphtype;
    // xmodname
    mmodel->xmodname[i+1] = mmodel->xmodname[i] + models[i].modname.size();
    for (std::size_t j = 0; j < models[i].modname.size(); ++j) {
      mmodel->modname[mmodel->xmodname[i] + j] = models[i].modname[j];
    }
    // nstate
    mmodel->nstate[i] = models[i].nstate;
    // nstick
    mmodel->nstick[i] = models[i].nstick;
    // nparam
    mmodel->nparam[i] = models[i].nparam;
    // xstatename
    for (std::size_t j = 0; j < models[i].statename.size(); ++j) {
      mmodel->xstatename[jstatename+1] = mmodel->xstatename[jstatename];
      // statename
      for (std::size_t c = 0; c < models[i].statename[j].size(); ++c) {
        mmodel->statename[mmodel->xstatename[jstatename]+c] = models[i].statename[j][c];
      }
      mmodel->statename[mmodel->xstatename[jstatename]+models[i].statename[j].size()] = '\0';
      // xstatename update
      mmodel->xstatename[jstatename+1] += models[i].statename[j].size() + 1;
      ++jstatename;
    }
    // xstatetype
    mmodel->xstatetype[i+1] = mmodel->xstatetype[i] + models[i].statetype.size();
    for (std::size_t j = 0; j < models[i].statetype.size(); ++j) {
      // statetype
      mmodel->statetype[mmodel->xstatetype[i]+j] = models[i].statetype[j];
      for (std::size_t s = 0; s < models[i].stateparam[j].size(); ++s) {
        mmodel->stateparam[jstateparam++] = models[i].stateparam[j][s];
      }
    }
    // xstickname
    for (std::size_t j = 0; j < models[i].stickname.size(); ++j) {
      mmodel->xstickname[jstickname+1] = mmodel->xstickname[jstickname];
      // stickname
      for (std::size_t c = 0; c < models[i].stickname[j].size(); ++c) {
        mmodel->stickname[mmodel->xstickname[jstickname]+c] = models[i].stickname[j][c];
      }
      mmodel->stickname[mmodel->xstickname[jstickname]+models[i].stickname[j].size()] = '\0';
      // xstickname update
      mmodel->xstickname[jstickname+1] += models[i].stickname[j].size() + 1;
      ++jstickname;
    }
    // xsticktype
    mmodel->xsticktype[i+1] = mmodel->xsticktype[i] + models[i].sticktype.size();
    for (std::size_t j = 0; j < models[i].sticktype.size(); ++j) {
      // sticktype
      mmodel->sticktype[mmodel->xsticktype[i]+j] = models[i].sticktype[j];
      for (std::size_t s = 0; s < models[i].stickparam[j].size(); ++s) {
        mmodel->stickparam[jstickparam++] = models[i].stickparam[j][s];
      }
    }
    // xparamname
    for (std::size_t j = 0; j < models[i].paramname.size(); ++j) {
      mmodel->xparamname[jparamname+1] = mmodel->xparamname[jparamname];
      // paramname
      for (std::size_t c = 0; c < models[i].paramname[j].size(); ++c) {
        mmodel->paramname[mmodel->xparamname[jparamname]+c] = models[i].paramname[j][c];
      }
      mmodel->paramname[mmodel->xparamname[jparamname]+models[i].paramname[j].size()] = '\0';
      // xparamname update
      mmodel->xparamname[jparamname+1] += models[i].paramname[j].size() + 1;
      ++jparamname;
    }
    // xparam
    mmodel->xparam[i+1] = mmodel->xparam[i] + models[i].param.size();
    for (std::size_t j = 0; j < models[i].param.size(); ++j) {
      // param
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
    // grpactive
    mmodel->grpactive[i] = models[i].grpactive;
    // grpmother
    mmodel->grpmother[i] = models[i].grpmother;
    // grpanchor
    mmodel->grpanchor[i] = models[i].grpanchor;
  }
  CkAssert(jstatename == nstatename);
  CkAssert(jstickname == nstickname);
  CkAssert(jparamname == nparamname);
  CkAssert(jstateparam == nstateparam);
  CkAssert(jstickparam == nstickparam);

  // Data files
  for (std::size_t i = 0; i < datafiles.size(); ++i) {
    // xdatafiles
    mmodel->xdatafiles[i+1] = mmodel->xdatafiles[i] + datafiles[i].size();
    for (std::size_t j = 0; j < datafiles[i].size(); ++j) {
      // datafiles
      mmodel->datafiles[mmodel->xdatafiles[i] + j] = datafiles[i][j];
    }
    // datatypes
    mmodel->datatypes[i] = datatypes[i];
  }

  // Return model
  return mmodel;
}

// Build graph information for network
//
mGraph* Main::BuildGraph() {
  /* Bookkeeping */
  idx_t nvtxparam;
  idx_t jvtxparam;
  idx_t nedgtarget;
  idx_t jedgtarget;
  idx_t nedgconntype;
  idx_t jedgconntype;
  idx_t nedgprobparam;
  idx_t jedgprobparam;
  idx_t nedgmaskparam;
  idx_t jedgmaskparam;

  // get total size of param
  nvtxparam = 0;
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    nvtxparam += vertices[i].param.size();
  }
  nedgtarget = 0;
  nedgconntype = 0;
  nedgprobparam = 0;
  nedgmaskparam = 0;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    nedgtarget += edges[i].target.size();
    nedgconntype += edges[i].conntype.size();
    for (std::size_t j = 0; j < edges[i].conntype.size(); ++j) {
      nedgprobparam += edges[i].probparam[j].size();
      nedgmaskparam += edges[i].maskparam[j].size();
    }
  }

  // Initialize graph message
  int msgSize[MSG_Graph];
  msgSize[0] = vertices.size();   // vtxmodidx
  msgSize[1] = vertices.size();   // vtxorder
  msgSize[2] = vertices.size();   // vtxshape
  msgSize[3] = vertices.size()+1; // xvtxparam
  msgSize[4] = nvtxparam;         // vtxparam
  msgSize[5] = vertices.size()*3; // vtxcoord
  msgSize[6] = edges.size();      // edgsource
  msgSize[7] = edges.size()+1;    // xedgtarget
  msgSize[8] = nedgtarget;        // edgtarget
  msgSize[9] = edges.size();      // edgmodidx
  msgSize[10] = edges.size();     // edgcutoff
  msgSize[11] = edges.size()+1;   // xedgconntype
  msgSize[12] = nedgconntype;     // edgconntype
  msgSize[13] = nedgconntype;     // medgprobparam
  msgSize[14] = nedgprobparam;    // edgprobparam
  msgSize[15] = nedgconntype;     // medgmaskparam
  msgSize[16] = nedgmaskparam;    // edgmaskparam
  mGraph *mgraph = new(msgSize, 0) mGraph;
  // Sizes
  mgraph->nvtx = vertices.size();
  mgraph->nvtxparam = nvtxparam;
  mgraph->nedg = edges.size();
  mgraph->nedgtarget = nedgtarget;
  mgraph->nedgconntype = nedgconntype;
  mgraph->nedgprobparam = nedgprobparam;
  mgraph->nedgmaskparam = nedgmaskparam;

  // prefixes start at zero
  mgraph->xvtxparam[0] = 0;
  
  // set up counters
  jvtxparam = 0;

  // Streams and Vertices
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    mgraph->vtxmodidx[i] = vertices[i].modidx;
    mgraph->vtxorder[i] = vertices[i].order;
    mgraph->vtxshape[i] = vertices[i].shape;
    mgraph->xvtxparam[i+1] = mgraph->xvtxparam[i] + vertices[i].param.size();
    for (std::size_t j = 0; j < vertices[i].param.size(); ++j) {
      mgraph->vtxparam[jvtxparam++] = vertices[i].param[j];
    }
    mgraph->vtxcoord[i*3+0] = vertices[i].coord[0];
    mgraph->vtxcoord[i*3+1] = vertices[i].coord[1];
    mgraph->vtxcoord[i*3+2] = vertices[i].coord[2];
  }
  // sanity check
  CkAssert(jvtxparam == nvtxparam);

  // prefixes start at zero
  mgraph->xedgtarget[0] = 0;
  mgraph->xedgconntype[0] = 0;

  // set up counters
  jedgtarget = 0;
  jedgconntype = 0;
  jedgprobparam = 0;
  jedgmaskparam = 0;

  // Edges
  for (std::size_t i = 0; i < edges.size(); ++i) {
    mgraph->edgsource[i] = edges[i].source;
    mgraph->edgmodidx[i] = edges[i].modidx;
    mgraph->edgcutoff[i] = edges[i].cutoff;

    mgraph->xedgtarget[i+1] = mgraph->xedgtarget[i] + edges[i].target.size();
    for (std::size_t j = 0; j < edges[i].target.size(); ++j) {
      mgraph->edgtarget[jedgtarget++] = edges[i].target[j];
    }

    mgraph->xedgconntype[i+1] = mgraph->xedgconntype[i] + edges[i].conntype.size();
    for (std::size_t j = 0; j < edges[i].conntype.size(); ++j) {
      mgraph->edgconntype[jedgconntype] = edges[i].conntype[j];
      mgraph->medgprobparam[jedgconntype] = edges[i].probparam[j].size();
      mgraph->medgmaskparam[jedgconntype++] = edges[i].maskparam[j].size();
      for (std::size_t k = 0; k < edges[i].probparam[j].size(); ++k) {
        mgraph->edgprobparam[jedgprobparam++] = edges[i].probparam[j][k];
      }
      for (std::size_t k = 0; k < edges[i].maskparam[j].size(); ++k) {
        mgraph->edgmaskparam[jedgmaskparam++] = edges[i].maskparam[j][k];
      }
    }
  }
  CkAssert(jedgtarget == nedgtarget);
  CkAssert(jedgconntype == nedgconntype);
  CkAssert(jedgprobparam == nedgprobparam);
  CkAssert(jedgmaskparam == nedgmaskparam);

  // return graph
  return mgraph;
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
  msgSize[0] = netparts+1; // vtxdist
  msgSize[1] = rpcport.size(); // rpcport
  mVtxs *mvtxs = new(msgSize, 0) mVtxs;

  // Get distribution info
  for (int i = 0; i < netparts+1; ++i) {
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

// Build RPC message (sync)
//
mRPC* Stream::BuildRPCSync(idx_t synciter) {
  // Initialize rpc message
  int msgSize[MSG_RPC];
  msgSize[0] = 0;    //rpcdata
  mRPC *mrpc = new(msgSize, 0) mRPC;
  // Sizes
  mrpc->nrpcdata = synciter;
  mrpc->command = RPCCOMMAND_PAUSED;

  // Return built message
  return mrpc;
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
