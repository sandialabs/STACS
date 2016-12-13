/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 *
 */

#include "stacs.h"
#include "network.h"
//#include "stream.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ idx_t npnet;


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
  msgSize[0] = models.size();     //modtype
  msgSize[1] = models.size();     //state
  msgSize[2] = models.size();     //stick
  msgSize[3] = models.size()+1;   //xparam
  msgSize[4] = nparam;            //param
  msgSize[5] = models.size()+1;   //xport
  msgSize[6] = nport;             //port
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

  }

  // Return model
  return mmodel;
}

#ifdef STACS_WITH_YARP
// Build graph adjacency information (just vertices)
//
mVtxDist* Main::BuildVtxDist() {
  // Initialize distribution message
  int msgSize[MSG_VtxDist];
  msgSize[0] = npnet+1; // vtxdist
  mVtxDist *mvtxdist = new(msgSize, 0) mVtxDist;

  // Get distribution info
  for (idx_t i = 0; i < npnet+1; ++i) {
    //vtxdist
    mvtxdist->vtxdist[i] = netdist[i].nvtx;
  }

  // Return distribution message
  return mvtxdist;
}
#endif
