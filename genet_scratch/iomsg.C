/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 *
 */

#include "genet.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ int netfiles;


/**************************************************************************
* Build Messages (Main)
**************************************************************************/

// Build models for network
//
mModel* Main::BuildModel() {
  /* Bookkeeping */
  idx_t nmodname;
  idx_t nstatetype;
  idx_t nstateparam;
  idx_t jstateparam;
  idx_t nsticktype;
  idx_t nstickparam;
  idx_t jstickparam;
  idx_t ndatafile;

  // Get total size of param
  nmodname = 0;
  nstatetype = 0;
  nstateparam = 0;
  nsticktype = 0;
  nstickparam = 0;
  for (std::size_t i = 0; i < models.size(); ++i) {
    nmodname += models[i].modname.size();
    nstatetype += models[i].statetype.size();
    nsticktype += models[i].sticktype.size();
    for (std::size_t j = 0; j < models[i].statetype.size(); ++j) {
      nstateparam += models[i].stateparam[j].size();
    }
    for (std::size_t j = 0; j < models[i].sticktype.size(); ++j) {
      nstickparam += models[i].stickparam[j].size();
    }
  }
  ndatafile = 0;
  for (std::size_t i = 0; i < datafiles.size(); ++i) {
    ndatafile += datafiles[i].size();
  }

  // Initialize model message
  int msgSize[MSG_Model];
  msgSize[0] = models.size();     // type
  msgSize[1] = models.size()+1;   // xmodname
  msgSize[2] = nmodname;          // modname
  msgSize[3] = models.size()+1;   // xstatetype
  msgSize[4] = models.size()+1;   // xsticktype
  msgSize[5] = nstatetype;        // statetype
  msgSize[6] = nsticktype;        // sticktype
  msgSize[7] = nstateparam;       // stateparam
  msgSize[8] = nstickparam;       // stickparam
  msgSize[9] = datafiles.size()+1;  // xdatafiles
  msgSize[10] = ndatafile;          // datafiles
  mModel *mmodel = new(msgSize, 0) mModel;
  // Sizes
  mmodel->nmodel = models.size();
  mmodel->nstateparam = nstateparam;
  mmodel->nstickparam = nstickparam;
  mmodel->ndatafiles = datafiles.size();

  // Prefixes starts with zero
  mmodel->xmodname[0] = 0;
  mmodel->xstatetype[0] = 0;
  mmodel->xsticktype[0] = 0;
  mmodel->xdatafiles[0] = 0;

  // Set up counters
  jstateparam = 0;
  jstickparam = 0;

  // Copy over model information
  for (std::size_t i = 0; i < models.size(); ++i) {
    // type
    mmodel->type[i] = models[i].type;
    // xmodname
    mmodel->xmodname[i+1] = mmodel->xmodname[i] + models[i].modname.size();
    for (std::size_t j = 0; j < models[i].modname.size(); ++j) {
      // modname
      mmodel->modname[mmodel->xmodname[i] + j] = models[i].modname[j];
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
    // xsticktype
    mmodel->xsticktype[i+1] = mmodel->xsticktype[i] + models[i].sticktype.size();
    for (std::size_t j = 0; j < models[i].sticktype.size(); ++j) {
      // sticktype
      mmodel->sticktype[mmodel->xsticktype[i]+j] = models[i].sticktype[j];
      for (std::size_t s = 0; s < models[i].stickparam[j].size(); ++s) {
        mmodel->stickparam[jstickparam++] = models[i].stickparam[j][s];
      }
    }
  }
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

// Build metis information for network
//
mMetis* Main::BuildMetis() {
  // Initialize metis message
  int msgSize[MSG_Metis];
  msgSize[0] = vtxdist.size();   // vtxdist
  msgSize[1] = edgdist.size();   // edgdist
  mMetis *mmetis = new(msgSize, 0) mMetis;

  // Sanity check
  CkAssert(vtxdist.size() == netfiles+1);
  CkAssert(edgdist.size() == netfiles+1);

  // copy over metis information
  for (int i = 0; i < netfiles+1; ++i) {
    mmetis->vtxdist[i] = vtxdist[i];
    mmetis->edgdist[i] = edgdist[i];
  }

  return mmetis;
}
