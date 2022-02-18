/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "genet.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ idx_t netparts;
extern /*readonly*/ int netfiles;


/**************************************************************************
* Build
**************************************************************************/

// Build Network
//
void GeNet::Build(mGraph *msg) {
  /* Bookkeeping */
  idx_t jvtxparam;
  idx_t jedgtarget;
  idx_t jedgconntype;
  idx_t jedgprobparam;
  idx_t jedgmaskparam;
  
  // initialize counters
  jvtxparam = 0;
  norder = 0;

  // Vertices
  vertices.resize(msg->nvtx);
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    vertices[i].modidx = msg->vtxmodidx[i];
    vertices[i].order = msg->vtxorder[i];
    norder += vertices[i].order;
    vertices[i].shape = msg->vtxshape[i];
    vertices[i].param.resize(msg->xvtxparam[i+1] - msg->xvtxparam[i]);
    for (std::size_t j = 0; j < vertices[i].param.size(); ++j) {
      vertices[i].param[j] = msg->vtxparam[jvtxparam++];
    }
    vertices[i].coord.resize(3);
    vertices[i].coord[0] = msg->vtxcoord[i*3+0];
    vertices[i].coord[1] = msg->vtxcoord[i*3+1];
    vertices[i].coord[2] = msg->vtxcoord[i*3+2];
  }
  // Sanity check
  CkAssert(jvtxparam == msg->nvtxparam);

  // initialize counters
  jedgtarget = 0;
  jedgconntype = 0;
  jedgprobparam = 0;
  jedgmaskparam = 0;

  // Edges
  edges.resize(msg->nedg);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    edges[i].source = msg->edgsource[i];
    edges[i].modidx = msg->edgmodidx[i];
    edges[i].cutoff = msg->edgcutoff[i];
    
    edges[i].target.resize(msg->xedgtarget[i+1] - msg->xedgtarget[i]);
    for (std::size_t j = 0; j < edges[i].target.size(); ++j) {
      edges[i].target[j] = msg->edgtarget[jedgtarget++];
    }
    
    edges[i].conntype.resize(msg->xedgconntype[i+1] - msg->xedgconntype[i]);
    edges[i].probparam.resize(msg->xedgconntype[i+1] - msg->xedgconntype[i]);
    edges[i].maskparam.resize(msg->xedgconntype[i+1] - msg->xedgconntype[i]);
    for (std::size_t j = 0; j < edges[i].conntype.size(); ++j) {
      edges[i].conntype[j] = msg->edgconntype[jedgconntype];
      edges[i].probparam[j].resize(msg->medgprobparam[jedgconntype]);
      edges[i].maskparam[j].resize(msg->medgmaskparam[jedgconntype++]);
      for (std::size_t k = 0; k < edges[i].probparam[j].size(); ++k) {
        edges[i].probparam[j][k] = msg->edgprobparam[jedgprobparam++];
      }
      for (std::size_t k = 0; k < edges[i].maskparam[j].size(); ++k) {
        edges[i].maskparam[j][k] = msg->edgmaskparam[jedgmaskparam++];
      }
    }
  }
  // Sanity checks
  CkAssert(jedgtarget == msg->nedgtarget);
  CkAssert(jedgconntype == msg->nedgconntype);
  CkAssert(jedgprobparam == msg->nedgprobparam);
  CkAssert(jedgmaskparam == msg->nedgmaskparam);

  // cleanup
  delete msg;

  // Print out some information
  if (datidx == 0) {
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      CkPrintf("Vertex: %" PRIidx "   Order: %" PRIidx"\n", vertices[i].modidx, vertices[i].order);
    }
    for (std::size_t i = 0; i < edges.size(); ++i) {
      std::string edgetargets;
      // collect edgetargets
      for (std::size_t j = 0; j < edges[i].target.size(); ++j) {
        std::ostringstream edgetarget;
        edgetarget << " " << edges[i].target[j];
        edgetargets.append(edgetarget.str());
      }
      CkPrintf("Edges:  %" PRIidx "   Source: %" PRIidx"   Targets:%s\n", edges[i].modidx, edges[i].source, edgetargets.c_str());
    }
  }

  // Bookkeeping to see how much each chare builds
  // Taking into account the different parts too
  // Initial distribution is simply contructing
  // vertices as evenly as possible across the parts
  // TODO: implement a mode to distribute vertices in chunks by population
  nordervtx.resize(nprt);
  xordervtx.resize(nprt);
  for (idx_t k = 0; k < nprt; ++k) {
    nordervtx[k].resize(vertices.size());
    xordervtx[k].resize(vertices.size());
  }
  idx_t xremvtx = 0;
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    idx_t ndivvtx = (vertices[i].order)/netparts;
    idx_t nremvtx = (vertices[i].order)%netparts;
    for (idx_t k = 0; k < nprt; ++k) {
      // looping integer magic
      nordervtx[k][i] = ndivvtx + (((xprt+k) >= xremvtx && (xprt+k) < nremvtx+xremvtx) ||
          (nremvtx+xremvtx >= netparts && (xprt+k) < xremvtx && (xprt+k) < (nremvtx+xremvtx)%netparts));
      xordervtx[k][i] = ndivvtx * (xprt+k);
      // TODO: see if you can get rid of the loop
      for (idx_t j = 0; j < xprt+k; ++j) {
        xordervtx[k][i] += ((j >= xremvtx && j < nremvtx+xremvtx) ||
            (nremvtx+xremvtx >= netparts && j < xremvtx && j < (nremvtx+xremvtx)%netparts));
      }
    }
    xremvtx = (xremvtx+nremvtx)%netparts;
  }
  
  // counting
  norderdat = 0;
  norderprt.resize(nprt);
  for (idx_t k = 0; k < nprt; ++k) {
    norderprt[k] = 0;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      norderprt[k] += nordervtx[k][i];
    }
    norderdat += norderprt[k];
  }

  // Print out vertex distribution information
  std::string orderprts;
  // collect order parts
  for (idx_t k = 0; k < nprt; ++k) {
    std::ostringstream orderprt;
    orderprt << " " << norderprt[k];
    orderprts.append(orderprt.str());
    orderprts.append(" [");
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      std::ostringstream ordervtx;
      ordervtx << " " << nordervtx[k][i] << "(" << xordervtx[k][i] << ")";
      orderprts.append(ordervtx.str());
    }
    orderprts.append(" ]");
  }
  CkPrintf("  Building File: %d   Vertices: %" PRIidx " {%s }\n", datidx, norderdat, orderprts.c_str());

  // From here, work in terms of the netfiles instead of by netparts
  // But keep in mind the per part artificial splitting

  // Create model indices
  vtxmodidx.resize(norderdat);
  vtxordidx.resize(norderdat);
  edgmodidx.resize(norderdat);
  xyz.resize(norderdat*3);
  idx_t jvtxidx = 0;
  for (idx_t k = 0; k < nprt; ++k) {
    // set with modidx
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      for (idx_t j = 0; j < nordervtx[k][i]; ++j) {
        // Set the model index
        vtxmodidx[jvtxidx] = vertices[i].modidx;
        vtxordidx[jvtxidx] = xordervtx[k][i] + j;
        edgmodidx[jvtxidx].clear();
        // Generate coordinates
        if (vertices[i].shape == VTXSHAPE_POINT) {
          // at a point
          xyz[jvtxidx*3+0] = vertices[i].coord[0];
          xyz[jvtxidx*3+1] = vertices[i].coord[1];
          xyz[jvtxidx*3+2] = vertices[i].coord[2];
        }
        else if (vertices[i].shape == VTXSHAPE_CIRCLE) {
          // uniformly inside circle
          real_t t = 2*M_PI*((*unifdist)(rngine));
          real_t r = vertices[i].param[0] * std::sqrt((*unifdist)(rngine));
          xyz[jvtxidx*3+0] = vertices[i].coord[0] + r*std::cos(t);
          xyz[jvtxidx*3+1] = vertices[i].coord[1] + r*std::sin(t);
          xyz[jvtxidx*3+2] = vertices[i].coord[2] + 0;
        }
        else if (vertices[i].shape == VTXSHAPE_SPHERE) {
          // uniformly inside sphere
          real_t u = ((*unifdist)(rngine));
          real_t x = ((*normdist)(rngine));
          real_t y = ((*normdist)(rngine));
          real_t z = ((*normdist)(rngine));
          real_t r = vertices[i].param[0] * std::cbrt(u) / std::sqrt(x*x+y*y+z*z);
          xyz[jvtxidx*3+0] = vertices[i].coord[0] + r*x;
          xyz[jvtxidx*3+1] = vertices[i].coord[1] + r*y;
          xyz[jvtxidx*3+2] = vertices[i].coord[2] + r*z;
        }
        else if (vertices[i].shape == VTXSHAPE_SPHERE) {
          // uniformly inside rectangle
        }
        // Increment for the next vertex
        ++jvtxidx;
      }
    }
  }
  CkAssert(jvtxidx == norderdat);

  // At this point vtxmodidx should have the modidx of all the vertices
  // TODO: enable edges connecting to edges at some point (by vertex id?)

  // Build vertices from model
  state.resize(norderdat);
  stick.resize(norderdat);
  event.resize(norderdat);
  for (idx_t i = 0; i < norderdat; ++i) {
    // Sanity check
    // 0 is reserved for 'none' edge type
    CkAssert(vtxmodidx[i] > 0);
    idx_t modidx = vtxmodidx[i] - 1;
    CkAssert(models[modidx].type == GRAPHTYPE_VTX || models[modidx].type == GRAPHTYPE_STR);
    // Allocate space for states
    std::vector<real_t> rngstate;
    std::vector<tick_t> rngstick;
    rngstate.resize(models[modidx].statetype.size());
    rngstick.resize(models[modidx].sticktype.size());
    // Randomly generate state
    for (std::size_t s = 0; s < models[modidx].statetype.size(); ++s) {
      if (models[modidx].statetype[s] == RNGTYPE_CONST) {
        rngstate[s] = rngconst(models[modidx].stateparam[s].data());
      }
      else if (models[modidx].statetype[s] == RNGTYPE_UNIF) {
        rngstate[s] = rngunif(models[modidx].stateparam[s].data());
      }
      else if (models[modidx].statetype[s] == RNGTYPE_UNINT) {
        rngstate[s] = rngunint(models[modidx].stateparam[s].data());
      }
      else if (models[modidx].statetype[s] == RNGTYPE_NORM) {
        rngstate[s] = rngnorm(models[modidx].stateparam[s].data());
      }
      else if (models[modidx].statetype[s] == RNGTYPE_BNORM) {
        rngstate[s] = rngbnorm(models[modidx].stateparam[s].data());
      }
      else if (models[modidx].statetype[s] == RNGTYPE_FILE) {
        rngstate[s] = rngfile(models[modidx].stateparam[s].data(), 0, vtxordidx[i]);
      }
      else {
        CkPrintf("  error: statetype %s is not valid for vertex\n", rngtype[models[modidx].statetype[s]].c_str());
        // TODO: cleaner error checking here?
        CkExit();
      }
    }
    // Randomly generate stick
    for (std::size_t s = 0; s < models[modidx].sticktype.size(); ++s) {
      if (models[modidx].sticktype[s] == RNGTYPE_CONST) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngconst(models[modidx].stickparam[s].data()));
      }
      else if (models[modidx].sticktype[s] == RNGTYPE_UNIF) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngunif(models[modidx].stickparam[s].data()));
      }
      else if (models[modidx].sticktype[s] == RNGTYPE_UNINT) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngunint(models[modidx].stickparam[s].data()));
      }
      else if (models[modidx].sticktype[s] == RNGTYPE_NORM) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngnorm(models[modidx].stickparam[s].data()));
      }
      else if (models[modidx].sticktype[s] == RNGTYPE_BNORM) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngbnorm(models[modidx].stickparam[s].data()));
      }
      else if (models[modidx].sticktype[s] == RNGTYPE_FILE) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngfile(models[modidx].stickparam[s].data(), 0, vtxordidx[i]));
      }
      else {
        CkPrintf("  error: statetype %s is not valid for vertex\n", rngtype[models[modidx].sticktype[s]].c_str());
        // TODO: cleaner error checking here?
        CkExit();
      }
    }
    // Add to state
    state[i].push_back(rngstate);
    stick[i].push_back(rngstick);
    // Empty events
    event[i].clear();
  }

  // Prepare for connection
  cpdat = 0;
  vtxdist.resize(netfiles+1);
  vtxdist[0] = 0;
  adjcy.clear();
  edgmodidx.clear();
  adjcy.resize(norderdat);
  edgmodidx.resize(norderdat);
  // Only need to worry about future edges
  adjcyconn.clear();
  adjcyconn.resize(netfiles);
  edgmodidxconn.clear();
  edgmodidxconn.resize(netfiles);

  // Request data from prev part
  if (cpdat < datidx) {
    thisProxy(cpdat).ConnRequest(datidx);
  }
  // Connect to curr part
  else if (cpdat == datidx) {
    mConn *mconn = BuildCurrConn();
    thisProxy(cpdat).Connect(mconn);
  }
  // Request data from next part
  else if (cpdat > datidx) {
    // this shouldn't happen
    thisProxy(cpdat).ConnRequest(datidx);
  }
}


/**************************************************************************
* Connect
**************************************************************************/

// Connect Network
//
void GeNet::Connect(mConn *msg) {
  // Sanity check
  CkAssert(msg->datidx == cpdat);
  // Some basic information on what's being connected
  CkPrintf("  Connecting %d to %d\n", datidx, msg->datidx);

  // Add to vtxdist
  vtxdist[cpdat+1] = vtxdist[cpdat] + msg->nvtx;

  // Perform connections
  adjcyconn[msg->datidx].resize(norderdat);
  edgmodidxconn[msg->datidx].resize(norderdat);

  // Prev
  //
  if (msg->datidx < datidx) {
    // initialize counters
    std::vector<idx_t> jadjcy(msg->nvtx, 0);
    // load prev adjacency (create connection states from previous)
    for (idx_t i = 0; i < norderdat; ++i) {
      for (idx_t j = 0; j < msg->nvtx; ++j) {
        if (jadjcy[j] == msg->xadj[j+1] - msg->xadj[j]) {
          // finished copying for this vertex
          continue;
        }
        while (jadjcy[j] < msg->xadj[j+1] - msg->xadj[j] &&
            i >= msg->adjcy[msg->xadj[j] + jadjcy[j]]) {
          if (msg->adjcy[msg->xadj[j]+jadjcy[j]] == i) {
            adjcy[i].push_back(vtxdist[msg->datidx]+j);
            idx_t modidx = msg->edgmodidx[msg->xadj[j]+jadjcy[j]];
            edgmodidx[i].push_back(modidx);
            // check if state needs to be built from j to i
            if (modidx) {
              // build state from j to i
              real_t distance = sqrt((xyz[i*3]-msg->xyz[j*3])*(xyz[i*3]-msg->xyz[j*3])+
                                (xyz[i*3+1]-msg->xyz[j*3+1])*(xyz[i*3+1]-msg->xyz[j*3+1])+
                                (xyz[i*3+2]-msg->xyz[j*3+2])*(xyz[i*3+2]-msg->xyz[j*3+2]));              
              state[i].push_back(BuildEdgState(modidx, distance, msg->vtxordidx[j], vtxordidx[i]));
              stick[i].push_back(BuildEdgStick(modidx, distance, msg->vtxordidx[j], vtxordidx[i]));
            }
            else {
              // build empty state
              state[i].push_back(std::vector<real_t>());
              stick[i].push_back(std::vector<tick_t>());
            }
          }
          ++jadjcy[j];
        }
      }
    }
  }
  
  // Curr
  //
  else if (msg->datidx == datidx) {
    CkAssert(msg->nvtx == norderdat);
    // initialize counters
    std::vector<idx_t> jadjcy(norderdat, 0);
    // perform self connection (create connection states)
    for (idx_t i = 0; i < norderdat; ++i) {
      for (idx_t j = 0; j < msg->nvtx; ++j) {
        if (i > j) {
          // copy over existing connections
          for (idx_t ji = 0; ji < norderdat; ++ji) {
            if (jadjcy[ji] == adjcyconn[datidx][ji].size()) {
              // don't copy more than once
              continue;
            }
            while (jadjcy[ji] < adjcyconn[datidx][ji].size() &&
                   i >= adjcyconn[datidx][ji][jadjcy[ji]]) {
              if (adjcyconn[datidx][ji][jadjcy[ji]] == i) {
                adjcy[i].push_back(vtxdist[datidx]+ji);
                idx_t modidx = edgmodidxconn[datidx][ji][jadjcy[ji]];
                edgmodidx[i].push_back(modidx);
                // check if state needs to be built from j to i
                if (modidx) {
                real_t distance = sqrt((xyz[i*3]-xyz[j*3])*(xyz[i*3]-xyz[j*3])+
                                  (xyz[i*3+1]-xyz[j*3+1])*(xyz[i*3+1]-xyz[j*3+1])+
                                  (xyz[i*3+2]-xyz[j*3+2])*(xyz[i*3+2]-xyz[j*3+2]));
                  // build state from j to i
                  state[i].push_back(BuildEdgState(modidx, distance, vtxordidx[j], vtxordidx[i]));
                  stick[i].push_back(BuildEdgStick(modidx, distance, vtxordidx[j], vtxordidx[i]));
                }
                else {
                  // build empty state
                  state[i].push_back(std::vector<real_t>());
                  stick[i].push_back(std::vector<tick_t>());
                }
              }
              ++jadjcy[ji];
            }
          }
        }
        else if (i < j) {
          real_t distance = sqrt((xyz[i*3]-xyz[j*3])*(xyz[i*3]-xyz[j*3])+
                            (xyz[i*3+1]-xyz[j*3+1])*(xyz[i*3+1]-xyz[j*3+1])+
                            (xyz[i*3+2]-xyz[j*3+2])*(xyz[i*3+2]-xyz[j*3+2]));
          idx_t modidx;
          // check possible connections from i to j
          if ((modidx = MakeConnection(vtxmodidx[i], vtxmodidx[j], vtxordidx[i], vtxordidx[j], distance))) {
            // add first connection
            adjcyconn[datidx][i].push_back(j);
            edgmodidxconn[datidx][i].push_back(modidx);
          }
          // check possible connections from j to i
          if ((modidx = MakeConnection(vtxmodidx[j], vtxmodidx[i], vtxordidx[j], vtxordidx[i], distance))) {
            // add second connection if it's not there
            if (adjcyconn[datidx][i].size()) {
              if (adjcyconn[datidx][i].back() != j) {
                adjcyconn[datidx][i].push_back(j);
                edgmodidxconn[datidx][i].push_back(0);
              }
            }
            else {
              adjcyconn[datidx][i].push_back(j);
              edgmodidxconn[datidx][i].push_back(0);
            }
          }
          // update adjacency with any new connections
          if (adjcyconn[datidx][i].size() && adjcyconn[datidx][i].back() == j) {
            adjcy[i].push_back(vtxdist[datidx]+j);
            edgmodidx[i].push_back(modidx);
            if (modidx) {
              // build state from j to i
              state[i].push_back(BuildEdgState(modidx, distance, vtxordidx[j], vtxordidx[i]));
              stick[i].push_back(BuildEdgStick(modidx, distance, vtxordidx[j], vtxordidx[i]));
            }
            else {
              // build empty state
              state[i].push_back(std::vector<real_t>());
              stick[i].push_back(std::vector<tick_t>());
            }
          }
        }
      }
    }
  }

  // Next
  //
  else if (msg->datidx > datidx) {
    // connect to later part (connections both ways)
    for (idx_t i = 0; i < norderdat; ++i) {
      for (idx_t j = 0; j < msg->nvtx; ++j) {
        real_t distance = sqrt((xyz[i*3]-msg->xyz[j*3])*(xyz[i*3]-msg->xyz[j*3])+
                          (xyz[i*3+1]-msg->xyz[j*3+1])*(xyz[i*3+1]-msg->xyz[j*3+1])+
                          (xyz[i*3+2]-msg->xyz[j*3+2])*(xyz[i*3+2]-msg->xyz[j*3+2]));
        idx_t modidx;
        // check possible connections from i to j
        if ((modidx = MakeConnection(vtxmodidx[i], msg->vtxmodidx[j], vtxordidx[i], msg->vtxordidx[j], distance))) {
          adjcyconn[msg->datidx][i].push_back(j);
          edgmodidxconn[msg->datidx][i].push_back(modidx);
        }
        // check possible connections from j to i
        if ((modidx = MakeConnection(msg->vtxmodidx[j], vtxmodidx[i], msg->vtxordidx[j], vtxordidx[i], distance))) {
          if (adjcyconn[msg->datidx][i].size()) {
            if (adjcyconn[msg->datidx][i].back() != j) {
              adjcyconn[msg->datidx][i].push_back(j);
              edgmodidxconn[msg->datidx][i].push_back(0);
            }
          }
          else {
            adjcyconn[msg->datidx][i].push_back(j);
            edgmodidxconn[msg->datidx][i].push_back(0);
          }
        }
        // update adjacency with any new connections
        if (adjcyconn[msg->datidx][i].size() && adjcyconn[msg->datidx][i].back() == j) {
          adjcy[i].push_back(vtxdist[msg->datidx]+j);
          edgmodidx[i].push_back(modidx);
          if (modidx) {
            // build state from j to i
            state[i].push_back(BuildEdgState(modidx, distance, msg->vtxordidx[j], vtxordidx[i]));
            stick[i].push_back(BuildEdgStick(modidx, distance, msg->vtxordidx[j], vtxordidx[i]));
          }
          else {
            // build empty state
            state[i].push_back(std::vector<real_t>());
            stick[i].push_back(std::vector<tick_t>());
          }
        }
      }
    }
  }

  // cleanup
  delete msg;

  // send any outstanding requests of built adjcy
  for (std::list<idx_t>::iterator ireqidx = adjcyreq.begin(); ireqidx != adjcyreq.end(); ++ireqidx) {
    if (*ireqidx == cpdat) {
      mConn *mconn = BuildPrevConn(*ireqidx);
      thisProxy(*ireqidx).Connect(mconn);
      ireqidx = adjcyreq.erase(ireqidx);
    }
  }

  // Move to next part
  ++cpdat;
  // return control to main when done
  if (cpdat == netfiles) {
    contribute(0, NULL, CkReduction::nop);
  }
  // Request data from next part
  else if (cpdat < datidx) {
    thisProxy(cpdat).ConnRequest(datidx);
  }
  // Connect to curr part
  else if (cpdat == datidx) {
    mConn *mconn = BuildCurrConn();
    thisProxy(cpdat).Connect(mconn);
  }
  // Request data from next part
  else if (cpdat > datidx) {
    thisProxy(cpdat).ConnRequest(datidx);
  }
}

// Connect Network Finished
//
void GeNet::ConnRequest(idx_t reqidx) {
  // send data adjacency and vertex info to requesting part
  if (reqidx > datidx) {
    // check if adjcy is built
    if (adjcyconn.size() && adjcyconn[reqidx].size()) {
      mConn *mconn = BuildPrevConn(reqidx);
      thisProxy(reqidx).Connect(mconn);
    }
    else {
      // record request for when adjcy is built
      adjcyreq.push_back(reqidx);
    }
  }
  // send data (vertex info) to requesting part
  else if (reqidx < datidx) {
    mConn *mconn = BuildNextConn();
    thisProxy(reqidx).Connect(mconn);
  }
}


/**************************************************************************
* Connection messages
**************************************************************************/

// Build Previous (includes vertices and adjacency)
//
mConn* GeNet::BuildPrevConn(idx_t reqidx) {
  /* Bookkeeping */
  idx_t nsizedat;
  idx_t jadjcyidx;

  // Sanity check
  CkAssert(adjcyconn[reqidx].size());

  // Count the sizes
  nsizedat = 0;
  for (idx_t i = 0; i < norderdat; ++i) {
    nsizedat += adjcyconn[reqidx][i].size();
  }

  // Initialize connection message
  int msgSize[MSG_Conn];
  msgSize[0] = norderdat;   // vtxmodidx
  msgSize[1] = norderdat;   // vtxordidx
  msgSize[2] = norderdat*3; // xyz
  msgSize[3] = norderdat+1; // xadj
  msgSize[4] = nsizedat;    // adjcy
  msgSize[5] = nsizedat;    // edgmodidx
  mConn *mconn = new(msgSize, 0) mConn;
  // Sizes
  mconn->datidx = datidx;
  mconn->nvtx = norderdat;

  // Prefix zero is zero
  mconn->xadj[0] = 0;
  // Initialize counter
  jadjcyidx = 0;
  
  // Build message
  for (idx_t i = 0; i < norderdat; ++i) {
    mconn->vtxmodidx[i] = vtxmodidx[i];
    mconn->vtxordidx[i] = vtxordidx[i];
    // xyz
    mconn->xyz[i*3+0] = xyz[i*3+0];
    mconn->xyz[i*3+1] = xyz[i*3+1];
    mconn->xyz[i*3+2] = xyz[i*3+2];
    // xadj
    mconn->xadj[i+1] = mconn->xadj[i] + adjcyconn[reqidx][i].size();
    for (std::size_t j = 0; j < adjcyconn[reqidx][i].size(); ++j) {
      // adjcy (of next parts)
      mconn->adjcy[jadjcyidx] = adjcyconn[reqidx][i][j];
      // edgmodidx (of next parts)
      mconn->edgmodidx[jadjcyidx++] = edgmodidxconn[reqidx][i][j];
    }
  }
  CkAssert(jadjcyidx == nsizedat);

  return mconn;
}

// Build Current (just includes sizes)
//
mConn* GeNet::BuildCurrConn() {
  // Initialize connection message
  int msgSize[MSG_Conn];
  msgSize[0] = 0;     // vtxmodidx
  msgSize[1] = 0;     // vtxordidx
  msgSize[2] = 0;     // xyz
  msgSize[3] = 0;     // xadj
  msgSize[4] = 0;     // adjcy
  msgSize[5] = 0;     // edgmodidx
  mConn *mconn = new(msgSize, 0) mConn;
  // Sizes
  mconn->datidx = datidx;
  mconn->nvtx = norderdat;

  return mconn;
}

// Build Next (includes vertices)
//
mConn* GeNet::BuildNextConn() {
  // Initialize connection message
  int msgSize[MSG_Conn];
  msgSize[0] = norderdat;   // vtxmodidx
  msgSize[1] = norderdat;   // vtxordidx
  msgSize[2] = norderdat*3; // xyz
  msgSize[3] = 0;           // xadj (we don't know these yet)
  msgSize[4] = 0;           // adjcy
  msgSize[5] = 0;           // edgmodidx
  mConn *mconn = new(msgSize, 0) mConn;
  // Sizes
  mconn->datidx = datidx;
  mconn->nvtx = norderdat;

  // Build message
  for (idx_t i = 0; i < norderdat; ++i) {
    mconn->vtxmodidx[i] = vtxmodidx[i];
    mconn->vtxordidx[i] = vtxordidx[i];
    mconn->xyz[i*3+0] = xyz[i*3+0];
    mconn->xyz[i*3+1] = xyz[i*3+1];
    mconn->xyz[i*3+2] = xyz[i*3+2];
  }

  return mconn;
}


/**************************************************************************
* Connection Construction
**************************************************************************/

idx_t GeNet::MakeConnection(idx_t source, idx_t target, idx_t sourceidx, idx_t targetidx, real_t dist) {
  // Go through edges to find the right connection
  // TODO: Use an unordered map to do the connection matching
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (source == edges[i].source) {
      for (std::size_t j = 0; j < edges[i].target.size(); ++j) {
        if (target == edges[i].target[j]) {
          // test for cutoff
          if (edges[i].cutoff != 0.0 && dist > edges[i].cutoff) {
            return 0;
          }
          // Connection computation
          real_t prob = 0.0;
          idx_t mask = 0;
          for (std::size_t k = 0; k < edges[i].conntype.size(); ++k) {
            if (edges[i].conntype[k] == CONNTYPE_UNIF) {
              prob += edges[i].probparam[k][0];
            }
            else if (edges[i].conntype[k] == CONNTYPE_SIG) {
              prob += sigmoid(dist, edges[i].probparam[k][0],
                        edges[i].probparam[k][1], edges[i].probparam[k][2]);
            }
            else if (edges[i].conntype[k] == CONNTYPE_IDX) {
              mask += (((sourceidx * edges[i].maskparam[k][2]) + edges[i].maskparam[k][3]) == targetidx);
            }
            else if (edges[i].conntype[k] == CONNTYPE_FILE) {
              // Check to see if it's in the file list
              // set mask to 1 if there is a non-zero entry
              // TODO: make sure file-based connections completely override
              //       other connection types (or make them mutually exclusive)
              if (sourceidx >= datafiles[(idx_t) (edges[i].probparam[k][0])].matrix.size()) {
                CkPrintf("  error: datafile %s does not have row for %" PRIidx "\n",
                         datafiles[(idx_t) (edges[i].probparam[k][0])].filename.c_str(), sourceidx);
              } else if (datafiles[(idx_t) (edges[i].probparam[k][0])].matrix[sourceidx].find(targetidx) ==
                         datafiles[(idx_t) (edges[i].probparam[k][0])].matrix[sourceidx].end()) {
                prob = 0.0;
                mask = 0;
              } else {
                mask = 1;
              }
            }
            else {
              // Shouldn't reach here due to prior error checking
              CkPrintf("  error: connection type %" PRIidx " undefined\n", edges[i].conntype[k]);
              CkExit();
            }
          }
          // Compute probability of connection
          if ((((*unifdist)(rngine)) < prob) || mask) {
            return edges[i].modidx;
          }
          else {
            return 0;
          }
        }
      }
    }
  }

  // no connection found, return 'none'
  return 0;
}


/**************************************************************************
* Edge State Building
**************************************************************************/

// States
//
std::vector<real_t> GeNet::BuildEdgState(idx_t modidx, real_t dist, idx_t sourceidx, idx_t targetidx) {
  // Sanity check
  // 0 is reserved for 'none' edge type
  CkAssert(modidx > 0);
  --modidx;
  CkAssert(models[modidx].type == GRAPHTYPE_EDG);
  // Allocate space for states
  std::vector<real_t> rngstate;
  rngstate.resize(models[modidx].statetype.size());
  // Randomly generate state
  for (std::size_t j = 0; j < rngstate.size(); ++j) {
    if (models[modidx].statetype[j] == RNGTYPE_CONST) {
      rngstate[j] = rngconst(models[modidx].stateparam[j].data());
    }
    else if (models[modidx].statetype[j] == RNGTYPE_UNIF) {
      rngstate[j] = rngunif(models[modidx].stateparam[j].data());
    }
    else if (models[modidx].statetype[j] == RNGTYPE_UNINT) {
      rngstate[j] = rngunint(models[modidx].stateparam[j].data());
    }
    else if (models[modidx].statetype[j] == RNGTYPE_NORM) {
      rngstate[j] = rngnorm(models[modidx].stateparam[j].data());
    }
    else if (models[modidx].statetype[j] == RNGTYPE_BNORM) {
      rngstate[j] = rngbnorm(models[modidx].stateparam[j].data());
    }
    else if (models[modidx].statetype[j] == RNGTYPE_LBNORM) {
      rngstate[j] = rnglbnorm(models[modidx].stateparam[j].data());
    }
    else if (models[modidx].statetype[j] == RNGTYPE_LIN) {
      rngstate[j] = rnglin(models[modidx].stateparam[j].data(), dist);
    }
    else if (models[modidx].statetype[j] == RNGTYPE_LBLIN) {
      rngstate[j] = rnglblin(models[modidx].stateparam[j].data(), dist);
    }
    else if (models[modidx].statetype[j] == RNGTYPE_BLIN) {
      rngstate[j] = rngblin(models[modidx].stateparam[j].data(), dist);
    }
    else if (models[modidx].statetype[j] == RNGTYPE_FILE) {
      rngstate[j] = rngfile(models[modidx].stateparam[j].data(), sourceidx, targetidx);
    }
    else {
      CkPrintf("  error: statetype %s is not valid for edge\n", rngtype[models[modidx].statetype[j]].c_str());
      // TODO: cleaner error checking here?
      CkExit();
    }
  }
  // return generated state
  return rngstate;
}

// Sticks
//
std::vector<tick_t> GeNet::BuildEdgStick(idx_t modidx, real_t dist, idx_t sourceidx, idx_t targetidx) {
  // Sanity check
  // 0 is reserved for 'none' edge type
  CkAssert(modidx > 0);
  --modidx;
  CkAssert(models[modidx].type == GRAPHTYPE_EDG);
  // Allocate space for sticks
  std::vector<tick_t> rngstick;
  rngstick.resize(models[modidx].sticktype.size());
  // Randomly generate stick
  for (std::size_t j = 0; j < rngstick.size(); ++j) {
    if (models[modidx].sticktype[j] == RNGTYPE_CONST) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngconst(models[modidx].stickparam[j].data()));
    }
    else if (models[modidx].sticktype[j] == RNGTYPE_UNIF) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngunif(models[modidx].stickparam[j].data()));
    }
    else if (models[modidx].sticktype[j] == RNGTYPE_UNINT) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngunint(models[modidx].stickparam[j].data()));
    }
    else if (models[modidx].sticktype[j] == RNGTYPE_NORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngnorm(models[modidx].stickparam[j].data()));
    }
    else if (models[modidx].sticktype[j] == RNGTYPE_BNORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngbnorm(models[modidx].stickparam[j].data()));
    }
    else if (models[modidx].sticktype[j] == RNGTYPE_LBNORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglbnorm(models[modidx].stickparam[j].data()));
    }
    else if (models[modidx].sticktype[j] == RNGTYPE_LIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglin(models[modidx].stickparam[j].data(), dist));
    }
    else if (models[modidx].sticktype[j] == RNGTYPE_LBLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglblin(models[modidx].stickparam[j].data(), dist));
    }
    else if (models[modidx].sticktype[j] == RNGTYPE_BLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngblin(models[modidx].stickparam[j].data(), dist));
    }
    else if (models[modidx].sticktype[j] == RNGTYPE_FILE) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngfile(models[modidx].stickparam[j].data(), sourceidx, targetidx));
    }
    else {
      CkPrintf("  error: statetype %s is not valid for edge\n", rngtype[models[modidx].sticktype[j]].c_str());
      // TODO: cleaner error checking here?
      CkExit();
    }
  }
  // return generated stick
  return rngstick;
}
