/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "stacs.h"
#include "network.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ unsigned randseed;
extern /*readonly*/ int netparts;
extern /*readonly*/ int netfiles;


/**************************************************************************
* Build
**************************************************************************/

// Build Network
//
void Netdata::Build(mGraph *msg) {
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
  
  // Bookkeeping across all populations
  xpopidxprt.resize(vertices.size());
  xvtxidxprt.resize(vertices.size());
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    xvtxidxprt[i].resize(netparts);
    xpopidxprt[i].resize(netparts+1);
  }
  xremvtx = 0;
  // Population vertex index
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    idx_t ndivvtx = (vertices[i].order)/netparts;
    idx_t nremvtx = (vertices[i].order)%netparts;
    for (int k = 0; k < netparts; ++k) {
      // looping integer magic
      xpopidxprt[i][k] = ndivvtx * k;
      // TODO: see if you can get rid of the loop
      for (int j = 0; j < k; ++j) {
        xpopidxprt[i][k] += ((j >= xremvtx && j < nremvtx+xremvtx) ||
            (nremvtx+xremvtx >= netparts && j < xremvtx && j < (nremvtx+xremvtx)%netparts));
      }
    }
    // Add order to final set
    xpopidxprt[i][netparts] = vertices[i].order;
    // update remainder
    xremvtx = (xremvtx+nremvtx)%netparts;
  }
  // Global population vertex index
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    for (int k = 0; k < netparts; ++k) {
      xvtxidxprt[i][k] = 0;
      for (std::size_t j = 0; j < i; ++j) {
        xvtxidxprt[i][k] += xpopidxprt[j][k+1];
      }
      for (std::size_t j = i; j < vertices.size(); ++j) {
        xvtxidxprt[i][k] += xpopidxprt[j][k];
      }
    }
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
      ordervtx << " " << nordervtx[k][i] << "(" << xordervtx[k][i] << ", " << xvtxidxprt[i][xprt+k] << ")";
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
        else if (vertices[i].shape == VTXSHAPE_LINE) {
          // uniformly along a line
          real_t l = vertices[i].param[0]*((*unifdist)(rngine));
          xyz[jvtxidx*3+0] = vertices[i].coord[0] + l;
          xyz[jvtxidx*3+1] = vertices[i].coord[1];
          xyz[jvtxidx*3+2] = vertices[i].coord[2];
        }
        else if (vertices[i].shape == VTXSHAPE_RECT) {
          // uniformly inside rectangle
          real_t w = vertices[i].param[0]*((*unifdist)(rngine));
          real_t h = vertices[i].param[1]*((*unifdist)(rngine));
          xyz[jvtxidx*3+0] = vertices[i].coord[0] + w;
          xyz[jvtxidx*3+1] = vertices[i].coord[1] + h;
          xyz[jvtxidx*3+2] = vertices[i].coord[2];
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
    CkAssert(modeldata[modidx].graphtype == GRAPHTYPE_VTX || modeldata[modidx].graphtype == GRAPHTYPE_STR);
    // Allocate space for states
    std::vector<real_t> rngstate;
    std::vector<tick_t> rngstick;
    rngstate.resize(modeldata[modidx].statetype.size());
    rngstick.resize(modeldata[modidx].sticktype.size());
    // Randomly generate state
    for (std::size_t s = 0; s < modeldata[modidx].statetype.size(); ++s) {
      if (modeldata[modidx].statetype[s] == RNGTYPE_CONST) {
        rngstate[s] = rngconst(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].statetype[s] == RNGTYPE_UNIF) {
        rngstate[s] = rngunif(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].statetype[s] == RNGTYPE_UNINT) {
        rngstate[s] = rngunint(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].statetype[s] == RNGTYPE_NORM) {
        rngstate[s] = rngnorm(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].statetype[s] == RNGTYPE_BNORM) {
        rngstate[s] = rngbnorm(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].statetype[s] == RNGTYPE_FILE) {
        rngstate[s] = rngfile(modeldata[modidx].stateparam[s].data(), 0, vtxordidx[i]);
      }
      else {
        CkPrintf("  error: statetype %s is not valid for vertex\n", rngtype[modeldata[modidx].statetype[s]].c_str());
        // TODO: cleaner error checking here?
        CkExit();
      }
    }
    // Randomly generate stick
    for (std::size_t s = 0; s < modeldata[modidx].sticktype.size(); ++s) {
      if (modeldata[modidx].sticktype[s] == RNGTYPE_CONST) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngconst(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].sticktype[s] == RNGTYPE_UNIF) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngunif(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].sticktype[s] == RNGTYPE_UNINT) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngunint(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].sticktype[s] == RNGTYPE_NORM) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngnorm(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].sticktype[s] == RNGTYPE_BNORM) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngbnorm(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].sticktype[s] == RNGTYPE_FILE) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngfile(modeldata[modidx].stickparam[s].data(), 0, vtxordidx[i]));
      }
      else {
        CkPrintf("  error: statetype %s is not valid for vertex\n", rngtype[modeldata[modidx].sticktype[s]].c_str());
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
  adjcylocalcount.resize(norderdat);
  // Only need to worry about future edges
  adjcyconn.clear();
  adjcyconn.resize(netfiles);
  edgmodidxconn.clear();
  edgmodidxconn.resize(netfiles);
  adjcyset.resize(norderdat);

  // Any index-based sample connectivity occurs first
  // for each vertex
  // TODO: this is only for incoming connections, also need outgoing connections (none models)
  for (idx_t i = 0; i < norderdat; ++i) {
    idx_t globalthisidx = 0;
    adjcyset[i].clear();
    // TODO: this is highly innefficient...
    for (int prt = 0; prt < netparts; ++prt) {
      if (vtxordidx[i] > xpopidxprt[vtxmodidx[i]-1][prt]) {
        globalthisidx = xvtxidxprt[vtxmodidx[i]-1][prt] + vtxordidx[i] - xpopidxprt[vtxmodidx[i]-1][prt];
        break;
      }
    }
    // for each edge population
    for (std::size_t e = 0; e < edges.size(); ++e) {
      // Sample types should only have one target at the moment
      if (vtxmodidx[i] == edges[e].target[0]) {
        // But they're still structured to have multiple connection types
        for (std::size_t k = 0; k < edges[e].conntype.size(); ++k) {
          if (edges[e].conntype[k] == CONNTYPE_SMPL) {
            // Sample sourceidx from the source vertex population order
            // generate the sample cache (once per target vertex per connection type)
            std::vector<idx_t> sourceorder(edges[e].maskparam[k][0]);
            std::iota(sourceorder.begin(), sourceorder.end(), 0);
            // pick the seed based on the targetidx so it is consistent across cores
            // The 32768 is 2^15, is just a reasonably large number to not get repeating seeds
            unsigned sampleseed = (randseed + (unsigned)(vtxordidx[i])) ^ ((unsigned)(e*32768));
            std::shuffle(sourceorder.begin(), sourceorder.end(), std::mt19937{sampleseed});
            // want to make sure sample number is less than source order
            CkAssert(edges[e].maskparam[k][0] >= edges[e].maskparam[k][1]);
            // copy over the shuffled indices for the sampling
            for (idx_t j = 0; j < edges[e].maskparam[k][1]; ++j) {
              // Convert from population index to global index
              idx_t globalsourceidx = 0;
              // TODO: this is highly innefficient...
              for (int prt = 0; prt < netparts; ++prt) {
                if (sourceorder[j] >= xpopidxprt[edges[e].source-1][prt] && sourceorder[j] < xpopidxprt[edges[e].source-1][prt+1]) {
                  globalsourceidx = xvtxidxprt[edges[e].source-1][prt] + (sourceorder[j] - xpopidxprt[edges[e].source-1][prt]);
                  break;
                }
              }
              // TODO: put in a check to not insert self-connections
              if (globalsourceidx == globalthisidx) {
                continue;
              }
              // this is just to check on memory usage
              adjcy[i].push_back(globalsourceidx);
              adjcyset[i].insert(globalsourceidx); // The set is useful for faster searching of edge existence
              edgmodidx[i].push_back(edges[e].modidx);
              // The state/stick will need to be reparameterized with correct distance information later
              // TODO: make a version that doesn't use dist, idxs, depending on the edge model
              state[i].push_back(BuildEdgState(edges[e].modidx, 0.0, sourceorder[j], vtxordidx[i]));
              stick[i].push_back(BuildEdgStick(edges[e].modidx, 0.0, sourceorder[j], vtxordidx[i]));
            }
          }
          else if (edges[e].conntype[k] == CONNTYPE_SMPL_NORM) {
            // Sample sourceidx from the source vertex population order
            // generate the sample cache (once per target vertex per connection type)
            std::vector<real_t> sourceweights(edges[e].maskparam[k][0]);
            for (idx_t j = 0; j < edges[e].maskparam[k][0]; ++j) {
              // (x_i - x_j)^2 / var(x_ij)
              real_t x_ij = ((((real_t) vtxordidx[i])/vertices[vtxmodidx[i]-1].order)-(((real_t) j)/edges[e].maskparam[k][0]));
              sourceweights[j] = std::exp(-(x_ij*x_ij)/(2*edges[e].probparam[k][0])); // Don't worry about normalizing
            }
            // pick the seed based on the targetidx so it is consistent across cores
            // The 32768 is 2^15, is just a reasonably large number to not get repeating seeds
            unsigned sampleseed = (randseed + (unsigned)(vtxordidx[i])) ^ ((unsigned)(e*32768));
            std::mt19937 rngsample(sampleseed);
            // want to make sure sample number is less than source order
            CkAssert(edges[e].maskparam[k][0] >= edges[e].maskparam[k][1]);
            adjcyset[i].clear();
            while (adjcyset[i].size() < edges[e].maskparam[k][1]) {
              std::discrete_distribution<idx_t> sampledist(sourceweights.begin(), sourceweights.end()); 
              idx_t sourceorder = sampledist(rngsample);
              // Convert from population index to global index
              idx_t globalsourceidx = 0;
              // TODO: this is highly innefficient...
              for (int prt = 0; prt < netparts; ++prt) {
                if (sourceorder >= xpopidxprt[edges[e].source-1][prt] && sourceorder < xpopidxprt[edges[e].source-1][prt+1]) {
                  globalsourceidx = xvtxidxprt[edges[e].source-1][prt] + (sourceorder - xpopidxprt[edges[e].source-1][prt]);
                  break;
                }
              }
              if (globalsourceidx == globalthisidx) {// || adjcyset[i].find(globalsourceidx) == adjcyset[i].end()) {
                continue;
              } else {
                sourceweights[sourceorder] = 0.0;
                adjcy[i].push_back(globalsourceidx);
                adjcyset[i].insert(globalsourceidx); // The set is useful for faster searching of edge existence
                edgmodidx[i].push_back(edges[e].modidx);
                // The state/stick will need to be reparameterized with correct distance information later
                state[i].push_back(BuildEdgState(edges[e].modidx, 0.0, sourceorder, vtxordidx[i]));
                stick[i].push_back(BuildEdgStick(edges[e].modidx, 0.0, sourceorder, vtxordidx[i]));
              }
            }
          }
          else if (edges[e].conntype[k] == CONNTYPE_SMPL_ANORM) {
            // Sample sourceidx from the source vertex population order
            // generate the sample cache (once per target vertex per connection type)
            std::vector<real_t> sourceweights(edges[e].maskparam[k][0]);
            for (idx_t j = 0; j < edges[e].maskparam[k][0]; ++j) {
              // ((x_i - x_j)^2 / var(x_ij)) - ((x_i - x_j)^2 / var(x_ij)/2)
              // Not using variance-y for this for now (needs additional information)
              real_t x_ij = ((((real_t) vtxordidx[i])/vertices[vtxmodidx[i]-1].order)-(((real_t) j)/edges[e].maskparam[k][0]));
              real_t var = edges[e].probparam[k][0];
              // Normalizing a bit more important here
              real_t wgt = std::exp(-(x_ij*x_ij)/(2*var))/std::sqrt(var) *
                (std::exp(-(x_ij*x_ij)/(2*var))/std::sqrt(var) - std::exp(-(x_ij*x_ij)/(var))/std::sqrt(var/2));
              sourceweights[j] = std::max(0.0, wgt);
            }
            // pick the seed based on the targetidx so it is consistent across cores
            // The 32768 is 2^15, is just a reasonably large number to not get repeating seeds
            unsigned sampleseed = (randseed + (unsigned)(vtxordidx[i])) ^ ((unsigned)(e*32768));
            std::mt19937 rngsample(sampleseed);
            std::uniform_real_distribution<real_t> sampleunifdist;
            // want to make sure sample number is less than source order
            CkAssert(edges[e].maskparam[k][0] >= edges[e].maskparam[k][1]);
            adjcyset[i].clear();
            // From: https://stackoverflow.com/questions/53632441/c-sampling-from-discrete-distribution-without-replacement
            std::vector<real_t> sourcevals;
            for (auto iter : sourceweights) {
              sourcevals.push_back(std::pow(sampleunifdist(rngsample), 1. / iter));
            }
            // Sorting vals, but retain the indices. 
            // There is unfortunately no easy way to do this with STL.
            std::vector<std::pair<idx_t, real_t>> valswithindices;
            for (std::size_t iter = 0; iter < sourcevals.size(); iter++) {
              valswithindices.emplace_back(iter, sourcevals[iter]);
            }
            std::sort(valswithindices.begin(), valswithindices.end(), [](std::pair<idx_t,real_t> x, std::pair<idx_t,real_t> y) {return x.second > y.second; });
            for (idx_t iter = 0; iter < edges[e].maskparam[k][1]; iter++) {
              idx_t sourceorder = valswithindices[iter].first;
              // Convert from population index to global index
              idx_t globalsourceidx = 0;
              // TODO: this is highly innefficient...
              for (int prt = 0; prt < netparts; ++prt) {
                if (sourceorder >= xpopidxprt[edges[e].source-1][prt] && sourceorder < xpopidxprt[edges[e].source-1][prt+1]) {
                  globalsourceidx = xvtxidxprt[edges[e].source-1][prt] + (sourceorder - xpopidxprt[edges[e].source-1][prt]);
                  break;
                }
              }
              if (globalsourceidx == globalthisidx) {// || adjcyset[i].find(globalsourceidx) == adjcyset[i].end()) {
                continue;
              } else {
                adjcy[i].push_back(globalsourceidx);
                adjcyset[i].insert(globalsourceidx); // The set is useful for faster searching of edge existence
                edgmodidx[i].push_back(edges[e].modidx);
                // The state/stick will need to be reparameterized with correct distance information later
                state[i].push_back(BuildEdgState(edges[e].modidx, 0.0, sourceorder, vtxordidx[i]));
                stick[i].push_back(BuildEdgStick(edges[e].modidx, 0.0, sourceorder, vtxordidx[i]));
              }
            }
          }
        }
      }
    }
    // This is just the size of the adjcy before adding the none models
    adjcylocalcount[i] = adjcy[i].size();
  }
  
  // Reorder edges vertex-by-vertex
  for (idx_t i = 0; i < norderdat; ++i) {
    edgorder.clear();
    for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
      edgorder.push_back(edgorder_t());
      edgorder.back().edgidx = adjcy[i][j];
      edgorder.back().modidx = edgmodidx[i][j];
      edgorder.back().state = state[i][j+1];
      edgorder.back().stick = stick[i][j+1];
    }
    // sort edge indices by global ordering
    std::sort(edgorder.begin(), edgorder.end());
    // add indices to data structures
    for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
      adjcy[i][j] = edgorder[j].edgidx;
      edgmodidx[i][j] = edgorder[j].modidx;
      state[i][j+1] = edgorder[j].state;
      stick[i][j+1] = edgorder[j].stick;
    }
  }

  //Build up the vtxdist from knowledge of how vertices will be distributed

  // Contribute for now
  //contribute(0, NULL, CkReduction::nop);

  // At this point, the partition should have all incoming connnections, but no knowledge of its outgoing connections
  // Request data from prev part
  if (cpdat < datidx) {
    thisProxy(cpdat).ConnNoneRequest(datidx);
  }
  // Connect to curr part
  else if (cpdat == datidx) {
    mConnNone *mconn = BuildConnNone(datidx);
    thisProxy(cpdat).ConnectNone(mconn);
  }
  // Request data from next part
  else if (cpdat > datidx) {
    // this shouldn't happen
    thisProxy(cpdat).ConnNoneRequest(datidx);
  }
  
  /*
  //TODO: Figure out a way to perform both sampling as above
  //      and the per-connection evaluation as below (hopefully faster)
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
  */
}


/**************************************************************************
* Connect
**************************************************************************/

// Some notes on Prev, Curr, and Next:
// The convention describes the sending partition with respect to the recieving partition
// Connecting to Prev means that the sending partition is at a lower index
// Connecting to Curr means that the sending partition is the same index
// Connecting to Next means that the sending partition is at a higher index
// Connections are built in a pipeline, in order, partition by partition which allows for
//   the indices to be written into the adjcy list in order (no need to sort)
//   however, this may be a bit slow when dealing with much larger networks
//   as it effectively has to evaluate the pair-wise connection probability between
//   each vertex (in order), but this also allows it to use proximity/distance information
//   there is an early cutoff where you may specify a maximum range to evaluate, and also
//   evaluating if there is a edge model that goes between the two vertex populations,
//   but these evaluations still need to happen for connections requiring local information
//   this is because the local vertex information is stored on a per-partition level and
//   thus needs to be communicated to the connecting partition when needed
// The work on this new branch for sample-based connections only works for networks that
//   connect based on ordinal (index-based) connectivity, where the assumption is that
//   we no longer need proximity information to perform connection existence, but we may
//   need it to determine the instantiation of the connection parameters (e.g. delay)
//   as a result, this modification will go through and create the adjcy lists first
//   then go through and instantiate the connection parameters via the pipeline

// Connect Network
//
void Netdata::Connect(mConn *msg) {
  // Sanity check
  CkAssert(msg->datidx == cpdat);
  // Some basic information on what's being connected
  CkPrintf("  Connecting %d to %d\n", datidx, msg->datidx);

  // Add to vtxdist
  vtxdist[cpdat+1] = vtxdist[cpdat] + msg->nvtx;

  // Perform connections
  adjcyconn[msg->datidx].resize(norderdat);
  edgmodidxconn[msg->datidx].resize(norderdat);

  // Sample-based cache
  samplecache.resize(edges.size());
  
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
    // We can reorder all the edges by global ordering now
    // Reorder edges vertex-by-vertex
    for (idx_t i = 0; i < norderdat; ++i) {
      edgorder.clear();
      for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
        edgorder.push_back(edgorder_t());
        edgorder.back().edgidx = adjcy[i][j];
        edgorder.back().modidx = edgmodidx[i][j];
        edgorder.back().state = state[i][j+1];
        edgorder.back().stick = stick[i][j+1];
      }
      // sort edge indices by global ordering
      std::sort(edgorder.begin(), edgorder.end());
      // add indices to data structures
      for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
        adjcy[i][j] = edgorder[j].edgidx;
        edgmodidx[i][j] = edgorder[j].modidx;
        state[i][j+1] = edgorder[j].state;
        stick[i][j+1] = edgorder[j].stick;
      }
    }

    // Print memory allocated
    int adjcysize = 0;
    int adjcycap = 0;
    int edgmodsize = 0;
    int edgmodcap = 0;
    for (size_t i = 0; i < adjcy.size(); ++i) {
      adjcy[i].shrink_to_fit();
      edgmodidx[i].shrink_to_fit();
      adjcysize += adjcy[i].size();
      adjcycap += adjcy[i].capacity();
      edgmodsize += edgmodidx[i].size();
      edgmodcap += edgmodidx[i].capacity();
    }
    CkPrintf("Part %d size/cap: adjcy: %d , %d edgmodidx: %d , %d\n", datidx, adjcysize, adjcycap, edgmodsize, edgmodcap);
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

  // cleanup
  adjcyconn[msg->datidx].clear();
  edgmodidxconn[msg->datidx].clear();
  samplecache.clear();
  delete msg;
}

// Connect Network Finished
//
void Netdata::ConnRequest(idx_t reqidx) {
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
mConn* Netdata::BuildPrevConn(idx_t reqidx) {
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
mConn* Netdata::BuildCurrConn() {
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
mConn* Netdata::BuildNextConn() {
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

idx_t Netdata::MakeConnection(idx_t source, idx_t target, idx_t sourceidx, idx_t targetidx, real_t dist) {
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
            // TODO: this is currently set to sample from source order
            // TODO: we want to enable a weighting of the indices w.r.t. distance
            else if (edges[i].conntype[k] == CONNTYPE_SMPL) {
              // Sample sourceidx from the source vertex population order
              if (samplecache[i].find(targetidx) == samplecache[i].end()) {
                // generate the sample cache (once per target vertex per connection type)
                std::vector<idx_t> sourceorder(edges[i].maskparam[k][0]);
                std::iota(sourceorder.begin(), sourceorder.end(), 0);
                // pick the seed based on the targetidx so it is consistent across cores
                unsigned sampleseed = (randseed + (unsigned)(targetidx)) ^ ((unsigned)(i*32768));
                std::shuffle(sourceorder.begin(), sourceorder.end(), std::mt19937{sampleseed});
                // want to make sure sample number is less than source order
                CkAssert(edges[i].maskparam[k][0] >= edges[i].maskparam[k][1]);
                // copy over the shuffled indices for the sampling
                samplecache[i][targetidx].resize(edges[i].maskparam[k][1]);
                std::copy(sourceorder.begin(), sourceorder.begin() + edges[i].maskparam[k][1], samplecache[i][targetidx].begin());
                //CkPrintf("%" PRIidx ", %" PRIidx" \n",targetidx, samplecache[i][targetidx][0]);
              }
              // Check if the sourceidx is within the sampled sources
              // TODO: make this into a binary search?
              if (std::find(samplecache[i][targetidx].begin(), samplecache[i][targetidx].end(), sourceidx) != samplecache[i][targetidx].end()) {
                mask = 1;
              }
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
std::vector<real_t> Netdata::BuildEdgState(idx_t modidx, real_t dist, idx_t sourceidx, idx_t targetidx) {
  // Sanity check
  // 0 is reserved for 'none' edge type
  CkAssert(modidx > 0);
  --modidx;
  CkAssert(modeldata[modidx].graphtype == GRAPHTYPE_EDG);
  // Allocate space for states
  std::vector<real_t> rngstate;
  rngstate.resize(modeldata[modidx].statetype.size());
  // Randomly generate state
  for (std::size_t j = 0; j < rngstate.size(); ++j) {
    if (modeldata[modidx].statetype[j] == RNGTYPE_CONST) {
      rngstate[j] = rngconst(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_UNIF) {
      rngstate[j] = rngunif(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_UNINT) {
      rngstate[j] = rngunint(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_NORM) {
      rngstate[j] = rngnorm(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_BNORM) {
      rngstate[j] = rngbnorm(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_LBNORM) {
      rngstate[j] = rnglbnorm(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_LIN) {
      rngstate[j] = rnglin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_LBLIN) {
      rngstate[j] = rnglblin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_BLIN) {
      rngstate[j] = rngblin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_FILE) {
      rngstate[j] = rngfile(modeldata[modidx].stateparam[j].data(), sourceidx, targetidx);
    }
    else {
      CkPrintf("  error: statetype %s is not valid for edge\n", rngtype[modeldata[modidx].statetype[j]].c_str());
      // TODO: cleaner error checking here?
      CkExit();
    }
  }
  // return generated state
  return rngstate;
}

// Sticks
//
std::vector<tick_t> Netdata::BuildEdgStick(idx_t modidx, real_t dist, idx_t sourceidx, idx_t targetidx) {
  // Sanity check
  // 0 is reserved for 'none' edge type
  CkAssert(modidx > 0);
  --modidx;
  CkAssert(modeldata[modidx].graphtype == GRAPHTYPE_EDG);
  // Allocate space for sticks
  std::vector<tick_t> rngstick;
  rngstick.resize(modeldata[modidx].sticktype.size());
  // Randomly generate stick
  for (std::size_t j = 0; j < rngstick.size(); ++j) {
    if (modeldata[modidx].sticktype[j] == RNGTYPE_CONST) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngconst(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_UNIF) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngunif(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_UNINT) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngunint(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_NORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngnorm(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_BNORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngbnorm(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_LBNORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglbnorm(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_LIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_LBLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglblin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_BLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngblin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_FILE) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngfile(modeldata[modidx].stickparam[j].data(), sourceidx, targetidx));
    }
    else {
      CkPrintf("  error: statetype %s is not valid for edge\n", rngtype[modeldata[modidx].sticktype[j]].c_str());
      // TODO: cleaner error checking here?
      CkExit();
    }
  }
  // return generated stick
  return rngstick;
}


// The following is W.I.P for sending none edges between fully build directed graph
// This enables the source partition to know where to send

// Connect Network
//
void Netdata::ConnectNone(mConnNone *msg) {
  // Sanity check
  CkAssert(msg->datidx == cpdat);
  // Some basic information on what's being connected
  CkPrintf("  Connecting %d to %d\n", datidx, msg->datidx);
  
  // Go through the vertices
  for (idx_t j = 0; j < msg->nvtx; ++j) {
    idx_t targetidx = msg->vtxidx[j];
    for (idx_t i = msg->xadj[j]; i < msg->xadj[j+1]; ++i) {
      idx_t sourceidx = msg->adjcy[i];
      //CkPrintf("    i: %" PRIidx ", j: %" PRIidx "\n", sourceidx, targetidx);
      // Make sure the index isn't already in adjcy (both outgoing and incoming edges)
      if (adjcyset[sourceidx].find(targetidx) == adjcyset[sourceidx].end()) {
        adjcy[sourceidx].push_back(targetidx);
        edgmodidx[sourceidx].push_back(0); // These are all 'none' models
        // build empty state
        state[sourceidx].push_back(std::vector<real_t>());
        stick[sourceidx].push_back(std::vector<tick_t>());
      }
    }
  }
  // Go through the incoming state and reparameterize the state/sticks if needed
  idx_t globalsourceidx_min = msg->vtxidx[0];
  idx_t globalsourceidx_max = msg->vtxidx[msg->nvtx-1];
  for (idx_t i = 0; i < norderdat; ++i) {
    for (std::size_t j = 0; j < adjcylocalcount[i]; ++j) {
      if (adjcy[i][j] >= globalsourceidx_min && adjcy[i][j] <= globalsourceidx_max) {
        // Reparameterize with information
        idx_t localsourceidx = adjcy[i][j] - globalsourceidx_min;
        real_t distance = sqrt((xyz[i*3]-msg->xyz[localsourceidx*3])*(xyz[i*3]-msg->xyz[localsourceidx*3])+
                          (xyz[i*3+1]-msg->xyz[localsourceidx*3+1])*(xyz[i*3+1]-msg->xyz[localsourceidx*3+1])+
                          (xyz[i*3+2]-msg->xyz[localsourceidx*3+2])*(xyz[i*3+2]-msg->xyz[localsourceidx*3+2]));              
        ReBuildEdgState(edgmodidx[i][j], distance, state[i][j+1]);
        ReBuildEdgStick(edgmodidx[i][j], distance, stick[i][j+1]);
      }
    }
  }

  // send any outstanding requests of built adjcy
  for (std::list<idx_t>::iterator ireqidx = adjcyreq.begin(); ireqidx != adjcyreq.end(); ++ireqidx) {
    if (*ireqidx == cpdat) {
      mConnNone *mconn = BuildConnNone(*ireqidx);
      thisProxy(*ireqidx).ConnectNone(mconn);
      ireqidx = adjcyreq.erase(ireqidx);
    }
  }

  delete msg;

  // Move to next part
  ++cpdat;
  // return control to main when done
  if (cpdat == netfiles) {
    // We can reorder all the edges by global ordering now
    // Reorder edges vertex-by-vertex
    for (idx_t i = 0; i < norderdat; ++i) {
      edgorder.clear();
      for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
        edgorder.push_back(edgorder_t());
        edgorder.back().edgidx = adjcy[i][j];
        edgorder.back().modidx = edgmodidx[i][j];
        edgorder.back().state = state[i][j+1];
        edgorder.back().stick = stick[i][j+1];
      }
      // sort edge indices by global ordering
      std::sort(edgorder.begin(), edgorder.end());
      // add indices to data structures
      for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
        adjcy[i][j] = edgorder[j].edgidx;
        edgmodidx[i][j] = edgorder[j].modidx;
        state[i][j+1] = edgorder[j].state;
        stick[i][j+1] = edgorder[j].stick;
      }
    }
    /*
    // Print memory allocated
    int adjcysize = 0;
    int adjcycap = 0;
    int edgmodsize = 0;
    int edgmodcap = 0;
    for (size_t i = 0; i < adjcy.size(); ++i) {
      adjcy[i].shrink_to_fit();
      edgmodidx[i].shrink_to_fit();
      adjcysize += adjcy[i].size();
      adjcycap += adjcy[i].capacity();
      edgmodsize += edgmodidx[i].size();
      edgmodcap += edgmodidx[i].capacity();
    }
    CkPrintf("Part %d size/cap: adjcy: %d , %d edgmodidx: %d , %d\n", datidx, adjcysize, adjcycap, edgmodsize, edgmodcap);
    */
    contribute(0, NULL, CkReduction::nop);
  }
  // Request data from next part
  else if (cpdat < datidx) {
    thisProxy(cpdat).ConnNoneRequest(datidx);
  }
  // Connect to curr part
  else if (cpdat == datidx) {
    mConnNone *mconn = BuildConnNone(datidx);
    thisProxy(cpdat).ConnectNone(mconn);
  }
  // Request data from next part
  else if (cpdat > datidx) {
    thisProxy(cpdat).ConnNoneRequest(datidx);
  }
}

// Connect Network Finished
//
void Netdata::ConnNoneRequest(idx_t reqidx) {
  // send data adjacency and vertex info to requesting part
  if (reqidx > datidx) {
    // check if vertex is built
    if (vtxmodidx.size()) {
      mConnNone *mconn = BuildConnNone(reqidx);
      thisProxy(reqidx).ConnectNone(mconn);
    }
    else {
      // record request for when adjcy is built
      adjcyreq.push_back(reqidx);
    }
  }
  // send data (vertex info) to requesting part
  else if (reqidx < datidx) {
    mConnNone *mconn = BuildConnNone(reqidx);
    thisProxy(reqidx).ConnectNone(mconn);
  }
}


/**************************************************************************
* Connection messages
**************************************************************************/

// Build Previous (includes vertices and adjacency)
//
mConnNone* Netdata::BuildConnNone(idx_t reqidx) {
  /* Bookkeeping */
  idx_t nsizedat;
  idx_t jadjcyidx;

  std::vector<std::vector<idx_t>> adjcyconnnone;
  adjcyconnnone.resize(norderdat);

  // Figure out the min and max global indices on partition reqidx
  idx_t globalsource_min = xvtxidxprt[0][reqidx];
  idx_t globalsource_max = 0;
  if (reqidx+1 == netparts) {
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      globalsource_max += vertices[i].order;
    }
  } else {
    globalsource_max = xvtxidxprt[0][reqidx+1];
  }

  // Count the sizes
  nsizedat = 0;
  for (idx_t i = 0; i < norderdat; ++i) {
    adjcyconnnone[i].clear();
    for (std::size_t j = 0; j < adjcylocalcount[i]; ++j) {
      if (globalsource_min <= adjcy[i][j] && adjcy[i][j] < globalsource_max) {
        // global source idx to local
        // should really have an xvtxidxdat
        adjcyconnnone[i].push_back(adjcy[i][j] - xvtxidxprt[0][reqidx]);
      }
    }
    // Add none connections to size
    nsizedat += adjcyconnnone[i].size();
  }
  CkPrintf("   reqidx: %" PRIidx ", min: %" PRIidx ", max: %" PRIidx ", adj: %" PRIidx "\n", reqidx, globalsource_min, globalsource_max, nsizedat);

  // Initialize connection message
  int msgSize[MSG_ConnNone];
  msgSize[0] = norderdat;   // vtxidx
  msgSize[1] = norderdat*3; // xyz
  msgSize[2] = norderdat+1; // xadj
  msgSize[3] = nsizedat;    // adjcy
  mConnNone *mconn = new(msgSize, 0) mConnNone;
  // Sizes
  mconn->datidx = datidx;
  mconn->nvtx = norderdat;

  // Prefix zero is zero
  mconn->xadj[0] = 0;
  // Initialize counter
  jadjcyidx = 0;
  
  // Build message
  for (idx_t i = 0; i < norderdat; ++i) {
    mconn->vtxidx[i] = xvtxidxprt[0][datidx]+i; // need to convert from local to global
    mconn->xyz[i*3+0] = xyz[i*3+0];
    mconn->xyz[i*3+1] = xyz[i*3+1];
    mconn->xyz[i*3+2] = xyz[i*3+2];
    // xadj
    mconn->xadj[i+1] = mconn->xadj[i] + adjcyconnnone[i].size();
    for (std::size_t j = 0; j < adjcyconnnone[i].size(); ++j) {
      // adjcy (of next parts)
      mconn->adjcy[jadjcyidx++] = adjcyconnnone[i][j];
    }
  }
  CkAssert(jadjcyidx == nsizedat);

  return mconn;
}

/**************************************************************************
* Edge State Building (reparameterize)
**************************************************************************/

// States
//
void Netdata::ReBuildEdgState(idx_t modidx, real_t dist, std::vector<real_t>& rngstate) {
  // Sanity check
  // 0 is reserved for 'none' edge type
  CkAssert(modidx > 0);
  --modidx;
  CkAssert(modeldata[modidx].graphtype == GRAPHTYPE_EDG);
  // Randomly generate state
  for (std::size_t j = 0; j < rngstate.size(); ++j) {
    if (modeldata[modidx].statetype[j] == RNGTYPE_LIN) {
      rngstate[j] = rnglin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_LBLIN) {
      rngstate[j] = rnglblin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].statetype[j] == RNGTYPE_BLIN) {
      rngstate[j] = rngblin(modeldata[modidx].stateparam[j].data(), dist);
    }
  }
}

// Sticks
//
void Netdata::ReBuildEdgStick(idx_t modidx, real_t dist, std::vector<tick_t>& rngstick) {
  // Sanity check
  // 0 is reserved for 'none' edge type
  CkAssert(modidx > 0);
  --modidx;
  CkAssert(modeldata[modidx].graphtype == GRAPHTYPE_EDG);
  // Randomly generate stick
  for (std::size_t j = 0; j < rngstick.size(); ++j) {
    if (modeldata[modidx].sticktype[j] == RNGTYPE_LIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_LBLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglblin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].sticktype[j] == RNGTYPE_BLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngblin(modeldata[modidx].stickparam[j].data(), dist));
    }
  }
}

