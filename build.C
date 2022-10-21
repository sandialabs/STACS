/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
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
  idx_t jedgdistparam;
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
  jedgdistparam = 0;
  jedgconntype = 0;
  jedgprobparam = 0;
  jedgmaskparam = 0;

  // Edges
  edges.resize(msg->nedg);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    edges[i].source = msg->edgsource[i];
    edges[i].modidx = msg->edgmodidx[i];
    edges[i].cutoff = msg->edgcutoff[i];
    edges[i].distype = msg->edgdistype[i];
    
    edges[i].target.resize(msg->xedgtarget[i+1] - msg->xedgtarget[i]);
    for (std::size_t j = 0; j < edges[i].target.size(); ++j) {
      edges[i].target[j] = msg->edgtarget[jedgtarget++];
    }

    edges[i].distparam.resize(msg->medgdistparam[i]);
    for (std::size_t j = 0; j < edges[i].distparam.size(); ++j) {
      edges[i].distparam[j] = msg->edgdistparam[jedgdistparam++];
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
  CkAssert(jedgdistparam == msg->nedgdistparam);
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

  // Create mapping from pairs of vertices (source, target) to edge index
  connmodmap.clear();
  // Create mapping from vertices to sample-based edge sources
  connsampleset.resize(vertices.size()+1);
  for (std::size_t t = 0; t < vertices.size(); ++t) {
    connsampleset[t].clear();
  }
  // populate connection maps
  for (std::size_t e = 0; e < edges.size(); ++e) {
    for (std::size_t t = 0; t < edges[e].target.size(); ++t) {
      // general connections
      connmodmap[edges[e].source*edges.size() + edges[e].target[t]] = (idx_t) e;
      // sample-based connections
      for (std::size_t k = 0; k < edges[e].conntype.size(); ++k) {
        if (edges[e].conntype[k] == CONNTYPE_SMPL ||
            edges[e].conntype[k] == CONNTYPE_SMPL_NORM ||
            edges[e].conntype[k] == CONNTYPE_SMPL_ANORM) {
          connsampleset[edges[e].target[t]].insert((idx_t) e);
        }
      }
    }
  }


  // Bookkeeping to see how much each chare builds
  // Taking into account the different parts too
  // Initial distribution is simply contructing
  // vertices as evenly as possible across the parts
  xpopvtxidxprt.resize(vertices.size()); // Prefix vertex index within population (per partition)
  xglbvtxidxprt.resize(vertices.size()); // Prefix vertex index globally (per partition)
  // Loop through vertex populations
  idx_t xremvtx = 0;
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    idx_t ndivvtx = (vertices[i].order)/netparts;
    idx_t nremvtx = (vertices[i].order)%netparts;
    xpopvtxidxprt[i].resize(netparts+1);
    // By population
    for (idx_t k = 0; k < netparts; ++k) {
      xpopvtxidxprt[i][k] = ndivvtx * k;
      for (idx_t j = 0; j < k; ++j) {
        xpopvtxidxprt[i][k] += ((j >= xremvtx && j < nremvtx+xremvtx) ||
            (nremvtx+xremvtx >= netparts && j < xremvtx && j < (nremvtx+xremvtx)%netparts));
      }
    }
    // Add order to final set
    xpopvtxidxprt[i][netparts] = vertices[i].order;
    // update remainder
    xremvtx = (xremvtx+nremvtx)%netparts;
  }
  // Global population vertex index
  for (std::size_t i = 0; i < vertices.size(); ++i) {
    xglbvtxidxprt[i].resize(netparts+1);
    for (int k = 0; k < netparts; ++k) {
      xglbvtxidxprt[i][k] = 0;
      for (std::size_t j = 0; j < i; ++j) {
        xglbvtxidxprt[i][k] += xpopvtxidxprt[j][k+1];
      }
      for (std::size_t j = i; j < vertices.size(); ++j) {
        xglbvtxidxprt[i][k] += xpopvtxidxprt[j][k];
      }
    }
    xglbvtxidxprt[i][netparts] = norder;
  }
  // Note: xglbvtxidxprt[0][xprt] gives the global vtx offset for datidx
  vtxdist.resize(netparts+1);
  for (int k = 0; k < netparts; ++k) {
    vtxdist[k] = xglbvtxidxprt[0][k];
  }
  vtxdist[netparts] = norder;
  xorderdat.resize(netfiles+1);
  for (int k = 0; k < netfiles; ++k) {
    int ndiv = netparts/netfiles;
    int nrem = netparts%netfiles;
    int xprt = k*ndiv + (k < nrem ? k : nrem);
    xorderdat[k] = xglbvtxidxprt[0][xprt];
  }
  xorderdat[netfiles] = norder;

  // Compute how much this chare builds
  norderdat = 0;
  norderprt.resize(nprt);
  // population specific sizes
  std::vector<std::vector<idx_t>> nordervtx;
  std::vector<std::vector<idx_t>> xordervtx;
  nordervtx.resize(nprt);
  xordervtx.resize(nprt);
  for (idx_t k = 0; k < nprt; ++k) {
    norderprt[k] = 0;
    nordervtx[k].resize(vertices.size());
    xordervtx[k].resize(vertices.size());
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      nordervtx[k][i] = xpopvtxidxprt[i][xprt+k+1] - xpopvtxidxprt[i][xprt+k];
      xordervtx[k][i] = xpopvtxidxprt[i][xprt+k];
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
      ordervtx << " " << nordervtx[k][i] << "(" << xordervtx[k][i] << ", " << xglbvtxidxprt[i][xprt+k] << ")";
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
  adjcy.resize(norderdat);
  adjcyset.resize(norderdat);
  nadjcysample.resize(norderdat);
  idx_t jvtxidx = 0;
  for (idx_t k = 0; k < nprt; ++k) {
    // set with modidx
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      for (idx_t j = 0; j < nordervtx[k][i]; ++j) {
        // Set the model index
        vtxmodidx[jvtxidx] = vertices[i].modidx;
        vtxordidx[jvtxidx] = xordervtx[k][i] + j;
        edgmodidx[jvtxidx].clear();
        adjcy[jvtxidx].clear();
        adjcyset[jvtxidx].clear();
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
        else if (vertices[i].shape == VTXSHAPE_SPHERE_SURFACE) {
          // uniformly on surface of sphere
          real_t x = ((*normdist)(rngine));
          real_t y = ((*normdist)(rngine));
          real_t z = ((*normdist)(rngine));
          real_t r = vertices[i].param[0] / std::sqrt(x*x+y*y+z*z);
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
    rngstate.resize(modeldata[modidx].stateinit.size());
    rngstick.resize(modeldata[modidx].stickinit.size());
    // Randomly generate state
    for (std::size_t s = 0; s < modeldata[modidx].stateinit.size(); ++s) {
      if (modeldata[modidx].stateinit[s] == RNGTYPE_CONST) {
        rngstate[s] = rngconst(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].stateinit[s] == RNGTYPE_UNIF) {
        rngstate[s] = rngunif(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].stateinit[s] == RNGTYPE_UNINT) {
        rngstate[s] = rngunint(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].stateinit[s] == RNGTYPE_NORM) {
        rngstate[s] = rngnorm(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].stateinit[s] == RNGTYPE_BNORM) {
        rngstate[s] = rngbnorm(modeldata[modidx].stateparam[s].data());
      }
      else if (modeldata[modidx].stateinit[s] == RNGTYPE_FILE) {
        rngstate[s] = rngfile(modeldata[modidx].stateparam[s].data(), 0, vtxordidx[i]);
      }
      else {
        CkPrintf("  error: stateinit %s is not valid for vertex\n", rngtype[modeldata[modidx].stateinit[s]].c_str());
        CkExit();
      }
    }
    // Randomly generate stick
    for (std::size_t s = 0; s < modeldata[modidx].stickinit.size(); ++s) {
      if (modeldata[modidx].stickinit[s] == RNGTYPE_CONST) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngconst(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].stickinit[s] == RNGTYPE_UNIF) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngunif(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].stickinit[s] == RNGTYPE_UNINT) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngunint(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].stickinit[s] == RNGTYPE_NORM) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngnorm(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].stickinit[s] == RNGTYPE_BNORM) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngbnorm(modeldata[modidx].stickparam[s].data()));
      }
      else if (modeldata[modidx].stickinit[s] == RNGTYPE_FILE) {
        rngstick[s] = (tick_t)(TICKS_PER_MS * rngfile(modeldata[modidx].stickparam[s].data(), 0, vtxordidx[i]));
      }
      else {
        CkPrintf("  error: stateinit %s is not valid for vertex\n", rngtype[modeldata[modidx].stickinit[s]].c_str());
        CkExit();
      }
    }
    // Add to state
    state[i].push_back(rngstate);
    stick[i].push_back(rngstick);
    // Empty events
    event[i].clear();
  }

  // Any index-based sample connectivity occurs first
  // for each vertex
  std::set<idx_t>::iterator jedg;
  for (idx_t i = 0; i < norderdat; ++i) {
    idx_t glbtargetidx = xorderdat[datidx] + i;
    for (jedg = connsampleset[vtxmodidx[i]].begin(); jedg != connsampleset[vtxmodidx[i]].end(); ++jedg) {
      idx_t edg = *jedg;
      for (std::size_t k = 0; k < edges[edg].conntype.size(); ++k) {
        // Uniform weights is easier to compute
        if (edges[edg].conntype[k] == CONNTYPE_SMPL) {
          // Compute samples w.r.t. source population order
          CkAssert(vertices[edges[edg].source-1].order == edges[edg].maskparam[k][0]);
          std::vector<idx_t> sourceordix(edges[edg].maskparam[k][0]);
          std::iota(sourceordix.begin(), sourceordix.end(), 0);
          // pick the seed based on the targetidx so it is consistent across cores
          // The 32768 is 2^15, is just a reasonably large number to not get repeating seeds
          unsigned sampleseed = (randseed + (unsigned)(vtxordidx[i])) ^ ((unsigned)(edg*32768));
          std::shuffle(sourceordix.begin(), sourceordix.end(), std::mt19937{sampleseed});
          // want to make sure sample number is less than source order
          CkAssert(edges[edg].maskparam[k][0] >= edges[edg].maskparam[k][1]);
          // copy over the shuffled indices for the sampling
          for (idx_t j = 0; j < edges[edg].maskparam[k][1]; ++j) {
            // Convert from population index to global index
            std::vector<idx_t>::iterator iprt;
            iprt = std::upper_bound(xpopvtxidxprt[edges[edg].source-1].begin(), xpopvtxidxprt[edges[edg].source-1].end(), sourceordix[j]);
            int prt = (iprt - xpopvtxidxprt[edges[edg].source-1].begin()) - 1;
            idx_t glbsourceidx = xglbvtxidxprt[edges[edg].source-1][prt] + (sourceordix[j] - xpopvtxidxprt[edges[edg].source-1][prt]);
            // Check for self connections
            if (glbsourceidx == glbtargetidx) {
              continue;
            } else {
              adjcy[i].push_back(glbsourceidx);
              adjcyset[i].insert(glbsourceidx); // The set is useful for faster searching of edge existence
              edgmodidx[i].push_back(edges[edg].modidx);
              // The state/stick will need to be reparameterized with correct distance information later
              state[i].push_back(BuildEdgState(edges[edg].modidx, 0.0, sourceordix[j], vtxordidx[i]));
              stick[i].push_back(BuildEdgStick(edges[edg].modidx, 0.0, sourceordix[j], vtxordidx[i]));
            }
          }
        }
        // non-uniform weights need some additional helper structures
        else if (edges[edg].conntype[k] == CONNTYPE_SMPL_NORM ||
                 edges[edg].conntype[k] == CONNTYPE_SMPL_ANORM) {
          // Compute samples w.r.t. source population order
          CkAssert(vertices[edges[edg].source-1].order == edges[edg].maskparam[k][0]);
          std::vector<real_t> samplewgt(edges[edg].maskparam[k][0]);
          if (edges[edg].conntype[k] == CONNTYPE_SMPL_NORM) {
            for (idx_t j = 0; j < edges[edg].maskparam[k][0]; ++j) {
              // (x_i - x_j)^2 / var(x_ij)
              real_t x_ij = ((((real_t) vtxordidx[i])*vertices[vtxmodidx[i]-1].param[0]/vertices[vtxmodidx[i]-1].order)
                  - (((real_t) j)*vertices[edges[edg].source-1].param[0]/vertices[edges[edg].source-1].order));
              samplewgt[j] = std::exp(-(x_ij*x_ij)/(2*edges[edg].probparam[k][0])); // Don't worry about normalizing
            }
          }
          else if (edges[edg].conntype[k] == CONNTYPE_SMPL_ANORM) {
            for (idx_t j = 0; j < edges[edg].maskparam[k][0]; ++j) {
              // ((x_i - x_j)^2 / var(x_ij)) - ((x_i - x_j)^2 / var(x_ij)/2)
              // Not using variance-y for this for now (needs additional information)
              real_t x_ij = ((((real_t) vtxordidx[i])*vertices[vtxmodidx[i]-1].param[0]/vertices[vtxmodidx[i]-1].order)
                  - (((real_t) j)*vertices[edges[edg].source-1].param[0]/vertices[edges[edg].source-1].order));
              real_t var = edges[edg].probparam[k][0];
              // Normalizing a bit more important here
              real_t wgt = std::exp(-(x_ij*x_ij)/(2*var))/std::sqrt(var) *
                (std::exp(-(x_ij*x_ij)/(2*var))/std::sqrt(var) - std::exp(-(x_ij*x_ij)/(var))/std::sqrt(var/2));
              samplewgt[j] = std::max(0.0, wgt);
            }
          }
          // pick the seed based on the targetidx so it is consistent across cores
          // The 32768 is 2^15, is just a reasonably large number to not get repeating seeds
          unsigned sampleseed = (randseed + (unsigned)(vtxordidx[i])) ^ ((unsigned)(edg*32768));
          std::mt19937 rngsample(sampleseed);
          std::uniform_real_distribution<real_t> sampleunifdist;
          // want to make sure sample number is less than source order
          CkAssert(edges[edg].maskparam[k][0] >= edges[edg].maskparam[k][1]);
          // From: https://stackoverflow.com/questions/53632441/c-sampling-from-discrete-distribution-without-replacement
          std::vector<real_t> sampleprb;
          std::vector<real_t>::iterator iwgt;
          for (iwgt = samplewgt.begin(); iwgt != samplewgt.end(); ++iwgt) {
            sampleprb.push_back(std::pow(sampleunifdist(rngsample), 1.0 / (*iwgt)));
          }
          // Sorting vals, but retain the indices. 
          // There is unfortunately no easy way to do this with STL.
          std::vector<std::pair<idx_t, real_t>> sampleprbidx;
          for (std::size_t iter = 0; iter < sampleprb.size(); iter++) {
            sampleprbidx.emplace_back(iter, sampleprb[iter]);
          }
          std::sort(sampleprbidx.begin(), sampleprbidx.end(), [](std::pair<idx_t,real_t> x, std::pair<idx_t,real_t> y) {return x.second > y.second; });
          // Place samples into adjcy after shuffling
          for (idx_t iter = 0; iter < edges[edg].maskparam[k][1]; iter++) {
            idx_t sourceordix = sampleprbidx[iter].first;
            // Convert from population index to global index
            std::vector<idx_t>::iterator iprt;
            iprt = std::upper_bound(xpopvtxidxprt[edges[edg].source-1].begin(), xpopvtxidxprt[edges[edg].source-1].end(), sourceordix);
            int prt = (iprt - xpopvtxidxprt[edges[edg].source-1].begin()) - 1;
            idx_t glbsourceidx = xglbvtxidxprt[edges[edg].source-1][prt] + (sourceordix - xpopvtxidxprt[edges[edg].source-1][prt]);
            if (glbsourceidx == glbtargetidx) {
              continue;
            } else {
              adjcy[i].push_back(glbsourceidx);
              adjcyset[i].insert(glbsourceidx); // The set is useful for faster searching of edge existence
              edgmodidx[i].push_back(edges[edg].modidx);
              // The state/stick will need to be reparameterized with correct distance information later
              state[i].push_back(BuildEdgState(edges[edg].modidx, 0.0, sourceordix, vtxordidx[i]));
              stick[i].push_back(BuildEdgStick(edges[edg].modidx, 0.0, sourceordix, vtxordidx[i]));
            }
          }
        }
      }
    }
    // Count of adjcy from just index-based connections
    nadjcysample[i] = adjcy[i].size();
  }
  
  // Prepare for connection
  cpdat = 0;
  cphnd = 0;

  // Sample-based connectivity w.r.t. the source (without distance information) is done
  // TODO: still need sample-based w.r.t. the target
  // Connect to this part
  if (cpdat == datidx) {
    mConn *mconn = BuildConnVtx(cpdat);
    thisProxy(cpdat).ConnectVtx(mconn);
  }
  // Request data from remote part
  else {
    thisProxy(cpdat).RequestConnVtx(datidx);
  }
}


/**************************************************************************
* Connect
**************************************************************************/

// Connect Network: distance-based connections
//
void Netdata::ConnectVtx(mConn *msg) {
  // Sanity check
  CkAssert(msg->datidx == cpdat);
  // Some basic information on what's being connected
  //CkPrintf("  Connecting %d to %d\n", datidx, msg->datidx);
  
  // Go through the incoming state and reparameterize the state/sticks if needed
  for (idx_t i = 0; i < norderdat; ++i) {
    for (std::size_t j = 0; j < nadjcysample[i]; ++j) {
      if (xorderdat[msg->datidx] <= adjcy[i][j] && adjcy[i][j] < xorderdat[msg->datidx+1]) {
        // Reparameterize with information
        idx_t locsourceidx = adjcy[i][j] - xorderdat[msg->datidx];
        idx_t edg = connmodmap[msg->vtxmodidx[locsourceidx]*edges.size()+vtxmodidx[i]];
        //CkPrintf("%" PRIidx ", %" PRIidx ": %" PRIidx ", %" PRIidx"\n", i, adjcy[i][j], locsourceidx, edg);
        real_t distance = distfunc(edges[edg].distype, xyz.data()+i*3, msg->xyz+locsourceidx*3, edges[edg].distparam.data());
        ReBuildEdgState(edgmodidx[i][j], distance, state[i][j+1]);
        ReBuildEdgStick(edgmodidx[i][j], distance, stick[i][j+1]);
      }
    }
  }

  // Build connections from j to i (only)
  //
  for (idx_t i = 0; i < norderdat; ++i) {
    for (idx_t j = 0; j < msg->nvtx; ++j) {
      // Skip unconnected populations
      if (connmodmap.find(msg->vtxmodidx[j]*edges.size()+vtxmodidx[i]) == connmodmap.end()) {
        continue;
      }
      // Skip same index i == j
      else if(xorderdat[datidx]+i == xorderdat[msg->datidx]+j) {
        continue;
      }
      // Evaluate connection between i and j
      else {
        idx_t edg = connmodmap[msg->vtxmodidx[j]*edges.size()+vtxmodidx[i]];
        real_t distance = distfunc(edges[edg].distype, xyz.data()+i*3, msg->xyz+j*3, edges[edg].distparam.data());
        CkAssert(distance >= 0.0);
        idx_t modidx;
        // check possible connections from j (source) to i (target)
        if ((modidx = MakeConnection(edg, msg->vtxordidx[j], vtxordidx[i], distance))) {
          adjcy[i].push_back(xorderdat[msg->datidx]+j);
          adjcyset[i].insert(xorderdat[msg->datidx]+j); // The set is useful for faster searching of edge existence
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

  // send any outstanding requests of built vertieces
  for (std::list<idx_t>::iterator ireqidx = connvtxreq.begin(); ireqidx != connvtxreq.end(); ++ireqidx) {
    mConn *mconn = BuildConnVtx(*ireqidx);
    thisProxy(*ireqidx).ConnectVtx(mconn);
    ireqidx = connvtxreq.erase(ireqidx);
  }

  // Move to next part
  ++cpdat;
  // return control to main when done
  if (cpdat == netfiles) {
    // Go to checkpoint now
    thisProxy.ConnectHandover();
  }
  // Connect to curr part
  else if (cpdat == datidx) {
    mConn *mconn = BuildConnVtx(datidx);
    thisProxy(cpdat).ConnectVtx(mconn);
  }
  // Request data from next part
  else {
    thisProxy(cpdat).RequestConnVtx(datidx);
  }

  // cleanup
  delete msg;
}

// Connect Network Request Data
//
void Netdata::RequestConnVtx(idx_t reqidx) {
  // send data adjacency and vertex info to requesting part
  // check if vertex is built
  if (vtxmodidx.size()) {
    mConn *mconn = BuildConnVtx(reqidx);
    thisProxy(reqidx).ConnectVtx(mconn);
  }
  else {
    // record request for when adjcy is built
    connvtxreq.push_back(reqidx);
  }
}

// Handover between sample and distance-based connectivity
//
void Netdata::ConnectHandover() {
  ++cphnd;
  // return control to main when done
  if (cphnd == netfiles) {
    cpdat = 0;
    cphnd = 0; // maybe use again?
    //CkPrintf("Build Handover\n");
    /*
    idx_t nadjcy = 0;
    for (int i = 0; i < norderdat; ++i) {
      nadjcy += adjcyset[i].size();
    }
    CkPrintf("File: %d, directed edges: %" PRIidx "\n", datidx, nadjcy);
    */
    
    if (cpdat == datidx) {
      mConn *mconn = BuildConnEdg(datidx);
      thisProxy(cpdat).ConnectEdg(mconn);
    }
    // Work on connecting none edges now
    else {
      thisProxy(cpdat).RequestConnEdg(datidx);
    }
  }
}

// Connect Network (directed to undirected edges)
//
void Netdata::ConnectEdg(mConn *msg) {
  // Sanity check
  CkAssert(msg->datidx == cpdat);
  // Some basic information on what's being connected
  //CkPrintf("  Connecting %d to %d\n", datidx, msg->datidx);
  
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

  delete msg;

  // Move to next part
  ++cpdat;
  // return control to main when done
  if (cpdat == netfiles) {
    // We can reorder all the edges by global ordering now
    // Reorder edges vertex-by-vertex
    for (idx_t i = 0; i < norderdat; ++i) {
      edgreord.clear();
      for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
        edgreord.push_back(edgreord_t());
        edgreord.back().edgidx = adjcy[i][j];
        edgreord.back().modidx = edgmodidx[i][j];
        edgreord.back().state = state[i][j+1];
        edgreord.back().stick = stick[i][j+1];
      }
      // sort edge indices by global ordering
      std::sort(edgreord.begin(), edgreord.end());
      // add indices to data structures
      for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
        adjcy[i][j] = edgreord[j].edgidx;
        edgmodidx[i][j] = edgreord[j].modidx;
        state[i][j+1] = edgreord[j].state;
        stick[i][j+1] = edgreord[j].stick;
      }
    }
    /*
    // Print memory allocated
    int adjcysize = 0;
    int adjcycap = 0;
    int edgmodsize = 0;
    int edgmodcap = 0;
    for (size_t i = 0; i < adjcy.size(); ++i) {
      //adjcy[i].shrink_to_fit();
      //edgmodidx[i].shrink_to_fit();
      adjcysize += adjcy[i].size();
      adjcycap += adjcy[i].capacity();
      edgmodsize += edgmodidx[i].size();
      edgmodcap += edgmodidx[i].capacity();
    }
    CkPrintf("Part %d size/cap: adjcy: %d , %d edgmodidx: %d , %d\n", datidx, adjcysize, adjcycap, edgmodsize, edgmodcap);
    */

    // Build parts from file-based network information
    BuildParts();
    
    // Done building all edges, return control to main
    contribute(0, NULL, CkReduction::nop);
  }
  // Connect to curr part
  else if (cpdat == datidx) {
    mConn *mconn = BuildConnEdg(datidx);
    thisProxy(cpdat).ConnectEdg(mconn);
  }
  // Request data from next part
  else {
    thisProxy(cpdat).RequestConnEdg(datidx);
  }
}

// Connect Network Finished
//
void Netdata::RequestConnEdg(idx_t reqidx) {
  // send data adjacency info to requesting part
  mConn *mconn = BuildConnEdg(reqidx);
  thisProxy(reqidx).ConnectEdg(mconn);
}


/**************************************************************************
* Connection messages
**************************************************************************/


// Build Vertices
//
mConn* Netdata::BuildConnVtx(idx_t reqidx) {
  // Initialize connection message
  int msgSize[MSG_Conn];
  msgSize[0] = norderdat;   // vtxmodidx
  msgSize[1] = norderdat;   // vtxordidx
  msgSize[2] = norderdat*3; // xyz
  msgSize[3] = 0;           // vtxidx
  msgSize[4] = 0;           // xadj
  msgSize[5] = 0;           // adjcy
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

// Build Edge (includes vertices and adjacency)
//
mConn* Netdata::BuildConnEdg(idx_t reqidx) {
  /* Bookkeeping */
  idx_t nsizedat;
  idx_t jadjcyidx;

  std::vector<std::vector<idx_t>> adjcyconn;
  adjcyconn.resize(norderdat);

  // Count the sizes
  nsizedat = 0;
  for (idx_t i = 0; i < norderdat; ++i) {
    adjcyconn[i].clear();
    std::set<idx_t>::iterator jadjcy;
    // Loop through the directed afferent edges
    for (jadjcy = adjcyset[i].begin(); jadjcy != adjcyset[i].end(); ++jadjcy) {
      if (xorderdat[reqidx] <= *jadjcy && *jadjcy < xorderdat[reqidx+1]) {
        // global source idx to local
        adjcyconn[i].push_back(*jadjcy - xorderdat[reqidx]);
      }
    }
    // Add none connections to size
    nsizedat += adjcyconn[i].size();
  }
  //CkPrintf("   reqidx: %" PRIidx ", min: %" PRIidx ", max: %" PRIidx ", adj: %" PRIidx "\n", reqidx, xorderdat[reqidx], xorderdat[reqidx+1], nsizedat);

  // Initialize connection message
  int msgSize[MSG_Conn];
  msgSize[0] = 0;           // vtxmodidx
  msgSize[1] = 0;           // vtxordidx
  msgSize[2] = 0;           // xyz
  msgSize[3] = norderdat;   // vtxidx
  msgSize[4] = norderdat+1; // xadj
  msgSize[5] = nsizedat;    // adjcy
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
    mconn->vtxidx[i] = xorderdat[datidx]+i; // need to convert from local to global
    // xadj
    mconn->xadj[i+1] = mconn->xadj[i] + adjcyconn[i].size();
    for (std::size_t j = 0; j < adjcyconn[i].size(); ++j) {
      // adjcy (of next parts)
      mconn->adjcy[jadjcyidx++] = adjcyconn[i][j];
    }
  }
  CkAssert(jadjcyidx == nsizedat);

  return mconn;
}


/**************************************************************************
* Connection Construction
**************************************************************************/

idx_t Netdata::MakeConnection(idx_t edg, idx_t sourceidx, idx_t targetidx, real_t dist) {
  // test for cutoff
  if (edges[edg].cutoff != 0.0 && dist > edges[edg].cutoff) {
    return 0;
  }
  // Connection computation
  real_t prob = 0.0;
  idx_t mask = 0;
  for (std::size_t k = 0; k < edges[edg].conntype.size(); ++k) {
    if (edges[edg].conntype[k] == CONNTYPE_UNIF) {
      prob += edges[edg].probparam[k][0];
    }
    else if (edges[edg].conntype[k] == CONNTYPE_SIG) {
      prob += sigmoid(dist, edges[edg].probparam[k][0],
          edges[edg].probparam[k][1], edges[edg].probparam[k][2]);
    }
    else if (edges[edg].conntype[k] == CONNTYPE_IDX) {
      mask += (((sourceidx * edges[edg].maskparam[k][2]) + edges[edg].maskparam[k][3]) == targetidx);
    }
    // TODO: we want to enable a weighting of the indices w.r.t. distance
    else if (edges[edg].conntype[k] == CONNTYPE_SMPL ||
        edges[edg].conntype[k] == CONNTYPE_SMPL_NORM ||
        edges[edg].conntype[k] == CONNTYPE_SMPL_ANORM) {
      // sample-based edges already computed
      return 0;
    }
    else if (edges[edg].conntype[k] == CONNTYPE_FILE) {
      // Check to see if it's in the file list
      // Dimensions are stored: targetdim x sourcedim
      // set mask to 1 if there is a non-zero entry
      // TODO: make sure file-based connections completely override
      //       other connection types (or make them mutually exclusive)
      // TODO: split up files in dCSR order as well
      if (targetidx >= datafiles[(idx_t) (edges[edg].probparam[k][0])].matrix.size()) {
        CkPrintf("  error: datafile %s does not have row for %" PRIidx "\n",
            datafiles[(idx_t) (edges[edg].probparam[k][0])].filename.c_str(), targetidx);
      } else if (datafiles[(idx_t) (edges[edg].probparam[k][0])].matrix[targetidx].find((real_t) sourceidx) ==
          datafiles[(idx_t) (edges[edg].probparam[k][0])].matrix[targetidx].end()) {
        prob = 0.0;
        mask = 0;
      } else {
        mask = 1;
      }
    }
    else {
      // Shouldn't reach here due to prior error checking
      CkPrintf("  error: connection type %" PRIidx " undefined\n", edges[edg].conntype[k]);
      CkExit();
    }
  }
  // Compute probability of connection
  if ((((*unifdist)(rngine)) < prob) || mask) {
    return edges[edg].modidx;
  }
  else {
    return 0;
  }
  return 0; // no connection found, return 'none'
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
  rngstate.resize(modeldata[modidx].stateinit.size());
  // Randomly generate state
  for (std::size_t j = 0; j < rngstate.size(); ++j) {
    if (modeldata[modidx].stateinit[j] == RNGTYPE_CONST) {
      rngstate[j] = rngconst(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_UNIF) {
      rngstate[j] = rngunif(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_UNINT) {
      rngstate[j] = rngunint(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_NORM) {
      rngstate[j] = rngnorm(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_BNORM) {
      rngstate[j] = rngbnorm(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_LBNORM) {
      rngstate[j] = rnglbnorm(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_LBLOGNORM) {
      rngstate[j] = rnglblognorm(modeldata[modidx].stateparam[j].data());
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_LIN) {
      rngstate[j] = rnglin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_LBLIN) {
      rngstate[j] = rnglblin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_BLIN) {
      rngstate[j] = rngblin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_FILE) {
      rngstate[j] = rngfile(modeldata[modidx].stateparam[j].data(), sourceidx, targetidx);
    }
    else {
      CkPrintf("  error: stateinit %s is not valid for edge\n", rngtype[modeldata[modidx].stateinit[j]].c_str());
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
  rngstick.resize(modeldata[modidx].stickinit.size());
  // Randomly generate stick
  for (std::size_t j = 0; j < rngstick.size(); ++j) {
    if (modeldata[modidx].stickinit[j] == RNGTYPE_CONST) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngconst(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_UNIF) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngunif(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_UNINT) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngunint(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_NORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngnorm(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_BNORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngbnorm(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_LBNORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglbnorm(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_LBLOGNORM) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglblognorm(modeldata[modidx].stickparam[j].data()));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_LIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_LBLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglblin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_BLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngblin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_FILE) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngfile(modeldata[modidx].stickparam[j].data(), sourceidx, targetidx));
    }
    else {
      CkPrintf("  error: stateinit %s is not valid for edge\n", rngtype[modeldata[modidx].stickinit[j]].c_str());
      CkExit();
    }
  }
  // return generated stick
  return rngstick;
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
    if (modeldata[modidx].stateinit[j] == RNGTYPE_LIN) {
      rngstate[j] = rnglin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_LBLIN) {
      rngstate[j] = rnglblin(modeldata[modidx].stateparam[j].data(), dist);
    }
    else if (modeldata[modidx].stateinit[j] == RNGTYPE_BLIN) {
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
    if (modeldata[modidx].stickinit[j] == RNGTYPE_LIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_LBLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rnglblin(modeldata[modidx].stickparam[j].data(), dist));
    }
    else if (modeldata[modidx].stickinit[j] == RNGTYPE_BLIN) {
      rngstick[j] = (tick_t)(TICKS_PER_MS * rngblin(modeldata[modidx].stickparam[j].data(), dist));
    }
  }
}

