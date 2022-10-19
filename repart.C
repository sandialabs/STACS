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
extern /*readonly*/ int netparts;
extern /*readonly*/ int netfiles;


/**************************************************************************
* Repartitioning
**************************************************************************/

// Load data from files into partitions
//
void Netdata::LoadPart(mDist *msg) {
  // Copy over metis distribution from dist
  vtxmetis.resize(netfiles+1);
  edgmetis.resize(netfiles+1);
  int ndiv = netparts/netfiles;
  int nrem = netparts%netfiles;
  for (int i = 0; i < netfiles ; ++i) {
    int j = i*ndiv + (i < nrem ? i : nrem);
    vtxmetis[i] = msg->vtxdist[j];
    edgmetis[i] = msg->edgdist[j];
  }
  vtxmetis[netfiles] = msg->vtxdist[netfiles];
  edgmetis[netfiles] = msg->edgdist[netfiles];
  // cleanup
  delete msg;

  // Read in files
  CkPrintf("Reading network data files %" PRIidx "\n", datidx);
  ReadPart();
  
  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// Load data from network into partitions
//
void Netdata::LoadRepart(mPart *msg) {
  // Stash part
  parts[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++cprt == nprt) {
    cprt = 0;

    // Convert parts to local data
    BuildRepart();

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete parts[i];
    }

    // Return control to main
    contribute(0, NULL, CkReduction::nop);
  }
}

// Scatter Partitions across Network
//
void Netdata::ScatterPart() {
  // Compute which part goes to which data
  for (int prtdatidx = 0; prtdatidx < netfiles; ++prtdatidx) {
    // Which parts on this data
    int ndiv = netparts/netfiles;
    int nrem = netparts%netfiles;
    int nprtidx = ndiv + (prtdatidx < nrem);
    int xprtidx = prtdatidx*ndiv + (prtdatidx < nrem ? prtdatidx : nrem);

    // Loop through parts (global prtidx)
    for (int prtidx = xprtidx; prtidx < xprtidx + nprtidx; ++prtidx) {
      // Count sizes
      idx_t nedgidx = 0;
      idx_t nstate = 0;
      idx_t nstick = 0;
      idx_t nevent = 0;
      for (std::size_t i = 0; i < vtxidxreprt[prtidx].size(); ++i) {
        nedgidx += adjcyreprt[prtidx][i].size();
      }
      for (std::size_t i = 0; i < statereprt[prtidx].size(); ++i) {
        nstate += statereprt[prtidx][i].size();
      }
      for (std::size_t i = 0; i < stickreprt[prtidx].size(); ++i) {
        nstick += stickreprt[prtidx][i].size();
      }
      for (std::size_t i = 0; i < eventreprt[prtidx].size(); ++i) {
        nevent += eventreprt[prtidx][i].size();
      }

      // Initialize connection message
      int msgSize[MSG_Part];
      msgSize[0] = vtxidxreprt[prtidx].size();     // vtxidx
      msgSize[1] = vtxidxreprt[prtidx].size();     // vtxmodidx
      msgSize[2] = vtxidxreprt[prtidx].size() * 3; // xyz
      msgSize[3] = vtxidxreprt[prtidx].size() + 1; // xadj
      msgSize[4] = nedgidx;                       // adjcy
      msgSize[5] = nedgidx;                       // edgmodidx
      msgSize[6] = nstate;                        // state
      msgSize[7] = nstick;                        // stick
      msgSize[8] = vtxidxreprt[prtidx].size() + 1; // xevent
      msgSize[9] = nevent;                        // diffuse
      msgSize[10] = nevent;                       // type
      msgSize[11] = nevent;                       // source
      msgSize[12] = nevent;                       // index
      msgSize[13] = nevent;                       // data
      mPart *mpart = new(msgSize, 0) mPart;
      // Sizes
      mpart->nvtx = vtxidxreprt[prtidx].size();
      mpart->nedg = nedgidx; // currently unused by gather
      mpart->nstate = nstate;
      mpart->nstick = nstick;
      mpart->nevent = nevent;
      mpart->prtidx = prtidx;

      // set up counters
      idx_t jedgidx = 0;
      idx_t jstate = 0;
      idx_t jstick = 0;
      idx_t jevent = 0;
      // prefixes start at zero
      mpart->xadj[0] = 0;
      mpart->xevent[0] = 0;

      for (std::size_t i = 0; i < vtxidxreprt[prtidx].size(); ++i) {
        // vtxidx (as vtxdist)
        mpart->vtxdist[i] = vtxidxreprt[prtidx][i];
        // vtxmodidx
        mpart->vtxmodidx[i] = vtxmodidxreprt[prtidx][i];
        // xyz
        mpart->xyz[i*3+0] = xyzreprt[prtidx][i*3+0];
        mpart->xyz[i*3+1] = xyzreprt[prtidx][i*3+1];
        mpart->xyz[i*3+2] = xyzreprt[prtidx][i*3+2];
        // xadj
        mpart->xadj[i+1] = mpart->xadj[i] + adjcyreprt[prtidx][i].size();
        for (std::size_t j = 0; j < adjcyreprt[prtidx][i].size(); ++j) {
          // adjncy
          mpart->adjcy[jedgidx] = adjcyreprt[prtidx][i][j];
          // edgmodidx
          mpart->edgmodidx[jedgidx++] = edgmodidxreprt[prtidx][i][j];
        }
        // xevent
        mpart->xevent[i+1] = mpart->xevent[i] + eventreprt[prtidx][i].size();
        for (std::size_t j = 0; j < eventreprt[prtidx][i].size(); ++j) {
          //event
          mpart->diffuse[jevent] = eventreprt[prtidx][i][j].diffuse;
          mpart->type[jevent] = eventreprt[prtidx][i][j].type;
          mpart->source[jevent] = eventreprt[prtidx][i][j].source;
          mpart->index[jevent] = eventreprt[prtidx][i][j].index;
          mpart->data[jevent++] = eventreprt[prtidx][i][j].data;
        }
      }
      CkAssert(jedgidx == nedgidx);
      CkAssert(jevent == nevent);
      for (std::size_t i = 0; i < statereprt[prtidx].size(); ++i) {
        for (std::size_t s = 0; s < statereprt[prtidx][i].size(); ++s) {
          mpart->state[jstate++] = statereprt[prtidx][i][s];
        }
        for (std::size_t s = 0; s < stickreprt[prtidx][i].size(); ++s) {
          mpart->stick[jstick++] = stickreprt[prtidx][i][s];
        }
      }
      CkAssert(jstate == nstate);
      CkAssert(jstick == nstick);

      // Send part
      thisProxy(prtdatidx).GatherPart(mpart);
    }
  }
}


// Gather Partitions and perform reordering
//
void Netdata::GatherPart(mPart *msg) {
  // Bookkeeping
  int prtidx = msg->prtidx - xprt; // global to local
  idx_t jstate = 0;
  idx_t jstick = 0;
  idx_t jevent = 0;
  idx_t xvtx = vtxprted[prtidx].size();
  norderdat += msg->nvtx;
  norderprt[prtidx] += msg->nvtx;

  // allocate data for repartitioning
  vtxprted[prtidx].resize(norderprt[prtidx]);
  xyzprted[prtidx].resize(norderprt[prtidx]*3);
  adjcyprted[prtidx].resize(norderprt[prtidx]);
  edgmodidxprted[prtidx].resize(norderprt[prtidx]);
  stateprted[prtidx].resize(norderprt[prtidx]);
  stickprted[prtidx].resize(norderprt[prtidx]);
  eventprted[prtidx].resize(norderprt[prtidx]);
  // also allocate for reordering
  adjcyreord[prtidx].resize(norderprt[prtidx]);
  edgmodidxreord[prtidx].resize(norderprt[prtidx]);
  statereord[prtidx].resize(norderprt[prtidx]);
  stickreord[prtidx].resize(norderprt[prtidx]);

  // copy part data (unsorted)
  for (idx_t i = 0; i < msg->nvtx; ++i) {
    // vtxidx
    vtxprted[prtidx][xvtx+i].vtxidx = msg->vtxdist[i];
    // vtxmodidx
    vtxprted[prtidx][xvtx+i].modidx = msg->vtxmodidx[i];
    // localidx
    vtxprted[prtidx][xvtx+i].vtxidxloc = xvtx+i;
    // xyz
    xyzprted[prtidx][(xvtx+i)*3+0] = msg->xyz[i*3+0];
    xyzprted[prtidx][(xvtx+i)*3+1] = msg->xyz[i*3+1];
    xyzprted[prtidx][(xvtx+i)*3+2] = msg->xyz[i*3+2];
    // vertex state
    stateprted[prtidx][xvtx+i].push_back(std::vector<real_t>());
    stickprted[prtidx][xvtx+i].push_back(std::vector<tick_t>());
    CkAssert(vtxprted[prtidx][xvtx+i].modidx > 0);
    stateprted[prtidx][xvtx+i][0].resize(modeldata[vtxprted[prtidx][xvtx+i].modidx-1].nstate);
    stickprted[prtidx][xvtx+i][0].resize(modeldata[vtxprted[prtidx][xvtx+i].modidx-1].nstick);
    for(std::size_t s = 0; s < stateprted[prtidx][xvtx+i][0].size(); ++s) {
      stateprted[prtidx][xvtx+i][0][s] = msg->state[jstate++];
    }
    for(std::size_t s = 0; s < stickprted[prtidx][xvtx+i][0].size(); ++s) {
      stickprted[prtidx][xvtx+i][0][s] = msg->stick[jstick++];
    }

    // handle edges
    idx_t xedg = adjcyprted[prtidx][xvtx+i].size();
    adjcyprted[prtidx][xvtx+i].resize(xedg + msg->xadj[i+1] - msg->xadj[i]);
    edgmodidxprted[prtidx][xvtx+i].resize(xedg + msg->xadj[i+1] - msg->xadj[i]);
    for (idx_t j = 0; j < msg->xadj[i+1] - msg->xadj[i]; ++j) {
      // adjcy
      adjcyprted[prtidx][xvtx+i][xedg+j] = msg->adjcy[msg->xadj[i] + j];
      // edgmodidx
      edgmodidxprted[prtidx][xvtx+i][xedg+j] = msg->edgmodidx[msg->xadj[i] + j];
      // state
      stateprted[prtidx][xvtx+i].push_back(std::vector<real_t>());
      stickprted[prtidx][xvtx+i].push_back(std::vector<tick_t>());
      // only push edge state if model and not 'none'
      if (edgmodidxprted[prtidx][xvtx+i][xedg+j] > 0) {
        stateprted[prtidx][xvtx+i][xedg+j+1].resize(modeldata[edgmodidxprted[prtidx][xvtx+i][xedg+j]-1].nstate);
        stickprted[prtidx][xvtx+i][xedg+j+1].resize(modeldata[edgmodidxprted[prtidx][xvtx+i][xedg+j]-1].nstick);
        for(std::size_t s = 0; s < stateprted[prtidx][xvtx+i][xedg+j+1].size(); ++s) {
          stateprted[prtidx][xvtx+i][xedg+j+1][s] = msg->state[jstate++];
        }
        for(std::size_t s = 0; s < stickprted[prtidx][xvtx+i][xedg+j+1].size(); ++s) {
          stickprted[prtidx][xvtx+i][xedg+j+1][s] = msg->stick[jstick++];
        }
      }
    }

    // events
    eventprted[prtidx][xvtx+i].resize(msg->xevent[i+1] - msg->xevent[i]);
    for (idx_t e = 0; e < msg->xevent[i+1] - msg->xevent[i]; ++e) {
      eventprted[prtidx][xvtx+i][e].diffuse = msg->diffuse[jevent];
      eventprted[prtidx][xvtx+i][e].type = msg->type[jevent];
      eventprted[prtidx][xvtx+i][e].source = msg->source[jevent];
      eventprted[prtidx][xvtx+i][e].index = msg->index[jevent];
      eventprted[prtidx][xvtx+i][e].data = msg->data[jevent++];
    }
  }
  CkAssert(jstate == msg->nstate);
  CkAssert(jstick == msg->nstick);
  CkAssert(jevent == msg->nevent);

  // cleanup
  delete msg;

  // When all parts are gathered from all other data,
  // Perform reordering of vertex indices and start reordering
  if (++cpprt == netfiles*nprt) {
    cpprt = 0;
    // cleanup finished data structures
    vtxmetis.clear();
    edgmetis.clear();
    vtxidxreprt.clear();
    vtxmodidxreprt.clear();
    xyzreprt.clear();
    adjcyreprt.clear();
    edgmodidxreprt.clear();
    statereprt.clear();
    stickreprt.clear();
    eventreprt.clear();
    vtxidxreprt.resize(netparts);
    vtxmodidxreprt.resize(netparts);
    xyzreprt.resize(netparts);
    edgmodidxreprt.resize(netparts);
    adjcyreprt.resize(netparts);
    statereprt.resize(netparts);
    stickreprt.resize(netparts);
    eventreprt.resize(netparts);

    // collect order parts
    std::string orderprts;
    for (idx_t jprt = 0; jprt < nprt; ++jprt) {
      std::ostringstream orderprt;
      orderprt << " " << norderprt[jprt];
      orderprts.append(orderprt.str());
    }
    CkPrintf("  Reparted File: %d   Vertices: %" PRIidx " {%s }\n",
             datidx, norderdat, orderprts.c_str());

    // set up containers
    vtxmodidx.resize(norderdat);
    xyz.resize(norderdat*3);
    adjcy.resize(norderdat);
    edgmodidx.resize(norderdat);
    state.resize(norderdat);
    stick.resize(norderdat);
    event.resize(norderdat);
    eventsourcereord.resize(norderdat);
    eventindexreord.resize(norderdat);

    // Go through part data and reorder
    idx_t xvtx = 0;
    for (idx_t jprt = 0; jprt < nprt; ++jprt) {
      //CkPrintf("  Reordering Part %" PRIidx "\n", xprt+jprt);
      // reorder based on modidx
      std::sort(vtxprted[jprt].begin(), vtxprted[jprt].end());

      // add to data structures
      for (idx_t i = 0; i < norderprt[jprt]; ++i) {
        // vtxmodidx
        vtxmodidx[xvtx+i] = vtxprted[jprt][i].modidx;
        // xyz
        xyz[(xvtx+i)*3+0] = xyzprted[jprt][(vtxprted[jprt][i].vtxidxloc)*3+0];
        xyz[(xvtx+i)*3+1] = xyzprted[jprt][(vtxprted[jprt][i].vtxidxloc)*3+1];
        xyz[(xvtx+i)*3+2] = xyzprted[jprt][(vtxprted[jprt][i].vtxidxloc)*3+2];
        // adjcy
        adjcyreord[jprt][i] = adjcyprted[jprt][vtxprted[jprt][i].vtxidxloc];
        edgmodidxreord[jprt][i] = edgmodidxprted[jprt][vtxprted[jprt][i].vtxidxloc];
        // state (for vertex)
        state[xvtx+i].push_back(stateprted[jprt][vtxprted[jprt][i].vtxidxloc][0]);
        stick[xvtx+i].push_back(stickprted[jprt][vtxprted[jprt][i].vtxidxloc][0]);
        // rest of state
        statereord[jprt][i] = stateprted[jprt][vtxprted[jprt][i].vtxidxloc];
        stickreord[jprt][i] = stickprted[jprt][vtxprted[jprt][i].vtxidxloc];
        // events
        event[xvtx+i] = eventprted[jprt][vtxprted[jprt][i].vtxidxloc];
        eventsourcereord[xvtx+i].resize(event[xvtx+i].size());
        for (std::size_t e = 0; e < event[xvtx+i].size(); ++e) {
          eventsourcereord[xvtx+i][e] = event[xvtx+i][e].source;
        }
        eventindexreord[xvtx+i].push_back(0);
      }

      // increment xvtx
      xvtx += norderprt[jprt];
    }

    // Take care of any ordering that may have come in
    // Go through ordering list
    for (std::list<mReorder *>::iterator iordidx = reordlist.begin(); iordidx != reordlist.end(); ++iordidx) {
      if ((*iordidx)->datidx == cpdat) {
        // Perform reordering
        ReordEdg((*iordidx));
        // Move to next part
        ++cpdat;
        // Erase element from list
        iordidx = reordlist.erase(iordidx);
      }
    }

    // Broadcast reordlist in order of file
    if (cpdat == datidx) {
      mReorder *morder = BuildReorder();
      thisProxy.Reorder(morder);
    }
  }
}


// Handle ordering messages
//
void Netdata::Reorder(mReorder *msg) {
  // Save message for processing
  reordlist.push_back(msg);

  if (adjcy.size()) {
    // Go through ordering list
    for (std::list<mReorder *>::iterator iordidx = reordlist.begin(); iordidx != reordlist.end(); ++iordidx) {
      if ((*iordidx)->datidx == cpdat) {
        // Perform reordering
        ReordEdg((*iordidx));
        // Move to next part
        ++cpdat;
        // Erase element from list
        iordidx = reordlist.erase(iordidx);
      }
    }
  }

  // Check if done reordering
  if (cpdat == netfiles) {
    cpdat = 0;
    // reindex events
    for (idx_t i = 0; i < norderdat; ++i) {
      CkAssert(eventindexreord[i].size() == adjcy[i].size()+1);
      // create map
      std::unordered_map<idx_t, idx_t> oldtonew;
      for (std::size_t j = 0; j < eventindexreord[i].size(); ++j) {
        oldtonew[eventindexreord[i][j]] = j;
      }
      // modify index
      for (std::size_t j = 0; j < event[i].size(); ++j) {
        event[i][j].index = oldtonew[event[i][j].index];
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

    // Build Parts
    BuildParts();
    
    // return control to main when done
    contribute(0, NULL, CkReduction::nop);
  }
  // Broadcast ordering in order of file
  else if (cpdat == datidx) {
    mReorder *morder = BuildReorder();
    thisProxy.Reorder(morder);
  }
}

// Reordering indices based on given ordering
//
void Netdata::ReordEdg(mReorder *msg) {
  // Add to vtxmetis
  vtxmetis[cpdat+1] = vtxmetis[cpdat] + msg->nvtx;
  int ndiv = netparts/netfiles;
  int nrem = netparts%netfiles;
  int xdatprt = cpdat*ndiv + (cpdat < nrem ? cpdat : nrem);
  // add to vtxdist
  for (idx_t k = 0; k < msg->nprt; ++k) {
    vtxdist[xdatprt+k+1] = vtxdist[xdatprt+k] + msg->vtxdist[k];
  }

  // create map
  std::unordered_map<idx_t, idx_t> oldtonew;
  for (idx_t i = 0; i < msg->nvtx; ++i) {
    oldtonew[msg->vtxidxold[i]] = vtxmetis[cpdat] + msg->vtxidxnew[i];
  }

  // cleanup
  delete msg;

  // Reorder edges
  idx_t xvtx = 0;
  for (idx_t jprt = 0; jprt < nprt; ++jprt) {
    for (idx_t i = 0; i < norderprt[jprt]; ++i) {
      edgreord.clear();
      for (std::size_t j = 0; j < adjcyreord[jprt][i].size(); ++j) {
        if (oldtonew.find(adjcyreord[jprt][i][j]) == oldtonew.end()) {
          continue;
        }
        else {
          edgreord.push_back(edgreord_t());
          edgreord.back().edgidx = oldtonew[adjcyreord[jprt][i][j]];
          edgreord.back().modidx = edgmodidxreord[jprt][i][j];
          edgreord.back().state = statereord[jprt][i][j+1];
          edgreord.back().stick = stickreord[jprt][i][j+1];
          edgreord.back().evtidx = j+1;
          // resource events
          for (std::size_t e = 0; e < eventsourcereord[xvtx+i].size(); ++e) {
            if (eventsourcereord[xvtx+i][e] == adjcyreord[jprt][i][j]) {
              event[xvtx+i][e].source = oldtonew[eventsourcereord[xvtx+i][e]];
            }
          }
        }
      }
      // sort newly added indices
      std::sort(edgreord.begin(), edgreord.end());
      // add indices to data structures
      for (std::size_t j = 0; j < edgreord.size(); ++j) {
        adjcy[xvtx+i].push_back(edgreord[j].edgidx);
        edgmodidx[xvtx+i].push_back(edgreord[j].modidx);
        state[xvtx+i].push_back(edgreord[j].state);
        stick[xvtx+i].push_back(edgreord[j].stick);
        eventindexreord[xvtx+i].push_back(edgreord[j].evtidx);
      }
    }
    xvtx += norderprt[jprt];
  }
}


/**************************************************************************
* Reordering messages
**************************************************************************/

// Build Reorder (from old vertices to new indices)
//
mReorder* Netdata::BuildReorder() {
  // Initialize connection message
  int msgSize[MSG_Reorder];
  msgSize[0] = nprt;        // vtxdist
  msgSize[1] = norderdat;   // vtxidxold
  msgSize[2] = norderdat;   // vtxidxnew
  mReorder *morder = new(msgSize, 0) mReorder;
  // sizes
  morder->datidx = datidx;
  morder->nvtx = norderdat;
  morder->nprt = nprt;

  // set up counters
  idx_t jvtxidx = 0;

  // load data
  for (idx_t jprt = 0; jprt < nprt; ++jprt) {
    morder->vtxdist[jprt] = norderprt[jprt];
    for (idx_t i = 0; i < norderprt[jprt]; ++i) {
      morder->vtxidxold[jvtxidx] = vtxprted[jprt][i].vtxidx;
      morder->vtxidxnew[jvtxidx] = jvtxidx;
      ++jvtxidx;
    }
  }
  CkAssert(jvtxidx == norderdat);

  return morder;
}

// Build part message
//
mPart* Network::BuildRepart() {
  /* Bookkeeping */
  idx_t nadjcy;
  idx_t nstate;
  idx_t nstick;
  idx_t nevent;
  
  // Get total size of adjcy
  // TODO: for structural plasticity, keep track of vertices/edges that
  //       have been created or need to be purged (e.g. reprtidx = -1)
  nadjcy = 0;
  nstate = 0;
  nstick = 0;
  nevent = 0;
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    nadjcy += adjcy[i].size();
    for (std::size_t j = 0; j < adjcy[i].size()+1; ++j) {
      nstate += state[i][j].size();
      nstick += stick[i][j].size();
    }
    for (std::size_t j = 0; j < evtcal[i].size(); ++j) {
      nevent += evtcal[i][j].size();
    }
    nevent += evtcol[i].size();
  }

  // Initialize partition data message
  int msgSize[MSG_Part];
  msgSize[0] = adjcy.size();  // reprtidx (as vtxdist)
  msgSize[1] = adjcy.size();  // vtxmodidx
  msgSize[2] = adjcy.size()*3;// xyz
  msgSize[3] = adjcy.size()+1;// xadj
  msgSize[4] = nadjcy;        // adjcy
  msgSize[5] = nadjcy;        // edgmodidx
  msgSize[6] = nstate;        // state
  msgSize[7] = nstick;        // stick
  msgSize[8] = adjcy.size()+1;// xevent
  msgSize[9] = nevent;        // diffuse
  msgSize[10] = nevent;       // type
  msgSize[11] = nevent;       // source
  msgSize[12] = nevent;       // index
  msgSize[13] = nevent;       // data
  mPart *mpart = new(msgSize, 0) mPart;

  // Data sizes
  mpart->nvtx = adjcy.size();
  mpart->nedg = nadjcy;
  mpart->nstate = nstate;
  mpart->nstick = nstick;
  mpart->nevent = nevent;
  mpart->prtidx = prtidx;

  // Vertex and Edge Information
  idx_t jstate = 0;
  idx_t jstick = 0;
  idx_t jevent = 0;
  mpart->xadj[0] = 0;
  mpart->xevent[0] = 0;
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    // reprtidx (for now keep all vertices on same partition)
    mpart->vtxdist[i] = prtidx;
    // vtxmodidx
    mpart->vtxmodidx[i] = vtxmodidx[i];
    // xyz
    mpart->xyz[i*3+0] = xyz[i*3+0];
    mpart->xyz[i*3+1] = xyz[i*3+1];
    mpart->xyz[i*3+2] = xyz[i*3+2];
    // vertex state
    for (std::size_t s = 0; s < state[i][0].size(); ++s) {
      mpart->state[jstate++] = state[i][0][s];
    }
    for (std::size_t s = 0; s < stick[i][0].size(); ++s) {
      mpart->stick[jstick++] = stick[i][0][s];
    }

    // xadj
    mpart->xadj[i+1] = mpart->xadj[i] + adjcy[i].size();
    for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
      // adjcy
      mpart->adjcy[mpart->xadj[i] + j] = adjcy[i][j];
      // edgmodidx
      mpart->edgmodidx[mpart->xadj[i] + j] = edgmodidx[i][j];
      // state
      for (std::size_t s = 0; s < state[i][j+1].size(); ++s) {
        mpart->state[jstate++] = state[i][j+1][s];
      }
      for (std::size_t s = 0; s < stick[i][j+1].size(); ++s) {
        mpart->stick[jstick++] = stick[i][j+1][s];
      }
    }

    // events
    for (std::size_t j = 0; j < evtcal[i].size(); ++j) {
      for (std::size_t e = 0; e < evtcal[i][j].size(); ++e) {
        mpart->diffuse[jevent] = evtcal[i][j][e].diffuse;// - tsim;
        mpart->type[jevent] = evtcal[i][j][e].type;
        mpart->source[jevent] = evtcal[i][j][e].source;
        mpart->index[jevent] = evtcal[i][j][e].index;
        mpart->data[jevent++] = evtcal[i][j][e].data;
      }
    }
    // events spillover
    for (std::size_t e = 0; e < evtcol[i].size(); ++e) {
      mpart->diffuse[jevent] = evtcol[i][e].diffuse;// - tsim;
      mpart->type[jevent] = evtcol[i][e].type;
      mpart->source[jevent] = evtcol[i][e].source;
      mpart->index[jevent] = evtcol[i][e].index;
      mpart->data[jevent++] = evtcol[i][e].data;
    }
    // xevent
    mpart->xevent[i+1] = jevent;
  }
  CkAssert(jstate == nstate);
  CkAssert(jstick == nstick);
  CkAssert(jevent == nevent);

  return mpart;
}
