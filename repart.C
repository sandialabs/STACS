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
extern /*readonly*/ CProxy_Main mainProxy;
extern /*readonly*/ int netparts;
extern /*readonly*/ int netfiles;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ idx_t nevtday;


/**************************************************************************
* Repartitioning
**************************************************************************/

// Load data from files into partitions
//
void Netdata::LoadRepart() {
  // Read in files
  CkPrintf("  Reading repartitioning files %d\n", datidx);
  ReadNetpart();
  
  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// Coordination with NetData chare array
//
void Network::LoadRepart() {
  // Request network part from input
  netdata(datidx).LoadNetpart(prtidx,
      CkCallback(CkIndex_Network::LoadNetpart(NULL), thisProxy(prtidx)));
}

void Netdata::LoadNetpart(int prtidx, const CkCallback &cbpart) {
  // Initialize distribution message
  int msgSize[MSG_Repart];
  msgSize[0] = vtxdist[prtidx+1] - vtxdist[prtidx]; // reprtidx
  mRepart *mrepart = new(msgSize, 0) mRepart;
  mrepart->nvtx = vtxdist[prtidx+1] - vtxdist[prtidx];
  mrepart->prtidx = prtidx;
  // Some sanity checks (no more vtxmetis)
  CkAssert(prtidx >= xprt);
  CkAssert(reprtidx.size() == (vtxdist[xprt+nprt] - vtxdist[xprt]));

  // Get distribution info
  for (int i = 0; i < vtxdist[prtidx+1] - vtxdist[prtidx]; ++i) {
    // reprtidx
    mrepart->reprtidx[i] = reprtidx[i + vtxdist[prtidx] - vtxdist[xprt]];
  }

  // Send part to network
  cbpart.send(mrepart);
}

void Network::LoadNetpart(mRepart *msg) {
  // Copy over repart data
  CkAssert(msg->prtidx == prtidx);
  
  reprtidx.resize(msg->nvtx);
  for (std::size_t i = 0; i < reprtidx.size(); ++i) {
    reprtidx[i] = msg->reprtidx[i];
  }
  
  // Clear out partitioning
  vtxprted.clear();
  xyzprted.clear();
  adjcyprted.clear();
  edgmodidxprted.clear();
  stateprted.clear();
  stickprted.clear();
  eventprted.clear();
  // Prepare for reordering partitions
  norderprt = 0;
  cpprt = 0;
  cphnd = 0;

  // cleanup
  delete msg;
  
  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

void Network::RebalNetwork() {
  // Compute Repartitioning
  reprtidx.resize(adjcy.size());
  evtflat.clear();
  evtflat.resize(adjcy.size());
  for (std::size_t i = 0; i < adjcy.size(); ++i) {
    reprtidx[i] = prtidx;
  }
  
  // Clear out partitioning
  vtxprted.clear();
  xyzprted.clear();
  adjcyprted.clear();
  edgmodidxprted.clear();
  stateprted.clear();
  stickprted.clear();
  eventprted.clear();
  // Prepare for reordering partitions
  norderprt = 0;
  cpprt = 0;
  cphnd = 0;

  //thisProxy(prtidx).Repart();
  contribute(0, NULL, CkReduction::nop, maincheck);
}

void Network::Repart() {
  // Repartitioning
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
  adjcyreprt.resize(netparts);
  edgmodidxreprt.resize(netparts);
  statereprt.resize(netparts);
  stickreprt.resize(netparts);
  eventreprt.resize(netparts);

  for (std::size_t i = 0; i < reprtidx.size(); ++i) {
    CkAssert(reprtidx[i] < netparts);
    // Preallocation
    adjcyreprt[reprtidx[i]].push_back(std::vector<idx_t>());
    edgmodidxreprt[reprtidx[i]].push_back(std::vector<idx_t>());
    statereprt[reprtidx[i]].push_back(std::vector<real_t>());
    stickreprt[reprtidx[i]].push_back(std::vector<tick_t>());
    eventreprt[reprtidx[i]].push_back(std::vector<event_t>());
    // Vertices
    vtxidxreprt[reprtidx[i]].push_back(vtxdist[prtidx] + i);
    vtxmodidxreprt[reprtidx[i]].push_back(vtxmodidx[i]);
    xyzreprt[reprtidx[i]].push_back(xyz[i*3+0]);
    xyzreprt[reprtidx[i]].push_back(xyz[i*3+1]);
    xyzreprt[reprtidx[i]].push_back(xyz[i*3+2]);
    // Vertex state
    for (std::size_t s = 0; s < state[i][0].size(); ++s) {
      statereprt[reprtidx[i]].back().push_back(state[i][0][s]);
    }
    for (std::size_t s = 0; s < stick[i][0].size(); ++s) {
      stickreprt[reprtidx[i]].back().push_back(stick[i][0][s]);
    }
    // Edges
    for (std::size_t j = 0; j < adjcy[i].size(); ++j) {
      adjcyreprt[reprtidx[i]].back().push_back(adjcy[i][j]);
      edgmodidxreprt[reprtidx[i]].back().push_back(edgmodidx[i][j]);
      for (std::size_t s = 0; s < state[i][j+1].size(); ++s) {
        statereprt[reprtidx[i]].back().push_back(state[i][j+1][s]);
      }
      for (std::size_t s = 0; s < stick[i][j+1].size(); ++s) {
        stickreprt[reprtidx[i]].back().push_back(stick[i][j+1][s]);
      }
    }
    // Events
    for (std::size_t j = 0; j < evtcal[i].size(); ++j) {
      for (std::size_t e = 0; e < evtcal[i][j].size(); ++e) {
        eventreprt[reprtidx[i]].back().push_back(evtcal[i][j][e]);
      }
    }
    for (std::size_t e = 0; e < evtcol[i].size(); ++e) {
      eventreprt[reprtidx[i]].back().push_back(evtcol[i][e]);
    }
  }
  // Clear data once copied
  vtxmodidx.clear();
  xyz.clear();
  adjcy.clear();
  edgmodidx.clear();
  state.clear();
  stick.clear();
  evtcal.clear();
  evtcol.clear();
  
  thisProxy(prtidx).ScatterPart();
}

// Scatter Partitions across Network
//
void Network::ScatterPart() {
  // Compute which part goes to which data
  for (int prtprtidx = 0; prtprtidx < netparts; ++prtprtidx) {
    // Count sizes
    idx_t nedgidx = 0;
    idx_t nstate = 0;
    idx_t nstick = 0;
    idx_t nevent = 0;
    for (std::size_t i = 0; i < vtxidxreprt[prtprtidx].size(); ++i) {
      nedgidx += adjcyreprt[prtprtidx][i].size();
    }
    for (std::size_t i = 0; i < statereprt[prtprtidx].size(); ++i) {
      nstate += statereprt[prtprtidx][i].size();
    }
    for (std::size_t i = 0; i < stickreprt[prtprtidx].size(); ++i) {
      nstick += stickreprt[prtprtidx][i].size();
    }
    for (std::size_t i = 0; i < eventreprt[prtprtidx].size(); ++i) {
      nevent += eventreprt[prtprtidx][i].size();
    }

    // Initialize connection message
    int msgSize[MSG_Part];
    msgSize[0] = 0;                              // vtxdist (to be recomputed)
    msgSize[1] = vtxidxreprt[prtprtidx].size();     // vtxidx
    msgSize[2] = vtxidxreprt[prtprtidx].size();     // vtxmodidx
    msgSize[3] = vtxidxreprt[prtprtidx].size() * 3; // xyz
    msgSize[4] = vtxidxreprt[prtprtidx].size() + 1; // xadj
    msgSize[5] = nedgidx;                       // adjcy
    msgSize[6] = nedgidx;                       // edgmodidx
    msgSize[7] = nstate;                        // state
    msgSize[8] = nstick;                        // stick
    msgSize[9] = vtxidxreprt[prtprtidx].size() + 1; // xevent
    msgSize[10] = nevent;                        // diffuse
    msgSize[11] = nevent;                       // type
    msgSize[12] = nevent;                       // source
    msgSize[13] = nevent;                       // index
    msgSize[14] = nevent;                       // data
    mPart *mpart = new(msgSize, 0) mPart;
    // Sizes
    mpart->nvtx = vtxidxreprt[prtprtidx].size();
    mpart->nedg = nedgidx; // currently unused by gather
    mpart->nstate = nstate;
    mpart->nstick = nstick;
    mpart->nevent = nevent;
    mpart->prtidx = prtprtidx;

    // set up counters
    idx_t jedgidx = 0;
    idx_t jstate = 0;
    idx_t jstick = 0;
    idx_t jevent = 0;
    // prefixes start at zero
    mpart->xadj[0] = 0;
    mpart->xevent[0] = 0;

    for (std::size_t i = 0; i < vtxidxreprt[prtprtidx].size(); ++i) {
      // vtxidx
      mpart->vtxidx[i] = vtxidxreprt[prtprtidx][i];
      // vtxmodidx
      mpart->vtxmodidx[i] = vtxmodidxreprt[prtprtidx][i];
      // xyz
      mpart->xyz[i*3+0] = xyzreprt[prtprtidx][i*3+0];
      mpart->xyz[i*3+1] = xyzreprt[prtprtidx][i*3+1];
      mpart->xyz[i*3+2] = xyzreprt[prtprtidx][i*3+2];
      // xadj
      mpart->xadj[i+1] = mpart->xadj[i] + adjcyreprt[prtprtidx][i].size();
      for (std::size_t j = 0; j < adjcyreprt[prtprtidx][i].size(); ++j) {
        // adjncy
        mpart->adjcy[jedgidx] = adjcyreprt[prtprtidx][i][j];
        // edgmodidx
        mpart->edgmodidx[jedgidx++] = edgmodidxreprt[prtprtidx][i][j];
      }
      // xevent
      mpart->xevent[i+1] = mpart->xevent[i] + eventreprt[prtprtidx][i].size();
      for (std::size_t j = 0; j < eventreprt[prtprtidx][i].size(); ++j) {
        //event
        mpart->diffuse[jevent] = eventreprt[prtprtidx][i][j].diffuse;
        mpart->type[jevent] = eventreprt[prtprtidx][i][j].type;
        mpart->source[jevent] = eventreprt[prtprtidx][i][j].source;
        mpart->index[jevent] = eventreprt[prtprtidx][i][j].index;
        mpart->data[jevent++] = eventreprt[prtprtidx][i][j].data;
      }
    }
    CkAssert(jedgidx == nedgidx);
    CkAssert(jevent == nevent);
    for (std::size_t i = 0; i < statereprt[prtprtidx].size(); ++i) {
      for (std::size_t s = 0; s < statereprt[prtprtidx][i].size(); ++s) {
        mpart->state[jstate++] = statereprt[prtprtidx][i][s];
      }
      for (std::size_t s = 0; s < stickreprt[prtprtidx][i].size(); ++s) {
        mpart->stick[jstick++] = stickreprt[prtprtidx][i][s];
      }
    }
    CkAssert(jstate == nstate);
    CkAssert(jstick == nstick);

    // Send part
    thisProxy(prtprtidx).GatherPart(mpart);
  }
}


// Gather Partitions and perform reordering
//
void Network::GatherPart(mPart *msg) {
  // Bookkeeping
  //int prtprtidx = msg->prtidx - xprt; // global to local
  idx_t jstate = 0;
  idx_t jstick = 0;
  idx_t jevent = 0;
  idx_t xvtx = vtxprted.size(); // current size
  norderprt += msg->nvtx;

  // increase data allocated for repartitioning
  vtxprted.resize(norderprt);
  xyzprted.resize(norderprt*3);
  adjcyprted.resize(norderprt);
  edgmodidxprted.resize(norderprt);
  stateprted.resize(norderprt);
  stickprted.resize(norderprt);
  eventprted.resize(norderprt);

  // copy part data (unsorted)
  for (idx_t i = 0; i < msg->nvtx; ++i) {
    // vtxidx
    vtxprted[xvtx+i].vtxidx = msg->vtxidx[i];
    // vtxmodidx
    vtxprted[xvtx+i].modidx = msg->vtxmodidx[i];
    // localidx
    vtxprted[xvtx+i].vtxidxloc = xvtx+i;
    // xyz
    xyzprted[(xvtx+i)*3+0] = msg->xyz[i*3+0];
    xyzprted[(xvtx+i)*3+1] = msg->xyz[i*3+1];
    xyzprted[(xvtx+i)*3+2] = msg->xyz[i*3+2];
    // vertex state
    stateprted[xvtx+i].push_back(std::vector<real_t>());
    stickprted[xvtx+i].push_back(std::vector<tick_t>());
    CkAssert(vtxprted[xvtx+i].modidx > 0);
    stateprted[xvtx+i][0].resize(modelconf[vtxprted[xvtx+i].modidx-1].nstate);
    stickprted[xvtx+i][0].resize(modelconf[vtxprted[xvtx+i].modidx-1].nstick);
    for(std::size_t s = 0; s < stateprted[xvtx+i][0].size(); ++s) {
      stateprted[xvtx+i][0][s] = msg->state[jstate++];
    }
    for(std::size_t s = 0; s < stickprted[xvtx+i][0].size(); ++s) {
      stickprted[xvtx+i][0][s] = msg->stick[jstick++];
    }

    // handle edges
    idx_t xedg = adjcyprted[xvtx+i].size();
    adjcyprted[xvtx+i].resize(xedg + msg->xadj[i+1] - msg->xadj[i]);
    edgmodidxprted[xvtx+i].resize(xedg + msg->xadj[i+1] - msg->xadj[i]);
    for (idx_t j = 0; j < msg->xadj[i+1] - msg->xadj[i]; ++j) {
      // adjcy
      adjcyprted[xvtx+i][xedg+j] = msg->adjcy[msg->xadj[i] + j];
      // edgmodidx
      edgmodidxprted[xvtx+i][xedg+j] = msg->edgmodidx[msg->xadj[i] + j];
      // state
      stateprted[xvtx+i].push_back(std::vector<real_t>());
      stickprted[xvtx+i].push_back(std::vector<tick_t>());
      // only push edge state if model and not 'none'
      if (edgmodidxprted[xvtx+i][xedg+j] > 0) {
        stateprted[xvtx+i][xedg+j+1].resize(modelconf[edgmodidxprted[xvtx+i][xedg+j]-1].nstate);
        stickprted[xvtx+i][xedg+j+1].resize(modelconf[edgmodidxprted[xvtx+i][xedg+j]-1].nstick);
        for(std::size_t s = 0; s < stateprted[xvtx+i][xedg+j+1].size(); ++s) {
          stateprted[xvtx+i][xedg+j+1][s] = msg->state[jstate++];
        }
        for(std::size_t s = 0; s < stickprted[xvtx+i][xedg+j+1].size(); ++s) {
          stickprted[xvtx+i][xedg+j+1][s] = msg->stick[jstick++];
        }
      }
    }

    // events
    eventprted[xvtx+i].resize(msg->xevent[i+1] - msg->xevent[i]);
    for (idx_t e = 0; e < msg->xevent[i+1] - msg->xevent[i]; ++e) {
      eventprted[xvtx+i][e].diffuse = msg->diffuse[jevent];
      eventprted[xvtx+i][e].type = msg->type[jevent];
      eventprted[xvtx+i][e].source = msg->source[jevent];
      eventprted[xvtx+i][e].index = msg->index[jevent];
      eventprted[xvtx+i][e].data = msg->data[jevent++];
    }
  }
  CkAssert(jstate == msg->nstate);
  CkAssert(jstick == msg->nstick);
  CkAssert(jevent == msg->nevent);

  // cleanup
  delete msg;

  // When all parts are gathered from all other data,
  // Perform reordering of vertex indices and start reordering
  if (++cphnd == netparts) {
    cphnd = 0;
    // cleanup finished data structures
    vtxidxreprt.clear();
    vtxmodidxreprt.clear();
    xyzreprt.clear();
    adjcyreprt.clear();
    edgmodidxreprt.clear();
    statereprt.clear();
    stickreprt.clear();
    eventreprt.clear();
    
    // Set up recomputation of vtxdist
    vtxdist.clear();
    vtxdist.resize(netparts+1);
    vtxdist[0] = 0;

    CkPrintf("  Repartition: %d (%d)   Vertices: %" PRIidx "\n",
             prtidx, datidx, norderprt);

    // set up containers
    vtxmodidx.resize(norderprt);
    xyz.resize(norderprt*3);
    adjcy.resize(norderprt);
    edgmodidx.resize(norderprt);
    state.resize(norderprt);
    stick.resize(norderprt);
    evtflat.resize(norderprt);
    // reordering
    adjcyreord.resize(norderprt);
    edgmodidxreord.resize(norderprt);
    statereord.resize(norderprt);
    stickreord.resize(norderprt);
    eventindexreord.clear();
    eventsourcereord.resize(norderprt);
    eventindexreord.resize(norderprt);

    // Go through part data and reorder based on modidx
    std::sort(vtxprted.begin(), vtxprted.end());

    // add to data structures
    for (idx_t i = 0; i < norderprt; ++i) {
      // vtxmodidx
      vtxmodidx[i] = vtxprted[i].modidx;
      // xyz
      xyz[i*3+0] = xyzprted[(vtxprted[i].vtxidxloc)*3+0];
      xyz[i*3+1] = xyzprted[(vtxprted[i].vtxidxloc)*3+1];
      xyz[i*3+2] = xyzprted[(vtxprted[i].vtxidxloc)*3+2];
      // adjcy
      adjcyreord[i] = adjcyprted[vtxprted[i].vtxidxloc];
      edgmodidxreord[i] = edgmodidxprted[vtxprted[i].vtxidxloc];
      // state (for vertex)
      state[i].push_back(stateprted[vtxprted[i].vtxidxloc][0]);
      stick[i].push_back(stickprted[vtxprted[i].vtxidxloc][0]);
      // rest of state
      statereord[i] = stateprted[vtxprted[i].vtxidxloc];
      stickreord[i] = stickprted[vtxprted[i].vtxidxloc];
      // events
      evtflat[i] = eventprted[vtxprted[i].vtxidxloc];
      eventsourcereord[i].resize(evtflat[i].size());
      for (std::size_t e = 0; e < evtflat[i].size(); ++e) {
        eventsourcereord[i][e] = evtflat[i][e].source;
      }
      eventindexreord[i].push_back(0);
    }
    // Clear completed structures
    xyzprted.clear();
    adjcyprted.clear();
    edgmodidxprted.clear();
    stateprted.clear();
    stickprted.clear();
    eventprted.clear();

    // Take care of any ordering that may have come in
    // Go through ordering list
    for (std::list<mReorder *>::iterator iordidx = reordlist.begin(); iordidx != reordlist.end(); ++iordidx) {
      if ((*iordidx)->prtidx == cpprt) {
        // Perform reordering
        ReordEdg((*iordidx));
        // Move to next part
        ++cpprt;
        // Erase element from list
        iordidx = reordlist.erase(iordidx);
      }
    }

    // Broadcast reordlist in order of file
    if (cpprt == prtidx) {
      mReorder *morder = BuildReorder();
      thisProxy.Reorder(morder);
    }
  }
}


// Handle ordering messages
//
void Network::Reorder(mReorder *msg) {
  // Save message for processing
  reordlist.push_back(msg);

  if (adjcy.size()) {
    // Go through ordering list
    for (std::list<mReorder *>::iterator iordidx = reordlist.begin(); iordidx != reordlist.end(); ++iordidx) {
      if ((*iordidx)->prtidx == cpprt) {
        // Perform reordering
        ReordEdg((*iordidx));
        // Move to next part
        ++cpprt;
        // Erase element from list
        iordidx = reordlist.erase(iordidx);
      }
    }
  }

  // Check if done reordering
  if (cpprt == netparts) {
    cpprt = 0;
    // reindex events
    for (idx_t i = 0; i < norderprt; ++i) {
      CkAssert(eventindexreord[i].size() == adjcy[i].size()+1);
      // create map
      std::unordered_map<idx_t, idx_t> oldtonew;
      for (std::size_t j = 0; j < eventindexreord[i].size(); ++j) {
        oldtonew[eventindexreord[i][j]] = j;
      }
      // modify index
      for (std::size_t j = 0; j < evtflat[i].size(); ++j) {
        evtflat[i][j].index = oldtonew[evtflat[i][j].index];
      }
    }
    // Copy events back on to queue
    evtcal.clear();
    evtcal.resize(norderprt);
    evtcol.clear();
    evtcol.resize(norderprt);
    for (idx_t i = 0; i < norderprt; ++i) {
      evtcal[i].resize(nevtday);
      idx_t arrival;
      for (std::size_t e = 0; e < evtflat[i].size(); ++e) {
        // Add to event queue or spillover
        arrival = (idx_t) (evtflat[i][e].diffuse/tstep);
        CkAssert(arrival >= iter);
        if (arrival - iter < nevtday) {
          evtcal[i][(arrival)%nevtday].push_back(evtflat[i][e]);
        }
        else {
          evtcol[i].push_back(evtflat[i][e]);
        }
      }
    }
    // Clear finished data structs
    vtxprted.clear();
    evtflat.clear();
    adjcyreord.clear();
    edgmodidxreord.clear();
    statereord.clear();
    stickreord.clear();
    eventindexreord.clear();
    eventsourcereord.clear();
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
    CkPrintf("Part %d size/cap: adjcy: %d , %d edgmodidx: %d , %d\n", prtidx, adjcysize, adjcycap, edgmodsize, edgmodcap);
    */

    if (tsim > 0) {
      contribute(0, NULL, CkReduction::nop, maincheck);
    }
    else {
      // return control to main when done
      contribute(0, NULL, CkReduction::nop);
    }
  }
  // Broadcast ordering in order of file
  else if (cpprt == prtidx) {
    mReorder *morder = BuildReorder();
    thisProxy.Reorder(morder);
  }
}

// Reordering indices based on given ordering
//
void Network::ReordEdg(mReorder *msg) {
  // Add to vtxdist
  vtxdist[cpprt+1] = vtxdist[cpprt] + msg->nvtx;

  // create map
  std::unordered_map<idx_t, idx_t> oldtonew;
  for (idx_t i = 0; i < msg->nvtx; ++i) {
    oldtonew[msg->vtxidxold[i]] = vtxdist[cpprt] + msg->vtxidxnew[i];
  }

  // cleanup
  delete msg;

  // Reorder edges
  for (idx_t i = 0; i < norderprt; ++i) {
    edgreord.clear();
    for (std::size_t j = 0; j < adjcyreord[i].size(); ++j) {
      if (oldtonew.find(adjcyreord[i][j]) == oldtonew.end()) {
        continue;
      }
      else {
        edgreord.push_back(edgreord_t());
        edgreord.back().edgidx = oldtonew[adjcyreord[i][j]];
        edgreord.back().modidx = edgmodidxreord[i][j];
        edgreord.back().state = statereord[i][j+1];
        edgreord.back().stick = stickreord[i][j+1];
        edgreord.back().evtidx = j+1;
        // resource events
        for (std::size_t e = 0; e < eventsourcereord[i].size(); ++e) {
          if (eventsourcereord[i][e] == adjcyreord[i][j]) {
            evtflat[i][e].source = oldtonew[eventsourcereord[i][e]];
          }
        }
      }
    }
    // sort newly added indices
    std::sort(edgreord.begin(), edgreord.end());
    // add indices to data structures
    for (std::size_t j = 0; j < edgreord.size(); ++j) {
      adjcy[i].push_back(edgreord[j].edgidx);
      edgmodidx[i].push_back(edgreord[j].modidx);
      state[i].push_back(edgreord[j].state);
      stick[i].push_back(edgreord[j].stick);
      eventindexreord[i].push_back(edgreord[j].evtidx);
    }
  }
}


/**************************************************************************
* Reordering messages
**************************************************************************/

// Build Reorder (from old vertices to new indices)
//
mReorder* Network::BuildReorder() {
  // Initialize connection message
  int msgSize[MSG_Reorder];
  msgSize[0] = norderprt;   // vtxidxold
  msgSize[1] = norderprt;   // vtxidxnew
  mReorder *morder = new(msgSize, 0) mReorder;
  // sizes
  morder->prtidx = prtidx;
  morder->nvtx = norderprt;
  //CkPrintf("build reorder %" PRIidx "\n",norderprt);

  // set up counters
  idx_t jvtxidx = 0;

  // load data
  for (idx_t i = 0; i < norderprt; ++i) {
    morder->vtxidxold[jvtxidx] = vtxprted[i].vtxidx;
    morder->vtxidxnew[jvtxidx] = jvtxidx;
    ++jvtxidx;
  }
  CkAssert(jvtxidx == norderprt);

  return morder;
}

