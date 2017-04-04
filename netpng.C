/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"


/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ idx_t npnet;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ idx_t nevtday;


/**************************************************************************
* Reduction for PNGs
**************************************************************************/

CkReduction::reducerType net_png;
/*initnode*/
void registerNetPNG(void) {
  net_png = CkReduction::addReducer(netPNG);
}

CkReductionMsg *netPNG(int nMsg, CkReductionMsg **msgs) {
  std::vector<stamp_t> ret;
  ret.clear();
  for (int i = 0; i < nMsg; i++) {
    for (std::size_t j = 0; j < msgs[i]->getSize()/sizeof(stamp_t); ++j) {
      // Extract data and reduce 
      ret.push_back(*((stamp_t *)msgs[i]->getData() + j));
    }
  }
  return CkReductionMsg::buildNew(ret.size()*sizeof(stamp_t), ret.data());
}


/**************************************************************************
* Polychronization Initialization
**************************************************************************/

// Coordination with NetData chare array
//
void Network::InitPNG(CProxy_Netdata cpdat) {
  // Set proxies
  netdata = cpdat;
  cbcycleprt = CkCallback(CkIndex_Network::CyclePNG(), prtidx, thisProxy);
  
  // Request network part from input
  CkCallback *cb = new CkCallback(CkIndex_Network::LoadNetwork(NULL), prtidx, thisProxy);
  netdata(datidx).LoadNetwork(prtidx, *cb);
}


/**************************************************************************
* Finding Polychronous Neuronal Groups (PNGs)
**************************************************************************/

// Find PNG (main control loop)
//
void Network::FindPNG() {
  // Loop through all vertices
  if (compidx < vtxdist[npnet]) {
    pngseeds.clear();
    // Only one vertex containing partition performs control
    std::unordered_map<idx_t, idx_t>::iterator mother = vtxmap.find(compidx);
    if (mother != vtxmap.end()) {
      if (netmodel[vtxmodidx[mother->second]]->getPNGMother()) {
        // Bookkeeping
        idx_t i = mother->second;
        pngs[i].clear();

        // Skip vertices with less than three
        // TODO: make the number of anchor vertices configurable
        if (edgmodidx[i].size() < 3) {
          // continue to next vertex
          thisProxy.FindPNG();
        }
        else {
          // Use only valid anchor edges
          std::vector<idx_t> anchor;
          anchor.clear();
          for (idx_t j = 0; j < edgmodidx[i].size(); ++j) {
            if (netmodel[edgmodidx[i][j]]->getPNGAnchor()) {
              // TODO: Based off of spiking property of the 
              //       netmodel instead of just anchor models
              if (netmodel[edgmodidx[i][j]]->getStickIdx("delay") == 0) {
                anchor.push_back(j);
              }
            }
          }

          // PNG combinatorics
          for (idx_t j0 = 0; j0 < anchor.size(); ++j0) {
            for (idx_t j1 = j0+1; j1 < anchor.size(); ++j1) {
              for (idx_t j2 = j1+1; j2 < anchor.size(); ++j2) {
                // Test for spiking of mother neuron
                // assuming perfect timing of anchor
                netmodel[vtxmodidx[i]]->Reset(state[i][0], stick[i][0]);
                event_t event;
                event.diffuse = 0;
                event.type = EVENT_SPIKE;
                event.source = i;
                event.data = 0.0;
                event.index = anchor[j0]+1;
                netmodel[edgmodidx[i][anchor[j0]]]->Hop(event, state[i], stick[i], edgaux[edgmodidx[i][anchor[j0]]][vtxmodidx[i]]);
                event.index = anchor[j1]+1;
                netmodel[edgmodidx[i][anchor[j1]]]->Hop(event, state[i], stick[i], edgaux[edgmodidx[i][anchor[j1]]][vtxmodidx[i]]);
                event.index = anchor[j2]+1;
                netmodel[edgmodidx[i][anchor[j2]]]->Hop(event, state[i], stick[i], edgaux[edgmodidx[i][anchor[j2]]][vtxmodidx[i]]);
                tick_t tdrift = 0;
                tick_t tstop = tstep * 5; // Strongly spiking triplets only
                while (tdrift < tstop) {
                  // Step through model drift (vertex)
                  tdrift += netmodel[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], evtlog);
                }
                // Potential PNG only if mother vertex spiked
                if (evtlog.size() && evtlog[0].type == EVENT_SPIKE) {
                  evtlog.clear();
                  std::vector<event_t> pngseed;
                  pngseed.clear();
                  event.type = EVENT_SPIKE;
                  event.data = 0.0;
                  event.diffuse = stick[i][anchor[j0]+1][0];
                  event.source = adjcy[i][anchor[j0]];
                  event.index = adjcy[i][anchor[j0]];
                  pngseed.push_back(event);
                  event.diffuse = stick[i][anchor[j1]+1][0];
                  event.source = adjcy[i][anchor[j1]];
                  event.index = adjcy[i][anchor[j1]];
                  pngseed.push_back(event);
                  event.diffuse = stick[i][anchor[j2]+1][0];
                  event.source = adjcy[i][anchor[j2]];
                  event.index = adjcy[i][anchor[j2]];
                  pngseed.push_back(event);
                  // correctly order the timing
                  std::sort(pngseed.begin(), pngseed.end());
                  pngseed[0].diffuse = pngseed[2].diffuse - pngseed[0].diffuse;
                  pngseed[1].diffuse = pngseed[2].diffuse - pngseed[1].diffuse;
                  pngseed[2].diffuse = pngseed[2].diffuse - pngseed[2].diffuse;
                  std::sort(pngseed.begin(), pngseed.end());
                  // Push to PNG seeds
                  pngseeds.push_back(pngseed);
                }
              }
            }
          }
          ncomp = pngseeds.size();
          // Display computation information
          CkPrintf("    Computing PNGs for vertex %" PRIidx " with %" PRIidx "\n", compidx, ncomp);
          thisProxy.ComputePNG(((ncomp > 10) ? 10 : ncomp), prtidx);
        }
      }
      else {
        // continue to next vertex
        thisProxy.FindPNG();
      }
    }
    ++compidx;
  }
  else {
    // return control to main
    contribute(0, NULL, CkReduction::nop);
  }
}

/**************************************************************************
* Computing PNGs
**************************************************************************/

// Initial setup for PNG computation
//
void Network::ComputePNG(idx_t nseeds, idx_t pngidx) {
  // Bookkeeping
  ccomp = 0;
  ncomp = nseeds;
  evalidx = pngidx;
  tcomp = tstep * 150;
  // TODO: make this max png length

  // Coordination after reset
  thisProxy(prtidx).ComputePNG();
}

// Compute PNG (vertex control loop)
//
void Network::ComputePNG() {
  // Loop through PNG seeds
  if (ccomp < ncomp) {
    pngcan.clear();
    pngaux.clear();
    if (!pngseeds.empty()) {
      // Initialize candidate PNG
      pngcan.resize(pngseeds[ccomp].size());
      for (std::size_t i = 0; i < pngseeds[ccomp].size(); ++i) {
        pngcan[i].diffuse = pngseeds[ccomp][i].diffuse;
        pngcan[i].source = pngseeds[ccomp][i].source;
        pngcan[i].origin = -1;
        pngcan[i].departure = 0;
        pngcan[i].arrival = 0;
      }
      // Seed spikes for simulation
      mEvent *mevent = BuildPNGSeed(pngseeds[ccomp]);
      thisProxy.SeedPNG(mevent);
    }
    ++ccomp;
  }
  else {
    std::unordered_map<idx_t, idx_t>::iterator mother = vtxmap.find(compidx-1);
    if (mother != vtxmap.end()) {
      idx_t i = mother->second;
      CkPrintf("  Found %d PNGs\n", pngs[i].size());
      // Write to file
      WritePNG(i);
    }
    // Return control to main loop
    thisProxy(prtidx).FindPNG();
  }
}

// Coordination for PNG computation
//
void Network::EvalPNG(CkReductionMsg *msg) {
  // Add to png candidate
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(stamp_t); ++i) {
    pngcan.push_back(*((stamp_t *)msg->getData()+i));
  }
  delete msg;

  //CkPrintf("Evaluating\n");
  // Sorting
  std::sort(pngcan.begin(), pngcan.end());

  // Max path of png should be longer than min path threshold (7)
  std::unordered_map<idx_t, int> pngpath;
  int maxpath = 0;
  for (std::size_t i = 0; i < pngcan.size(); ++i) {
    pngpath[pngcan[i].source] = std::max(pngpath[pngcan[i].source], 1+pngpath[pngcan[i].origin]);
    maxpath = std::max(maxpath, pngpath[pngcan[i].source]);
  }
  if (maxpath > 7) {
    // Anchors should contribute to more than just the mother neuron
    bool alluseful = true;
    for (std::size_t j = 0; j < pngseeds[ccomp-1].size(); ++j) {
      int useful = 0;
      for (std::size_t i = 0; i < pngcan.size(); ++i) {
        if (pngcan[i].origin == pngseeds[ccomp-1][j].source) {
          if (++useful >= 2) { break; }
        }
      }
      if (useful < 2) {
        alluseful = false;
        break;
      }
    }
    if (alluseful) {
      std::unordered_map<idx_t, idx_t>::iterator mother = vtxmap.find(compidx-1);
      idx_t i = mother->second;
      pngs[i].push_back(pngcan);
    }
  }

  // Compute next PNG
  thisProxy.ComputePNG();
}

/**************************************************************************
* Computing PNGs (simulation)
**************************************************************************/

// Simulation loop for PNG computation
//
void Network::CyclePNG() {
  // Check if computation is complete
  if (tsim >= tcomp) {
    // Reset network
    ResetNetwork();
    tsim = 0;
    iter = 0;
    cadjprt[0] = 0;
    cadjprt[1] = 0;
    prtiter = 0;
    
    // Coordination after reset
    // Reduce PNG information
    CkCallback *cb = new CkCallback(CkIndex_Network::EvalPNG(NULL), evalidx, thisProxy);
    contribute(pngaux.size()*sizeof(stamp_t), pngaux.data(), net_png, *cb);
  }
#ifdef STACS_WITH_YARP
  // Synchronization from RPC
  else if (iter == synciter) {
    // Bookkkeeping
    synciter = IDX_T_MAX;

    // Display synchronization information
    if (prtidx == 0) {
      CkPrintf("  Synchronizing at iteration %" PRIidx "\n", iter);
    }

    // move control to sychronization callback
    contribute(0, NULL, CkReduction::nop);
  }
#endif
  else {
    // Bookkeeping
    idx_t evtday = iter%nevtday;
    tick_t tstop = tsim + tstep;

    // Clear event buffer
    evtext.clear();
    idx_t nevent = 0;
    // Redistribute any events (on new year)
    if (evtday == 0) {
      MarkEvent();
    }
    
    // Check for repeating events
    if (tsim >= trep) {
      std::vector<event_t>::iterator event = repevt.begin();
      // Compute periodic events
      while (event != repevt.end() && event->diffuse <= tsim) {
        // Set temporary model index
        idx_t modidx = event->index;
        // Loop through all models
        for (std::size_t i = 0; i < repidx[modidx].size(); ++i) {
          event->index = repidx[modidx][i][1];
          if (event->index) {
            netmodel[modidx]->Hop(*event, state[repidx[modidx][i][0]], stick[repidx[modidx][i][0]], edgaux[modidx][vtxmodidx[repidx[modidx][i][0]]]);
          }
          else {
            netmodel[modidx]->Hop(*event, state[repidx[modidx][i][0]], stick[repidx[modidx][i][0]], vtxaux[repidx[modidx][i][0]]);
          }
        }
        // Return model index
        event->index = modidx;
        // Update timing
        event->diffuse += ((tick_t) event->data)* TICKS_PER_MS;
        ++event;
      }
      std::sort(repevt.begin(), repevt.end());
      trep = repevt[0].diffuse;
    }
    
    // Perform computation
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
      if (netmodel[vtxmodidx[i]]->getPNGActive() == false) {
        evtcal[i][evtday].clear();
        continue;
      }
      // Timing
      tick_t tdrift = tsim;

      // Sort events
      std::sort(evtcal[i][evtday].begin(), evtcal[i][evtday].end());
      nevent += evtcal[i][evtday].size();

      // Perform events starting at beginning of step
      std::vector<event_t>::iterator event = evtcal[i][evtday].begin();
      while (event != evtcal[i][evtday].end() && event->diffuse <= tdrift) {
        // edge events
        if (event->index) {
          netmodel[edgmodidx[i][event->index-1]]->Hop(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
          // Add to PNG log
          if (event->type == EVENT_SPIKE) {
            route_t pngpre;
            pngpre.origin = adjcy[i][event->index-1];
            pngpre.departure = event->diffuse - stick[i][event->index][0];
            pngpre.arrival = event->diffuse;
            pnglog[i].push_back(pngpre);
          }
        }
        // vertex events
        else {
          netmodel[vtxmodidx[i]]->Hop(*event, state[i], stick[i], vtxaux[i]);
        }
        ++event;
      }
      
      // Move sliding window of contributing routes forward
      // TODO: make this value user adjustable
      while (!pnglog[i].empty()) {
        if (pnglog[i].front().arrival + ((tick_t)10.0)*TICKS_PER_MS <= tdrift) {
          pnglog[i].pop_front();
        }
        else {
          break;
        }
      }

      // Computation
      while (tdrift < tstop) {
        // Step through model drift (vertex)
        tdrift += netmodel[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], evtlog);

        // Handle generated events (if any)
        // TODO: Conversion from edge indices to global (for individual output)
        if (evtlog.size()) {
          for (std::size_t e = 0; e < evtlog.size(); ++e) {
            // Polychronization information
            if (evtlog[e].type == EVENT_SPIKE) {
              stamp_t pngpre;
              pngpre.diffuse = evtlog[e].diffuse;
              pngpre.source = vtxidx[i];
              // go through contribution log
              for (std::size_t j = 0; j < pnglog[i].size(); ++j) {
                pngpre.origin = pnglog[i][j].origin;
                pngpre.departure = pnglog[i][j].departure;
                pngpre.arrival = pnglog[i][j].arrival;
                pngaux.push_back(pngpre);
              }
            }
            // Get information
            idx_t target = evtlog[e].source;
            idx_t index = evtlog[e].index;
            // Remote events (multicast to edges)
            if (target & REMOTE_EDGES) {
              // reindex to global
              evtlog[e].source = vtxidx[i];
              evtlog[e].index = vtxidx[i];
              // push to communication
              evtext.push_back(evtlog[e]);
            }
            // Remote event (singlecast to edge)
            else if (target & REMOTE_EDGE) {
              // reindex to global
              evtlog[e].source = vtxidx[i];
              // TODO: get this value from the target mapping
              evtlog[e].index = adjcy[i][index];
              // push to communication
              evtext.push_back(evtlog[e]);
            }
            // Remote event (singlecast to vertex)
            else if (target & REMOTE_VERTEX) {
              // reindex to global
              evtlog[e].source = vtxidx[i];
              // TODO: get this value from the target mapping
              evtlog[e].index = -adjcy[i][index]-1; // negative index indicates vertex
              // push to communication
              evtext.push_back(evtlog[e]);
            }
            // Local events (multicast to edges)
            if (target & LOCAL_EDGES) {
              evtlog[e].source = -vtxidx[i]-1; // negative source indicates local event
              // Jump loops
              if ((evtlog[e].diffuse - tsim - tstep)/tstep < nevtday) {
                for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
                  if (edgmodidx[i][j]) {
                    evtlog[e].index = j+1;
                    evtcal[i][(evtlog[e].diffuse/tstep)%nevtday].push_back(evtlog[e]);
                  }
                }
              }
              else if (evtlog[e].diffuse < tsim + tstep) {
                for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
                  if (edgmodidx[i][j]) {
                    evtlog[e].index = j+1;
                    // Jump now
                    netmodel[edgmodidx[i][j]]->Hop(evtlog[e], state[i], stick[i], edgaux[edgmodidx[i][j]][vtxmodidx[i]]);
                  }
                }
              }
              else {
                for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
                  if (edgmodidx[i][j]) {
                    evtlog[e].index = j+1;
                    evtcol[i].push_back(evtlog[e]);
                  }
                }
              }
            }
            // Local event (singlecast to vertex)
            if (target & LOCAL_VERTEX) {
              // vertex to itself
              evtlog[e].source = -vtxidx[i]-1; // negative source indicates local event
              evtlog[e].index = 0;
              if ((evtlog[e].diffuse - tsim - tstep)/tstep < nevtday) {
                evtcal[i][(evtlog[e].diffuse/tstep)%nevtday].push_back(evtlog[e]);
              }
              else if (evtlog[e].diffuse < tsim + tstep) {
                // Jump now
                netmodel[vtxmodidx[i]]->Hop(evtlog[e], state[i], stick[i], vtxaux[i]);
              }
              else {
                evtcol[i].push_back(evtlog[e]);
              }
            }
          }
          // clear log for next time
          evtlog.clear();
        }
        
        // Perform events up to tdrift
        while (event != evtcal[i][evtday].end() && event->diffuse <= tdrift) {
          // edge events
          if (event->index) {
            netmodel[edgmodidx[i][event->index-1]]->Hop(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
          }
          // vertex events
          else {
            netmodel[vtxmodidx[i]]->Hop(*event, state[i], stick[i], vtxaux[i]);
          }
          ++event;
        }
      }

      // Clear event queue
      evtcal[i][evtday].clear();
    }

    // Send messages to neighbors
    mEvent *mevent = BuildEvent();
    netgroup.CommEvent(mevent);
    
    // Increment simulated time
    tsim += tstep;
  }
}

/**************************************************************************
* PNG Events
**************************************************************************/

// Seed events for PNG computation
//
void Network::SeedPNG(mEvent *msg) {
  // Event prototype
  event_t event;
  tick_t departure;

  // Distribute events
  for (std::size_t i = 0; i < msg->nevent; ++i) {
    // Fill in prototype
    departure = msg->diffuse[i];
    event.type = msg->type[i];
    event.source = msg->source[i];
    event.data = msg->data[i];
    // Determine local event target(s)
    // If index == source (multicast to edges)
    // Find target mapping from source
    std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>>::iterator targets = adjmap.find(msg->source[i]);
    if (targets != adjmap.end()) {
      for (std::vector<std::array<idx_t, 2>>::iterator target = targets->second.begin(); target != targets->second.end(); ++target) {
        event.diffuse = departure + stick[(*target)[0]][(*target)[1]][0]; // delay always first stick of edge
        event.index = (*target)[1];
        // Add to event queue or spillover
        if (event.diffuse/tstep < nevtday) {
          evtcal[(*target)[0]][(event.diffuse/tstep)%nevtday].push_back(event);
        }
        else {
          evtcol[(*target)[0]].push_back(event);
        }
      }
    }
  }
  delete msg;

  // Start cycle after seeding events
  thisProxy(prtidx).CyclePNG();
}


/**************************************************************************
* Build Messages
**************************************************************************/

// Build event seed for PNG computation
//
mEvent* Network::BuildPNGSeed(std::vector<event_t>& pngseed) {
  // Initialize distribution message
  int msgSize[MSG_Event];
  msgSize[0] = pngseed.size();     // diffuse
  msgSize[1] = pngseed.size();     // type
  msgSize[2] = pngseed.size();     // source
  msgSize[3] = pngseed.size();     // index
  msgSize[4] = pngseed.size();     // data
  mEvent *mevent = new(msgSize, 0) mEvent;
  mevent->nevent = pngseed.size();
  mevent->iter = 0;
  
  // Pack event information
  for (std::size_t i = 0; i < pngseed.size(); ++i) {
    // Add event to message
    mevent->diffuse[i] = pngseed[i].diffuse;
    mevent->type[i] = pngseed[i].type;
    mevent->source[i] = pngseed[i].source;
    mevent->index[i] = pngseed[i].index;
    mevent->data[i] = pngseed[i].data;
  }

  return mevent;
}

