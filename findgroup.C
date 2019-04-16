/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"


/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ idx_t nevtday;
extern /*readonly*/ int grpminlen;
extern /*readonly*/ tick_t grpmaxdur;
extern /*readonly*/ idx_t grpvtxmin;
extern /*readonly*/ idx_t grpvtxmax;
extern CkReduction::reducerType max_idx;


/**************************************************************************
* Reduction for Groups
**************************************************************************/

CkReduction::reducerType net_group;
/*initnode*/
void registerNetGroup(void) {
  net_group = CkReduction::addReducer(netGroup);
}

CkReductionMsg *netGroup(int nMsg, CkReductionMsg **msgs) {
  std::vector<route_t> ret;
  ret.clear();
  for (int i = 0; i < nMsg; i++) {
    for (std::size_t j = 0; j < msgs[i]->getSize()/sizeof(route_t); ++j) {
      // Extract data and reduce 
      ret.push_back(*((route_t *)msgs[i]->getData() + j));
    }
  }
  return CkReductionMsg::buildNew(ret.size()*sizeof(route_t), ret.data());
}


/**************************************************************************
* Polychronization Initialization
**************************************************************************/

// Coordination with NetData chare array
//
void Network::InitGroup(CProxy_Netdata cpdata) {
  // Set proxies
  netdata = cpdata;
  cyclepart = CkCallback(CkIndex_Network::CycleGroup(), partidx, thisProxy);

  // Initialization
  tcomp = 0;
  compidx = grpvtxmin;
  ccomp = 0;
  ncomp = 0;
  compart = 0;
  
  // Request network part from input
  netdata(fileidx).LoadNetwork(partidx, 
      CkCallback(CkIndex_Network::LoadNetwork(NULL), partidx, thisProxy));
}


/**************************************************************************
* Finding Polychronous Neuronal Groups
**************************************************************************/

// Find Group (main control loop)
//
void Network::FindGroup() {
  // Loop through range of vertices
  if (compidx < grpvtxmax) {
    grpseeds.clear();
    // Only one vertex containing partition performs control
    std::unordered_map<idx_t, idx_t>::iterator mother = vtxmap.find(compidx);
    if (mother != vtxmap.end()) {
      if (model[vtxmodidx[mother->second]]->getMother()) {
        // Bookkeeping
        idx_t i = mother->second;
        grpstamps[i].clear();
        grproutes.clear();

        // Skip vertices with less than three
        // TODO: make the number of anchor vertices configurable
        if (edgmodidx[i].size() < 3) {
          // continue to next vertex
          thisProxy.FindGroup();
        }
        else {
          // Use only valid anchor edges
          std::vector<idx_t> anchor;
          anchor.clear();
          for (idx_t j = 0; j < edgmodidx[i].size(); ++j) {
            if (model[edgmodidx[i][j]]->getAnchor()) {
              // TODO: Based off of spiking property of the 
              //       model instead of just anchor models
              if (model[edgmodidx[i][j]]->getStickIdx("delay") == 0) {
                anchor.push_back(j);
              }
            }
          }

          // Group combinatorics
          for (idx_t j0 = 0; j0 < anchor.size(); ++j0) {
            for (idx_t j1 = j0+1; j1 < anchor.size(); ++j1) {
              for (idx_t j2 = j1+1; j2 < anchor.size(); ++j2) {
                // Test for spiking of mother neuron
                // assuming perfect timing of anchor
                model[vtxmodidx[i]]->Reset(state[i][0], stick[i][0]);
                event_t event;
                event.diffuse = 0;
                event.type = EVENT_SPIKE;
                event.source = i;
                event.data = 0.0;
                event.index = anchor[j0]+1;
                model[edgmodidx[i][anchor[j0]]]->Jump(event, state[i], stick[i], edgaux[edgmodidx[i][anchor[j0]]][vtxmodidx[i]]);
                event.index = anchor[j1]+1;
                model[edgmodidx[i][anchor[j1]]]->Jump(event, state[i], stick[i], edgaux[edgmodidx[i][anchor[j1]]][vtxmodidx[i]]);
                event.index = anchor[j2]+1;
                model[edgmodidx[i][anchor[j2]]]->Jump(event, state[i], stick[i], edgaux[edgmodidx[i][anchor[j2]]][vtxmodidx[i]]);
                tick_t tdrift = 0;
                tick_t tstop = tstep * 5; // Strongly spiking triplets only
                while (tdrift < tstop) {
                  // Step through model drift (vertex)
                  tdrift += model[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], events);
                }
                // Potential Group only if mother vertex spiked
                if (events.size() && events[0].type == EVENT_SPIKE) {
                  events.clear();
                  std::vector<event_t> grpseed;
                  grpseed.clear();
                  event.type = EVENT_SPIKE;
                  event.data = 0.0;
                  event.diffuse = stick[i][anchor[j0]+1][0];
                  event.source = adjcy[i][anchor[j0]];
                  event.index = adjcy[i][anchor[j0]];
                  grpseed.push_back(event);
                  event.diffuse = stick[i][anchor[j1]+1][0];
                  event.source = adjcy[i][anchor[j1]];
                  event.index = adjcy[i][anchor[j1]];
                  grpseed.push_back(event);
                  event.diffuse = stick[i][anchor[j2]+1][0];
                  event.source = adjcy[i][anchor[j2]];
                  event.index = adjcy[i][anchor[j2]];
                  grpseed.push_back(event);
                  // correctly order the timing
                  std::sort(grpseed.begin(), grpseed.end());
                  grpseed[0].diffuse = grpseed[2].diffuse - grpseed[0].diffuse;
                  grpseed[1].diffuse = grpseed[2].diffuse - grpseed[1].diffuse;
                  grpseed[2].diffuse = grpseed[2].diffuse - grpseed[2].diffuse;
                  std::sort(grpseed.begin(), grpseed.end());
                  // Push to Group seeds
                  grpseeds.push_back(grpseed);
                }
              }
            }
          }
          ncomp = grpseeds.size();
          // Display computation information
          CkPrintf("  Computing vertex %" PRIidx " groups %" PRIidx "\n", compidx, ncomp);
          thisProxy.ComputeGroup(ncomp, partidx);
        }
      }
      else {
        // continue to next vertex
        thisProxy.FindGroup();
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
* Computing Groups
**************************************************************************/

// Initial setup for Group computation
//
void Network::ComputeGroup(idx_t nseeds, int grpart) {
  // Bookkeeping
  ccomp = 0;
  ncomp = nseeds;
  compart = grpart;
  tcomp = grpmaxdur;

  thisProxy(partidx).ComputeGroup();
}

// Compute Group (vertex control loop)
//
void Network::ComputeGroup() {
  // Loop through Group seeds
  if (ccomp < ncomp) {
    grpleg.clear();
    grproute.clear();
    if (!grpseeds.empty()) {
      // Initialize candidate group
      grproute.resize(grpseeds[ccomp].size());
      for (std::size_t i = 0; i < grpseeds[ccomp].size(); ++i) {
        grproute[i].diffuse = grpseeds[ccomp][i].diffuse;
        grproute[i].source = grpseeds[ccomp][i].source;
        grproute[i].origin = -1;
        grproute[i].departure = 0;
        grproute[i].arrival = 0;
      }
  
      // Seed spikes for simulation
      mEvent *mevent = BuildGroupSeed(grpseeds[ccomp]);
      thisProxy.SeedGroup(mevent);
    }
    ++ccomp;
  }
  else {
    std::unordered_map<idx_t, idx_t>::iterator mother = vtxmap.find(compidx-1);
    if (!grpseeds.empty() && mother != vtxmap.end()) {
      idx_t groupidx = mother->second;
      CkPrintf("  Groups found %d\n", grpstamps[groupidx].size());
      if (grpstamps[groupidx].size()) {
        // Write to file
        WriteGroup(groupidx);
      }
      // Clear found groups after writing
      grpstamps[groupidx].clear();
      grproutes.clear();
    }
    // Return control to main loop
    thisProxy(partidx).FindGroup();
  }
}

// Coordination for Group computation
//
void Network::EvalGroup(CkReductionMsg *msg) {
  // Add to group candidate
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(route_t); ++i) {
    grproute.push_back(*((route_t *)msg->getData()+i));
  }
  delete msg;

  //CkPrintf("Evaluating\n");
  // Sorting
  std::sort(grproute.begin(), grproute.end());

  // Max path of group should be longer than min path length
  std::unordered_map<idx_t, int> grpath;
  int maxlen = 0;
  for (std::size_t i = 0; i < grproute.size(); ++i) {
    grpath[grproute[i].source] = std::max(grpath[grproute[i].source], 1+grpath[grproute[i].origin]);
    maxlen = std::max(maxlen, grpath[grproute[i].source]);
  }
  if (maxlen >= grpminlen) {
    // Anchors should contribute to more than just the mother neuron
    bool alluseful = true;
    for (std::size_t j = 0; j < grpseeds[ccomp-1].size(); ++j) {
      int useful = 0;
      for (std::size_t i = 0; i < grproute.size(); ++i) {
        if (grproute[i].origin == grpseeds[ccomp-1][j].source) {
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
      idx_t groupidx = mother->second;
      std::set<stamp_t> grpset;
      for (std::size_t i = 0; i < grproute.size(); ++i) {
        grpset.insert((stamp_t){grproute[i].diffuse, grproute[i].source});
      }
      std::vector<stamp_t> grpvec;
      grpvec.assign(grpset.begin(), grpset.end());
      grpstamps[groupidx].push_back(grpvec);
      grproutes.push_back(grproute);
    }
  }

  // Compute next Group
  thisProxy.ComputeGroup();
}

/**************************************************************************
* Computing Groups (simulation)
**************************************************************************/

// Simulation loop for Group computation
//
void Network::CycleGroup() {
  // Check if computation is complete
  if (tsim >= tcomp) {
    // Reset network for next group computation
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
      // Clear events
      for (idx_t j = 0; j < nevtday; ++j) {
        evtcal[i][j].clear();
      }
      evtcol[i].clear();
      // Reset vertices
      model[vtxmodidx[i]]->Reset(state[i][0], stick[i][0]);
    }
    // Reset timing
    tsim = 0;
    iter = 0;
    // Reset coordination
    commiter = 0;
    cadjpart[0] = 0;
    cadjpart[1] = 0;
    partiter = 0;

    // Remove excess grptraces
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
      grptraces[i].clear();
    }
    
    // Reduce Group information
    contribute(grpleg.size()*sizeof(route_t), grpleg.data(), net_group, 
        CkCallback(CkIndex_Network::EvalGroup(NULL), compart, thisProxy));
  }
#ifdef STACS_WITH_YARP
  // Synchronization from RPC
  else if (syncing && synciter == IDX_T_MAX) {
    // nop
  }
  else if (iter == synciter) {
    if (!syncing) {
      // Bookkkeeping
      synciter = IDX_T_MAX;
      syncing = true;
      
      idx_t contiter = iter;
      // move control to sychronization callback
      contribute(sizeof(idx_t), &contiter, max_idx);
    }
    else {
      // Bookkkeeping
      synciter = IDX_T_MAX;
      
      // Display synchronization information
      if (partidx == 0) {
        CkPrintf("  Synchronizing at iteration %" PRIidx "\n", iter);
      }

      // move control to sychronization callback
      contribute(0, NULL, CkReduction::nop);
    }
  }
#endif
  else {
    // Bookkeeping
    idx_t evtday = iter%nevtday;
    tick_t tstop = tsim + tstep;

    // Clear event buffer
    evtext.clear();
    // Redistribute any events (on new year)
    if (evtday == 0) {
      SortEventCalendar();
    }
    
    // Check for periodic events
    if (tsim >= tskip) {
      SkipEvent();
    }
    
    // Perform computation
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
      if (model[vtxmodidx[i]]->getActive() == false) {
        evtcal[i][evtday].clear();
        continue;
      }
      // Timing
      tick_t tdrift = tsim;

      // Sort events
      std::sort(evtcal[i][evtday].begin(), evtcal[i][evtday].end());

      // Perform events starting at beginning of step
      std::vector<event_t>::iterator event = evtcal[i][evtday].begin();
      while (event != evtcal[i][evtday].end() && event->diffuse <= tdrift) {
        // edge events
        if (event->index) {
          model[edgmodidx[i][event->index-1]]->Jump(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
          // Add to contribution log
          if (event->type == EVENT_SPIKE) {
            trace_t trace;
            trace.origin = adjcy[i][event->index-1];
            trace.departure = event->diffuse - stick[i][event->index][0];
            trace.arrival = event->diffuse;
            grptraces[i].push_back(trace);
          }
        }
        // vertex events
        else {
          model[vtxmodidx[i]]->Jump(*event, state[i], stick[i], vtxaux[i]);
        }
        ++event;
      }
      
      // Move sliding window of contributing routes forward
      // TODO: make this value user adjustable
      while (!grptraces[i].empty()) {
        if (grptraces[i].front().arrival + ((tick_t)10.0)*TICKS_PER_MS <= tdrift) {
          grptraces[i].pop_front();
        }
        else {
          break;
        }
      }

      // Computation
      while (tdrift < tstop) {
        // Step through model drift (vertex)
        tdrift += model[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], events);

        // Handle generated events (if any)
        if (events.size()) {
          for (std::size_t e = 0; e < events.size(); ++e) {
            // Polychronization information
            if (events[e].type == EVENT_SPIKE) {
              route_t route;
              route.diffuse = events[e].diffuse;
              route.source = vtxidx[i];
              // go through contribution log
              for (std::size_t j = 0; j < grptraces[i].size(); ++j) {
                route.origin = grptraces[i][j].origin;
                route.departure = grptraces[i][j].departure;
                route.arrival = grptraces[i][j].arrival;
                grpleg.push_back(route);
              }
            }
            // TODO: Conversion from edge indices to global (for individual output)
            // Get information
            idx_t target = events[e].source;
            idx_t index = events[e].index;
            // Reindex to global
            events[e].source = vtxidx[i];
            // Remote events (multicast to edges)
            if (target & REMOTE_EDGES) {
              // reindex to global
              events[e].index = vtxidx[i];
              // push to communication
              evtext.push_back(events[e]);
            }
            // Remote event (singlecast to edge)
            else if (target & REMOTE_EDGE) {
              // reindex to global
              // TODO: get this value from the target mapping
              events[e].index = adjcy[i][index];
              // push to communication
              evtext.push_back(events[e]);
            }
            // Remote event (singlecast to vertex)
            else if (target & REMOTE_VERTEX) {
              // reindex to global
              // TODO: get this value from the target mapping
              events[e].index = -adjcy[i][index]-1; // negative index indicates vertex
              // push to communication
              evtext.push_back(events[e]);
            }
            // Local events (multicast to edges)
            if (target & LOCAL_EDGES) {
              events[e].source = -1; // negative source indicates local event
              // Jump loops
              if ((events[e].diffuse - tsim - tstep)/tstep < nevtday) {
                for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
                  if (edgmodidx[i][j]) {
                    events[e].index = j+1;
                    evtcal[i][(events[e].diffuse/tstep)%nevtday].push_back(events[e]);
                  }
                }
              }
              else if (events[e].diffuse < tsim + tstep) {
                for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
                  if (edgmodidx[i][j]) {
                    events[e].index = j+1;
                    // Jump now
                    model[edgmodidx[i][j]]->Jump(events[e], state[i], stick[i], edgaux[edgmodidx[i][j]][vtxmodidx[i]]);
                  }
                }
              }
              else {
                for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
                  if (edgmodidx[i][j]) {
                    events[e].index = j+1;
                    evtcol[i].push_back(events[e]);
                  }
                }
              }
            }
            // Local event (singlecast to vertex)
            if (target & LOCAL_VERTEX) {
              // vertex to itself
              events[e].source = -1; // negative source indicates local event
              events[e].index = 0;
              if ((events[e].diffuse - tsim - tstep)/tstep < nevtday) {
                evtcal[i][(events[e].diffuse/tstep)%nevtday].push_back(events[e]);
              }
              else if (events[e].diffuse < tsim + tstep) {
                // Jump now
                model[vtxmodidx[i]]->Jump(events[e], state[i], stick[i], vtxaux[i]);
              }
              else {
                evtcol[i].push_back(events[e]);
              }
            }
          }
          // clear log for next time
          events.clear();
        }
        
        // Perform events up to tdrift
        while (event != evtcal[i][evtday].end() && event->diffuse <= tdrift) {
          // edge events
          if (event->index) {
            model[edgmodidx[i][event->index-1]]->Jump(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
          }
          // vertex events
          else {
            model[vtxmodidx[i]]->Jump(*event, state[i], stick[i], vtxaux[i]);
          }
          ++event;
        }
      }

      // Clear event queue
      evtcal[i][evtday].clear();
    }

    // Send messages to neighbors
    mEvent *mevent = BuildEvent();
    netcomm.CommEvent(mevent);
    
    // Increment simulated time
    tsim += tstep;

    // Increment iteration
    ++iter;
  }
}

/**************************************************************************
* Group Events
**************************************************************************/

// Seed events for Group computation
//
void Network::SeedGroup(mEvent *msg) {
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
  thisProxy(partidx).CycleGroup();
}


/**************************************************************************
* Build Messages
**************************************************************************/

// Build event seed for Group computation
//
mEvent* Network::BuildGroupSeed(std::vector<event_t>& grpseed) {
  // Initialize distribution message
  int msgSize[MSG_Event];
  msgSize[0] = grpseed.size();     // diffuse
  msgSize[1] = grpseed.size();     // type
  msgSize[2] = grpseed.size();     // source
  msgSize[3] = grpseed.size();     // index
  msgSize[4] = grpseed.size();     // data
  mEvent *mevent = new(msgSize, 0) mEvent;
  mevent->nevent = grpseed.size();
  mevent->iter = 0;
  
  // Pack event information
  for (std::size_t i = 0; i < grpseed.size(); ++i) {
    // Add event to message
    mevent->diffuse[i] = grpseed[i].diffuse;
    mevent->type[i] = grpseed[i].type;
    mevent->source[i] = grpseed[i].source;
    mevent->index[i] = grpseed[i].index;
    mevent->data[i] = grpseed[i].data;
  }

  return mevent;
}

