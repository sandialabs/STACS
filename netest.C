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
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ tick_t tcheck;
extern /*readonly*/ tick_t tdisplay;
extern /*readonly*/ idx_t nevtday;
extern /*readonly*/ idx_t ntrials;


/**************************************************************************
* Reduction for Events
**************************************************************************/

CkReduction::reducerType net_event;
/*initnode*/
void registerNetEvent(void) {
  net_event = CkReduction::addReducer(netEvent);
}

CkReductionMsg *netEvent(int nMsg, CkReductionMsg **msgs) {
  std::vector<event_t> ret;
  ret.clear();
  for (int i = 0; i < nMsg; i++) {
    for (std::size_t j = 0; j < msgs[i]->getSize()/sizeof(event_t); ++j) {
      // Extract data and reduce 
      ret.push_back(*((event_t *)msgs[i]->getData() + j));
    }
  }
  return CkReductionMsg::buildNew(ret.size()*sizeof(event_t), ret.data());
}


/**************************************************************************
* Network Estimation Initialization
**************************************************************************/

// Coordination with NetData chare array
//
void Network::InitEstStatic(CProxy_Netdata cpdat) {
  // Set proxies
  netdata = cpdat;
  cbcycleprt = CkCallback(CkIndex_Network::CycleEstStatic(), prtidx, thisProxy);
  
  // Request network part from input
  CkCallback *cb = new CkCallback(CkIndex_Network::LoadNetwork(NULL), prtidx, thisProxy);
  netdata(datidx).LoadNetwork(prtidx, *cb);
}


/**************************************************************************
* Estimation Recording
**************************************************************************/

// Send Estimates for writing
//
void Network::SaveEstimate(CkReductionMsg *msg) {
  // Add to png list
  pnglist.clear();
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(route_t); ++i) {
    pnglist.push_back(*((event_t *)msg->getData()+i));
  }
  delete msg;

  // Sorting
  std::sort(pnglist.begin(), pnglist.end());

  // Write estimate
  WriteEstimate(evalidx-1);
  
  // Start a new cycle (checked data sent)
  thisProxy.CycleEstStatic();
}


/**************************************************************************
* Network Estimation Cycle (no plasticity)
**************************************************************************/

// Main control flow
//
void Network::CycleEstStatic() {
  // Check if computation is complete
  if (tsim >= tcomp) {
    // Reset network
    ResetNetwork();
    tsim = 0;
    tleap = 0;
    iter = 0;
    cadjprt[0] = 0;
    cadjprt[1] = 0;
    prtiter = 0;
    
    // One second
    tcomp = tstep * 1000;
    
    for (std::size_t i = 0; i < pngwin.size(); ++i) {
      for (std::size_t p = 0; p < pngwin[i].size(); ++p) {
        pngwin[i][p].clear();
      }
    }
    
    // Coordination after reset
    if (evalidx < ntrials) {
      // Reduce PNG information
      CkCallback *cb = new CkCallback(CkIndex_Network::SaveEstimate(NULL), 0, thisProxy);
      contribute(pnglog.size()*sizeof(event_t), pnglog.data(), net_event, *cb);

      pnglog.clear();
      ++evalidx;
    }
    else {
      // return control to main
      contribute(0, NULL, CkReduction::nop);
    }
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
  // Simulate next cycle
  else {
    // Display iteration information
    if (tsim == 0 && prtidx == 0) {
      CkPrintf("  Estimating iteration %" PRIidx "\n", evalidx);
      //CkPrintf("    Simulating time %" PRIrealsec " seconds\n", ((real_t) tsim)/(TICKS_PER_MS*1000));
    }
    
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
    
    // Check for periodic events
    if (tsim >= tleap) {
      std::vector<event_t>::iterator event = evtleap.begin();
      // Compute periodic events
      while (event != evtleap.end() && event->diffuse <= tsim) {
        // Set netmodel index
        idx_t n = event->source;
        // Loop through all models
        for (std::size_t m = 0; m < leapidx[n].size(); ++m) {
          event->index = leapidx[n][m][1];
          if (event->index) {
            netmodel[n]->Hop(*event, state[leapidx[n][m][0]], stick[leapidx[n][m][0]], edgaux[n][vtxmodidx[leapidx[n][m][0]]]);
          }
          else {
            netmodel[n]->Hop(*event, state[leapidx[n][m][0]], stick[leapidx[n][m][0]], vtxaux[leapidx[n][m][0]]);
          }
        }
        // Update timing
        event->diffuse += (tick_t)(event->data*TICKS_PER_MS);
        ++event;
      }
      std::sort(evtleap.begin(), evtleap.end());
      tleap = evtleap.front().diffuse;
    }
    
    // Perform computation
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
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
        }
        // vertex events
        else {
          netmodel[vtxmodidx[i]]->Hop(*event, state[i], stick[i], vtxaux[i]);
        }
        ++event;
      }

      // Polychronization
      for (std::size_t p = 0; p < pngs[i].size(); ++p) {
        // Pop out any old stamps
        while (!pngwin[i][p].empty()) {
          if (pngwin[i][p].front().diffuse + pnglen[i][p] <= tdrift) {
            pngwin[i][p].pop_front();
          }
          else {
            break;
          }
        }
        // Template event
        event_t pngevent;
        pngevent.diffuse = tdrift;
        pngevent.type = EVENT_GROUP;
        pngevent.data = 0.0;
        // Check for threshold number of stamps
        // TODO: Threshold based off of excitatory neurons only?
        if (pngwin[i][p].size() > (pngs[i][p].size() / 2)) {
          // Compute group activation
          int nactive = 0;
          std::deque<stamp_t>::iterator pngcan = pngwin[i][p].begin();
          for (std::size_t t = 0; t < pngs[i][p].size(); ++t) {
            while (pngcan != pngwin[i][p].end()) {
              // TODO: Refine the acceptable jitter
              if ((pngcan->diffuse + pnglen[i][p] - tdrift + 5 * TICKS_PER_MS == pngs[i][p][t].diffuse && pngcan->source == pngs[i][p][t].source) ||
                  (pngcan->diffuse + pnglen[i][p] - tdrift + 4 * TICKS_PER_MS == pngs[i][p][t].diffuse && pngcan->source == pngs[i][p][t].source) ||
                  (pngcan->diffuse + pnglen[i][p] - tdrift + 3 * TICKS_PER_MS == pngs[i][p][t].diffuse && pngcan->source == pngs[i][p][t].source) ||
                  (pngcan->diffuse + pnglen[i][p] - tdrift + 2 * TICKS_PER_MS == pngs[i][p][t].diffuse && pngcan->source == pngs[i][p][t].source) ||
                  (pngcan->diffuse + pnglen[i][p] - tdrift + TICKS_PER_MS == pngs[i][p][t].diffuse && pngcan->source == pngs[i][p][t].source) ||
                  (pngcan->diffuse + pnglen[i][p] - tdrift == pngs[i][p][t].diffuse && pngcan->source == pngs[i][p][t].source)) {
                ++nactive;
              }
              else if (pngcan->diffuse + pnglen[i][p] - tdrift > pngs[i][p][t].diffuse) {
                break;
              }
              ++pngcan;
            }
          }
          if (nactive > (pngs[i][p].size() / 2)) {
            CkPrintf("PNG %d, %d activated\n", vtxidx[i], p);
            // Record group activation
            pngevent.source = vtxidx[i];
            pngevent.index = p;
            pngevent.data = ((real_t)(pnglen[i][p]/TICKS_PER_MS));
            pnglog.push_back(pngevent);
            // Clear window for repeats
            pngwin[i][p].clear();
          }
        }
      }

      // Computation
      while (tdrift < tstop) {
        // Step through model drift (vertex)
        tdrift += netmodel[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], events);

        // Handle generated events (if any)
        // TODO: Conversion from edge indices to global (for individual output)
        if (events.size()) {
          for (std::size_t e = 0; e < events.size(); ++e) {
            // Get information
            idx_t target = events[e].source;
            idx_t index = events[e].index;
            // reindex to global
            events[e].source = vtxidx[i];
            // Record listed event
            if (evtloglist[events[e].type]) {
              evtlog.push_back(events[e]);
            }
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
                    netmodel[edgmodidx[i][j]]->Hop(events[e], state[i], stick[i], edgaux[edgmodidx[i][j]][vtxmodidx[i]]);
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
                netmodel[vtxmodidx[i]]->Hop(events[e], state[i], stick[i], vtxaux[i]);
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
      //CkAssert(event == event[i][evtday].end());
      evtcal[i][evtday].clear();
    }
    //CkPrintf("    Events on %d: %d\n", prtidx, nevent);

    // Send messages to entire network
    // TODO: Reduce communication due to monitoring
    mEvent *mevent = BuildEvent();
    thisProxy.CommStamp(mevent);

    // Increment simulated time
    tsim += tstep;
    
    // Store new records
    StoreRecord();
  }
}
