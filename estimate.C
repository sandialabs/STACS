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
extern /*readonly*/ idx_t intdisp;
extern /*readonly*/ idx_t intrec;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tepisode;
extern /*readonly*/ idx_t episodes;


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
void Network::InitEstCont(CProxy_Netdata cpdata) {
  // Set proxies
  netdata = cpdata;
  cyclepart = CkCallback(CkIndex_Network::CycleEstCont(), partidx, thisProxy);
  
  // Request network part from input
  netdata(fileidx).LoadNetwork(partidx,
      CkCallback(CkIndex_Network::LoadNetwork(NULL), partidx, thisProxy));
}

// Coordination with NetData chare array
//
void Network::InitEstEpis(CProxy_Netdata cpdata) {
  // Set proxies
  netdata = cpdata;
  cyclepart = CkCallback(CkIndex_Network::CycleEstEpis(), partidx, thisProxy);

  // Initialization
  teps = 0;
  epsidx = -1;
  
  // Request network part from input
  netdata(fileidx).LoadNetwork(partidx,
      CkCallback(CkIndex_Network::LoadNetwork(NULL), partidx, thisProxy));
}


/**************************************************************************
* Estimation Recording
**************************************************************************/

// Send Estimates for writing
//
void Network::SaveEstimate() {
  // Build record message for saving
  mRecord* mrecord = BuildRecord();
  netdata(fileidx).SaveRecord(mrecord);

  if (partidx == 0) {
    // Recording information
    event_t event;
    event.diffuse = 0;
    event.type = EVENT_GROUP;
    event.source = -1;
    event.index = iter;
    event.data = 0.0;
    grplog.push_back(event);
  }
  // Reduce group information
  contribute(grplog.size()*sizeof(event_t), grplog.data(), net_event, 
      CkCallback(CkIndex_Netdata::SaveEstimate(NULL), 0, netdata));
  grplog.clear();
  
  // Start a new cycle (data sent)
  cyclepart.send();
}

// Send Estimates for writing (final)
//
void Network::SaveFinalEstimate() {
  // Build record message for saving
  mRecord* mrecord = BuildRecord();
  netdata(fileidx).SaveFinalRecord(mrecord);
  
  if (partidx == 0) {
    // Recording information
    event_t event;
    event.diffuse = 0;
    event.type = EVENT_GROUP;
    event.source = -1;
    event.index = iter;
    event.data = 0.0;
    grplog.push_back(event);
  }
  // Reduce group information
  contribute(grplog.size()*sizeof(event_t), grplog.data(), net_event, 
      CkCallback(CkIndex_Netdata::SaveFinalEstimate(NULL), 0, netdata));
  grplog.clear();
}


/**************************************************************************
* Network Estimation Cycle (continuous)
**************************************************************************/

// Main control flow
//
void Network::CycleEstCont() {
  // Check if simulation time is complete
  if (tsim >= tmax) {
    // return control to main
    contribute(0, NULL, CkReduction::nop);
  }
  // Recording
  else if (iter == reciter) {
    // Bookkeeping
    reciter += intrec;

    // Send records
    thisProxy(partidx).SaveEstimate();
  }
#ifdef STACS_WITH_YARP
  // Synchronization from RPC
  else if (iter == synciter) {
    // Bookkkeeping
    synciter = IDX_T_MAX;

    // Display synchronization information
    if (partidx == 0) {
      CkPrintf("  Synchronizing at iteration %" PRIidx "\n", iter);
    }

    // move control to sychronization callback
    contribute(0, NULL, CkReduction::nop);
  }
#endif
  // Simulate next cycle
  else {
    // Display iteration information
    if (iter >= dispiter && partidx == 0) {
      dispiter += intdisp;
      CkPrintf("  Estimating iteration %" PRIidx "\n", iter);
      //CkPrintf("    Simulating time %" PRIrealsec " seconds\n", ((real_t) tsim)/(TICKS_PER_MS*1000));
    }
    
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
        }
        // vertex events
        else {
          model[vtxmodidx[i]]->Jump(*event, state[i], stick[i], vtxaux[i]);
        }
        ++event;
      }

      // Polychronization
      EstimateGroup(i);

      // Computation
      while (tdrift < tstop) {
        // Step through model drift (vertex)
        tdrift += model[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], events);

        // Handle generated events (if any)
        if (events.size()) {
          for (std::size_t e = 0; e < events.size(); ++e) {
            HandleEvent(events[e], i);
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

    // Send messages to entire network
    // TODO: Reduce communication due to monitoring
    mEvent *mevent = BuildEvent();
    thisProxy.CommStamp(mevent);

    // Increment simulated time
    tsim += tstep;

    // Add new records
    AddRecord();
  }
}


/**************************************************************************
* Network Estimation Cycle (episodic)
**************************************************************************/

// Main control flow
//
void Network::CycleEstEpis() {
  // Check if episode is complete
  if (tsim >= teps) {
    // Check if all episodes are complete
    if (++epsidx >= episodes) {
      // return control to main
      contribute(0, NULL, CkReduction::nop);
    }
    else {
      teps = tsim + tepisode;

      // Rerun any episodic models
      for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
        model[vtxmodidx[i]]->Rerun(state[i][0], stick[i][0]);
      }

      // Clear out groups from current episode
      for (std::size_t i = 0; i < grpwindow.size(); ++i) {
        for (std::size_t p = 0; p < grpwindow[i].size(); ++p) {
          grpwindow[i][p].clear();
        }
      }

      // Send records and estimates
      thisProxy(partidx).SaveEstimate();
    }
  }
#ifdef STACS_WITH_YARP
  // Synchronization from RPC
  else if (iter == synciter) {
    // Bookkkeeping
    synciter = IDX_T_MAX;

    // Display synchronization information
    if (partidx == 0) {
      CkPrintf("  Synchronizing at iteration %" PRIidx "\n", iter);
    }

    // move control to sychronization callback
    contribute(0, NULL, CkReduction::nop);
  }
#endif
  // Simulate next cycle
  else {
    // Display iteration information
    if (iter == dispiter && partidx == 0) {
      dispiter += intdisp;
      CkPrintf("  Estimating episode %" PRIidx "\n", epsidx);
    }

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
        }
        // vertex events
        else {
          model[vtxmodidx[i]]->Jump(*event, state[i], stick[i], vtxaux[i]);
        }
        ++event;
      }

      // Polychronization
      EstimateGroup(i);

      // Computation
      while (tdrift < tstop) {
        // Step through model drift (vertex)
        tdrift += model[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], events);

        // Handle generated events (if any)
        if (events.size()) {
          for (std::size_t e = 0; e < events.size(); ++e) {
            HandleEvent(events[e], i);
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

    // Send messages to entire network
    // TODO: Reduce communication due to monitoring
    //       Make two communication channels
    //       One for events, and the other for stamps
    mEvent *mevent = BuildEvent();
    thisProxy.CommStamp(mevent);

    // Increment simulated time
    tsim += tstep;
    
    // Add new records
    AddRecord();
  }
}


/**************************************************************************
* Network Estimation Helpers
**************************************************************************/

// Estimate any group activations
//
void Network::EstimateGroup(idx_t i) {
  for (std::size_t p = 0; p < grpstamps[i].size(); ++p) {
    // Pop out any old stamps
    while (!grpwindow[i][p].empty()) {
      if (grpwindow[i][p].front().diffuse + grpdur[i][p] <= tsim) {
        grpwindow[i][p].pop_front();
      }
      else {
        break;
      }
    }
    // Template event
    event_t event;
    event.diffuse = tsim;
    event.type = EVENT_GROUP;
    event.data = 0.0;
    // Check for threshold number of stamps
    // TODO: Threshold based off of excitatory neurons only?
    if (grpwindow[i][p].size() > (grpstamps[i][p].size() / 2)) {
      // Compute group activation
      int nactive = 0;
      std::deque<stamp_t>::iterator stamp = grpwindow[i][p].begin();
      for (std::size_t t = 0; t < grpstamps[i][p].size(); ++t) {
        while (stamp != grpwindow[i][p].end()) {
          // TODO: Refine the acceptable jitter
          if ((stamp->diffuse + grpdur[i][p] - tsim + 5 * TICKS_PER_MS == grpstamps[i][p][t].diffuse && stamp->source == grpstamps[i][p][t].source) ||
              (stamp->diffuse + grpdur[i][p] - tsim + 4 * TICKS_PER_MS == grpstamps[i][p][t].diffuse && stamp->source == grpstamps[i][p][t].source) ||
              (stamp->diffuse + grpdur[i][p] - tsim + 3 * TICKS_PER_MS == grpstamps[i][p][t].diffuse && stamp->source == grpstamps[i][p][t].source) ||
              (stamp->diffuse + grpdur[i][p] - tsim + 2 * TICKS_PER_MS == grpstamps[i][p][t].diffuse && stamp->source == grpstamps[i][p][t].source) ||
              (stamp->diffuse + grpdur[i][p] - tsim + TICKS_PER_MS == grpstamps[i][p][t].diffuse && stamp->source == grpstamps[i][p][t].source) ||
              (stamp->diffuse + grpdur[i][p] - tsim == grpstamps[i][p][t].diffuse && stamp->source == grpstamps[i][p][t].source)) {
            ++nactive;
          }
          else if (stamp->diffuse + grpdur[i][p] - tsim > grpstamps[i][p][t].diffuse) {
            break;
          }
          ++stamp;
        }
      }
      if (nactive > (grpstamps[i][p].size() / 2)) {
        CkPrintf("PNG %d, %d activated\n", vtxidx[i], p);
        // Record group activation
        event.source = vtxidx[i];
        event.index = p;
        event.data = ((real_t)(grpdur[i][p]/TICKS_PER_MS));
        grplog.push_back(event);
        // Clear window for repeats
        grpwindow[i][p].clear();
      }
    }
  }
}
