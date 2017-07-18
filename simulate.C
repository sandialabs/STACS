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
extern /*readonly*/ idx_t intsave;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tepisode;
extern /*readonly*/ idx_t episodes;


/**************************************************************************
* Network Simulation Initialization
**************************************************************************/

// Coordination with NetData chare array
//
void Network::InitSimContPlas(CProxy_Netdata cpdata) {
  // Set proxies
  netdata = cpdata;
  cyclepart = CkCallback(CkIndex_Network::CycleSimContPlas(), partidx, thisProxy);

  // Request network part from input
  netdata(fileidx).LoadNetwork(partidx,
      CkCallback(CkIndex_Network::LoadNetwork(NULL), partidx, thisProxy));
}

// Coordination with NetData chare array
//
void Network::InitSimCont(CProxy_Netdata cpdata) {
  // Set proxies
  netdata = cpdata;
  cyclepart = CkCallback(CkIndex_Network::CycleSimCont(), partidx, thisProxy);

  // Request network part from input
  netdata(fileidx).LoadNetwork(partidx,
      CkCallback(CkIndex_Network::LoadNetwork(NULL), partidx, thisProxy));
}

// Coordination with NetData chare array
//
void Network::InitSimEpisPlas(CProxy_Netdata cpdata) {
  // Set proxies
  netdata = cpdata;
  cyclepart = CkCallback(CkIndex_Network::CycleSimEpisPlas(), partidx, thisProxy);

  // Initialization
  teps = 0;
  epsidx = -1;
  
  // Request network part from input
  netdata(fileidx).LoadNetwork(partidx,
      CkCallback(CkIndex_Network::LoadNetwork(NULL), partidx, thisProxy));
}

// Coordination with NetData chare array
//
void Network::InitSimEpis(CProxy_Netdata cpdata) {
  // Set proxies
  netdata = cpdata;
  cyclepart = CkCallback(CkIndex_Network::CycleSimEpis(), partidx, thisProxy);

  // Initialization
  teps = 0;
  epsidx = -1;
  
  // Request network part from input
  netdata(fileidx).LoadNetwork(partidx,
      CkCallback(CkIndex_Network::LoadNetwork(NULL), partidx, thisProxy));
}


/**************************************************************************
* Network Simulation Cycle (with plasticity)
**************************************************************************/

// Main control flow
//
void Network::CycleSimContPlas() {
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
    thisProxy(partidx).SaveRecord();
  }
  // Saving
  else if (iter == saveiter) {
    // Bookkeeping
    saveiter += intsave;
    
    // Display checkpointing information
    if (partidx == 0) {
      CkPrintf("  Saving network at iteration %" PRIidx "\n", iter);
    }

    // Checkpoint
    thisProxy(partidx).SaveNetwork();
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
      CkPrintf("  Simulating iteration %" PRIidx "\n", iter);
      //CkPrintf("    Simulating time %" PRIrealsec " seconds\n", ((real_t) tsim)/(TICKS_PER_MS*1000));
    }
    
    // Bookkeeping
    idx_t evtday = iter%nevtday;
    tick_t tstop = tsim + tstep;

    // Clear event buffer
    evtext.clear();
    //idx_t nevent = 0;
    // Redistribute any events (on new year)
    if (evtday == 0) {
      SortEventCalendar();
    }
    
    // Check for periodic events
    if (tsim >= tskip) {
      SkipEventPlas();
    }
    
    // Perform computation
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
      // Timing
      tick_t tdrift = tsim;

      // Sort events
      std::sort(evtcal[i][evtday].begin(), evtcal[i][evtday].end());
      //nevent += evtcal[i][evtday].size();

      // Perform events starting at beginning of step
      std::vector<event_t>::iterator event = evtcal[i][evtday].begin();
      while (event != evtcal[i][evtday].end() && event->diffuse <= tdrift) {
        // edge events
        if (event->index) {
          model[edgmodidx[i][event->index-1]]->Leap(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
        }
        // vertex events
        else {
          model[vtxmodidx[i]]->Leap(*event, state[i], stick[i], vtxaux[i]);
        }
        ++event;
      }

      // Computation
      while (tdrift < tstop) {
        // Step through model drift (vertex)
        tdrift += model[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], events);

        // Handle generated events (if any)
        if (events.size()) {
          for (std::size_t e = 0; e < events.size(); ++e) {
            HandleEventPlas(events[e], i);
          }
          // clear log for next time
          events.clear();
        }
        
        // Perform events up to tdrift
        while (event != evtcal[i][evtday].end() && event->diffuse <= tdrift) {
          // edge events
          if (event->index) {
            model[edgmodidx[i][event->index-1]]->Leap(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
          }
          // vertex events
          else {
            model[vtxmodidx[i]]->Leap(*event, state[i], stick[i], vtxaux[i]);
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

    // Add new records
    AddRecord();
  }
}


/**************************************************************************
* Network Simulation Cycle (no plasticity)
**************************************************************************/

// Main control flow
//
void Network::CycleSimCont() {
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
    thisProxy(partidx).SaveRecord();
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
      CkPrintf("  Simulating iteration %" PRIidx "\n", iter);
      //CkPrintf("    Simulating time %" PRIrealsec " seconds\n", ((real_t) tsim)/(TICKS_PER_MS*1000));
    }
    
    // Bookkeeping
    idx_t evtday = iter%nevtday;
    tick_t tstop = tsim + tstep;

    // Clear event buffer
    evtext.clear();
    //idx_t nevent = 0;
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
      //nevent += evtcal[i][evtday].size();

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
      //CkAssert(event == event[i][evtday].end());
      evtcal[i][evtday].clear();
    }
    //CkPrintf("    Events on %d: %d\n", partidx, nevent);

    // Send messages to neighbors
    mEvent *mevent = BuildEvent();
    netcomm.CommEvent(mevent);

    // Increment simulated time
    tsim += tstep;

    // Add new records
    AddRecord();
  }
}


/**************************************************************************
* Network Training Episodes Cycle (plasticity)
**************************************************************************/

// Main control flow
//
void Network::CycleSimEpisPlas() {
  // Check if episode is complete
  if (tsim >= teps) {
    // Check if all episodes are complete
    if (++epsidx >= episodes) {
      // return control to main
      contribute(0, NULL, CkReduction::nop);
    }
    else {
      teps = tsim + tepisode;
    
      // Start a new cycle (after checked data sent)
      thisProxy(partidx).SaveRecord();
    }
  }
  // Saving
  else if (iter == saveiter) {
    // Bookkeeping
    saveiter += intsave;
    
    // Display checkpointing information
    if (partidx == 0) {
      CkPrintf("  Saving network at iteration %" PRIidx "\n", iter);
    }

    // Checkpoint
    thisProxy(partidx).SaveNetwork();
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
      CkPrintf("  Simulating episode %" PRIidx "\n", epsidx);
    }
    
    // Bookkeeping
    idx_t evtday = iter%nevtday;
    tick_t tstop = tsim + tstep;

    // Clear event buffer
    evtext.clear();
    //idx_t nevent = 0;
    // Redistribute any events (on new year)
    if (evtday == 0) {
      SortEventCalendar();
    }
    
    // Check for periodic events
    if (tsim >= tskip) {
      SkipEventPlas();
    }
    
    // Perform computation
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
      // Timing
      tick_t tdrift = tsim;

      // Sort events
      std::sort(evtcal[i][evtday].begin(), evtcal[i][evtday].end());
      //nevent += evtcal[i][evtday].size();

      // Perform events starting at beginning of step
      std::vector<event_t>::iterator event = evtcal[i][evtday].begin();
      while (event != evtcal[i][evtday].end() && event->diffuse <= tdrift) {
        // edge events
        if (event->index) {
          model[edgmodidx[i][event->index-1]]->Leap(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
        }
        // vertex events
        else {
          model[vtxmodidx[i]]->Leap(*event, state[i], stick[i], vtxaux[i]);
        }
        ++event;
      }

      // Computation
      while (tdrift < tstop) {
        // Step through model drift (vertex)
        tdrift += model[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], events);

        // Handle generated events (if any)
        if (events.size()) {
          for (std::size_t e = 0; e < events.size(); ++e) {
            HandleEventPlas(events[e], i);
          }
          // clear log for next time
          events.clear();
        }
        
        // Perform events up to tdrift
        while (event != evtcal[i][evtday].end() && event->diffuse <= tdrift) {
          // edge events
          if (event->index) {
            model[edgmodidx[i][event->index-1]]->Leap(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
          }
          // vertex events
          else {
            model[vtxmodidx[i]]->Leap(*event, state[i], stick[i], vtxaux[i]);
          }
          ++event;
        }
      }

      // Clear event queue
      //CkAssert(event == event[i][evtday].end());
      evtcal[i][evtday].clear();
    }
    //CkPrintf("    Events on %d: %d\n", partidx, nevent);

    // Send messages to neighbors
    mEvent *mevent = BuildEvent();
    netcomm.CommEvent(mevent);

    // Increment simulated time
    tsim += tstep;
    
    // Add new records
    AddRecord();
  }
}


/**************************************************************************
* Network Training Episodes Cycle (no plasticity)
**************************************************************************/

// Main control flow
//
void Network::CycleSimEpis() {
  // Check if episode is complete
  if (tsim >= teps) {
    // Check if simulation episodes are complete
    if (++epsidx >= episodes) {
      // return control to main
      contribute(0, NULL, CkReduction::nop);
    }
    else {
      teps = tsim + tepisode;
    
      // Start a new cycle (after checked data sent)
      thisProxy(partidx).SaveRecord();
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
      CkPrintf("  Simulating episode %" PRIidx "\n", epsidx);
    }

    // Bookkeeping
    idx_t evtday = iter%nevtday;
    tick_t tstop = tsim + tstep;

    // Clear event buffer
    evtext.clear();
    //idx_t nevent = 0;
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
      //nevent += evtcal[i][evtday].size();

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
      //CkAssert(event == event[i][evtday].end());
      evtcal[i][evtday].clear();
    }
    //CkPrintf("    Events on %d: %d\n", partidx, nevent);

    // Send messages to neighbors
    mEvent *mevent = BuildEvent();
    netcomm.CommEvent(mevent);

    // Increment simulated time
    tsim += tstep;
    
    // Add new records
    AddRecord();
  }
}
