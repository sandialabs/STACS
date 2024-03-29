/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
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
* Charm++ Reduction
**************************************************************************/

// Reduction for type idx_t
//
extern CkReduction::reducerType max_idx;


/**************************************************************************
* Network Simulation Cycle
**************************************************************************/

// Main control flow
//
void Network::CycleSimGPU() {
  // Check if simulation time is complete
  if (tsim >= tmax && !episodic) {
    // return control to main
    contribute(0, NULL, CkReduction::nop);
  }
  // Recording
  else if (iter == reciter && !episodic) {
    // Bookkeeping
    reciter += intrec;

    // Send records
    thisProxy(prtidx).SaveRecord();
  }
  // Check if episode is complete
  else if (tsim >= teps && episodic) {
    // Check if all episodes are complete
    if (++epsidx >= episodes) {
      // return control to main
      contribute(0, NULL, CkReduction::nop);
    }
    else {
      teps = tsim + tepisode;
      
      // Renew any episodic models
      for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
        model[vtxmodidx[i]]->Renew(state[i][0], stick[i][0]);
      }
    
      // Start a new cycle (after checked data sent)
      thisProxy(prtidx).SaveRecord();
    }
  }
  // Saving
  else if (iter == saveiter) {
    // Bookkeeping
    saveiter += intsave;
    
    // Display checkpointing information
    if (prtidx == 0) {
      CkPrintf("  Saving network at iteration %" PRIidx "\n", iter);
    }

    // Checkpoint
    thisProxy(prtidx).SaveNetwork();
  }
#ifdef STACS_WITH_YARP
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
      if (prtidx == 0) {
        CkPrintf("  Synchronized at iteration %" PRIidx "\n", iter);
      }

      // move control to sychronization callback
      contribute(0, NULL, CkReduction::nop);
    }
  }
#endif
  // Simulate next cycle
  else {
    // Display iteration information
    if (iter >= dispiter && prtidx == 0) {
      dispiter += intdisp;
      if (episodic) {
        CkPrintf("  Simulating episode %" PRIidx "\n", epsidx);
      }
      else {
        CkPrintf("  Simulating iteration %" PRIidx "\n", iter);
        //CkPrintf("    Simulating time %" PRIrealsec " seconds\n", ((real_t) tsim)/(TICKS_PER_MS*1000));
      }
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
    if (tsim >= tleap) {
      LeapEvent();
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

    // Send messages to neighbors
    mEvent *mevent = BuildEvent();
    netcomm.CommEvent(mevent);

    // Increment simulated time
    tsim += tstep;

    // Add new records
    AddRecord();
    
    // Increment iteration
    ++iter;
  }
}

