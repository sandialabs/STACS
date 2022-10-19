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
extern /*readonly*/ idx_t intbal;
extern /*readonly*/ idx_t intsave;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tepisode;
extern /*readonly*/ idx_t episodes;


/**************************************************************************
* Charm++ Reduction
**************************************************************************/

// Reduction for type idx_t
//
CkReduction::reducerType max_idx;
// Initnode
void registerMaxIdx(void) {
  max_idx = CkReduction::addReducer(maxIdx);
}
// Reduction function
CkReductionMsg *maxIdx(int nMsg, CkReductionMsg **msgs) {
  // Initialize to 0
  idx_t ret = 0;
  for (int i = 0; i < nMsg; ++i) {
    // Sanity check
    CkAssert(msgs[i]->getSize() == sizeof(idx_t));
    // Extract data and reduce 
    idx_t m = *(idx_t *)msgs[i]->getData();
    ret = std::max(ret,m);
  }
  // Return maximum idx_t
  return CkReductionMsg::buildNew(sizeof(idx_t),&ret);
}


/**************************************************************************
* Network Simulation Initialization
**************************************************************************/

// Coordination with NetData chare array
//
void Network::InitSim(CProxy_Netdata cpdata) {
  // Set proxies
  netdata = cpdata;
  cyclepart = CkCallback(CkIndex_Network::CycleSim(), thisProxy(prtidx));

  // Request network part from input
  netdata(datidx).LoadNetwork(prtidx,
      CkCallback(CkIndex_Network::LoadNetwork(NULL), thisProxy(prtidx)));
}

// Return after repartitioning
//
void Network::ContSim() {
  netdata(datidx).LoadNetwork(prtidx,
      CkCallback(CkIndex_Network::ReloadNetwork(NULL), thisProxy(prtidx)));
}


/**************************************************************************
* Network Simulation Cycle
**************************************************************************/

// Main control flow
//
void Network::CycleSim() {
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
    
      if (iter >= reciter) {
        reciter += intrec;
        // Start a new cycle (after checked data sent)
        thisProxy(prtidx).SaveRecord();
      }
      else {
        // Start a new cycle (basically a continue)
        cyclepart.send();
      }
    }
  }
  // Repartitioning / migrating vertices
  else if (iter == baliter && loadbal) {
    // Bookkeeping
    baliter += intbal;

    // Send part
    thisProxy(prtidx).RebalNetwork();
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
    //idx_t nevent = 0;
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
    //CkPrintf("    Events on %d: %d\n", prtidx, nevent);

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

