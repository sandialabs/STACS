/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"


/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ tick_t tcheck;
extern /*readonly*/ tick_t trecord;
extern /*readonly*/ tick_t tdisplay;
extern /*readonly*/ idx_t nevtday;


/**************************************************************************
* Network Simulation Initialization
**************************************************************************/

// Coordination with NetData chare array
//
void Network::InitSimPlastic(CProxy_Netdata cpdat) {
  // Set proxies
  netdata = cpdat;
  cbcycleprt = CkCallback(CkIndex_Network::CycleSimPlastic(), prtidx, thisProxy);
  
  // Request network part from input
  CkCallback *cb = new CkCallback(CkIndex_Network::LoadNetwork(NULL), prtidx, thisProxy);
  netdata(datidx).LoadNetwork(prtidx, *cb);
}

// Coordination with NetData chare array
//
void Network::InitSimStatic(CProxy_Netdata cpdat) {
  // Set proxies
  netdata = cpdat;
  cbcycleprt = CkCallback(CkIndex_Network::CycleSimStatic(), prtidx, thisProxy);

  // Request network part from input
  CkCallback *cb = new CkCallback(CkIndex_Network::LoadNetwork(NULL), prtidx, thisProxy);
  netdata(datidx).LoadNetwork(prtidx, *cb);
}


/**************************************************************************
* Network Simulation Cycle (with plasticity)
**************************************************************************/

// Main control flow
//
void Network::CycleSimPlastic() {
  // Check if simulation time is complete
  if (tsim >= tmax) {
    // return control to main
    contribute(0, NULL, CkReduction::nop);
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
  // Checkpointing
  else if (iter == checkiter) {
    // Bookkeeping
    checkiter = checkiter + (idx_t) (tcheck/tstep);
    
    // Display checkpointing information
    if (prtidx == 0) {
      CkPrintf("  Checkpointing at iteration %" PRIidx "\n", iter);
    }

    // Checkpoint
    thisProxy(prtidx).SaveNetwork();
  }
  // Recording
  else if (iter == reciter) {
    // Bookkeeping
    reciter = reciter + (idx_t) (trecord/tstep);

    // Send records
    thisProxy(prtidx).SaveRecord();
  }
  // Simulate next cycle
  else {
    // Display iteration information
    if (tsim >= tdisp && prtidx == 0) {
      tdisp = tsim + tdisplay;
      CkPrintf("  Simulating iteration %" PRIidx "\n", iter);
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
            netmodel[modidx]->Jump(*event, state[repidx[modidx][i][0]], stick[repidx[modidx][i][0]], edgaux[modidx][vtxmodidx[repidx[modidx][i][0]]]);
          }
          else {
            netmodel[modidx]->Jump(*event, state[repidx[modidx][i][0]], stick[repidx[modidx][i][0]], vtxaux[repidx[modidx][i][0]]);
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
          netmodel[edgmodidx[i][event->index-1]]->Jump(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
        }
        // vertex events
        else {
          netmodel[vtxmodidx[i]]->Jump(*event, state[i], stick[i], vtxaux[i]);
        }
        ++event;
      }

      // Computation
      while (tdrift < tstop) {
        // Step through model drift (vertex)
        tdrift += netmodel[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], evtlog);

        // Handle generated events (if any)
        // TODO: Conversion from edge indices to global (for individual output)
        if (evtlog.size()) {
          for (std::size_t e = 0; e < evtlog.size(); ++e) {
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
                    netmodel[edgmodidx[i][j]]->Jump(evtlog[e], state[i], stick[i], edgaux[edgmodidx[i][j]][vtxmodidx[i]]);
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
                netmodel[vtxmodidx[i]]->Jump(evtlog[e], state[i], stick[i], vtxaux[i]);
              }
              else {
                evtcol[i].push_back(evtlog[e]);
              }
            }
            // Record listed event
            if (recevtlist[evtlog[e].type]) {
              // reindex to global
              evtlog[e].source = vtxidx[i];
              evtlog[e].index = index;
              recevt.push_back(evtlog[e]);
            }
          }
          // clear log for next time
          evtlog.clear();
        }
        
        // Perform events up to tdrift
        while (event != evtcal[i][evtday].end() && event->diffuse <= tdrift) {
          // edge events
          if (event->index) {
            netmodel[edgmodidx[i][event->index-1]]->Jump(*event, state[i], stick[i], edgaux[edgmodidx[i][event->index-1]][vtxmodidx[i]]);
          }
          // vertex events
          else {
            netmodel[vtxmodidx[i]]->Jump(*event, state[i], stick[i], vtxaux[i]);
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
    netgroup.CommEvent(mevent);

    // Increment simulated time
    tsim += tstep;

    // Store new records
    StoreRecord();
  }
}


/**************************************************************************
* Network Simulation Cycle (no plasticity)
**************************************************************************/

// Main control flow
//
void Network::CycleSimStatic() {
  // Check if simulation time is complete
  if (tsim >= tmax) {
    // return control to main
    contribute(0, NULL, CkReduction::nop);
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
  // Recording
  else if (iter == reciter) {
    // Bookkeeping
    reciter = reciter + (idx_t) (trecord/tstep);

    // Send records
    thisProxy(prtidx).SaveRecord();
  }
  // Simulate next cycle
  else {
    // Display iteration information
    if (tsim >= tdisp && prtidx == 0) {
      tdisp = tsim + tdisplay;
      CkPrintf("  Simulating iteration %" PRIidx "\n", iter);
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

      // Computation
      while (tdrift < tstop) {
        // Step through model drift (vertex)
        tdrift += netmodel[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], evtlog);

        // Handle generated events (if any)
        // TODO: Conversion from edge indices to global (for individual output)
        if (evtlog.size()) {
          for (std::size_t e = 0; e < evtlog.size(); ++e) {
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
            // Record listed event
            if (recevtlist[evtlog[e].type]) {
              // reindex to global
              evtlog[e].source = vtxidx[i];
              evtlog[e].index = index;
              recevt.push_back(evtlog[e]);
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
      //CkAssert(event == event[i][evtday].end());
      evtcal[i][evtday].clear();
    }
    //CkPrintf("    Events on %d: %d\n", prtidx, nevent);

    // Send messages to neighbors
    mEvent *mevent = BuildEvent();
    netgroup.CommEvent(mevent);

    // Increment simulated time
    tsim += tstep;

    // Store new records
    StoreRecord();
  }
}
