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
extern /*readonly*/ idx_t equeue;


/**************************************************************************
* Reduction for PNGs
**************************************************************************/

CkReduction::reducerType net_png;
/*initnode*/
void registerNetPNG(void) {
  net_png = CkReduction::addReducer(netPNG);
}

CkReductionMsg *netPNG(int nMsg, CkReductionMsg **msgs) {
  std::vector<png_t> ret;
  ret.clear();
  for (int i = 0; i < nMsg; i++) {
    for (std::size_t j = 0; j < msgs[i]->getSize()/sizeof(png_t); ++j) {
      // Extract data and reduce 
      ret.push_back(*((png_t *)msgs[i]->getData() + j));
    }
  }
  return CkReductionMsg::buildNew(ret.size()*sizeof(png_t), ret.data());
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
      if (netmodel[vtxmodidx[mother->second]]->getPNGMod()) {
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
            if (netmodel[edgmodidx[i][j]]->getPNGMod()) {
              // TODO: Based off of spiking property of the 
              //       netmodel instead of just active models
              if (stick[i][j+1].size()) {
                anchor.push_back(j);
              }
            }
          }

          // PNG combinatorics
          for (idx_t j0 = 0; j0 < anchor.size(); ++j0) {
            for (idx_t j1 = j0; j1 < anchor.size(); ++j1) {
              for (idx_t j2 = j1; j2 < anchor.size(); ++j2) {
                // Test for spiking of mother neuron
                // assuming perfect timing of anchor
                netmodel[vtxmodidx[i]]->Reset(state[i][0], stick[i][0]);
                event_t evtpre;
                evtpre.diffuse = 0;
                evtpre.type = EVENT_SPIKE;
                evtpre.source = i;
                evtpre.data = 0.0;
                evtpre.index = anchor[j0]+1;
                netmodel[edgmodidx[i][anchor[j0]]]->Hop(evtpre, state[i], stick[i], edgaux[edgmodidx[i][anchor[j0]]][vtxmodidx[i]]);
                evtpre.index = anchor[j1]+1;
                netmodel[edgmodidx[i][anchor[j1]]]->Hop(evtpre, state[i], stick[i], edgaux[edgmodidx[i][anchor[j1]]][vtxmodidx[i]]);
                evtpre.index = anchor[j2]+1;
                netmodel[edgmodidx[i][anchor[j2]]]->Hop(evtpre, state[i], stick[i], edgaux[edgmodidx[i][anchor[j2]]][vtxmodidx[i]]);
                tick_t tdrift = 0;
                tick_t tstop = tstep * 4; // Strongly spiking triplets only
                while (tdrift < tstop) {
                  // Step through model drift (vertex)
                  tdrift += netmodel[vtxmodidx[i]]->Step(tdrift, tstop - tdrift, state[i][0], stick[i][0], evtlog);
                }
                // Potential PNG only if mother vertex spiked
                if (evtlog.size() && evtlog[0].type == EVENT_SPIKE) {
                  evtlog.clear();
                  std::vector<event_t> pngseed;
                  pngseed.clear();
                  evtpre.type = EVENT_SPIKE;
                  evtpre.data = 0.0;
                  evtpre.diffuse = stick[i][anchor[j0]+1][0];
                  evtpre.source = adjcy[i][anchor[j0]];
                  evtpre.index = adjcy[i][anchor[j0]];
                  pngseed.push_back(evtpre);
                  evtpre.diffuse = stick[i][anchor[j1]+1][0];
                  evtpre.source = adjcy[i][anchor[j1]];
                  evtpre.index = adjcy[i][anchor[j1]];
                  pngseed.push_back(evtpre);
                  evtpre.diffuse = stick[i][anchor[j2]+1][0];
                  evtpre.source = adjcy[i][anchor[j2]];
                  evtpre.index = adjcy[i][anchor[j2]];
                  pngseed.push_back(evtpre);
                  // correctly order the timing
                  std::sort(pngseed.begin(), pngseed.end());
                  pngseed[2].diffuse = pngseed[2].diffuse - pngseed[2].diffuse;
                  pngseed[1].diffuse = pngseed[2].diffuse - pngseed[1].diffuse;
                  pngseed[0].diffuse = pngseed[2].diffuse - pngseed[0].diffuse;
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

// Coordination for PNG computation
//
void Network::EvalPNG(CkReductionMsg *msg) {
  // Add to png candidate
  for (std::size_t i = 0; i < (msg->getSize())/sizeof(png_t); ++i) {
    pngcan.push_back(*((png_t *)msg->getData()+i));
  }
  delete msg;

  //if (pngcan.size() > minpngsize) {
  //  pngs[i].push_back(pngcan);
  //}
    
  // Compute next PNG
  thisProxy.ComputePNG();
}

// Compute PNG (vertex control loop)
//
void Network::ComputePNG() {
  // Loop through PNG seeds
  if (ccomp < ncomp) {
    pngcan.clear();
    pngaux.clear();
    if (pngseeds.size()) {
      // Initialize candidate PNG
      pngcan.resize(pngseeds[ccomp].size());
      for (std::size_t i = 0; i < pngseeds[ccomp].size(); ++i) {
        pngcan[i].diffuse = pngseeds[ccomp][i].diffuse;
        pngcan[i].source = pngseeds[ccomp][i].source;
      }
      // Seed spikes for simulation
      mEvent *mevent = BuildPNGSeed(pngseeds[ccomp]);
      thisProxy.SeedPNG(mevent);
    }
    ++ccomp;
  }
  else {
    // Return control to main loop
    thisProxy(prtidx).FindPNG();
  }
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
    
    // Coordination after reset
    // Reduce PNG information
    CkCallback *cb = new CkCallback(CkIndex_Network::EvalPNG(NULL), evalidx, thisProxy);
    //contribute(0, NULL, CkReduction::nop, *cb);
    contribute(sizeof(png_t), pngaux.data(), net_png, *cb);
  }
#ifdef STACS_WITH_YARP
  // Synchronization from RPC
  else if (iter == synciter) {
    // Bookkkeeping
    synciter = IDX_T_MAX;

    // Display synchronization information
    if (prtidx == 0) {
      CkPrintf("  Synchronized at iteration %" PRIidx "\n", iter);
    }

    // move control to sychronization callback
    contribute(0, NULL, CkReduction::nop);
  }
#endif
  else {
    // Bookkeeping
    idx_t evtiter = iter%equeue;
    tick_t tstop = tsim + tstep;

    // Clear event buffer
    evtext.clear();
    idx_t nevent = 0;
    // Redistribute any events (on new year)
    if (evtiter == 0) {
      RedisEvent();
    }
    
    // Check for repeating events
    if (tsim >= trep) {
      std::vector<event_t>::iterator evt = repevt.begin();
      // Compute periodic events
      while (evt != repevt.end() && evt->diffuse <= tsim) {
        // Set temporary model index
        idx_t modidx = evt->index;
        // Loop through all models
        for (std::size_t i = 0; i < repidx[modidx].size(); ++i) {
          evt->index = repidx[modidx][i][1];
          if (evt->index) {
            netmodel[modidx]->Hop(*evt, state[repidx[modidx][i][0]], stick[repidx[modidx][i][0]], edgaux[modidx][vtxmodidx[repidx[modidx][i][0]]]);
          }
          else {
            netmodel[modidx]->Hop(*evt, state[repidx[modidx][i][0]], stick[repidx[modidx][i][0]], vtxaux[repidx[modidx][i][0]]);
          }
        }
        // Return model index
        evt->index = modidx;
        // Update timing
        evt->diffuse += ((tick_t) evt->data)* TICKS_PER_MS;
        ++evt;
      }
      std::sort(repevt.begin(), repevt.end());
      trep = repevt[0].diffuse;
    }
    
    // Perform computation
    for (std::size_t i = 0; i < vtxmodidx.size(); ++i) {
      if (netmodel[vtxmodidx[i]]->getActive() == false) {
        event[i][evtiter].clear();
        continue;
      }
      // Timing
      tick_t tdrift = tsim;

      // Sort events
      std::sort(event[i][evtiter].begin(), event[i][evtiter].end());
      nevent += event[i][evtiter].size();

      // Perform events starting at beginning of step
      std::vector<event_t>::iterator evt = event[i][evtiter].begin();
      while (evt != event[i][evtiter].end() && evt->diffuse <= tdrift) {
        // edge events
        if (evt->index) {
          netmodel[edgmodidx[i][evt->index-1]]->Hop(*evt, state[i], stick[i], edgaux[edgmodidx[i][evt->index-1]][vtxmodidx[i]]);
        }
        // vertex events
        else {
          netmodel[vtxmodidx[i]]->Hop(*evt, state[i], stick[i], vtxaux[i]);
        }
        ++evt;
      }
      
      // TODO: Method to grow PNGs (asynchronously?)
      //       Also a spike-timing relationship graph

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
              if ((evtlog[e].diffuse - tsim - tstep)/tstep < equeue) {
                for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
                  if (edgmodidx[i][j]) {
                    evtlog[e].index = j+1;
                    event[i][(evtlog[e].diffuse/tstep)%equeue].push_back(evtlog[e]);
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
                    evtaux[i].push_back(evtlog[e]);
                  }
                }
              }
            }
            // Local event (singlecast to vertex)
            if (target & LOCAL_VERTEX) {
              // vertex to itself
              evtlog[e].source = -vtxidx[i]-1; // negative source indicates local event
              evtlog[e].index = 0;
              if ((evtlog[e].diffuse - tsim - tstep)/tstep < equeue) {
                event[i][(evtlog[e].diffuse/tstep)%equeue].push_back(evtlog[e]);
              }
              else if (evtlog[e].diffuse < tsim + tstep) {
                // Jump now
                netmodel[vtxmodidx[i]]->Hop(evtlog[e], state[i], stick[i], vtxaux[i]);
              }
              else {
                evtaux[i].push_back(evtlog[e]);
              }
            }
          }
          // clear log for next time
          evtlog.clear();
        }
        
        // Perform events up to tdrift
        while (evt != event[i][evtiter].end() && evt->diffuse <= tdrift) {
          // edge events
          if (evt->index) {
            netmodel[edgmodidx[i][evt->index-1]]->Hop(*evt, state[i], stick[i], edgaux[edgmodidx[i][evt->index-1]][vtxmodidx[i]]);
          }
          // vertex events
          else {
            netmodel[vtxmodidx[i]]->Hop(*evt, state[i], stick[i], vtxaux[i]);
          }
          ++evt;
        }
      }

      // Clear event queue
      event[i][evtiter].clear();
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

// Seed events for PNG computation
//
void Network::SeedPNG(mEvent *msg) {
  // Event prototype
  event_t evtpre;
  tick_t evtdif;

  // Distribute events
  for (std::size_t i = 0; i < msg->nevent; ++i) {
    // Fill in prototype
    evtdif = msg->diffuse[i];
    evtpre.type = msg->type[i];
    evtpre.source = msg->source[i];
    evtpre.data = msg->data[i];
    // Determine local event target(s)
    // If index == source (multicast to edges)
    // Find target mapping from source
    std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>>::iterator targets = adjmap.find(msg->source[i]);
    if (targets != adjmap.end()) {
      for (std::vector<std::array<idx_t, 2>>::iterator target = targets->second.begin(); target != targets->second.end(); ++target) {
        evtpre.diffuse = evtdif + stick[(*target)[0]][(*target)[1]][0]; // delay always first stick of edge
        evtpre.index = (*target)[1];
        // Add to event queue or spillover
        if (evtpre.diffuse/tstep < equeue) {
          event[(*target)[0]][(evtpre.diffuse/tstep)%equeue].push_back(evtpre);
        }
        else {
          evtaux[(*target)[0]].push_back(evtpre);
        }
      }
    }
  }
  delete msg;

  // Start cycle after seeding events
  thisProxy(prtidx).CyclePNG();
}
