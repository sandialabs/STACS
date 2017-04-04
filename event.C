/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"


/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ idx_t nevtday;


/**************************************************************************
* Network Event Communication
**************************************************************************/

// Simple go-ahead message
//
void Network::GoAhead(mGo *msg) {
  // Increment coordination
  ++cadjprt[(prtiter + (msg->iter - iter))%2];
  delete msg;

  // Start next cycle
  if (cadjprt[prtiter] == nadjprt) {
    // Bookkeepping
    cadjprt[prtiter] = 0;
    prtiter = (prtiter+1)%2;

    // Increment iteration
    ++iter;

    // Start a new cycle
    //thisProxy(prtidx).CycleNetwork();
    cbcycleprt.send();
  }
}

// Multicast communication of events
//
void Network::CommEvent(mEvent *msg) {
  // Increment coordination
  ++cadjprt[(prtiter + (msg->iter - iter))%2];

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
    // If index != source (singlecast to edge)
    // If index < 0 (singlecast to vertex)
    if (msg->index[i] == msg->source[i]) {
      // Find target mapping from source
      std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>>::iterator targets = adjmap.find(msg->source[i]);
      if (targets != adjmap.end()) {
        for (std::vector<std::array<idx_t, 2>>::iterator target = targets->second.begin(); target != targets->second.end(); ++target) {
          event.diffuse = departure + stick[(*target)[0]][(*target)[1]][0]; // delay always first stick of edge
          event.index = (*target)[1];
          // Add to event queue or spillover
          if ((event.diffuse/tstep - msg->iter) < nevtday) {
            evtcal[(*target)[0]][(event.diffuse/tstep)%nevtday].push_back(event);
          }
          else if (event.diffuse/tstep < msg->iter) {
            evtcal[(*target)[0]][(msg->iter+1)%nevtday].push_back(event);
          }
          else {
            evtcol[(*target)[0]].push_back(event);
          }
        }
      }
    }
    else if (msg->index[i] < 0) {
      // Find local target
      std::unordered_map<idx_t, idx_t>::iterator target = vtxmap.find(-(msg->index[i]+1));
      if (target != vtxmap.end()) {
        event.diffuse = departure; // direct events to vertices have no edge delay
        event.index = 0;
        // Add to event queue or spillover
        if ((event.diffuse/tstep - msg->iter) < nevtday) {
          evtcal[target->second][(event.diffuse/tstep)%nevtday].push_back(event);
        }
        else if (event.diffuse/tstep < msg->iter) {
          evtcal[target->second][(msg->iter+1)%nevtday].push_back(event);
        }
        else {
          evtcol[target->second].push_back(event);
        }
      }
    }
    else {
      // Find local target
      std::unordered_map<idx_t, idx_t>::iterator target = vtxmap.find(msg->index[i]);
      if (target != vtxmap.end()) {
        // Find target mapping from source
        for (std::size_t j = 0; j < adjcy[target->second].size(); ++j) {
          if (adjcy[target->second][j] == msg->source[i]) {
            event.diffuse = departure + stick[target->second][j+1][0]; // delay always first stick of edge
            event.index = j+1; // 0'th entry is vertex
            // Add to event queue or spillover
            if ((event.diffuse/tstep - msg->iter) < nevtday) {
              evtcal[target->second][(event.diffuse/tstep)%nevtday].push_back(event);
            }
            else if (event.diffuse/tstep < msg->iter) {
              evtcal[target->second][(msg->iter+1)%nevtday].push_back(event);
            }
            else {
              evtcol[target->second].push_back(event);
            }
          }
        }
      }
    }
  }

  delete msg;

  // Start next cycle
  if (cadjprt[prtiter] == nadjprt) {
    // Bookkeepping
    cadjprt[prtiter] = 0;
    prtiter = (prtiter+1)%2;

    // Increment iteration
    ++iter;

    // Start a new cycle
    //thisProxy(prtidx).CycleNetwork();
    cbcycleprt.send();
  }
}


/**************************************************************************
* Network Event Helpers
**************************************************************************/

// Move from Event Collection to Calendar (on new year)
//
void Network::MarkEvent() {
  for (std::size_t i = 0; i < evtcol.size(); ++i) {
    for (std::size_t j = 0; j < evtcol[i].size(); ++j) {
      // Add to event queue or back onto spillover
      if ((evtcol[i][j].diffuse - tsim)/tstep < nevtday) {
        evtcal[i][(evtcol[i][j].diffuse/tstep)%nevtday].push_back(evtcol[i][j]);
      }
      else {
        evtext.push_back(evtcol[i][j]);
      }
    }
    // Copy back spillover
    if (evtext.size()) {
      evtcol[i] = evtext;
      evtext.clear();
    }
    else {
      evtcol[i].clear();
    }
  }
}


/**************************************************************************
* Build Messages
**************************************************************************/

// Build event messages (to multicast)
//
mEvent* Network::BuildEvent() {
  // Initialize distribution message
  int msgSize[MSG_Event];
  msgSize[0] = evtext.size();     // diffuse
  msgSize[1] = evtext.size();     // type
  msgSize[2] = evtext.size();     // source
  msgSize[3] = evtext.size();     // index
  msgSize[4] = evtext.size();     // data
  mEvent *mevent = new(msgSize, 0) mEvent;
  mevent->nevent = evtext.size();
  mevent->iter = iter;
  
  // Pack event information
  for (std::size_t i = 0; i < evtext.size(); ++i) {
    // Add event to message
    mevent->diffuse[i] = evtext[i].diffuse;
    mevent->type[i] = evtext[i].type;
    mevent->source[i] = evtext[i].source;
    mevent->index[i] = evtext[i].index;
    mevent->data[i] = evtext[i].data;
  }

  return mevent;
}
