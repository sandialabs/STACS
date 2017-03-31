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
extern /*readonly*/ idx_t equeue;


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
    // If index != source (singlecast to edge)
    // If index < 0 (singlecast to vertex)
    if (msg->index[i] == msg->source[i]) {
      // Find target mapping from source
      std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>>::iterator targets = adjmap.find(msg->source[i]);
      if (targets != adjmap.end()) {
        for (std::vector<std::array<idx_t, 2>>::iterator target = targets->second.begin(); target != targets->second.end(); ++target) {
          evtpre.diffuse = evtdif + stick[(*target)[0]][(*target)[1]][0]; // delay always first stick of edge
          evtpre.index = (*target)[1];
          // Add to event queue or spillover
          if ((evtpre.diffuse/tstep - msg->iter) < equeue) {
            event[(*target)[0]][(evtpre.diffuse/tstep)%equeue].push_back(evtpre);
          }
          else if (evtpre.diffuse/tstep < msg->iter) {
            event[(*target)[0]][(msg->iter+1)%equeue].push_back(evtpre);
          }
          else {
            evtaux[(*target)[0]].push_back(evtpre);
          }
        }
      }
    }
    else if (msg->index[i] < 0) {
      // Find local target
      std::unordered_map<idx_t, idx_t>::iterator target = vtxmap.find(-(msg->index[i]+1));
      if (target != vtxmap.end()) {
        evtpre.diffuse = evtdif; // direct events to vertices have no edge delay
        evtpre.index = 0;
        // Add to event queue or spillover
        if ((evtpre.diffuse/tstep - msg->iter) < equeue) {
          event[target->second][(evtpre.diffuse/tstep)%equeue].push_back(evtpre);
        }
        else if (evtpre.diffuse/tstep < msg->iter) {
          event[target->second][(msg->iter+1)%equeue].push_back(evtpre);
        }
        else {
          evtaux[target->second].push_back(evtpre);
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
            evtpre.diffuse = evtdif + stick[target->second][j+1][0]; // delay always first stick of edge
            evtpre.index = j+1; // 0'th entry is vertex
            // Add to event queue or spillover
            if ((evtpre.diffuse/tstep - msg->iter) < equeue) {
              event[target->second][(evtpre.diffuse/tstep)%equeue].push_back(evtpre);
            }
            else if (evtpre.diffuse/tstep < msg->iter) {
              event[target->second][(msg->iter+1)%equeue].push_back(evtpre);
            }
            else {
              evtaux[target->second].push_back(evtpre);
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

// Redistribute Event Spillover (on new year)
//
void Network::RedisEvent() {
  for (std::size_t i = 0; i < evtaux.size(); ++i) {
    for (std::size_t j = 0; j < evtaux[i].size(); ++j) {
      // Add to event queue or back onto spillover
      if ((evtaux[i][j].diffuse - tsim)/tstep < equeue) {
        event[i][(evtaux[i][j].diffuse/tstep)%equeue].push_back(evtaux[i][j]);
      }
      else {
        evtext.push_back(evtaux[i][j]);
      }
    }
    // Copy back spillover
    if (evtext.size()) {
      evtaux[i] = evtext;
      evtext.clear();
    }
    else {
      evtaux[i].clear();
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
