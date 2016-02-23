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
extern /*readonly*/ idx_t evtcal;


/**************************************************************************
* Network Event Communication
**************************************************************************/

// Build event messages (to multicast)
//
mEvent* Network::BuildEvent() {
  // Initialize distribution message
  int msgSize[MSG_Event];
  msgSize[0] = evtext.size();     // diffuse
  msgSize[1] = evtext.size();     // index
  msgSize[2] = evtext.size();     // type
  msgSize[3] = evtext.size();     // data
  mEvent *mevent = new(msgSize, 0) mEvent;
  mevent->nevent = evtext.size();
  mevent->iter = iter;
  //mevent->prtidx = prtidx;
  
  // Pack event information (0'th is event template)
  for (std::size_t i = 0; i < evtext.size(); ++i) {
    // Add event to message
    mevent->diffuse[i] = evtext[i].diffuse;
    mevent->index[i] = evtext[i].index;
    mevent->type[i] = evtext[i].type;
    mevent->data[i] = evtext[i].data;
  }

  return mevent;
}


// Multicast communication of event
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
    evtpre.data = msg->data[i];
    // Find target mapping from source
    std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>>::iterator targets = adjmap.find(msg->index[i]);
    if (targets != adjmap.end()) {
      for (std::vector<std::array<idx_t, 2>>::iterator target = targets->second.begin(); target != targets->second.end(); ++target) {
        evtpre.index = (*target)[1];
        evtpre.diffuse = evtdif + stick[(*target)[0]][(*target)[1]][0]; // delay always first stick of edge
        // Add to event queue or spillover
        if ((evtpre.diffuse/tstep - msg->iter) < evtcal) {
          event[(*target)[0]][(evtpre.diffuse/tstep)%evtcal].push_back(evtpre);
        }
        else if (evtpre.diffuse/tstep < msg->iter) {
          event[(*target)[0]][(msg->iter+1)%evtcal].push_back(evtpre);
        }
        else {
          evtaux[(*target)[0]].push_back(evtpre);
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
    thisProxy(prtidx).Cycle();
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
      if ((evtaux[i][j].diffuse - tsim)/tstep < evtcal) {
        event[i][(evtaux[i][j].diffuse/tstep)%evtcal].push_back(evtaux[i][j]);
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

