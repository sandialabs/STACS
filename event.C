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
  msgSize[0] = evtlog.size()-1;     // diffuse
  msgSize[1] = evtlog.size()-1;     // index
  msgSize[2] = evtlog.size()-1;     // type
  msgSize[3] = evtlog.size()-1;     // data
  mEvent *mevent = new(msgSize, 0) mEvent;
  mevent->nevent = evtlog.size()-1;
  mevent->iter = iter;
  //mevent->prtidx = prtidx;
  
  // Pack event information (0'th is event template)
  // Check if recording events
  if (recevtlist) {
    for (std::size_t i = 1; i < evtlog.size(); ++i) {
      // Add event to message
      mevent->diffuse[i-1] = evtlog[i].diffuse;
      mevent->index[i-1] = evtlog[i].index;
      mevent->type[i-1] = evtlog[i].type;
      mevent->data[i-1] = evtlog[i].data;
      // Add to record if listed
      if (evtlog[i].type & recevtlist) {
        recevt.push_back(evtlog[i]);
      }
    }
  }
  else {
    for (std::size_t i = 1; i < evtlog.size(); ++i) {
      // Add event to message
      mevent->diffuse[i-1] = evtlog[i].diffuse;
      mevent->index[i-1] = evtlog[i].index;
      mevent->type[i-1] = evtlog[i].type;
      mevent->data[i-1] = evtlog[i].data;
    }
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
        evtpre.diffuse = evtdif + stick[(*target)[0]][(*target)[1]][0];
        // Add to event queue or spillover
        if ((evtpre.diffuse - tsim)/tstep < evtcal) {
          event[(*target)[0]][(evtpre.diffuse/tstep)%evtcal].push_back(evtpre);
        }
        else if (evtpre.diffuse < tsim) {
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
        evtreaux.push_back(evtaux[i][j]);
      }
    }
    // Copy back spillover
    if (evtreaux.size()) {
      evtaux[i] = evtreaux;
      evtreaux.clear();
    }
    else {
      evtaux[i].clear();
    }
  }
}

