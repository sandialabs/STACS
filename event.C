/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"


/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ int netparts;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ idx_t nevtday;


/**************************************************************************
* Network Event Communication
**************************************************************************/

// Simple go-ahead message
//
void Network::GoAhead(mGo *msg) {
  // Increment coordination
  ++cadjpart[(partiter + (msg->iter - commiter))%2];
  delete msg;

  // Start next cycle
  if (cadjpart[partiter] == nadjpart) {
    // Bookkeepping
    cadjpart[partiter] = 0;
    partiter = (partiter+1)%2;

    // Increment iteration
    ++commiter;

    // Start a new cycle
    //thisProxy(partidx).CycleNetwork();
    cyclepart.send();
  }
}

// Multicast communication of events
//
void Network::CommEvent(mEvent *msg) {
  // Increment coordination
  ++cadjpart[(partiter + (msg->iter - commiter))%2];

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
  if (cadjpart[partiter] == nadjpart) {
    // Bookkeepping
    cadjpart[partiter] = 0;
    partiter = (partiter+1)%2;

    // Increment iteration
    ++commiter;

    // Start a new cycle
    //thisProxy(partidx).CycleNetwork();
    cyclepart.send();
  }
}


// Multicast communication of events (for monitoring)
//
void Network::CommStamp(mEvent *msg) {
  // Increment coordination
  ++cadjpart[(partiter + (msg->iter - commiter))%2];

  // Event prototype
  event_t event;
  tick_t departure;
  stamp_t stamp;

  // Distribute events
  for (std::size_t i = 0; i < msg->nevent; ++i) {
    // Fill in prototype stamp
    if (msg->type[i] == EVENT_SPIKE && msg->index[i] == msg->source[i]) {
      stamp.diffuse = msg->diffuse[i];
      stamp.source = msg->source[i];
      // distribute to polychronous groups
      std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>>::iterator targets = grpmap.find(msg->source[i]);
      if (targets != grpmap.end()) {
        for (std::vector<std::array<idx_t, 2>>::iterator target = targets->second.begin(); target != targets->second.end(); ++target) {
          grpwindow[(*target)[0]][(*target)[1]].push_back(stamp);
        }
      }
    }
    // Fill in prototype event
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
  if (cadjpart[partiter] == netparts) {
    // Bookkeepping
    cadjpart[partiter] = 0;
    partiter = (partiter+1)%2;

    // Increment iteration
    ++commiter;

    // Start a new cycle
    //thisProxy(partidx).CycleNetwork();
    cyclepart.send();
  }
}


/**************************************************************************
* Network Event Helpers
**************************************************************************/

// Move from Event Collection to Calendar (on new year)
//
void Network::SortEventCalendar() {
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
      
// Perform periodic events
//
void Network::SkipEventPlas() {
  std::vector<event_t>::iterator event = evtskip.begin();
  // Compute periodic events
  while (event != evtskip.end() && event->diffuse <= tsim) {
    // Set model index
    idx_t n = event->source;
    // Loop through all models
    for (std::size_t m = 0; m < skipidx[n].size(); ++m) {
      event->index = skipidx[n][m][1];
      if (event->index) {
        model[n]->Leap(*event, state[skipidx[n][m][0]], stick[skipidx[n][m][0]], edgaux[n][vtxmodidx[skipidx[n][m][0]]]);
      }
      else {
        model[n]->Leap(*event, state[skipidx[n][m][0]], stick[skipidx[n][m][0]], vtxaux[skipidx[n][m][0]]);
      }
    }
    // Update timing
    event->diffuse += (tick_t)(event->data*TICKS_PER_MS);
    ++event;
  }
  std::sort(evtskip.begin(), evtskip.end());
  tskip = evtskip.front().diffuse;
}

// Perform periodic events
//
void Network::SkipEvent() {
  std::vector<event_t>::iterator event = evtskip.begin();
  // Compute periodic events
  while (event != evtskip.end() && event->diffuse <= tsim) {
    // Set model index
    idx_t n = event->source;
    // Loop through all models
    for (std::size_t m = 0; m < skipidx[n].size(); ++m) {
      event->index = skipidx[n][m][1];
      if (event->index) {
        model[n]->Jump(*event, state[skipidx[n][m][0]], stick[skipidx[n][m][0]], edgaux[n][vtxmodidx[skipidx[n][m][0]]]);
      }
      else {
        model[n]->Jump(*event, state[skipidx[n][m][0]], stick[skipidx[n][m][0]], vtxaux[skipidx[n][m][0]]);
      }
    }
    // Update timing
    event->diffuse += (tick_t)(event->data*TICKS_PER_MS);
    ++event;
  }
  std::sort(evtskip.begin(), evtskip.end());
  tskip = evtskip.front().diffuse;
}

// Handle generated events (for communication, plastic)
//
void Network::HandleEventPlas(event_t& event, const idx_t i) {
  // TODO: Conversion from edge indices to global (for individual output)
  // Get information
  idx_t target = event.source;
  idx_t index = event.index;
  // Reindex to global
  event.source = vtxidx[i];
  // Record listed event
  if (evtloglist[event.type]) {
    evtlog.push_back(event);
  }
  // Remote events (multicast to edges)
  if (target & REMOTE_EDGES) {
    // reindex to global
    event.index = vtxidx[i];
    // push to communication
    evtext.push_back(event);
  }
  // Remote event (singlecast to edge)
  else if (target & REMOTE_EDGE) {
    // reindex to global
    // TODO: get this value from the target mapping
    event.index = adjcy[i][index];
    // push to communication
    evtext.push_back(event);
  }
  // Remote event (singlecast to vertex)
  else if (target & REMOTE_VERTEX) {
    // reindex to global
    // TODO: get this value from the target mapping
    event.index = -adjcy[i][index]-1; // negative index indicates vertex
    // push to communication
    evtext.push_back(event);
  }
  // Local events (multicast to edges)
  if (target & LOCAL_EDGES) {
    event.source = -1; // negative source indicates local event
    // Jump loops
    if ((event.diffuse - tsim - tstep)/tstep < nevtday) {
      for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
        if (edgmodidx[i][j]) {
          event.index = j+1;
          evtcal[i][(event.diffuse/tstep)%nevtday].push_back(event);
        }
      }
    }
    else if (event.diffuse < tsim + tstep) {
      for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
        if (edgmodidx[i][j]) {
          event.index = j+1;
          // Jump now
          model[edgmodidx[i][j]]->Leap(event, state[i], stick[i], edgaux[edgmodidx[i][j]][vtxmodidx[i]]);
        }
      }
    }
    else {
      for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
        if (edgmodidx[i][j]) {
          event.index = j+1;
          evtcol[i].push_back(event);
        }
      }
    }
  }
  // Local event (singlecast to vertex)
  if (target & LOCAL_VERTEX) {
    // vertex to itself
    event.source = -1; // negative source indicates local event
    event.index = 0;
    if ((event.diffuse - tsim - tstep)/tstep < nevtday) {
      evtcal[i][(event.diffuse/tstep)%nevtday].push_back(event);
    }
    else if (event.diffuse < tsim + tstep) {
      // Jump now
      model[vtxmodidx[i]]->Leap(event, state[i], stick[i], vtxaux[i]);
    }
    else {
      evtcol[i].push_back(event);
    }
  }
}

// Handle generated events (for communication)
//
void Network::HandleEvent(event_t& event, const idx_t i) {
  // TODO: Conversion from edge indices to global (for individual output)
  // Get information
  idx_t target = event.source;
  idx_t index = event.index;
  // Reindex to global
  event.source = vtxidx[i];
  // Record listed event
  if (evtloglist[event.type]) {
    evtlog.push_back(event);
  }
  // Remote events (multicast to edges)
  if (target & REMOTE_EDGES) {
    // reindex to global
    event.index = vtxidx[i];
    // push to communication
    evtext.push_back(event);
  }
  // Remote event (singlecast to edge)
  else if (target & REMOTE_EDGE) {
    // reindex to global
    // TODO: get this value from the target mapping
    event.index = adjcy[i][index];
    // push to communication
    evtext.push_back(event);
  }
  // Remote event (singlecast to vertex)
  else if (target & REMOTE_VERTEX) {
    // reindex to global
    // TODO: get this value from the target mapping
    event.index = -adjcy[i][index]-1; // negative index indicates vertex
    // push to communication
    evtext.push_back(event);
  }
  // Local events (multicast to edges)
  if (target & LOCAL_EDGES) {
    event.source = -1; // negative source indicates local event
    // Jump loops
    if ((event.diffuse - tsim - tstep)/tstep < nevtday) {
      for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
        if (edgmodidx[i][j]) {
          event.index = j+1;
          evtcal[i][(event.diffuse/tstep)%nevtday].push_back(event);
        }
      }
    }
    else if (event.diffuse < tsim + tstep) {
      for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
        if (edgmodidx[i][j]) {
          event.index = j+1;
          // Jump now
          model[edgmodidx[i][j]]->Jump(event, state[i], stick[i], edgaux[edgmodidx[i][j]][vtxmodidx[i]]);
        }
      }
    }
    else {
      for (std::size_t j = 0; j < edgmodidx[i].size(); ++j) {
        if (edgmodidx[i][j]) {
          event.index = j+1;
          evtcol[i].push_back(event);
        }
      }
    }
  }
  // Local event (singlecast to vertex)
  if (target & LOCAL_VERTEX) {
    // vertex to itself
    event.source = -1; // negative source indicates local event
    event.index = 0;
    if ((event.diffuse - tsim - tstep)/tstep < nevtday) {
      evtcal[i][(event.diffuse/tstep)%nevtday].push_back(event);
    }
    else if (event.diffuse < tsim + tstep) {
      // Jump now
      model[vtxmodidx[i]]->Jump(event, state[i], stick[i], vtxaux[i]);
    }
    else {
      evtcol[i].push_back(event);
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
