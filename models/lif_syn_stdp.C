/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class LIFSynSTDP : public ModelTmpl < 36, LIFSynSTDP > {
  public:
    /* Constructor */
    LIFSynSTDP() {
      // parameters
      paramlist.resize(3);
      paramlist[0] = "wmax";
      paramlist[1] = "tau";
      paramlist[2] = "update";
      // states
      statelist.resize(4);
      statelist[0] = "weight";
      statelist[1] = "wdelta";
      statelist[2] = "ptrace";
      statelist[3] = "ntrace";
      // sticks
      sticklist.resize(3);
      sticklist[0] = "delay";
      sticklist[1] = "ptlast";
      sticklist[2] = "ntlast";
      // auxiliary states
      auxstate.resize(1);
      auxstate[0] = "I";
      // auxiliary sticks
      auxstick.resize(0);
      // ports
      portlist.resize(0);
    }
    
    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx);
    
    /* Periodic Events */
    void Leap(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx);
    void getLeap(std::vector<event_t>& events);
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t LIFSynSTDP::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  return tdiff;
}

// Simulation jump
//
void LIFSynSTDP::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (!plastic) {
    if (event.type == EVENT_SPIKE && event.source >= 0) {
      // Apply effect to neuron (vertex)
      state[0][auxidx[0].stateidx[0]] += state[event.index][0];
    }
  }
  else if (event.type == EVENT_SPIKE) {
    // External spike event
    if (event.source >= 0) {
      idx_t e = event.index;
      // Apply effect to neuron (vertex)
      state[0][auxidx[0].stateidx[0]] += state[e][0];
      // Compute trace value
      state[e][3] = state[e][3]*exp(-((real_t) (event.diffuse - stick[e][2])/TICKS_PER_MS)/param[1]);
      // Depress weight
      state[e][1] += state[e][3];
      // Reset trace to 0.1 when pre-synaptic neuron fires
      state[e][2] = 0.1;
      // Update last event value for calculating trace
      stick[e][1] = event.diffuse;
    }
    // Internal spike event
    else if (event.source < 0) {
      idx_t e = event.index;
      // Compute trace value
      state[e][2] = state[e][2]*exp(-((real_t) (event.diffuse - stick[e][1])/TICKS_PER_MS)/param[1]);
      // Facilitate weight
      state[e][1] += state[e][2];
      // Reset negative trace to -0.12 when post-synaptic neuron fires
      state[e][3] = -0.12;
      // Update last event value for calculating trace
      stick[e][2] = event.diffuse;
    }
  }
}

// Periodic events
//
void LIFSynSTDP::Leap(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_SYNUP) {
    idx_t e = event.index;
    // Update weight only every second
    state[e][0] += 0.01 + state[e][1];
    if (state[e][0] < 0) {
      state[e][0] = 0;
    }
    else if (state[e][0] > param[0]) {
      state[e][0] = param[0];
    }
    // Filter the change in weight
    state[e][1] *= 0.9;
  }
}

void LIFSynSTDP::getLeap(std::vector<event_t>& events) {
  if (plastic) {
    event_t event;
    event.diffuse = 0;
    event.type = EVENT_SYNUP;
    event.source = -1;
    event.index = 0;
    event.data = param[2];
    events.push_back(event);
  }
}

