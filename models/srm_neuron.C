/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class SRMNeuron : public ModelTmpl < 20, SRMNeuron > {
  public:
    /* Constructor */
    SRMNeuron() {
      // parameters
      paramlist.resize(4);
      paramlist[0] = "vrest";
      paramlist[1] = "vthresh";
      paramlist[2] = "tau";
      paramlist[3] = "tabs";
      // states
      statelist.resize(2);
      statelist[0] = "v";
      statelist[1] = "v_app";
      // sticks
      sticklist.resize(1);
      sticklist[0] = "tlast";
      // auxiliary states
      auxstate.resize(0);
      // auxiliary sticks
      auxstick.resize(0);
      // ports
      portlist.resize(0);
    }

    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t diff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx);
    
    /* Protocol */
    void Reset(std::vector<real_t>& state, std::vector<tick_t>& stick);
};


/**************************************************************************
* Class methods
**************************************************************************/

// Reset model
//
void SRMNeuron::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = param[0];
    state[1] = 0.0;
    stick[0] = 0;
}

// Simulation step
//
tick_t SRMNeuron::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  real_t halfstep = ((real_t ) tdiff)/(TICKS_PER_MS * 2.0);
  // Refractory period
  if ((tdrift - stick[0]) > (tick_t)(param[3]*TICKS_PER_MS)) {
    // Add spikes
    state[0] += state[1];
    state[1] = 0.0;
    if (state[0] > param[1]) {
      // spike and reset
      state[0] = param[0];
      state[1] = 0.0;
      stick[0] = tdrift + tdiff;

      // generate events
      event_t event;
      event.diffuse = tdrift + tdiff;
      event.type = EVENT_SPIKE;
      event.source = REMOTE_EDGES | LOCAL_EDGES;
      event.index = 0;
      event.data = 0.0;
      events.push_back(event);
    }
    // Leak voltage
    state[0] = state[0] - halfstep * (state[0] - param[0]) / param[2];
    state[0] = state[0] - halfstep * (state[0] - param[0]) / param[2];
  }
  else {
    // Leak spikes
    state[1] = state[1] - halfstep * state[1] / param[2];
    state[1] = state[1] - halfstep * state[1] / param[2];
  }
  return tdiff;
}

// Simulation jump
//
void SRMNeuron::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_STIM) {
    // Add stim to applied current
    state[0][0] += event.data;
  }
}
