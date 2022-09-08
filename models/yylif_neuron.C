/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class YYLIFNeuron : public ModelTmpl < 39, YYLIFNeuron > {
  public:
    /* Constructor */
    YYLIFNeuron() {
      // parameters
      paramlist.resize(5);
      paramlist[0] = "v_thresh";
      paramlist[1] = "tau_m";
      paramlist[2] = "C";
      paramlist[3] = "tau_syn";
      paramlist[4] = "psn_rate";
      // states
      statelist.resize(4);
      statelist[0] = "v";
      statelist[1] = "I";
      statelist[2] = "I_stim";
      // sticks
      sticklist.resize(0);
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
void YYLIFNeuron::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = 0;
    state[1] = 0;
    state[2] = 0;
    state[3] = 0;
}

// Simulation step
//
tick_t YYLIFNeuron::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;
  // update state (tstep will mostly just be 1)
  //state[0] = state[0]*exp(-(tstep/param[1])) + (tstep/param[2])*(state[1] + state[2]);
  state[0] = state[0]*exp(-(tstep/param[1])) + state[1];//(tstep/param[2])*(state[1] + state[2]);

  // Update the current based on tau_syn
  //state[1] = state[1]*exp(-(tstep/param[3]));
  state[1] = 0.0;
  // Make sure membrane voltage always non-negative
  if (state[0] < 0.0) {
    state[0] = 0.0;
  }
  
  // shortcircuit the spiking activity with I_stim
  if (state[2] > 0.0 || (tstep*param[4]/1000.0) > (*unifdist)(*rngine)) {
    // reset
    state[0] = 0.0;
    
    // generate events
    event_t event;
    event.diffuse = tdrift + tickstep;
    event.type = EVENT_SPIKE;
    event.source = REMOTE_EDGES | LOCAL_EDGES;
    event.index = 0;
    event.data = 0.0;
    events.push_back(event);
  }
  else if (state[2] == 0.0) {
    // Random spiking, or regular spiking event
    if (state[0] >= param[0] || (tstep*param[4]/1000.0) > (*unifdist)(*rngine)) {
      // reset
      state[0] = 0.0;

      // generate events
      event_t event;
      event.diffuse = tdrift + tickstep;
      event.type = EVENT_SPIKE;
      event.source = REMOTE_EDGES | LOCAL_EDGES;
      event.index = 0;
      event.data = 0.0;
      events.push_back(event);
    }
  }
  else {
    // No spiking if supressed
    state[0] = 0.0;
  }

  return tickstep;
}

// Simulation jump
//
void YYLIFNeuron::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_STIM) {
    // Add stim to applied current
    state[0][2] += event.data;
  }
}
