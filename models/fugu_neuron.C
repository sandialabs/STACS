/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class FuguNeuron : public ModelTmpl < 80, FuguNeuron > {
  public:
    /* Constructor */
    FuguNeuron() {
      // parameters
      paramlist.resize(0);
      // states
      statelist.resize(8);
      statelist[0] = "v";
      statelist[1] = "v_thresh";
      statelist[2] = "v_reset";
      statelist[3] = "v_bias";
      statelist[4] = "v_leak";
      statelist[5] = "p_spike";
      statelist[6] = "I_syn";
      statelist[7] = "I_clamp";
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
void FuguNeuron::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = state[2];
    state[6] = 0.0;
    state[7] = 0.0;
}

// Simulation step
//
tick_t FuguNeuron::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;

  // update state (tstep will mostly just be 1)
  //state[0] = (state[0] * state[4]) + state[3] + state[6];
  state[0] = state[0] + state[3] + state[6];
  
  // Clear transient current for next time
  state[6] = 0;
  
  // shortcircuit the spiking activity with I_clamp
  if (state[7] > 0.0) {
    // reset
    state[0] = state[2];
    // return to default behavior next iteration
    state[7] = 0.0;
    
    // generate events
    event_t event;
    //event.diffuse = tdrift + tickstep;
    event.diffuse = tdrift; // same timestep
    event.type = EVENT_SPIKE;
    event.source = REMOTE_EDGES;
    event.index = 0;
    event.data = 0.0;
    events.push_back(event);
  }
  else if (state[7] == 0.0) {
    // Regular spiking event (with stochasticity)
    if (state[0] > state[1] && state[5] >= (*unifdist)(*rngine)) {
      // reset
      state[0] = state[2];

      // generate events
      event_t event;
      //event.diffuse = tdrift + tickstep;
      event.diffuse = tdrift; // same timestep
      event.type = EVENT_SPIKE;
      event.source = REMOTE_EDGES;
      event.index = 0;
      event.data = 0.0;
      events.push_back(event);
    }
    else {
      // decay the membrane
      state[0] = state[0] * (1 - state[4]);
    }
  }
  else { // I_clamp < 0.0
    // No spiking if supressed
    state[0] = state[2];
  }

  return tickstep;
}

// Simulation jump
//
void FuguNeuron::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_CLAMP) {
    // Add stim to applied current
    state[0][7] = event.data;
  }
}

