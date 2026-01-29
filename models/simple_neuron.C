/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class SimpleNeuron : public ModelTmpl < ModelHash("simple_neuron"), SimpleNeuron > {
  public:
    /* Constructor */
    SimpleNeuron() {
      // parameters
      paramlist.resize(0);
      // states
      statelist.resize(6);
      statelist[0] = "v";
      statelist[1] = "v_thresh";
      statelist[2] = "v_reset";
      statelist[3] = "v_bias";
      statelist[4] = "v_leak";
      statelist[5] = "I_syn";
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
void SimpleNeuron::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = state[2];
    state[5] = 0.0;
}

// Simulation step
//
tick_t SimpleNeuron::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;

  // update state (tstep will mostly just be 1)
  state[0] = state[0] + state[3] + state[5];
  
  // Clear transient current for next time
  state[5] = 0;
  
  // Regular spiking event
  if (state[0] > state[1]) {
    // reset
    state[0] = state[2];

    // generate events
    event_t event;
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

  return tickstep;
}

// Simulation jump
//
void SimpleNeuron::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_STIM) {
    // Add stim to applied current
    state[0][5] = event.data;
  }
}

