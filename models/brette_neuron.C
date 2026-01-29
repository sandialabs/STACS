/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class BretteNeuron : public ModelTmpl < ModelHash("brette_neuron"), BretteNeuron > {
  public:
    /* Constructor */
    BretteNeuron() {
      // parameters
      paramlist.resize(6);
      paramlist[0] = "v_thresh";
      paramlist[1] = "v_reset";
      paramlist[2] = "E_L";
      paramlist[3] = "tau_m";
      paramlist[4] = "tau_ge";
      paramlist[5] = "tau_gi";
      // states
      statelist.resize(3);
      statelist[0] = "v";
      statelist[1] = "ge";
      statelist[2] = "gi";
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
void BretteNeuron::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = param[1];
    state[1] = 0;
    state[2] = 0;
}

// Simulation step
//
tick_t BretteNeuron::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;
  // update state
  state[0] = state[0] + tstep * (state[1] + state[2] - (state[0] - param[2])) / param[3];
  state[1] = state[1] - tstep * (state[1]) / param[4];
  state[2] = state[2] - tstep * (state[2]) / param[5];

  // if spike occured, generate event
  if (state[0] >= param[0]) {
    // reset
    state[0] = param[1];
    state[1] = 0;
    state[2] = 0;

    // generate events
    event_t event;
    event.diffuse = tdrift + tickstep;
    event.type = EVENT_SPIKE;
    event.source = REMOTE_EDGES | LOCAL_EDGES;
    event.index = 0;
    event.data = 0.0;
    events.push_back(event);
  }

  return tickstep;
}

// Simulation jump
//
void BretteNeuron::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  // pass
}
