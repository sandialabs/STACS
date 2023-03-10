/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class IzhiNeuron : public ModelTmpl < 10, IzhiNeuron > {
  public:
    /* Constructor */
    IzhiNeuron() {
      // parameters
      paramlist.resize(4);
      paramlist[0] = "a";
      paramlist[1] = "b";
      paramlist[2] = "c";
      paramlist[3] = "d";
      // states
      statelist.resize(4);
      statelist[0] = "v";
      statelist[1] = "u";
      statelist[2] = "I";
      statelist[3] = "I_app";
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
void IzhiNeuron::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = -70.0;
    state[1] = param[1] * -70.0;
    state[2] = 0;
    state[3] = 0;
}

// Simulation step
//
tick_t IzhiNeuron::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;
  // update state
  state[0] = state[0] + (tstep/2)*((0.04*state[0]+5)*state[0] + 140 - state[1] + state[2] + state[3]);
  state[0] = state[0] + (tstep-tstep/2)*((0.04*state[0]+5)*state[0] + 140 - state[1] + state[2] + state[3]);
  state[1] = state[1] + tstep*param[0]*(0.2*state[0] - state[1]);

  // Clear transient current for next time
  state[2] = 0;
  
  // if spike occured, generate event
  if (state[0] >= 30) {
    // reset
    state[0] = param[2];
    state[1] = state[1] + param[3];

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
void IzhiNeuron::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_STIM) {
    // Add stim to applied current
    state[0][3] += event.data;
  }
}
