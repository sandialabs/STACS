/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class IzhiNeuronClamp : public ModelTmpl < 14, IzhiNeuronClamp > {
  public:
    /* Constructor */
    IzhiNeuronClamp() {
      // parameters
      paramlist.resize(7);
      paramlist[0] = "a";
      paramlist[1] = "b";
      paramlist[2] = "c";
      paramlist[3] = "d";
      paramlist[4] = "thal_rate";
      paramlist[5] = "thal_ampl";
      paramlist[6] = "psn_rate";
      // states
      statelist.resize(5);
      statelist[0] = "v";
      statelist[1] = "u";
      statelist[2] = "I";
      statelist[3] = "I_app";
      statelist[4] = "I_clamp";
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
void IzhiNeuronClamp::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = param[2];
    state[1] = param[1] * param[2];
    state[2] = 0;
    state[3] = 0;
    state[4] = 0;
}

// Simulation step
//
tick_t IzhiNeuronClamp::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;

  // random thalamic input
  if ((tstep*param[4]/1000.0) > (*unifdist)(*rngine)) {
    state[2] += param[5];
  }

  // update state
  state[0] = state[0] + (tstep/2)*((0.04*state[0]+5)*state[0] + 140 - state[1] + state[2] + state[3]);
  state[0] = state[0] + (tstep-tstep/2)*((0.04*state[0]+5)*state[0] + 140 - state[1] + state[2] + state[3]);
  state[1] = state[1] + tstep*param[0]*(0.2*state[0] - state[1]);

  // Clear transient current for next time
  state[2] = 0;
  
  // shortcircuit the spiking activity with I_clamp or poisson background rate
  if (state[4] > 0.0 || (tstep*param[6]/1000.0) > (*unifdist)(*rngine)) {
    // reset to default state
    state[0] = param[2];
    state[1] = param[1] * param[2];
    // return to default behavior next iteration
    state[4] = 0.0;
    
    // generate events
    event_t event;
    event.diffuse = tdrift + tickstep;
    event.type = EVENT_SPIKE;
    event.source = REMOTE_EDGES | LOCAL_EDGES;
    event.index = 0;
    event.data = 0.0;
    events.push_back(event);
  }
  // if spike occured, generate event (default behavior)
  else if (state[4] == 0.0 && state[0] >= 30.0) {
    // reset to dynamic state
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
  // suppress any spiking
  else {
    // keep at default state
    state[0] = param[2];
    state[1] = param[1] * param[2];
  }

  return tickstep;
}

// Simulation jump
//
void IzhiNeuronClamp::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_STIM) {
    // Add stim to applied current
    state[0][3] += event.data;
  }
  else if (event.type == EVENT_CLAMP) {
    // Add stim to applied current
    state[0][4] = event.data;
  }
}
