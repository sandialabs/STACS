/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DGIzhiNeuron : public ModelTmpl < 60, DGIzhiNeuron > {
  public:
    /* Constructor */
    DGIzhiNeuron() {
      // parameters
      paramlist.resize(13);
      paramlist[0] = "vt";
      paramlist[1] = "vr";
      paramlist[2] = "C";
      paramlist[3] = "a";
      paramlist[4] = "b";
      paramlist[5] = "c";
      paramlist[6] = "d";
      paramlist[7] = "k1";
      paramlist[8] = "k2";
      paramlist[9] = "tau_glu";
      paramlist[10] = "tau_gaba";
      paramlist[11] = "E_glu";
      paramlist[12] = "E_gaba";
      // states
      statelist.resize(4);
      statelist[0] = "v";
      statelist[1] = "u";
      statelist[2] = "g_act_glu";
      statelist[3] = "g_act_gaba";
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
void DGIzhiNeuron::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = param[1];
    state[1] = param[4] * param[1];
    state[2] = 0;
    state[3] = 0;
}

// Simulation step
//
tick_t DGIzhiNeuron::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;
  // update state
  real_t k = param[7] + param[8]*std::tanh(state[0] - param[0]);
  real_t I_glu = state[2] * (param[11] - state[0]);
  real_t I_gaba = state[3] * (param[12] - state[0]);
  state[0] = state[0] + (tstep/2)*(k*(state[0] - param[1])*(state[0] - param[0]) - state[1] + I_glu + I_gaba)/param[2];
  state[0] = state[0] + (tstep-tstep/2)*(k*(state[0] - param[1])*(state[0] - param[0]) - state[1] + I_glu + I_gaba)/param[2];
  state[1] = state[1] + tstep*param[3]*(param[4]*(state[0] - param[1]) - state[1]);

  // Update applied current
  state[2] = state[2]*exp(-(tstep/param[9]));
  state[3] = state[3]*exp(-(tstep/param[10]));

  if (std::isnan(state[0])) {
    state[0] = param[1];
    state[1] = 0.0;
  }

  // Make sure voltage > E_gaba
  if (state[0] < param[12]) {
    state[0] = param[12];
  }
  
  // if spike occured, generate event
  if (state[0] >= param[0]) {
    // reset
    state[0] = param[5];
    state[1] = state[1] + param[6];

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
void DGIzhiNeuron::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_STIM) {
    // Add stim to applied current
    //state[0][4] += event.data;
  }
}
