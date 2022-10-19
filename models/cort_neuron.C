/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class CortNeuron : public ModelTmpl < 50, CortNeuron > {
  public:
    /* Constructor */
    CortNeuron() {
      // parameters
      paramlist.resize(8);
      paramlist[0] = "v_thresh";
      paramlist[1] = "tau_ref";
      paramlist[2] = "tau_m";
      paramlist[3] = "C_m";
      paramlist[4] = "tau_syn";
      paramlist[5] = "psn_order";
      paramlist[6] = "psn_rate";
      paramlist[7] = "psn_weight";
      // states
      statelist.resize(2);
      statelist[0] = "v";
      statelist[1] = "I_syn";
      // sticks
      sticklist.resize(1);
      sticklist[0] = "t_last";
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
void CortNeuron::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = 0;
    state[1] = 0;
    stick[0] = ((tick_t) (-param[1]*TICKS_PER_MS));
}

// Simulation step
//
tick_t CortNeuron::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS; // we're using tstep = 0.1ms

  // poisson background input
  idx_t num_spikes = 0;
  for (idx_t i = 0; i < (idx_t) param[5]; ++i){
    if ((tstep*param[6]/1000.0) > (*unifdist)(*rngine)) {
      ++num_spikes;
    }
  }
  state[1] += param[7] * num_spikes;

  // update state (account for refractory period)
  if (tdrift - stick[0] > param[1]*TICKS_PER_MS) {
    state[0] = state[0]*exp(-(tstep/param[2])) + (tstep/param[3])*(state[1]);
  } else {
    state[0] = 0;
  }
  state[1] = state[1]*exp(-(tstep/param[4]));
  
  // if spike occured, generate event
  if (state[0] >= param[0]) {
    // reset
    state[0] = 0;
    stick[0] = tdrift + tickstep;

    // generate events
    event_t event;
    event.diffuse = tdrift + tickstep;
    event.type = EVENT_SPIKE;
    event.source = REMOTE_EDGES; // No need for local edges (no stdp)
    event.index = 0;
    event.data = 0.0;
    events.push_back(event);
  }

  return tickstep;
}

// Simulation jump
//
void CortNeuron::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  // This model has no jump processes
}
