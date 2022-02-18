/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class CortLIF : public ModelTmpl < 50, CortLIF > {
  public:
    /* Constructor */
    CortLIF() {
      // parameters
      paramlist.resize(5);
      paramlist[0] = "V_thresh"; // Should be computed as V_thresh - V_reset (-50mV - -65mV = 15mV)
      paramlist[1] = "tau_ref";
      paramlist[2] = "tau_m";
      paramlist[3] = "C_m";
      paramlist[4] = "tau_syn";
      // states
      statelist.resize(2);
      statelist[0] = "V";
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
void CortLIF::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = 0;
    state[1] = 0;
    stick[0] = ((tick_t) (-param[1]*TICKS_PER_MS));
}

// Simulation step
//
tick_t CortLIF::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS; // we're using tstep = 0.1ms
  // update state (is this the right update equation?)
  //if (stick[0] + param[1]*TICKS_PER_MS < tdrift) {
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
void CortLIF::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  // This model has no jump processes
}
