/**
 * Copyright (C) 2023 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class LIFSynCobaAmpa : public ModelTmpl < 71, LIFSynCobaAmpa > {
  public:
    /* Constructor */
    LIFSynCobaAmpa() {
      // parameters
      paramlist.resize(1);
      paramlist[0] = "R_ss";
      // states
      statelist.resize(1);
      statelist[0] = "g_ampa";
      // sticks
      sticklist.resize(1);
      sticklist[0] = "delay";
      // auxiliary states
      auxstate.resize(2);
      auxstate[0] = "g_ampa_rise";
      auxstate[1] = "g_ampa_fall";
      // auxiliary sticks
      auxstick.resize(0);
      // ports
      portlist.resize(0);
    }
    
    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx);
    
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t LIFSynCobaAmpa::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  return tdiff;
}

// Simulation jump
//
void LIFSynCobaAmpa::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_SPIKE && event.source >= 0) {
    // TODO: Implement short-term plasticity?
    // Apply effect to neuron (vertex)
    state[0][auxidx[0].stateidx[0]] += param[0] * state[event.index][0]; // ampa_rise
    state[0][auxidx[0].stateidx[1]] += param[0] * state[event.index][0]; // ampa_fall
  }
}
