/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class CortSyn : public ModelTmpl < 51, CortSyn > {
  public:
    /* Constructor */
    CortSyn() {
      // parameters
      paramlist.resize(0);
      // states
      statelist.resize(1);
      statelist[0] = "weight";
      // sticks
      sticklist.resize(1);
      sticklist[0] = "delay";
      // auxiliary states
      auxstate.resize(1);
      auxstate[0] = "I_syn";
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
tick_t CortSyn::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  return tdiff;
}

// Simulation jump
//
void CortSyn::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  // External spike event
  if (event.type == EVENT_SPIKE && event.source >= 0) {
    // Apply effect to neuron (vertex)
    state[0][auxidx[0].stateidx[0]] += state[event.index][0];
  }
  else if (event.type == EVENT_COUNT && event.source >= 0) {
    // Apply effect to neuron (vertex) multiple times
    state[0][auxidx[0].stateidx[0]] += state[event.index][0] * event.data;
  }
}
