/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DGSynChangeStim : public ModelTmpl < 68, DGSynChangeStim > {
  public:
    /* Constructor */
    DGSynChangeStim() {
      // parameters
      paramlist.resize(0);
      // states
      statelist.resize(0);
      // sticks
      sticklist.resize(1);
      sticklist[0] = "delay";
      // auxiliary states
      auxstate.resize(1);
      auxstate[0] = "I_app";
      // auxiliary sticks
      auxstick.resize(0);
      // ports
      portlist.resize(0);
    }

    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t diff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx);
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t DGSynChangeStim::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  return tdiff;
}

// Simulation jump
//
void DGSynChangeStim::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  // External change baseline firing rate event
  if (event.type == EVENT_STIM && event.source >= 0) {
    // Change the rate on target neuron
    state[0][auxidx[0].stateidx[0]] = event.data;
  }
}

