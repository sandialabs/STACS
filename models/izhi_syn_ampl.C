/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class IzhiSynAmpl : public NetModelTmpl < 11, IzhiSynAmpl > {
  public:
    /* Constructor */
    IzhiSynAmpl() {
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
    tick_t Step(tick_t tdrift, tick_t diff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog);
    void Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux);
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t IzhiSynAmpl::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) {
  return tdiff;
}

// Simulation jump
//
void IzhiSynAmpl::Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) {
  // External stim event
  if (evt.type == EVENT_STIM && evt.source >= 0) {
    // Add stim to applied current
    state[0][aux[0].stateidx[0]] += evt.data;
  }
}
