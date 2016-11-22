/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class IzhiSynNoSTDP : public NetModelTmpl < 13, IzhiSynNoSTDP > {
  public:
    /* Constructor */
    IzhiSynNoSTDP() {
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
      auxstate[0] = "I";
      // auxiliary sticks
      auxstick.resize(0);
      // ports
      portlist.resize(0);
    }
    
    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog);
    void Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux);
};

/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t IzhiSynNoSTDP::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) {
  CkPrintf("Stepping Edg\n");
  return tdiff;
}

// Simulation jump
//
void IzhiSynNoSTDP::Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) {
  // External spike event
  if (evt.type == EVTYPE_SPIKE && evt.index > 0) {
    // Apply effect to neuron (vertex)
    state[0][aux[0].stateidx[0]] += state[evt.index][0];
  }
}


