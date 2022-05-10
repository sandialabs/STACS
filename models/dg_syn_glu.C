/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DGIzhiSynGlu : public ModelTmpl < 62, DGIzhiSynGlu > {
  public:
    /* Constructor */
    DGIzhiSynGlu() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "n_sites";
      paramlist[1] = "p_rel";
      // states
      statelist.resize(1);
      statelist[0] = "g_act";
      // sticks
      sticklist.resize(1);
      sticklist[0] = "delay";
      // auxiliary states
      auxstate.resize(1);
      auxstate[0] = "g_act_glu";
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
tick_t DGIzhiSynGlu::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // This model has no step processes
  return tdiff;
}

// Simulation jump
//
void DGIzhiSynGlu::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  // External spike event
  if (event.type == EVENT_SPIKE && event.source >= 0) {
    // Apply effect to neuron (vertex)
    real_t g_act = 0.0;
    for (int i = 0; i < int(std::floor(n_sites)); ++i) {
      if ((*unifdist)(*rngine) < param[1]) {
        g_act += state[event.index][0];
      }
    }
    state[0][auxidx[0].stateidx[0]] += g_act;
  }
}
