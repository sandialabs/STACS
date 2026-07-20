/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class SwitchSynapse : public ModelTmpl < ModelHash("switch_synapse"), SwitchSynapse > {
  public:
    /* Constructor */
    SwitchSynapse() {
      // parameters
      paramlist.resize(0);
      // states
      statelist.resize(2);
      statelist[0] = "weight";
      statelist[1] = "switch";
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
tick_t SwitchSynapse::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  return tdiff;
}

// Simulation jump
//
void SwitchSynapse::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  // External spike event
  // Affected state is flexibly indexed by the switch (e.g. multiple input accumulators)
  // Note: be careful switch value doesn't exceed input accumulators
  if (event.type == EVENT_SPIKE && event.source >= 0) {
    // Apply effect to neuron (vertex)
    state[0][auxidx[0].stateidx[0] + ((std::size_t) state[event.index][1])] += state[event.index][0];
  }
  else if (event.type == EVENT_CURRENT && event.source >= 0) {
    // Current (graded spike) multiply by weight
    state[0][auxidx[0].stateidx[0] + ((std::size_t) state[event.index][1])] += state[event.index][0] * event.data;
  }
  else if (event.type == EVENT_STIM && event.source >= 0) {
    // Direct stimulation (just add event payload)
    state[0][auxidx[0].stateidx[0] + ((std::size_t) state[event.index][1])] += event.data;
  }
}
