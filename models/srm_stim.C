/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class SRMStim : public ModelTmpl < 200, SRMStim > {
  public:
    /* Constructor */
    SRMStim() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "points";
      paramlist[1] = "ampl";
      // states
      statelist.resize(0);
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
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) { }
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t SRMStim::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // Random thalamic input (for initialization)
  if (param[0] > 0 && tdrift < 300 * TICKS_PER_MS) {
    // generate events
    event_t event;
    event.diffuse = tdrift;
    event.type = EVENT_STIM;
    event.source = REMOTE_EDGE;
    event.index = std::floor(param[0]*(*unifdist)(*rngine));
    event.data = param[1];
    events.push_back(event);
  }
  return tdiff;
}
