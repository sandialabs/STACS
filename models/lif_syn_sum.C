/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class LIFSynSum : public ModelTmpl < 32, LIFSynSum > {
  public:
    /* Constructor */
    LIFSynSum() {
      // parameters
      paramlist.resize(1);
      paramlist[0] = "tau";
      // states
      statelist.resize(1);
      statelist[0] = "I";
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
void LIFSynSum::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = 0;
}

// Simulation step
//
tick_t LIFSynSum::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;
  // update state (tstep will mostly just be 1)
  state[0] = state[0]*exp(-(tstep/param[0]));
  
  // generate events
  event_t event;
  event.diffuse = tdrift + tickstep;
  event.type = EVENT_CURRENT;
  event.source = REMOTE_EDGES;
  event.index = 0;
  event.data = state[0];
  events.push_back(event);

  return tickstep;
}

// Simulation jump
//
void LIFSynSum::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
}
