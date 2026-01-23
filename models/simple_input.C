/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class SimpleInput : public ModelTmpl < 7, SimpleInput > {
  public:
    /* Constructor */
    SimpleInput() {
      // parameters
      paramlist.resize(0);
      // states
      statelist.resize(1);
      statelist[0] = "I_clamp";
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
void SimpleInput::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = 0.0;
}

// Simulation step
//
tick_t SimpleInput::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;

  // Spike when the clamped value is active
  if (state[0] > 0.0) {
    // generate events
    event_t event;
    event.diffuse = tdrift; // same timestep
    event.type = EVENT_SPIKE;
    event.source = REMOTE_EDGES;
    event.index = 0;
    event.data = 0.0;
    events.push_back(event);
  }
  // Clear clamp value
  state[0] = 0.0;

  return tickstep;
}

// Simulation jump
//
void SimpleInput::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_CLAMP) {
    // Add stim to applied current
    state[0][0] = event.data;
  }
}

