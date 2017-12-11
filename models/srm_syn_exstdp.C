/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class SRMSynExSTDP : public ModelTmpl < 21, SRMSynExSTDP > {
  public:
    /* Constructor */
    SRMSynExSTDP() {
      // parameters
      paramlist.resize(7);
      paramlist[0] = "wmax";
      paramlist[1] = "umax";
      paramlist[2] = "alpha";
      paramlist[3] = "pdw";
      paramlist[4] = "ptau";
      paramlist[5] = "ndw";
      paramlist[6] = "ntau";
      // states
      statelist.resize(1);
      statelist[0] = "weight";
      // sticks
      sticklist.resize(2);
      sticklist[0] = "delay";
      sticklist[1] = "tlast";
      // auxiliary states
      auxstate.resize(1);
      auxstate[0] = "v_app";
      // auxiliary sticks
      auxstick.resize(1);
      auxstick[0] = "tlast";
      // ports
      portlist.resize(0);
    }

    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t diff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx);
    void Leap(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx);
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t SRMSynExSTDP::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  return tdiff;
}

// Simulation jump
//
void SRMSynExSTDP::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  // External spike event
  if (event.type == EVENT_SPIKE && event.source >= 0) {
    // Apply effect to neuron (vertex)
    state[0][auxidx[0].stateidx[0]] += state[event.index][0];
  }
}

// Simulation leap
//
void SRMSynExSTDP::Leap(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_SPIKE) {
    // External spike event
    if (event.source >= 0) {
      idx_t e = event.index;
      // Apply effect to neuron (vertex)
      state[0][auxidx[0].stateidx[0]] += state[e][0] * param[1];

      // Depress weight
      real_t dt = ((real_t)(event.diffuse - stick[0][auxidx[0].stickidx[0]]))/TICKS_PER_MS;
      if (dt < param[6]) {
        real_t dw = param[5] * ((param[6] - dt)/param[6]);
        state[e][0] = state[e][0] - param[2] * (state[e][0]) * dw;
      }
      stick[e][1] = event.diffuse;
    }
    // Internal spike event
    else if (event.source < 0) {
      idx_t e = event.index;
      // Facilitate weight
      real_t dt = ((real_t)(event.diffuse - stick[e][1]))/TICKS_PER_MS;
      if (dt < param[4]) {
        real_t pthalf = param[4]/2.0;
        if (dt <= pthalf) {
          real_t dw = param[3] * (dt/pthalf);
          state[e][0] = state[e][0] + param[2] * (param[0] - state[e][0]) * dw;
        }
        else {
          real_t dw = param[3] * ((param[4] - dt)/pthalf);
          state[e][0] = state[e][0] + param[2] * (param[0] - state[e][0]) * dw;
        }
      }
    }
  }
}
