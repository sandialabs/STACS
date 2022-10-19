/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class LIFPoissonMulti : public ModelTmpl < 35, LIFPoissonMulti > {
  public:
    /* Constructor */
    LIFPoissonMulti() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "rate";
      paramlist[1] = "order";
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
tick_t LIFPoissonMulti::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // Poisson approximation of an input neuron
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;
  // Estimate spikes per ms
  idx_t num_spikes = 0;
  for (idx_t i = 0; i < (idx_t) param[1]; ++i) {
    if ((tstep*param[0]/1000.0) > (*unifdist)(*rngine)) {
      ++num_spikes;
    }
  }
  // generate events
  event_t event;
  event.diffuse = tdrift + tickstep;
  event.type = EVENT_COUNT;
  event.source = REMOTE_EDGE;
  event.index = 0;
  event.data = (real_t) num_spikes;
  events.push_back(event);

  return tickstep;
}

