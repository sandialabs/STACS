/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DGIzhiSynGluStdp : public ModelTmpl < 63, DGIzhiSynGluStdp > {
  public:
    /* Constructor */
    DGIzhiSynGluStdp() {
      // parameters
      paramlist.resize(7);
      paramlist[0] = "n_sites";
      paramlist[1] = "p_rel";
      paramlist[2] = "A_ltp";
      paramlist[3] = "A_ltd";
      paramlist[4] = "tau";
      paramlist[5] = "g_max";
      paramlist[6] = "update";
      // states
      statelist.resize(4);
      statelist[0] = "g_act";
      statelist[1] = "gdelta";
      statelist[2] = "ltptrace";
      statelist[3] = "ltdtrace";
      // sticks
      sticklist.resize(3);
      sticklist[0] = "delay";
      sticklist[1] = "ltptlast";
      sticklist[2] = "ltdtlast";
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
    /* Periodic Events */
    void Leap(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx);
    void getLeap(std::vector<event_t>& events);
};

/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t DGIzhiSynGluStdp::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // This model has no step processes
  return tdiff;
}

// Simulation jump
//
void DGIzhiSynGluStdp::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (!plastic) {
    // External spike event
    if (event.type == EVENT_SPIKE && event.source >= 0) {
      // Apply effect to neuron (vertex)
      real_t g_act = 0.0;
      for (int i = 0; i < int(std::floor(param[0])); ++i) {
        if ((*unifdist)(*rngine) < param[1]) {
          g_act += state[event.index][0];
        }
      }
      state[0][auxidx[0].stateidx[0]] += g_act;
    }
  }
  else if (event.type == EVENT_SPIKE) {
    // External spike event
    if (event.source >= 0) {
      idx_t e = event.index;
      // Apply effect to neuron (vertex)
      real_t g_act = 0.0;
      for (int i = 0; i < int(std::floor(param[0])); ++i) {
        if ((*unifdist)(*rngine) < param[1]) {
          g_act += state[e][0];
        }
      }
      // Apply effect to neuron (vertex)
      state[0][auxidx[0].stateidx[0]] += g_act;
      // Compute trace value
      state[e][3] = state[e][3]*exp(-((real_t) (event.diffuse - stick[e][2])/TICKS_PER_MS)/param[4]);
      // Depress weight
      state[e][1] += state[e][3];
      // Reset trace to 0.1 when pre-synaptic neuron fires
      state[e][2] = param[2];
      // Update last event value for calculating trace
      stick[e][1] = event.diffuse;
    }
    // Internal spike event
    else if (event.source < 0) {
      idx_t e = event.index;
      // Compute trace value
      state[e][2] = state[e][2]*exp(-((real_t) (event.diffuse - stick[e][1])/TICKS_PER_MS)/param[4]);
      // Facilitate weight
      state[e][1] += state[e][2];
      // Reset negative trace to -0.12 when post-synaptic neuron fires
      state[e][3] = -param[3];
      // Update last event value for calculating trace
      stick[e][2] = event.diffuse;
    }
  }
}

// Periodic events
//
void DGIzhiSynGluStdp::Leap(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_SYNUP) {
    idx_t e = event.index;
    // Update weight only every second
    state[e][0] += 0.000001 + state[e][1];
    if (state[e][0] < 0) {
      state[e][0] = 0;
    }
    else if (state[e][0] > param[5]) {
      state[e][0] = param[5];
    }
    // Filter the change in weight
    state[e][1] *= 0.9;
  }
}

void DGIzhiSynGluStdp::getLeap(std::vector<event_t>& events) {
  if (plastic) {
    event_t event;
    event.diffuse = 0;
    event.type = EVENT_SYNUP;
    event.source = -1;
    event.index = 0;
    event.data = param[6];
    events.push_back(event);
  }
}

