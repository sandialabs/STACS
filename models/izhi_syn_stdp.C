/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class IzhiSynSTDP : public NetModelTmpl < 12, IzhiSynSTDP > {
  public:
    /* Constructor */
    IzhiSynSTDP() {
      // parameters
      paramlist.resize(3);
      paramlist[0] = "wmax";
      paramlist[1] = "tau";
      paramlist[2] = "update";
      // states
      statelist.resize(4);
      statelist[0] = "weight";
      statelist[1] = "wdelta";
      statelist[2] = "ptrace";
      statelist[3] = "ntrace";
      // sticks
      sticklist.resize(3);
      sticklist[0] = "delay";
      sticklist[1] = "ptlast";
      sticklist[2] = "ntlast";
      // auxiliary states
      auxstate.resize(1);
      auxstate[0] = "I";
      // auxiliary sticks
      auxstick.resize(0);
      // ports
      portlist.resize(0);
    }
    
    /* Periodic events */
    void addRepeat(idx_t modidx, std::vector<event_t>& repevt);

    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog);
    void Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux);
    void Hop(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux);
};

/**************************************************************************
* Class methods
**************************************************************************/

// Periodic events
//
void IzhiSynSTDP::addRepeat(idx_t modidx, std::vector<event_t>& repevt) {
  event_t evtpre;
  evtpre.diffuse = 0;
  evtpre.source = 0;
  evtpre.index = modidx;
  evtpre.type = EVENT_EDGUP;
  evtpre.data = param[2];
  repevt.push_back(evtpre);
}


// Simulation step
//
tick_t IzhiSynSTDP::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) {
  return tdiff;
}

// Simulation jump
//
void IzhiSynSTDP::Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) {
  if (evt.type == EVENT_SPIKE) {
    // External spike event
    if (evt.source >= 0) {
      idx_t e = evt.index;
      // Apply effect to neuron (vertex)
      state[0][aux[0].stateidx[0]] += state[e][0];
      // Compute trace value
      state[e][3] = state[e][3]*exp(-((real_t) (evt.diffuse - stick[e][2])/TICKS_PER_MS)/param[1]);
      // Depress weight
      state[e][1] += state[e][3];
      // Reset trace to 0.1 when pre-synaptic neuron fires
      state[e][2] = 0.1;
      // Update last event value for calculating trace
      stick[e][1] = evt.diffuse;
    }
    // Internal spike event
    else if (evt.source < 0) {
      idx_t e = evt.index;
      // Compute trace value
      state[e][2] = state[e][2]*exp(-((real_t) (evt.diffuse - stick[e][1])/TICKS_PER_MS)/param[1]);
      // Facilitate weight
      state[e][1] += state[e][2];
      // Reset negative trace to -0.12 when post-synaptic neuron fires
      state[e][3] = -0.12;
      // Update last event value for calculating trace
      stick[e][2] = evt.diffuse;
    }
  }
  else if (evt.type == EVENT_EDGUP) {
    idx_t e = evt.index;
    // Update weight only every second
    state[e][0] += 0.01 + state[e][1];
    if (state[e][0] < 0) {
      state[e][0] = 0;
    }
    else if (state[e][0] > param[0]) {
      state[e][0] = param[0];
    }
    // Filter the change in weight
    state[e][1] *= 0.9;
  }
}

// Simulation hop
//
void IzhiSynSTDP::Hop(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) {
  if (evt.type == EVENT_SPIKE && evt.source >= 0) {
    // Apply effect to neuron (vertex)
    state[0][aux[0].stateidx[0]] += state[evt.index][0];
  }
}


