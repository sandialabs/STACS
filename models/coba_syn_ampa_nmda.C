/**
 * Copyright (C) 2023 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class LIFSynCobaAmpaNmda : public ModelTmpl < 74, LIFSynCobaAmpaNmda > {
  public:
    /* Constructor */
    LIFSynCobaAmpaNmda() {
      // parameters
      paramlist.resize(5);
      paramlist[0] = "R_ss";
      paramlist[1] = "ampa_nmda_ratio";
      paramlist[2] = "U_stp";
      paramlist[3] = "tau_facil";
      paramlist[4] = "tau_rec";
      // states
      statelist.resize(3);
      statelist[0] = "g_ampa_nmda";
      statelist[1] = "R_stp";
      statelist[2] = "u_stp";
      // sticks
      sticklist.resize(2);
      sticklist[0] = "delay";
      sticklist[1] = "t_last";
      // auxiliary states
      auxstate.resize(5);
      auxstate[0] = "g_ampa_rise";
      auxstate[1] = "g_ampa_fall";
      auxstate[2] = "g_nmda_rise";
      auxstate[3] = "g_nmda_fall";
      auxstate[4] = "t_ref_min";
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
tick_t LIFSynCobaAmpaNmda::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  return tdiff;
}

// Simulation jump
//
void LIFSynCobaAmpaNmda::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_SPIKE && event.source >= 0) {
    if (plastic) {
      // ex to ex
      if (state[0][auxidx[0].stateidx[4]] > 4.0) { //(*unifdist)(*rngine) < (1.0 - exp(-12.3*state[event.index][0]))) {
        // Short-term plasticity
        if (stick[event.index][1] == 0) {
          stick[event.index][1] = event.diffuse;
          state[event.index][2] = param[2];
          state[event.index][1] = 1.0 - param[2];
        }
        else {
          real_t t_delta = ((real_t) (event.diffuse - stick[event.index][1]))/ TICKS_PER_MS;
          stick[event.index][1] = event.diffuse;
          state[event.index][2] = state[event.index][2]*exp(-(t_delta/param[3])) + param[2] * (1.0 - state[event.index][2]*exp(-(t_delta/param[3])));
          state[event.index][1] = state[event.index][1]*(1.0 - state[event.index][2])*exp(-(t_delta/param[4])) + 1.0 - exp(-(t_delta/param[4]));
        }
        
        // variable conductance
        real_t g_act = state[event.index][0] + ((*normdist)(*rngine)) / (0.22 * sqrt(state[event.index][0]));

        // probability of transmission (otherwise, do nothing)
        if ((*unifdist)(*rngine) < (1.0 - exp(-12.3*state[event.index][0]))) {
          // Apply effect to neuron (vertex)
          //state[0][auxidx[0].stateidx[0]] += state[event.index][1] * state[event.index][2] * state[event.index][0]; // ampa_rise
          //state[0][auxidx[0].stateidx[1]] += state[event.index][1] * state[event.index][2] * state[event.index][0]; // ampa_fall
          //state[0][auxidx[0].stateidx[2]] += state[event.index][1] * state[event.index][2] * param[1] * state[event.index][0]; // nmda_rise
          //state[0][auxidx[0].stateidx[3]] += state[event.index][1] * state[event.index][2] * param[1] * state[event.index][0]; // nmda_fall
          state[0][auxidx[0].stateidx[0]] += state[event.index][1] * state[event.index][2] * g_act; // ampa_rise
          state[0][auxidx[0].stateidx[1]] += state[event.index][1] * state[event.index][2] * g_act; // ampa_fall
          state[0][auxidx[0].stateidx[2]] += state[event.index][1] * state[event.index][2] * param[1] * g_act; // nmda_rise
          state[0][auxidx[0].stateidx[3]] += state[event.index][1] * state[event.index][2] * param[1] * g_act; // nmda_fall
        }
      }
      // ex to in
      else {  
        // Short-term plasticity
        if (stick[event.index][1] == 0) {
          stick[event.index][1] = event.diffuse;
          state[event.index][2] = param[2];
          state[event.index][1] = 1.0 - param[2];
        }
        else {
          real_t t_delta = ((real_t) (event.diffuse - stick[event.index][1]))/ TICKS_PER_MS;
          stick[event.index][1] = event.diffuse;
          state[event.index][2] = state[event.index][2]*exp(-(t_delta/param[3])) + param[2] * (1.0 - state[event.index][2]*exp(-(t_delta/param[3])));
          state[event.index][1] = state[event.index][1]*(1.0 - state[event.index][2])*exp(-(t_delta/param[4])) + 1.0 - exp(-(t_delta/param[4]));
        }

        // Apply effect to neuron (vertex)
        state[0][auxidx[0].stateidx[0]] += state[event.index][1] * state[event.index][2] * state[event.index][0]; // ampa_rise
        state[0][auxidx[0].stateidx[1]] += state[event.index][1] * state[event.index][2] * state[event.index][0]; // ampa_fall
        state[0][auxidx[0].stateidx[2]] += state[event.index][1] * state[event.index][2] * param[1] * state[event.index][0]; // nmda_rise
        state[0][auxidx[0].stateidx[3]] += state[event.index][1] * state[event.index][2] * param[1] * state[event.index][0]; // nmda_fall
      }
    }
    else {
      // Apply effect to neuron (vertex)
      state[0][auxidx[0].stateidx[0]] += param[0] * state[event.index][0]; // ampa_rise
      state[0][auxidx[0].stateidx[1]] += param[0] * state[event.index][0]; // ampa_fall
      state[0][auxidx[0].stateidx[2]] += param[0] * param[1] * state[event.index][0]; // nmda_rise
      state[0][auxidx[0].stateidx[3]] += param[0] * param[1] * state[event.index][0]; // nmda_fall
    }
  }
}
