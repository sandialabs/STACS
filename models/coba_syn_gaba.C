/**
 * Copyright (C) 2023 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class LIFSynCobaGaba : public ModelTmpl < 73, LIFSynCobaGaba > {
  public:
    /* Constructor */
    LIFSynCobaGaba() {
      // parameters
      paramlist.resize(4);
      paramlist[0] = "R_ss";
      paramlist[1] = "U_stp";
      paramlist[2] = "tau_facil";
      paramlist[3] = "tau_rec";
      // states
      statelist.resize(3);
      statelist[0] = "g_gaba";
      statelist[1] = "R_stp";
      statelist[2] = "u_stp";
      // sticks
      sticklist.resize(2);
      sticklist[0] = "delay";
      sticklist[1] = "t_last";
      // auxiliary states
      auxstate.resize(2);
      auxstate[0] = "g_gaba_rise";
      auxstate[1] = "g_gaba_fall";
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
tick_t LIFSynCobaGaba::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  return tdiff;
}

// Simulation jump
//
void LIFSynCobaGaba::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  if (event.type == EVENT_SPIKE && event.source >= 0) {
    if (plastic) {
      // Short-term plasticity
      // R_1 = 1.0 - U
      if (stick[event.index][1] == 0) {
        stick[event.index][1] = event.diffuse;
        state[event.index][2] = param[1];
        state[event.index][1] = 1.0 - param[1];
      }
      else {
        // Get t_delta and t_last
        real_t t_delta = ((real_t) (event.diffuse - stick[event.index][1]))/ TICKS_PER_MS;
        stick[event.index][1] = event.diffuse;
        // Update R and u
        // u_{n+1} = u_n * exp( - t_delta / tau_facil) + U * ( 1 - u_n * exp( - t_delta / tau_facil))
        state[event.index][2] = state[event.index][2]*exp(-(t_delta/param[2])) + param[1] * (1.0 - state[event.index][2]*exp(-(t_delta/param[2])));
        // R_{n+1} = R_n * ( 1 - u_{n+1} ) * exp( - t_delta / tau_rec) + 1 - exp( - t_delta / tau_rec)
        state[event.index][1] = state[event.index][1]*(1.0 - state[event.index][2])*exp(-(t_delta/param[3])) + 1.0 - exp(-(t_delta/param[3]));
      }
      
      // Apply effect to neuron (vertex)
      state[0][auxidx[0].stateidx[0]] += state[event.index][1] * state[event.index][2] * state[event.index][0]; // gaba_rise
      state[0][auxidx[0].stateidx[1]] += state[event.index][1] * state[event.index][2] * state[event.index][0]; // gaba_fall
    }
    else {
      // Apply effect to neuron (vertex)
      state[0][auxidx[0].stateidx[0]] += param[0] * state[event.index][0]; // gaba_rise
      state[0][auxidx[0].stateidx[1]] += param[0] * state[event.index][0]; // gaba_fall
    }
  }
}
