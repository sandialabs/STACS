/**
 * Copyright (C) 2023 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class LIFNeuronCoba : public ModelTmpl < 70, LIFNeuronCoba > {
  public:
    /* Constructor */
    LIFNeuronCoba() {
      // parameters
      paramlist.resize(16);
      paramlist[0] = "v_reset";
      paramlist[1] = "v_thresh";
      paramlist[2] = "g_leak";
      paramlist[3] = "C";
      paramlist[4] = "E_rev_ampa";
      paramlist[5] = "tau_ampa_rise";
      paramlist[6] = "tau_ampa_fall";
      paramlist[7] = "E_rev_nmda";
      paramlist[8] = "tau_nmda_rise";
      paramlist[9] = "tau_nmda_fall";
      paramlist[10] = "mg_block_nmda";
      paramlist[11] = "E_rev_gaba";
      paramlist[12] = "tau_gaba_rise";
      paramlist[13] = "tau_gaba_fall";
      paramlist[14] = "t_ref_min";
      paramlist[15] = "t_ref_max";
      // states
      statelist.resize(7);
      statelist[0] = "v";
      statelist[1] = "g_ampa_rise";
      statelist[2] = "g_ampa_fall";
      statelist[3] = "g_nmda_rise";
      statelist[4] = "g_nmda_fall";
      statelist[5] = "g_gaba_rise";
      statelist[6] = "g_gaba_fall";
      // sticks
      sticklist.resize(2);
      sticklist[0] = "t_last";
      sticklist[1] = "t_refract";
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
void LIFNeuronCoba::Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) {
    state[0] = param[0];
    state[1] = 0;
    state[2] = 0;
    state[3] = 0;
    state[4] = 0;
    state[5] = 0;
    state[6] = 0;
    stick[0] = 0; // TODO: should be -inf
}

// Simulation step
//
tick_t LIFNeuronCoba::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;

  // random initial event
  if ((*unifdist)(*rngine) < 0.00001) {
    // reset
    state[0] = param[0];
    // capture time of spike
    stick[0] = tdrift + tdiff;
    // Reset conductances too?
    // Update refractory period randomly
    //real_t trefract = (param[14] + round((*unifdist)(*rngine) * (param[15] - param[14])));
    //stick[1] = ((tick_t) trefract * TICKS_PER_MS);

    // generate events
    event_t event;
    event.diffuse = tdrift + tickstep;
    event.type = EVENT_SPIKE;
    event.source = REMOTE_EDGES | LOCAL_EDGES;
    event.index = 0;
    event.data = 0.0;
    events.push_back(event);
  }
  
  // Normal current integration
  if ((tdrift - stick[0]) > stick[1] || stick[0] == 0) {
    // compute the various currents
    real_t I_ampa = (state[1] - state[2]) * (state[0] - param[4]); // ampa
    real_t mg_block = 1.0 / (1.0 + (param[10] / 3.57) * exp(-0.062 * state[0]));
    real_t I_nmda = (state[3] - state[4]) * (state[0] - param[7]); // nmda
    I_nmda = I_nmda * (1.0 - mg_block);
    real_t I_gaba = (state[5] - state[6]) * (state[0] - param[11]); // gaba
    real_t I_leak = param[2] * (state[0] - param[0]); // leak
    // update voltages (with leak)
    state[0] = state[0] - (tstep/param[3]) * (I_leak + I_ampa + I_nmda + I_gaba) / 100.0;
    // update the various conductances
    //state[1] = state[1]*exp(-(tstep/param[5])); // ampa_rise
    //state[2] = state[2]*exp(-(tstep/param[6])); // ampa_fall
    //state[3] = state[3]*exp(-(tstep/param[8])); // nmda_rise
    //state[4] = state[4]*exp(-(tstep/param[9])); // nmda_fall
    //state[5] = state[5]*exp(-(tstep/param[12])); // gaba_rise
    //state[6] = state[6]*exp(-(tstep/param[13])); // gaba_fall
    state[1] -= state[1]*(tstep/param[5]); // ampa_rise
    state[2] -= state[2]*(tstep/param[6]); // ampa_fall
    state[3] -= state[3]*(tstep/param[8]); // nmda_rise
    state[4] -= state[4]*(tstep/param[9]); // nmda_fall
    state[5] -= state[5]*(tstep/param[12]); // gaba_rise
    state[6] -= state[6]*(tstep/param[13]); // gaba_fall
    
    // if spike occured, generate event
    if (state[0] >= param[1]) {
      // reset
      state[0] = param[0];
      // capture time of spike
      stick[0] = tdrift + tdiff;
      // Update refractory period randomly
      //real_t trefract = (param[14] + round((*unifdist)(*rngine) * (param[15] - param[14])));
      //stick[1] = ((tick_t) trefract * TICKS_PER_MS);

      // generate events
      event_t event;
      event.diffuse = tdrift + tickstep;
      event.type = EVENT_SPIKE;
      event.source = REMOTE_EDGES | LOCAL_EDGES;
      event.index = 0;
      event.data = 0.0;
      events.push_back(event);
    }
  }
  // Refractory period
  else {
    state[0] = param[0];
    // TODO: Conductances set to zero or continually decay?
    state[1] = 0;
    state[2] = 0;
    state[3] = 0;
    state[4] = 0;
    state[5] = 0;
    state[6] = 0;
    //state[1] = state[1]*exp(-(tstep/param[5])); // ampa_rise
    //state[2] = state[2]*exp(-(tstep/param[6])); // ampa_fall
    //state[3] = state[3]*exp(-(tstep/param[8])); // nmda_rise
    //state[4] = state[4]*exp(-(tstep/param[9])); // nmda_fall
    //state[5] = state[5]*exp(-(tstep/param[12])); // gaba_rise
    //state[6] = state[6]*exp(-(tstep/param[13])); // gaba_fall
  }

  return tickstep;
}

// Simulation jump
//
void LIFNeuronCoba::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  // pass
}

