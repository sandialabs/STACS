/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class IzhiAudio : public NetModelTmpl < 15, IzhiAudio > {
  public:
    /* Constructor */
    IzhiAudio() {
      // parameters
      paramlist.resize(4);
      paramlist[0] = "a";
      paramlist[1] = "b";
      paramlist[2] = "c";
      paramlist[3] = "d";
      // states
      statelist.resize(3);
      statelist[0] = "v";
      statelist[1] = "u";
      statelist[2] = "I";
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
    tick_t Step(tick_t tdrift, tick_t diff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog);
    void Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux);
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t IzhiAudio::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) {
  // if spike occured, generate event
  if (state[0] > 30) {
    // reset
    state[0] = -65;
    state[1] = state[1] + param[3];

    // generate events
    event_t evtpre;
    evtpre.diffuse = tdrift;
    evtpre.index = EVENT_EXTERNAL | EVENT_LOCALEDG;
    evtpre.type = EVTYPE_SPIKE;
    evtlog.push_back(evtpre);
  }
  
  // Input from audio
  state[2] = 0;

  // for numerical stability, use timestep (at most) = 1ms
  tick_t tickstep = (tdiff > TICKS_PER_MS ? TICKS_PER_MS : tdiff);
  real_t tstep = ((real_t) tickstep)/TICKS_PER_MS;
  // update state
  state[0] = state[0] + (tstep/2)*((0.04*state[0]+5)*state[0] + 140 - state[1] + state[2] + state[3]);
  state[0] = state[0] + (tstep-tstep/2)*((0.04*state[0]+5)*state[0] + 140 - state[1] + state[2] + state[3]);
  state[1] = state[1] + tstep*param[0]*(0.2*state[0] - state[1]);

  return tickstep;
}

// Simulation jump
//
void IzhiAudio::Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) {
  //CkPrintf("Jumping Vtx\n");
}
