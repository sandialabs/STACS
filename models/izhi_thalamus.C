/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class IzhiThalamus : public NetModelTmpl < 100, IzhiThalamus > {
  public:
    /* Constructor */
    IzhiThalamus() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "rng";
      paramlist[1] = "ampl";
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
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog);
    void Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux);
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t IzhiThalamus::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) {
  // Random thalamic input (to transient current)
  if (param[0] > 0) {
    // generate events
    event_t evtpre;
    evtpre.diffuse = tdrift;
    evtpre.type = EVENT_STIM;
    evtpre.source = REMOTE_EDGE;
    evtpre.index = std::floor(param[0]*(*unifdist)(*rngine));
    evtpre.data = param[1];
    evtlog.push_back(evtpre);
    evtpre.diffuse = tdrift + TICKS_PER_MS;
    evtpre.data = -param[1];
    evtlog.push_back(evtpre);
  }
  return tdiff;
}

// Simulation jump
//
void IzhiThalamus::Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) {
  //CkPrintf("Jumping Vtx\n");
}

