/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DummyEdg : public NetModelTmpl < 2, DummyEdg > {
  public:
    /* Constructor */
    DummyEdg() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "dp0";
      paramlist[1] = "dp1";
      // states
      statelist.resize(2);
      statelist[0] = "weight";
      statelist[1] = "trace";
      // sticks
      sticklist.resize(1);
      sticklist[0] = "delay";
      // auxiliary states
      auxstate.resize(1);
      auxstate[0] = "v";
      // auxiliary sticks
      auxstick.resize(0);
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
tick_t DummyEdg::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) {
  CkPrintf("Stepping Edg\n");
  return tdiff;
}

// Simulation jump
//
void DummyEdg::Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) {
  //CkPrintf("Jumping Edg\n");
  std::this_thread::sleep_for(std::chrono::microseconds(10));
}

