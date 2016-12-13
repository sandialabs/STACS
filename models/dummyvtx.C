/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DummyVtx : public NetModelTmpl < 1, DummyVtx > {
  public:
    /* Constructor */
    DummyVtx() {
      // parameters
      paramlist.resize(1);
      paramlist[0] = "dp0";
      // states
      statelist.resize(2);
      statelist[0] = "v";
      statelist[1] = "u";
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
tick_t DummyVtx::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) {
  idx_t tcomp = ((idx_t) tdiff)*100/TICKS_PER_MS; // 1ms sleeps for 100us
  if (tcomp > 1000000) { tcomp = 1000000; }
  else if (tcomp == 0) { tcomp = 1; }
  std::this_thread::sleep_for(std::chrono::microseconds(tcomp));

  // generate events
  event_t evtpre;
  evtpre.diffuse = tdrift + tdiff/2;
  evtpre.index = EVENT_EXTERNAL | EVENT_LOCALEDG;
  evtpre.type = EVTYPE_SPIKE;
  evtlog.push_back(evtpre);

  return tdiff;
}

// Simulation jump
//
void DummyVtx::Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) {
  CkPrintf("Jumping Vtx\n");
}

