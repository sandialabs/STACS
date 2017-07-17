/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DummyVtx : public ModelTmpl < 1, DummyVtx > {
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
    tick_t Step(tick_t tdrift, tick_t diff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx);
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t DummyVtx::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  idx_t tcomp = ((idx_t) tdiff)*100/TICKS_PER_MS; // 1ms sleeps for 100us
  if (tcomp > 1000000) { tcomp = 1000000; }
  else if (tcomp == 0) { tcomp = 1; }
  std::this_thread::sleep_for(std::chrono::microseconds(tcomp));

  // generate events
  event_t event;
  event.diffuse = tdrift + tdiff/2;
  event.type = EVENT_SPIKE;
  event.source = REMOTE_EDGES | LOCAL_EDGES;
  event.index = 0;
  event.data = 0.0;
  events.push_back(event);

  return tdiff;
}

// Simulation jump
//
void DummyVtx::Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {
  CkPrintf("Jumping Vtx\n");
}

