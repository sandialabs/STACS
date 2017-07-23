/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class IzhiStimIzhiEpis : public ModelTmpl < 106, IzhiStimIzhiEpis > {
  private:
    bool newepisode;
  public:
    /* Constructor */
    IzhiStimIzhiEpis() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "points";
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

      newepisode = true;
    }
    
    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) { }
    void Rerun(std::vector<real_t>& state, std::vector<tick_t>& stick);
};


/**************************************************************************
* Class methods
**************************************************************************/

void IzhiStimIzhiEpis::Rerun(std::vector<real_t>& state, std::vector<tick_t>& stick) {
  newepisode = true;
}

// Simulation step
//
tick_t IzhiStimIzhiEpis::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  if (newepisode) {
    bool forward = ((*unifdist)(*rngine) > 0.5);
    if (forward) {
      for (idx_t i = 0; i < (idx_t) param[0]; ++i) {
        // generate events
        event_t event;
        event.diffuse = tdrift + i * tdiff;
        event.type = EVENT_STIM;
        event.source = REMOTE_EDGE;
        event.index = i;
        event.data = param[1];
        events.push_back(event);
        event.diffuse = tdrift + (i+1) * tdiff;
        event.data = -param[1];
        events.push_back(event);
      }
    }
    else {
      for (idx_t i = 0; i < (idx_t) param[0]; ++i) {
        // generate events
        event_t event;
        event.diffuse = tdrift + (((idx_t) param[0]) - i) * tdiff;
        event.type = EVENT_STIM;
        event.source = REMOTE_EDGE;
        event.index = i;
        event.data = param[1];
        events.push_back(event);
        event.diffuse = tdrift + (((idx_t) param[0]) - i + 1) * tdiff;
        event.data = -param[1];
        events.push_back(event);
      }
    }
    newepisode = false;
  }
  return tdiff;
}
