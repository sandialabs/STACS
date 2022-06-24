/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

// Using yaml-cpp (specification version 1.2)
#include "yaml-cpp/yaml.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DGInputLocation : public ModelTmpl < 65, DGInputLocation > {
  public:
    /* Constructor */
    DGInputLocation() {
      // parameters
      paramlist.resize(0);
      // states
      statelist.resize(0);
      // sticks
      sticklist.resize(0);
      // auxiliary states
      auxstate.resize(0);
      // auxiliary sticks
      auxstick.resize(0);
      // ports
      portlist.resize(1);
      portlist[0] = "input";
    }

    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t diff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) { }
    
    /* Protocol */
    void OpenPorts();
    void ClosePorts();
  
  private:
    YAML::Node input;
};


/**************************************************************************
* Class methods
**************************************************************************/


// Simulation step
//
tick_t DGInputLocation::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // Work with the input to compute spike rates
  // dependent on location
  /*
  for (idx_t i = 0; i < targets.size(); ++i) {
    // generate events
    event_t event;
    event.diffuse = tdrift + tdiff;
    event.type = EVENT_CHGRATE;
    event.source = REMOTE_EDGES;
    event.index = i;
    if (rates[i] != old_rates[i]) {
      old_rates[i] = rates[i];
      event.data = rates[i];
      events.push_back(event);
    }
  }
  */

  return tdiff;
}

// Open ports
//
void DGInputLocation::OpenPorts() {
  // Load input file
  CkPrintf("Reading input from %s\n", portname[0].c_str());
  try {
    input = YAML::LoadFile(portname[0].c_str());
  } catch (YAML::BadFile& e) {
    CkPrintf("  %s\n", e.what());
  }
  // Print some information to display
  // Error checking for file formatting
}

// Close ports
//
void DGInputLocation::ClosePorts() {
  // File was loaded into input YAML::Node and closed with LoadFile
}
