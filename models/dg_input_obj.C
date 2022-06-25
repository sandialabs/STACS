/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

// Using yaml-cpp (specification version 1.2)
#include "yaml-cpp/yaml.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ std::string netwkdir;

/**************************************************************************
* Class declaration
**************************************************************************/
class DGInputObject : public ModelTmpl < 66, DGInputObject > {
  public:
    /* Constructor */
    DGInputObject() {
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
    int traj_index;
    tick_t tinterval;
    tick_t tupdate;
};


/**************************************************************************
* Class methods
**************************************************************************/


// Simulation step
//
tick_t DGInputObject::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // Work with the input to compute spike rates
  // dependent on object salience
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
void DGInputObject::OpenPorts() {
  char ymlfile[100];
  sprintf(ymlfile, "%s/%s", netwkdir.c_str(), portname[0].c_str());
  // Load input file
  CkPrintf("Reading input from %s\n", ymlfile);
  try {
    input = YAML::LoadFile(ymlfile);
  } catch (YAML::BadFile& e) {
    CkPrintf("  %s\n", e.what());
  }
  // Print some information to display
  // Error checking for file formatting
  real_t tupdate_r;
  try {
    tupdate_r = input["update"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    tupdate_r = 100.0; // default of 100ms
  }
  CkPrintf("  object update set at: %.2g ms\n", tupdate_r);
  tinterval = (tick_t)(tupdate_r*TICKS_PER_MS);
  tupdate = 0;
  traj_index = 0;
}

// Close ports
//
void DGInputObject::ClosePorts() {
  // File was loaded into input YAML::Node and closed with LoadFile
}
