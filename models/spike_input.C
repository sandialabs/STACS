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
class SpikeInput : public ModelTmpl < 5, SpikeInput > {
  public:
    /* Constructor */
    SpikeInput() {
      // parameters
      paramlist.resize(1);
      paramlist[0] = "n";
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
    // Event list information
    std::vector<std::vector<real_t>> spike_list;
    // Control timing
    std::vector<tick_t> next_spike_tick;
    std::vector<int> next_spike_index;
    tick_t next_tick_update;
};


/**************************************************************************
* Class methods
**************************************************************************/


// Simulation step
//
tick_t SpikeInput::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // Generate spike times based on input
  event_t event;
  event.type = EVENT_CLAMP;
  event.source = REMOTE_EDGE;
  if (tdrift >= next_tick_update) {
    // TODO: there is probably a more efficient way of doing this
    next_tick_update = TICK_T_MAX;
    for (std::size_t i = 0; i < next_spike_tick.size(); ++i) {
      if (tdrift >= next_spike_tick[i]) {
        // Add event to queue
        event.index = i;
        event.data = 1.0;
        event.diffuse = tdrift + TICKS_PER_MS;
        events.push_back(event);
        // Update control timing
        next_spike_index[i] += 1;
        if (next_spike_index[i] >= spike_list[i].size()) {
          next_spike_tick[i] = TICK_T_MAX;
        }
        else {
          next_spike_tick[i] = (tick_t)(spike_list[i][next_spike_index[i]] * TICKS_PER_MS);
        }
        // Find min for next tick update
        if (next_spike_tick[i] < next_tick_update) {
          next_tick_update = next_spike_tick[i];
        }
      }
    }
    //next_tick_update = *std::min_element(next_spike_tick.begin(), next_spike_tick.end());
  }
  return tdiff;
}


// Open ports
//
void SpikeInput::OpenPorts() {
  char ymlfile[100];
  sprintf(ymlfile, "%s/%s", netwkdir.c_str(), portname[0].c_str());
  // Load input file
  CkPrintf("Reading input from %s\n", ymlfile);
  try {
    input = YAML::LoadFile(ymlfile);
  } catch (YAML::BadFile& e) {
    CkPrintf("  %s\n", e.what());
  }
  // Load spike input
  spike_list.clear();
  try {
    spike_list = input["spike_list"].as<std::vector<std::vector<real_t>>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: spike list not defined\n");
  }
  CkPrintf("  number of connected vertices: %d\n", spike_list.size());
  // Make sure the number of connected vertices matches the yaml file input
  CkAssert(param[0] == spike_list.size());
  // Set up spiking updates
  next_spike_tick.resize(spike_list.size());
  next_spike_index.resize(spike_list.size());
  for (std::size_t i = 0; i < spike_list.size(); ++i) {
    if (spike_list[i].size()) {
      next_spike_tick[i] = (tick_t)(spike_list[i][0] * TICKS_PER_MS);
      next_spike_index[i] = 0;
    }
    else {
      // for vertices with no input
      next_spike_tick[i] = TICK_T_MAX;
    }
  }
}

// Close ports
//
void SpikeInput::ClosePorts() {
  // File was loaded into input YAML::Node and closed with LoadFile
}
