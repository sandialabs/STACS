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
extern /*readonly*/ tick_t tepisode;

/**************************************************************************
* Class declaration
**************************************************************************/
class YinYangEpis : public ModelTmpl < 130, YinYangEpis > {
  public:
    /* Constructor */
    YinYangEpis() {
      // parameters
      paramlist.resize(3);
      paramlist[0] = "t_max";
      paramlist[1] = "ampl";
      paramlist[2] = "train";
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
    
    /* Control Flow */
    void Renew(std::vector<real_t>& state, std::vector<tick_t>& stick);

  private:
    YAML::Node input;
    // Sample information
    std::vector<real_t> sample_x;
    std::vector<real_t> sample_y;
    std::vector<int> sample_c; // class
    // Control flow
    tick_t tsample;
    int sample_idx;
};


/**************************************************************************
* Class methods
**************************************************************************/


// Renew model
//
void YinYangEpis::Renew(std::vector<real_t>& state, std::vector<tick_t>& stick) {
  // Start of new sample
  tsample = 0;
}


// Simulation step
//
tick_t YinYangEpis::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // Push new events at the start of the sample
  // Otherwise do nothing (e.g. variable sample length)
  if (tsample == 0) {
    // Get the new values
    real_t x = sample_x[sample_idx];
    real_t y = sample_y[sample_idx];
    int c = sample_c[sample_idx];

    // Generate the times
    // Align it to unit millisecond intervals
    // The addition of 1 is manage min delay
    int tx  = ((int) (round(x * param[0]))) + 1;
    int tx1 = ((int) (round((1-x) * param[0]))) + 1;
    int ty  = ((int) (round(y * param[0]))) + 1;
    int ty1 = ((int) (round((1-y) * param[0]))) + 1;
    //int tc  = ((int) (round(param[0]/2.0))) + 1;
    // Find the fastest time the network has all of its information
    int tc = std::max(std::min(tx,tx1),std::min(ty,ty1)) + 4;

    // Generate event times based on the sample information
    event_t event;
    event.type = EVENT_STIM;
    event.source = REMOTE_EDGE;
    // enough current in one time step to force a spike
    // Bias spike
    event.index = 0;
    event.data = param[1];
    event.diffuse = tdrift + TICKS_PER_MS;
    events.push_back(event);
    event.data = -param[1];
    event.diffuse = tdrift + 2 * TICKS_PER_MS;
    events.push_back(event);
    // x spike
    event.index = 1;
    event.data = param[1];
    event.diffuse = tdrift + (tx) * TICKS_PER_MS;
    events.push_back(event);
    event.data = -param[1];
    event.diffuse = tdrift + (tx+1) * TICKS_PER_MS;
    events.push_back(event);
    // 1-x spike
    event.index = 2;
    event.data = param[1];
    event.diffuse = tdrift + (tx1) * TICKS_PER_MS;
    events.push_back(event);
    event.data = -param[1];
    event.diffuse = tdrift + (tx1+1) * TICKS_PER_MS;
    events.push_back(event);
    // y spike
    event.index = 3;
    event.data = param[1];
    event.diffuse = tdrift + (ty) * TICKS_PER_MS;
    events.push_back(event);
    event.data = -param[1];
    event.diffuse = tdrift + (ty+1) * TICKS_PER_MS;
    events.push_back(event);
    // 1-y spike
    event.index = 4;
    event.data = param[1];
    event.diffuse = tdrift + (ty1) * TICKS_PER_MS;
    events.push_back(event);
    event.data = -param[1];
    event.diffuse = tdrift + (ty1+1) * TICKS_PER_MS;
    events.push_back(event);
    // Class label
    if (param[2] > 0.0) {
      event.index = 5;
      event.data = -param[1];
      event.diffuse = tdrift + TICKS_PER_MS;
      events.push_back(event);
      event.data = param[1];
      event.diffuse = tdrift + tepisode + TICKS_PER_MS;
      events.push_back(event);
      event.index = 6;
      event.data = -param[1];
      event.diffuse = tdrift +  TICKS_PER_MS;
      events.push_back(event);
      event.data = param[1];
      event.diffuse = tdrift + tepisode + TICKS_PER_MS;
      events.push_back(event);
      event.index = 7;
      event.data = -param[1];
      event.diffuse = tdrift +  TICKS_PER_MS;
      events.push_back(event);
      event.data = param[1];
      event.diffuse = tdrift + tepisode + TICKS_PER_MS;
      events.push_back(event);
      event.index = 5 + c;
      event.data = 2*param[1];
      event.diffuse = tdrift + tc * TICKS_PER_MS;
      events.push_back(event);
      event.data = -2*param[1];
      event.diffuse = tdrift + (tc+1) * TICKS_PER_MS;
      events.push_back(event);
    }
    
    // Update sample index for next time
    ++sample_idx;
    if (sample_idx >= sample_c.size()) { sample_idx = 0; }
    // Set tsample to something to not push spikes until next time
    tsample = TICKS_PER_MS;
  }

  return tdiff;
}


// Open ports
//
void YinYangEpis::OpenPorts() {
  char ymlfile[100];
  sprintf(ymlfile, "%s/%s", netwkdir.c_str(), portname[0].c_str());
  // Load input file
  CkPrintf("Reading input from %s\n", ymlfile);
  try {
    input = YAML::LoadFile(ymlfile);
  } catch (YAML::BadFile& e) {
    CkPrintf("  %s\n", e.what());
  }
  // Load samples
  sample_idx = 0;
  try {
    sample_x = input["sample_x"].as<std::vector<real_t>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: yin yang samples not defined\n");
  }
  try {
    sample_y = input["sample_y"].as<std::vector<real_t>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: yin yang samples not defined\n");
  }
  try {
    sample_c = input["sample_c"].as<std::vector<int>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: yin yang samples not defined\n");
  }
  CkPrintf("  number of yin-yang samples: %d\n", sample_c.size());
}

// Close ports
//
void YinYangEpis::ClosePorts() {
  // File was loaded into input YAML::Node and closed with LoadFile
}
