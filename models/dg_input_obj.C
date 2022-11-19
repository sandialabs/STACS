/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
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
      paramlist.resize(2);
      paramlist[0] = "I_rheo";
      paramlist[1] = "obj_mult";
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
    // Trajectory variables
    std::vector<std::vector<real_t>> trajectory;
    std::size_t traj_index;
    tick_t tinterval;
    tick_t tupdate;
    // Obj cell variables
    std::vector<std::vector<real_t>> object_location;
    std::vector<std::vector<real_t>> obj_act;
};


/**************************************************************************
* Class methods
**************************************************************************/


// Simulation step
//
tick_t DGInputObject::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // Work with the input to compute spike rates
  // dependent on object salience
  if (tdrift >= tupdate) {
    //real_t x = trajectory[traj_index][0];
    //real_t y = trajectory[traj_index][1];
    //CkPrintf("  updated location: %" PRIreal ", %" PRIreal "\n", x, y);
    event_t event;
    event.diffuse = tdrift + tdiff;
    event.type = EVENT_RATE;
    event.source = REMOTE_EDGE;

    // Figure out which object is going to be attended to (nearest object)
    int attended_obj = 0;
    real_t attended_salience = 0;
    for (int i = 0; (std::size_t) i < obj_act.size(); ++i) {
      real_t dist = std::sqrt((trajectory[traj_index][0] - object_location[i][0]) * 
                              (trajectory[traj_index][0] - object_location[i][0]) +
                              (trajectory[traj_index][1] - object_location[i][1]) *
                              (trajectory[traj_index][1] - object_location[i][1]));
      real_t theta = std::atan2(trajectory[traj_index][1] - object_location[i][1],
                                trajectory[traj_index][0] - object_location[i][0]) - trajectory[traj_index][2];
      real_t salience = exp(-dist*theta*theta);
      // TODO: attend probabilistically based on salience values
      if (salience < attended_salience) {
        attended_obj = i;
      }
    }
    // Compute salience cell activation
    for (std::size_t j = 0; j < obj_act[attended_obj].size(); ++j) {
      // Figure out the rate without discontinuities
      real_t I_obj = param[0]*param[1]*obj_act[attended_obj][j];
      
      // generate events
      event.index = (idx_t) (j);
      event.data = I_obj;
      events.push_back(event);
    }
    
    // Update time for next eval
    tupdate = tdrift + tinterval;
    ++traj_index;
    if (traj_index >= trajectory.size()) { traj_index = 0; }
  }

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
  try {
    trajectory = input["trajectory"].as<std::vector<std::vector<real_t>>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: trajectory not defined\n");
    trajectory = {{0.0,0.0,0.0}};
  }
  CkPrintf("  trajectory length: %zu\n", trajectory.size());
  try {
    object_location = input["object_location"].as<std::vector<std::vector<real_t>>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: objects not defined\n");
    object_location = {{0.0,0.0}};
  }
  CkPrintf("  number of objects: %zu\n", object_location.size());
  try {
    obj_act = input["obj_act"].as<std::vector<std::vector<real_t>>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: object cell activity distribution not defined\n");
  }
  CkPrintf("  object neurons: %zu\n", obj_act.size());
}

// Close ports
//
void DGInputObject::ClosePorts() {
  // File was loaded into input YAML::Node and closed with LoadFile
}
