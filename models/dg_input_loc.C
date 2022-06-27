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
class DGInputLocation : public ModelTmpl < 65, DGInputLocation > {
  public:
    /* Constructor */
    DGInputLocation() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "I_rheo";
      paramlist[1] = "loc_mult";
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
    int traj_index;
    tick_t tinterval;
    tick_t tupdate;
    // Grid cell variables
    std::vector<real_t> loc_act;
    std::vector<real_t> lambda;
    std::vector<real_t> theta;
    std::vector<real_t> psi_x;
    std::vector<real_t> psi_y;
};


/**************************************************************************
* Class methods
**************************************************************************/


// Simulation step
//
tick_t DGInputLocation::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
  // Work with the input to compute spike rates
  // dependent on location
  if (tdrift >= tupdate) {
    real_t x = trajectory[traj_index][0];
    real_t y = trajectory[traj_index][1];
    //CkPrintf("  updated location: %" PRIreal ", %" PRIreal "\n", x, y);
    event_t event;
    event.diffuse = tdrift + tdiff;
    event.type = EVENT_STIM;
    event.source = REMOTE_EDGES;

    // Compute grid cell activation
    for (std::size_t i = 0; i < loc_act.size(); ++i) {
      k1 = (4*M_PI*lambda[i]/std::sqrt(6))*((std::cos(theta[i]+M_PI/12) + std::sin(theta[i]+M_PI/12))*(x - psi_x[i]) +
                                            (std::cos(theta[i]+M_PI/12) - std::sin(theta[i]+M_PI/12))*(y - psi_y[i]));
      k2 = (4*M_PI*lambda[i]/std::sqrt(6))*((std::cos(theta[i]+M_PI*5/12) + std::sin(theta[i]+M_PI*5/12))*(x - psi_x[i]) +
                                            (std::cos(theta[i]+M_PI*5/12) - std::sin(theta[i]+M_PI*5/12))*(y - psi_y[i]));
      k3 = (4*M_PI*lambda[i]/std::sqrt(6))*((std::cos(theta[i]+M_PI*3/4) + std::sin(theta[i]+M_PI*3/4))*(x - psi_x[i]) +
                                            (std::cos(theta[i]+M_PI*3/4) - std::sin(theta[i]+M_PI*3/4))*(y - psi_y[i]));
      grid_act = 2/3*((std::cos(k1) + std::cos(k2) + std::cos(k3))/3+0.5);

      I_loc = param[0]*param[1]*loc_act[i]*grid_act;
      
      // generate events
      event.index = i;
      event.data = I_loc;
      events.push_back(event);
    }

    // Update time for next eval
    tupdate = tdrift + tinterval;
    ++traj_index;
    if (traj_index >= trajectory.size()) { traj_index = 0 }
  }

  return tdiff;
}

// Open ports
//
void DGInputLocation::OpenPorts() {
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
  CkPrintf("  location update set at: %.2g ms\n", tupdate_r);
  tinterval = (tick_t)(tupdate_r*TICKS_PER_MS);
  tupdate = 0;
  traj_index = 0;
  try {
    trajectory = input["trajectory"].as<std::vector<std::vector<real_t>>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: trajectory not defined\n");
    trajectory = {{0.0,0.0,0.0}};
  }
  CkPrintf("  trajectory length: %d\n", trajectory.size());
  try {
    loc_act = input["loc_act"].as<std::vector<real_t>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: grid cell activity distribution not defined\n");
  }
  try {
    lambda = input["lambda"].as<std::vector<real_t>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: grid cell spatial scales not defined\n");
  }
  try {
    theta = input["theta"].as<std::vector<real_t>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: grid cell orientation not defined\n");
  }
  try {
    psi_x = input["psi_x"].as<std::vector<real_t>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: grid cell offset x not defined\n");
  }
  try {
    psi_y = input["psi_y"].as<std::vector<real_t>>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  warning: grid cell offset y not defined\n");
  }
  CkPrintf("  grid neurons: %d\n", loc_act.size());
}

// Close ports
//
void DGInputLocation::ClosePorts() {
  // File was loaded into input YAML::Node and closed with LoadFile
}
