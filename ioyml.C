/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 *
 * ioyml.C
 * Handles YAML Ain't Markup Language format
 */

#include "stacs.h"
#include "network.h"
//#include "stream.h"

// Using yaml-cpp (specification version 1.2)
#include "yaml-cpp/yaml.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ std::string filebase;
extern /*readonly*/ std::string filemod;
extern /*readonly*/ idx_t npdat;
extern /*readonly*/ idx_t npnet;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ tick_t tcheck;
extern /*readonly*/ std::string rpcport;


/**************************************************************************
* Main Configuration
**************************************************************************/

// Parse configuration file
//
int Main::ParseConfig(std::string configfile) {
  // Load config file
  CkPrintf("Loading config from %s\n", configfile.c_str());
  YAML::Node config;
  try {
    config = YAML::LoadFile(configfile);
  } catch (YAML::BadFile& e) {
    CkPrintf("  %s\n", e.what());
    return 1;
  }

  // Get configuration
  // Network data file
  try {
    filebase = config["filebase"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  filebase: %s\n", e.what());
    return 1;
  }
  // output filename modifications
  try {
    filemod = config["filemod"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  filemod not defined, defaulting to: \".o\"\n");
    filemod = std::string(".o");
  }
  // Number of data files
  try {
    npdat = config["npdat"].as<idx_t>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  npdat: %s\n", e.what());
    return 1;
  }
  // Number of network parts
  try {
    npnet = config["npnet"].as<idx_t>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  npnet: %s\n", e.what());
    return 1;
  }
  // Maximum simulation time (in ms)
  real_t treal;
  try {
     treal = config["tmax"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = TMAX_DEFAULT;
    CkPrintf("  tmax not defined, defaulting to: %.2g ms\n", treal);
  }
  tmax = (tick_t)(treal*TICKS_PER_MS);
  // Time of a simulation step (in ms)
  try {
     treal = config["tstep"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = TSTEP_DEFAULT;
    CkPrintf("  tstep not defined, defaulting to: %.2g ms\n", treal);
  }
  tstep = (tick_t)(treal*TICKS_PER_MS);
  // Time between checkpoints (in ms)
  try {
     treal = config["tcheck"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = TCHECK_DEFAULT;
    CkPrintf("  tcheck not defined, defaulting to: %.2g ms\n", treal);
  }
  tcheck = (tick_t)(treal*TICKS_PER_MS);
  // RPC port
  try {
     rpcport = config["rpcport"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    rpcport = RPCPORT_DEFAULT;
    CkPrintf("  rpcport not defined, defaulting to: %s\n", rpcport.c_str());
  }

  // Return success
  return 0;
}


/**************************************************************************
* Read in Models
**************************************************************************/

// Parse model file (multiple yaml docs in one file)
//
int Main::ReadModel() {
  // Load model file
  CkPrintf("Loading models from %s.model\n", filebase.c_str());
  YAML::Node modfile;
  try {
    modfile = YAML::LoadAllFromFile(filebase + ".model");
  } catch (YAML::BadFile& e) {
    CkPrintf("  %s\n", e.what());
    return 1;
  }

  // Setup model data
  models.resize(modfile.size());

  // Get model data
  for (std::size_t i = 0; i < modfile.size(); ++i) {
    try {
      // modname
      models[i].modname = modfile[i]["modname"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  modname: %s\n", e.what());
      return 1;
    }
    try {
      // modtype
      models[i].modtype = modfile[i]["modtype"].as<idx_t>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  modtype: %s\n", e.what());
      return 1;
    }
    try {
      // param
      models[i].param = modfile[i]["param"].as<std::vector<real_t>>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  param: %s\n", e.what());
      return 1;
    }
    
    // States are their own 'node'
    YAML::Node state = modfile[i]["state"];
    if (state.size() == 0) {
      CkPrintf("  warning: %s has no state\n", models[i].modname.c_str());
    }
    // Count states and sticks
    models[i].nstate = 0;
    models[i].nstick = 0;
    for (std::size_t j = 0; j < state.size(); ++j) {
      std::string reptype;
      try {
        // reptype
        reptype = state[j]["rep"].as<std::string>();
      } catch (YAML::RepresentationException& e) {
        reptype = std::string("real");
      }
      if (reptype == "tick") {
        ++models[i].nstick;
      }
      else {
        ++models[i].nstate;
      }
    }
  }

  // Return success
  return 0;
}
