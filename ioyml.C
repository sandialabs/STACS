/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 *
 * ioyml.C
 * Handles YAML Ain't Markup Language format
 * Configuration and Model information
 */

#include <random>
#include "stacs.h"

// Using yaml-cpp (specification version 1.2)
#include "yaml-cpp/yaml.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ idx_t npdat;
extern /*readonly*/ idx_t npnet;
extern /*readonly*/ std::string filedir;
extern /*readonly*/ std::string filebase;
extern /*readonly*/ std::string filemod;
extern /*readonly*/ std::string recordir;
extern /*readonly*/ idx_t rngseed;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ tick_t tqueue;
extern /*readonly*/ tick_t tcheck;
extern /*readonly*/ tick_t trecord;
extern /*readonly*/ tick_t tdisplay;
extern /*readonly*/ idx_t equeue;


/**************************************************************************
* Main Configuration
**************************************************************************/

// Parse configuration file
//
int Main::ReadConfig(std::string configfile) {
  // Load configuration file
  CkPrintf("Reading config from %s\n", configfile.c_str());
  YAML::Node config;
  try {
    config = YAML::LoadFile(configfile);
  } catch (YAML::BadFile& e) {
    CkPrintf("  %s\n", e.what());
    return 1;
  }

  // Process configuration information

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
  
  // Network data directory
  try {
    filedir = config["filedir"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  filedir: %s\n", e.what());
    return 1;
  }
  // Records output directory
  try {
    recordir = config["recordir"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  recordir not defined, defaulting to: \"%s\"\n", filedir.c_str());
    recordir = filedir;
  }
  // Network file base name
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
  
  // Random number seed
  try {
     rngseed = config["rngseed"].as<idx_t>();
  } catch (YAML::RepresentationException& e) {
    std::random_device rd;
    rngseed = rd();
    CkPrintf("  rngseed not defined, seeding with: %" PRIidx "\n", rngseed);
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
  // Time of standard event queue (in ms)
  try {
    treal = config["tqueue"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = TQUEUE_DEFAULT;
    CkPrintf("  tqueue not defined, defaulting to: %.2g ms\n", treal);
  }
  tqueue = (tick_t)(treal*TICKS_PER_MS);
  // Time between checkpoints (in ms)
  try {
    treal = config["tcheck"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = TCHECK_DEFAULT;
    CkPrintf("  tcheck not defined, defaulting to: %.2g ms\n", treal);
  }
  tcheck = (tick_t)(treal*TICKS_PER_MS);
  // Time between recording points (in ms)
  try {
    treal = config["trecord"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = TRECORD_DEFAULT;
    CkPrintf("  trecord not defined, defaulting to: %.2g ms\n", treal);
  }
  trecord = (tick_t)(treal*TICKS_PER_MS);
  // How often to display the simulation time (in ms)
  try {
    treal = config["tdisplay"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = TDISPLAY_DEFAULT;
    CkPrintf("  tdisplay not defined, defaulting to: %.2g ms\n", treal);
  }
  tdisplay = (tick_t)(treal*TICKS_PER_MS);
  // Slots in event queue (calendar days per year)
  equeue = (tqueue / tstep) + 1;
  
  // Run mode
  try {
    runmode = config["runmode"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    runmode = std::string(RUNMODE_DEFAULT);
    CkPrintf("  runmode not defined, defaulting to: %s\n", runmode.c_str());
  }
  if (runmode != RUNMODE_SIM && runmode != RUNMODE_PNG) {
    runmode = std::string(RUNMODE_DEFAULT);
    CkPrintf("  runmode is invalid, defaulting to: %s\n", runmode.c_str());
  }
  // Network plasticity
  try {
    plasticity = config["plasticity"].as<bool>();
  } catch (YAML::RepresentationException& e) {
    plasticity = PLASTICITY_DEFAULT;
    CkPrintf("  plasticity not defined, defaulting to: %s\n", (plasticity ? "true" : "false"));
  }
#ifdef STACS_WITH_YARP
  // RPC port
  try {
    rpcport = config["rpcportname"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    rpcport = RPCPORTNAME_DEFAULT;
    CkPrintf("  rpcportname not defined, defaulting to: %s\n", rpcport.c_str());
  }
  // Start simulation paused
  try {
    startpaused = config["startpaused"].as<bool>();
  } catch (YAML::RepresentationException& e) {
    startpaused = STARTPAUSED_DEFAULT;
    CkPrintf("  startpaused not defined, defaulting to: %s\n", (startpaused ? "true" : "false"));
  }
#endif

  // Polychronization (active models)
  pngactives.clear();
  // Identifiers are their own 'node'
  YAML::Node pngactive = config["pngactive"];
  pngactives.resize(pngactive.size());
  for (std::size_t i = 0; i < pngactive.size(); ++i) {
    try {
      pngactives[i] = pngactive[i].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  pngactive: %s\n", e.what());
      return 1;
    }
  }
  // Polychronization (mothers)
  pngmothers.clear();
  // Identifiers are their own 'node'
  YAML::Node pngmother = config["pngmother"];
  pngmothers.resize(pngmother.size());
  for (std::size_t i = 0; i < pngmother.size(); ++i) {
    try {
      pngmothers[i] = pngmother[i].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  pngmother: %s\n", e.what());
      return 1;
    }
  }
  // Polychronization (anchors)
  pnganchors.clear();
  // Identifiers are their own 'node'
  YAML::Node pnganchor = config["pnganchor"];
  pnganchors.resize(pnganchor.size());
  for (std::size_t i = 0; i < pnganchor.size(); ++i) {
    try {
      pnganchors[i] = pnganchor[i].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  pnganchor: %s\n", e.what());
      return 1;
    }
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
  CkPrintf("Reading model information\n");// from %s/%s.model\n", filedir.c_str(), filebase.c_str());
  YAML::Node modfile;
  try {
    modfile = YAML::LoadAllFromFile(filedir + "/" + filebase + ".model");
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
    
    // Params are their own 'node'
    YAML::Node param = modfile[i]["param"];
    if (param.size() == 0) {
      CkPrintf("  warning: %s has no parameters\n", models[i].modname.c_str());
    }
    models[i].param.resize(param.size());
    for (std::size_t j = 0; j < param.size(); ++j) {
      try {
        models[i].param[j] = param[j]["value"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  param: %s\n", e.what());
        return 1;
      }
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

    // Ports are their own 'node'
    YAML::Node port = modfile[i]["port"];
    if (port.size() == 0) {
      CkPrintf("  warning: %s has no ports\n", models[i].modname.c_str());
    }
    models[i].port.resize(port.size());
    for (std::size_t j = 0; j < port.size(); ++j) {
      try {
        models[i].port[j] = port[j]["value"].as<std::string>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  port: %s\n", e.what());
        return 1;
      }
    }

    // Polychronization (active models)
    models[i].pngactive = false;
    for (std::size_t j = 0; j < pngactives.size(); ++j) {
      if (models[i].modname == pngactives[j]) {
        models[i].pngactive = true;
        break;
      }
    }
    // Polychronization (mothers)
    models[i].pngmother = false;
    for (std::size_t j = 0; j < pngmothers.size(); ++j) {
      if (models[i].modname == pngmothers[j]) {
        models[i].pngmother = true;
        break;
      }
    }
    // Polychronization (anchors)
    models[i].pnganchor = false;
    for (std::size_t j = 0; j < pnganchors.size(); ++j) {
      if (models[i].modname == pnganchors[j]) {
        models[i].pnganchor = true;
        break;
      }
    }
  }

  // Return success
  return 0;
}
