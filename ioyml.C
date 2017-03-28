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
extern /*readonly*/ idx_t npdat;
extern /*readonly*/ idx_t npnet;
extern /*readonly*/ std::string filebase;
extern /*readonly*/ std::string filein;
extern /*readonly*/ std::string fileout;
extern /*readonly*/ std::string recbase;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ tick_t tcheck;
extern /*readonly*/ tick_t trecord;
extern /*readonly*/ tick_t tdisplay;
extern /*readonly*/ idx_t evtcal;
extern /*readonly*/ bool simpause;
extern /*readonly*/ idx_t rngseed;
extern /*readonly*/ std::string rpcport;
extern /*readonly*/ idx_t runmode;


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
  // Network data file
  try {
    filebase = config["filebase"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  filebase: %s\n", e.what());
    return 1;
  }
  // input filename modifications
  try {
    filein = config["filein"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  filein not defined, defaulting to: \"\"\n");
    filein = std::string("");
  }
  // output filename modifications
  try {
    fileout = config["fileout"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  fileout not defined, defaulting to: \".o\"\n");
    fileout = std::string(".o");
  }
  // record output base
  try {
    recbase = config["recbase"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  recbase not defined, defaulting to: \"%s\"\n", filebase.c_str());
    recbase = filebase;
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
    treal = TDISP_DEFAULT;
    CkPrintf("  tdisplay not defined, defaulting to: %.2g ms\n", treal);
  }
  tdisplay = (tick_t)(treal*TICKS_PER_MS);
  // Event queue calendar days per year
  try {
    evtcal = config["evtcal"].as<idx_t>();
  } catch (YAML::RepresentationException& e) {
    evtcal = EVTCAL_DEFAULT;
    CkPrintf("  evtcal not defined, defaulting to: %d iters\n", evtcal);
  }
  // Simulation startup
  try {
    simpause = config["simpause"].as<bool>();
  } catch (YAML::RepresentationException& e) {
    simpause = SIMPAUSE_DEFAULT;
    CkPrintf("  simpause not defined, defaulting to: %s\n", (simpause ? "true" : "false"));
  }
  // Random number seed
  try {
     rngseed = config["rngseed"].as<idx_t>();
  } catch (YAML::RepresentationException& e) {
    std::random_device rd;
    rngseed = rd();
    CkPrintf("  rngseed not defined, seeding with: %" PRIidx "\n", rngseed);
  }
  // RPC port
  try {
    rpcport = config["rpcport"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    rpcport = RPCPORT_DEFAULT;
    CkPrintf("  rpcport not defined, defaulting to: %s\n", rpcport.c_str());
  }
  // Run mode
  std::string modestring;
  try {
    modestring = config["runmode"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    modestring = std::string("sim");
    CkPrintf("  runmode not defined, defaulting to: %s\n", modestring.c_str());
  }
  if (modestring == "sim") {
    runmode = RUNMODE_SIM;
  } else if (modestring == "png") {
    runmode = RUNMODE_PNG;
  } else {
    CkPrintf("  runmode is invalid, defaulting to: sim\n");
    runmode = RUNMODE_SIM;
  }
  // Active models
  actives.clear();
  // Identifiers are their own 'node'
  YAML::Node active = config["active"];
  if (active.size() == 0) {
    CkPrintf("  warning: network has no active models\n");
  }
  actives.resize(active.size());
  for (std::size_t i = 0; i < active.size(); ++i) {
    try {
      actives[i] = active[i].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  active: %s\n", e.what());
      return 1;
    }
  }
  // PNG models
  pngmods.clear();
  // Identifiers are their own 'node'
  YAML::Node pngmod = config["pngmod"];
  if (pngmod.size() == 0) {
    CkPrintf("  warning: network has no png models\n");
  }
  pngmods.resize(pngmod.size());
  for (std::size_t i = 0; i < pngmod.size(); ++i) {
    try {
      pngmods[i] = pngmod[i].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  pngmod: %s\n", e.what());
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

    // Flags
    models[i].modflag = MODFLAG_DEFAULT;

    // Active models
    if (actives.size()) {
      for (std::size_t j = 0; j < actives.size(); ++j) {
        if (models[i].modname == actives[j]) {
          models[i].modflag |= MODFLAG_ACTIVE;
          break;
        }
      }
    } else {
      models[i].modflag |= MODFLAG_ACTIVE;
    }
    
    // PNG models
    if (pngmods.size()) {
      for (std::size_t j = 0; j < pngmods.size(); ++j) {
        if (models[i].modname == pngmods[j]) {
          models[i].modflag |= MODFLAG_PNGMOD;
          break;
        }
      }
    }
  }

  // Return success
  return 0;
}
