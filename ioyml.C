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
extern /*readonly*/ unsigned randseed;
extern /*readonly*/ std::string netwkdir;
extern /*readonly*/ int netparts;
extern /*readonly*/ int netfiles;
extern /*readonly*/ std::string filebase;
extern /*readonly*/ std::string fileload;
extern /*readonly*/ std::string filesave;
extern /*readonly*/ std::string recordir;
extern /*readonly*/ std::string groupdir;
extern /*readonly*/ tick_t tstep;
extern /*readonly*/ idx_t nevtday;
extern /*readonly*/ idx_t intdisp;
extern /*readonly*/ idx_t intrec;
extern /*readonly*/ idx_t intsave;
extern /*readonly*/ tick_t tmax;
extern /*readonly*/ tick_t tepisode;
extern /*readonly*/ idx_t episodes;
extern /*readonly*/ int grpminlen;
extern /*readonly*/ tick_t grpmaxdur;


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
  // Simulation
  // Run mode
  try {
    runmode = config["runmode"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    runmode = std::string(RUNMODE_DEFAULT);
    CkPrintf("  runmode not defined, defaulting to: %s\n", runmode.c_str());
  }
  if (runmode != std::string(RUNMODE_SIMULATE) && 
      runmode != std::string(RUNMODE_FINDGROUP) && 
      runmode != std::string(RUNMODE_ESTIMATE)) {
    runmode = std::string(RUNMODE_DEFAULT);
    CkPrintf("  runmode is invalid, defaulting to: %s\n", runmode.c_str());
  }
  // Random number seed
  try {
     randseed = config["randseed"].as<unsigned>();
  } catch (YAML::RepresentationException& e) {
    std::random_device rd;
    randseed = rd();
    CkPrintf("  randseed not defined, seeding with: %u\n", randseed);
  }
  // Network plasticity
  try {
    plastic = config["plastic"].as<bool>();
  } catch (YAML::RepresentationException& e) {
    plastic = PLASTIC_DEFAULT;
    CkPrintf("  plastic not defined, defaulting to: %s\n", (plastic ? "true" : "false"));
  }
  // Episodic simulation
  try {
    episodic = config["episodic"].as<bool>();
  } catch (YAML::RepresentationException& e) {
    episodic = EPISODIC_DEFAULT;
    CkPrintf("  episodic not defined, defaulting to: %s\n", (episodic ? "true" : "false"));
  }
#ifdef STACS_WITH_YARP
  // RPC port
  try {
    rpcport = config["rpcport"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    rpcport = std::string(RPCPORT_DEFAULT);
    CkPrintf("  rpcportname not defined, defaulting to: %s\n", rpcport.c_str());
  }
  // Start simulation paused
  try {
    rpcpause = config["rpcpause"].as<bool>();
  } catch (YAML::RepresentationException& e) {
    rpcpause = RPCPAUSE_DEFAULT;
    CkPrintf("  rpcpause not defined, defaulting to: %s\n", (rpcpause ? "true" : "false"));
  }
#endif

  // Network
  // Network data directory
  try {
    netwkdir = config["netwkdir"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  netwkdir: %s\n", e.what());
    return 1;
  }
  // Number of network parts
  try {
    netparts = config["netparts"].as<int>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  netparts: %s\n", e.what());
    return 1;
  }
  // Number of data files
  try {
    netfiles = config["netfiles"].as<int>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  netfiles: %s\n", e.what());
    return 1;
  }
  // Network file base name
  try {
    filebase = config["filebase"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  filebase: %s\n", e.what());
    return 1;
  }
  // input filename modifications
  try {
    fileload = config["fileload"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    fileload = std::string(FILELOAD_DEFAULT);
    CkPrintf("  fileload not defined, defaulting to: \"%s\"\n", fileload.c_str());
  }
  // output filename modifications
  try {
    filesave = config["filesave"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    filesave = std::string(FILESAVE_DEFAULT);
    CkPrintf("  filesave not defined, defaulting to: \"%s\"\n", filesave.c_str());
  }
  // Records output directory
  try {
    recordir = config["recordir"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    recordir = std::string(RECORDIR_DEFAULT);
    CkPrintf("  recordir not defined, defaulting to: \"%s\"\n", recordir.c_str());
  }
  // Groups output directory
  try {
    groupdir = config["groupdir"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    groupdir = std::string(GROUPDIR_DEFAULT);
    CkPrintf("  groupdir not defined, defaulting to: \"%s\"\n", groupdir.c_str());
  }
  
  // Timing
  real_t treal;
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
    teventq = config["teventq"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    teventq = TEVENTQ_DEFAULT;
    CkPrintf("  teventq not defined, defaulting to: %.2g ms\n", teventq);
  }
  nevtday = (idx_t)(((tick_t)(teventq*TICKS_PER_MS))/tstep) + 1;
  // How often to display the simulation time (in ms)
  try {
    tdisplay = config["tdisplay"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    tdisplay = TDISPLAY_DEFAULT;
    CkPrintf("  tdisplay not defined, defaulting to: %.2g ms\n", tdisplay);
  }
  intdisp = (idx_t)(((tick_t)(tdisplay*TICKS_PER_MS))/tstep);
  // Time between recording points (in ms)
  try {
    trecord = config["trecord"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    trecord = TRECORD_DEFAULT;
    CkPrintf("  trecord not defined, defaulting to: %.2g ms\n", trecord);
  }
  intrec = (idx_t)(((tick_t)(trecord*TICKS_PER_MS))/tstep);
  // Time between checkpoints (in ms)
  try {
    tsave = config["tsave"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    tsave = TSAVE_DEFAULT;
    CkPrintf("  tsave not defined, defaulting to: %.2g ms\n", tsave);
  }
  intsave = (idx_t)(((tick_t)(tsave*TICKS_PER_MS))/tstep);
  // Maximum simulation time (in ms)
  try {
    treal = config["tmax"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = TMAX_DEFAULT;
    if (!episodic) {
      CkPrintf("  tmax not defined, defaulting to: %.2g ms\n", treal);
    }
  }
  tmax = (tick_t)(treal*TICKS_PER_MS);
  // Time per episode
  try {
    treal = config["tepisode"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = TEPISODE_DEFAULT;
    if (episodic) {
      CkPrintf("  tepisode not defined, defaulting to: %.2g ms\n", treal);
    }
  }
  tepisode = (tick_t)(treal*TICKS_PER_MS);
  // Number of episodes
  try {
    episodes = config["episodes"].as<idx_t>();
  } catch (YAML::RepresentationException& e) {
    episodes = EPISODES_DEFAULT;
    if (episodic) {
      CkPrintf("  episodes not defined, defaulting to: %" PRIidx " eps\n", episodes);
    }
  }
  // Modifications to episodic simulation
  if (episodic) {
    // Display iteration with episodes
    intdisp = (idx_t)(tepisode/tstep);
    // Save network according to episode boundaries
    idx_t savediv = intsave/intdisp;
    intsave = savediv*intdisp;
  }

  // Polychronous groups
  // Active models
  grpactives.clear();
  // Identifiers are their own 'node'
  YAML::Node grpactive = config["grpactive"];
  grpactives.resize(grpactive.size());
  for (std::size_t i = 0; i < grpactive.size(); ++i) {
    try {
      grpactives[i] = grpactive[i].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  grpactive: %s\n", e.what());
      return 1;
    }
  }
  // Mother vertices
  grpmothers.clear();
  // Identifiers are their own 'node'
  YAML::Node grpmother = config["grpmother"];
  grpmothers.resize(grpmother.size());
  for (std::size_t i = 0; i < grpmother.size(); ++i) {
    try {
      grpmothers[i] = grpmother[i].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  grpmother: %s\n", e.what());
      return 1;
    }
  }
  // Anchor edges
  grpanchors.clear();
  // Identifiers are their own 'node'
  YAML::Node grpanchor = config["grpanchor"];
  grpanchors.resize(grpanchor.size());
  for (std::size_t i = 0; i < grpanchor.size(); ++i) {
    try {
      grpanchors[i] = grpanchor[i].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  grpanchor: %s\n", e.what());
      return 1;
    }
  }
  // Minimum group path length
  try {
    grpminlen = config["grpminlen"].as<int>();
  } catch (YAML::RepresentationException& e) {
    grpminlen = GRPMINLEN_DEFAULT;
  }
  // Maximum group duration
  try {
    treal = config["grpmaxdur"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    treal = GRPMAXDUR_DEFAULT;
  }
  grpmaxdur = (tick_t)(treal*TICKS_PER_MS);
  // Minimum evaluated vertex (fraction of network)
  try {
    grpvtxminreal = config["grpvtxmin"].as<realidx_t>();
  } catch (YAML::RepresentationException& e) {
    grpvtxminreal = 0.0;
  }
  // Maximum evaluated vertex (fraction of network)
  try {
    grpvtxmaxreal = config["grpvtxmax"].as<realidx_t>();
  } catch (YAML::RepresentationException& e) {
    grpvtxmaxreal = 1.0;
  }
  if (grpvtxminreal < 0.0 || grpvtxminreal > 1.0 || grpvtxmaxreal < 0.0 || grpvtxmaxreal > 1.0) {
    CkPrintf("  grpvtxmin/max out of bounds (0.0 to 1.0)\n");
    return 1;
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
  CkPrintf("Reading model information\n");// from %s/%s.model\n", netwkdir.c_str(), filebase.c_str());
  YAML::Node modfile;
  try {
    modfile = YAML::LoadAllFromFile(netwkdir + "/" + filebase + ".model");
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
    models[i].grpactive = false;
    for (std::size_t j = 0; j < grpactives.size(); ++j) {
      if (models[i].modname == grpactives[j]) {
        models[i].grpactive = true;
        break;
      }
    }
    // Polychronization (mother vertices)
    models[i].grpmother = false;
    for (std::size_t j = 0; j < grpmothers.size(); ++j) {
      if (models[i].modname == grpmothers[j]) {
        models[i].grpmother = true;
        break;
      }
    }
    // Polychronization (anchor edges)
    models[i].grpanchor = false;
    for (std::size_t j = 0; j < grpanchors.size(); ++j) {
      if (models[i].modname == grpanchors[j]) {
        models[i].grpanchor = true;
        break;
      }
    }
  
    // Print out model information
    std::string modelports;
    // collect ports
    for (idx_t j = 0; j < models[i].port.size(); ++j) {
      std::ostringstream port;
      port << " " << models[i].port[j];
      modelports.append(port.str());
    }
    // TODO: modtype to name for base model
    CkPrintf("  Model: %d   Name: %s   Type: %d   States: %u   Params: %u   Ports:%s\n",
        i+1, models[i].modname.c_str(), models[i].modtype, models[i].nstate + models[i].nstick,
        models[i].param.size(), (modelports == "") ? " None" : modelports.c_str());
  }

  // Return success
  return 0;
}
