/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 *
 * ioyml.C
 * Handles YAML Ain't Markup Language format
 * Configuration and Model information
 */

#include <random>
#include "stacs.h"
#include "build.h"
#include "event.h"
#include "record.h"

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
extern /*readonly*/ idx_t intbal;
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
  // Only read runmode from config if it wasn't passed in as a command line argument
  if (runmode == std::string(RUNMODE_EMPTY)) {
    try {
      runmode = config["runmode"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      runmode = std::string(RUNMODE_DEFAULT);
      CkPrintf("  runmode not defined, defaulting to: %s\n", runmode.c_str());
    }
  }
  // Make sure it's a valid runmode
  if (runmode != std::string(RUNMODE_SIMULATE) && 
      runmode != std::string(RUNMODE_SIMGPU) &&
      runmode != std::string(RUNMODE_BUILDSIM) &&
      runmode != std::string(RUNMODE_BUILD) &&
      runmode != std::string(RUNMODE_REPART) &&
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
  // Load balancing
  try {
    loadbal = config["loadbal"].as<bool>();
  } catch (YAML::RepresentationException& e) {
    loadbal = LOADBAL_DEFAULT;
    CkPrintf("  loadbal not defined, defaulting to: %s\n", (loadbal ? "true" : "false"));
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
  if (runmode == "build") {
    filesave = std::string("");
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
  // Time between repartitioning (in ms)
  try {
    tbalance = config["tbalance"].as<real_t>();
  } catch (YAML::RepresentationException& e) {
    tbalance = TBALANCE_DEFAULT;
    CkPrintf("  tbalance not defined, defaulting to: %.2g ms\n", tbalance);
  }
  intbal = (idx_t)(((tick_t)(tbalance*TICKS_PER_MS))/tstep);
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
  // Count number of graph models (not metadata files)
  idx_t nmodels = 0;
  for (std::size_t i = 0; i < modfile.size(); ++i) {
    // graph type
    std::string graphtype;
    try {
      // type
      graphtype = modfile[i]["type"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  type: %s\n", e.what());
      return 1;
    }
    if (graphtype == "stream" ||
        graphtype == "vertex" ||
        graphtype == "edge") {
      ++nmodels;
    }
  }

  // Setup model data
  //models.resize(modfile.size());
  modelconf.resize(nmodels);
  // Model map
  //std::unordered_map<std::string, std::size_t> modmap;
  modmap.clear();
  modmap[std::string("none")] = 0;
  // Data file names
  datafiles.clear();
  datatypes.clear();

  // Get model data
  idx_t jmod = 0;
  for (std::size_t i = 0; i < modfile.size(); ++i) {
    // graph type
    std::string graphtype;
    try {
      // type
      graphtype = modfile[i]["type"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  type: %s\n", e.what());
      return 1;
    }
    // Set type
    if (graphtype == "stream") {
      modelconf[jmod].graphtype = GRAPHTYPE_STR;
    }
    else if (graphtype == "vertex") {
      modelconf[jmod].graphtype = GRAPHTYPE_VTX;
    }
    else if (graphtype == "edge") {
      modelconf[jmod].graphtype = GRAPHTYPE_EDG;
    }
    else if (graphtype == "record") {
      // don't add records to models
      continue;
    }
    else {
      CkPrintf("  type: '%s' unknown type\n", graphtype.c_str());
      return 1;
    }

    // Graph vertices information
    try {
      // modname
      modelconf[jmod].modname = modfile[i]["modname"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  modname: %s\n", e.what());
      return 1;
    }
    try {
      // modtype
      modelconf[jmod].modtype = modfile[i]["modtype"].as<idx_t>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  modtype: %s\n", e.what());
      return 1;
    }

    // Composite model (should be for vertices only?)
    // TODO: These would ideally allow for the dynamics of multiple
    // models to be combined, but its not actually implemented in build
    if (modelconf[jmod].modtype == 0) {
      // Parts are their own 'node'
      YAML::Node part = modfile[i]["part"];
      // Use port list for model names here
      modelconf[jmod].port.resize(part.size());
      for (std::size_t j = 0; j < part.size(); ++j) {
        try {
          // source model for part
          modelconf[jmod].port[j] = part[j]["source"].as<std::string>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  part source: %s\n", e.what());
          return 1;
        }
        std::vector<std::string> targets;
        try {
          // target list for part
          targets = part[j]["target"].as<std::vector<std::string>>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  part target: %s\n", e.what());
          return 1;
        }
        // Add targets to end of modnames
        for (std::size_t s = 0; s < targets.size(); ++s) {
          std::ostringstream target;
          target << " " << targets[s];
          modelconf[jmod].port[j].append(target.str());
        }
      }
      // Composite models don't have their own parameters
      modelconf[jmod].param.resize(0);
      // To be filled in after all models read in
      modelconf[jmod].nstate = 0;
      modelconf[jmod].nstick = 0;
    }
    // Concrete model
    else {
      // Create list of model indices (for composite models)
      modmap[modelconf[jmod].modname] = jmod+1;

      // Params are their own 'node'
      YAML::Node param = modfile[i]["param"];
      modelconf[jmod].nparam = param.size();
      modelconf[jmod].paramname.resize(param.size());
      modelconf[jmod].param.resize(param.size());
      for (std::size_t j = 0; j < param.size(); ++j) {
        try {
          modelconf[jmod].paramname[j] = param[j]["name"].as<std::string>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  param name: %s\n", e.what());
          return 1;
        }
        try {
          modelconf[jmod].param[j] = param[j]["value"].as<real_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  param value: %s\n", e.what());
          return 1;
        }
      }

      // States are their own 'node'
      YAML::Node state = modfile[i]["state"];
      // Count states and sticks
      modelconf[jmod].nstate = 0;
      modelconf[jmod].nstick = 0;
      for (std::size_t j = 0; j < state.size(); ++j) {
        std::string reptype;
        try {
          // reptype
          reptype = state[j]["rep"].as<std::string>();
        } catch (YAML::RepresentationException& e) {
          reptype = std::string("real");
        }
        if (reptype == "tick") {
          ++modelconf[jmod].nstick;
        }
        else {
          ++modelconf[jmod].nstate;
        }
      }
    
      // preallocate space for non-default parameters
      modelconf[jmod].statename.resize(modelconf[jmod].nstate);
      modelconf[jmod].stateinit.resize(modelconf[jmod].nstate);
      modelconf[jmod].stateparam.resize(modelconf[jmod].nstate);
      modelconf[jmod].stickname.resize(modelconf[jmod].nstick);
      modelconf[jmod].stickinit.resize(modelconf[jmod].nstick);
      modelconf[jmod].stickparam.resize(modelconf[jmod].nstick);
      idx_t jstate = 0;
      idx_t jstick = 0;
    
      // loop through the states
      for (std::size_t j = 0; j < state.size(); ++j) {
        std::string rngtype;
        std::string reptype;
        try {
          // rngtype
          rngtype = state[j]["init"].as<std::string>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  state init: %s\n", e.what());
          return 1;
        }

        // Representation type defaults to real unless specifically set as "tick"
        try {
          // reptype
          reptype = state[j]["rep"].as<std::string>();
        } catch (YAML::RepresentationException& e) {
          reptype = std::string("real");
        }
        if (reptype != "tick") {
          reptype = std::string("real");
        }

        // state names (for underspecified models, relying on model defaults)
        if (reptype == "tick") {
          try {
            modelconf[jmod].stickname[jstick] = state[j]["name"].as<std::string>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state name: %s\n", e.what());
            return 1;
          }
        }
        else {
          try {
            modelconf[jmod].statename[jstate] = state[j]["name"].as<std::string>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state name: %s\n", e.what());
            return 1;
          }
        }

        // based on rng type, get params
        if (rngtype == "constant") {
          if (reptype == "tick") {
            // Constant value
            modelconf[jmod].stickinit[jstick] = RNGTYPE_CONST;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_CONST);
            try {
              // value
              modelconf[jmod].stickparam[jstick][0] = state[j]["value"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state constant value: %s\n", e.what());
              return 1;
            }
            ++jstick;
          }
          else {
            // Constant value
            modelconf[jmod].stateinit[jstate] = RNGTYPE_CONST;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_CONST);
            try {
              // value
              modelconf[jmod].stateparam[jstate][0] = state[j]["value"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state constant value: %s\n", e.what());
              return 1;
            }
            ++jstate;
          }
        }
        else if (rngtype == "uniform") {
          if (reptype == "tick") {
            // Uniform distribution
            modelconf[jmod].stickinit[jstick] = RNGTYPE_UNIF;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_UNIF);
            try {
              // min value
              modelconf[jmod].stickparam[jstick][0] = state[j]["min"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform min: %s\n", e.what());
              return 1;
            }
            try {
              // max value
              modelconf[jmod].stickparam[jstick][1] = state[j]["max"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform max: %s\n", e.what());
              return 1;
            }
            ++jstick;
          }
          else {
            // Uniform distribution
            modelconf[jmod].stateinit[jstate] = RNGTYPE_UNIF;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_UNIF);
            try {
              // min value
              modelconf[jmod].stateparam[jstate][0] = state[j]["min"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform min: %s\n", e.what());
              return 1;
            }
            try {
              // max value
              modelconf[jmod].stateparam[jstate][1] = state[j]["max"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform max: %s\n", e.what());
              return 1;
            }
            ++jstate;
          }
        }
        else if (rngtype == "uniform interval") {
          if (reptype == "tick") {
            // Uniform distribution (intervalled)
            modelconf[jmod].stickinit[jstick] = RNGTYPE_UNINT;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_UNINT);
            try {
              // min value
              modelconf[jmod].stickparam[jstick][0] = state[j]["min"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform min: %s\n", e.what());
              return 1;
            }
            try {
              // max value
              modelconf[jmod].stickparam[jstick][1] = state[j]["max"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform max: %s\n", e.what());
              return 1;
            }
            try {
              // int value
              modelconf[jmod].stickparam[jstick][2] = state[j]["int"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform int: %s\n", e.what());
              return 1;
            }
            ++jstick;
          }
          else {
            // Uniform distribution (intervalled)
            modelconf[jmod].stateinit[jstate] = RNGTYPE_UNINT;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_UNINT);
            try {
              // min value
              modelconf[jmod].stateparam[jstate][0] = state[j]["min"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform min: %s\n", e.what());
              return 1;
            }
            try {
              // max value
              modelconf[jmod].stateparam[jstate][1] = state[j]["max"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform max: %s\n", e.what());
              return 1;
            }
            try {
              // int value
              modelconf[jmod].stateparam[jstate][2] = state[j]["int"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state uniform int: %s\n", e.what());
              return 1;
            }
            ++jstate;
          }
        }
        else if (rngtype == "normal") {
          if (reptype == "tick") {
            // Normal distribution
            modelconf[jmod].stickinit[jstick] = RNGTYPE_NORM;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_NORM);
            try {
              // mean
              modelconf[jmod].stickparam[jstick][0] = state[j]["mean"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state normal mean: %s\n", e.what());
              return 1;
            }
            try {
              // standard deviation
              modelconf[jmod].stickparam[jstick][1] = state[j]["std"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state normal std: %s\n", e.what());
              return 1;
            }
            ++jstick;
          }
          else {
            // Normal distribution
            modelconf[jmod].stateinit[jstate] = RNGTYPE_NORM;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_NORM);
            try {
              // mean
              modelconf[jmod].stateparam[jstate][0] = state[j]["mean"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state normal mean: %s\n", e.what());
              return 1;
            }
            try {
              // standard deviation
              modelconf[jmod].stateparam[jstate][1] = state[j]["std"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state normal std: %s\n", e.what());
              return 1;
            }
            ++jstate;
          }
        }
        else if (rngtype == "bounded normal") {
          if (reptype == "tick") {
            // Normal distribution
            modelconf[jmod].stickinit[jstick] = RNGTYPE_BNORM;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_BNORM);
            try {
              // mean
              modelconf[jmod].stickparam[jstick][0] = state[j]["mean"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded normal mean: %s\n", e.what());
              return 1;
            }
            try {
              // standard deviation
              modelconf[jmod].stickparam[jstick][1] = state[j]["std"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded normal std: %s\n", e.what());
              return 1;
            }
            try {
              // bounds
              modelconf[jmod].stickparam[jstick][2] = state[j]["bound"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded normal bound: %s\n", e.what());
              return 1;
            }
            ++jstick;
          }
          else {
            // Normal distribution
            modelconf[jmod].stateinit[jstate] = RNGTYPE_BNORM;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_BNORM);
            try {
              // mean
              modelconf[jmod].stateparam[jstate][0] = state[j]["mean"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded normal mean: %s\n", e.what());
              return 1;
            }
            try {
              // standard deviation
              modelconf[jmod].stateparam[jstate][1] = state[j]["std"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded normal std: %s\n", e.what());
              return 1;
            }
            try {
              // bounds
              modelconf[jmod].stateparam[jstate][2] = state[j]["bound"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded normal bound: %s\n", e.what());
              return 1;
            }
            ++jstate;
          }
        }
        else if (rngtype == "lower bounded normal") {
          if (reptype == "tick") {
            // Normal distribution
            modelconf[jmod].stickinit[jstick] = RNGTYPE_LBNORM;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_LBNORM);
            try {
              // mean
              modelconf[jmod].stickparam[jstick][0] = state[j]["mean"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded normal mean: %s\n", e.what());
              return 1;
            }
            try {
              // standard deviation
              modelconf[jmod].stickparam[jstick][1] = state[j]["std"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded normal std: %s\n", e.what());
              return 1;
            }
            try {
              // bounds
              modelconf[jmod].stickparam[jstick][2] = state[j]["bound"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded normal bound: %s\n", e.what());
              return 1;
            }
            ++jstick;
          }
          else {
            // Normal distribution
            modelconf[jmod].stateinit[jstate] = RNGTYPE_LBNORM;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_LBNORM);
            try {
              // mean
              modelconf[jmod].stateparam[jstate][0] = state[j]["mean"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded normal mean: %s\n", e.what());
              return 1;
            }
            try {
              // standard deviation
              modelconf[jmod].stateparam[jstate][1] = state[j]["std"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded normal std: %s\n", e.what());
              return 1;
            }
            try {
              // bounds
              modelconf[jmod].stateparam[jstate][2] = state[j]["bound"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded normal bound: %s\n", e.what());
              return 1;
            }
            ++jstate;
          }
        }
        else if (rngtype == "lower bounded lognorm") {
          if (reptype == "tick") {
            // Normal distribution
            modelconf[jmod].stickinit[jstick] = RNGTYPE_LBLOGNORM;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_LBLOGNORM);
            try {
              // mean
              modelconf[jmod].stickparam[jstick][0] = state[j]["mean"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded lognorm mean: %s\n", e.what());
              return 1;
            }
            try {
              // standard deviation
              modelconf[jmod].stickparam[jstick][1] = state[j]["std"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded lognorm std: %s\n", e.what());
              return 1;
            }
            try {
              // bounds
              modelconf[jmod].stickparam[jstick][2] = state[j]["bound"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded lognorm bound: %s\n", e.what());
              return 1;
            }
            try {
              // scale
              modelconf[jmod].stickparam[jstick][3] = state[j]["scale"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded lognorm scale: %s\n", e.what());
              return 1;
            }
            ++jstick;
          }
          else {
            // Normal distribution
            modelconf[jmod].stateinit[jstate] = RNGTYPE_LBLOGNORM;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_LBLOGNORM);
            try {
              // mean
              modelconf[jmod].stateparam[jstate][0] = state[j]["mean"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded lognorm mean: %s\n", e.what());
              return 1;
            }
            try {
              // standard deviation
              modelconf[jmod].stateparam[jstate][1] = state[j]["std"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded lognorm std: %s\n", e.what());
              return 1;
            }
            try {
              // bounds
              modelconf[jmod].stateparam[jstate][2] = state[j]["bound"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded lognorm bound: %s\n", e.what());
              return 1;
            }
            try {
              // scale
              modelconf[jmod].stateparam[jstate][3] = state[j]["scale"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state lower bounded lognorm scale: %s\n", e.what());
              return 1;
            }
            ++jstate;
          }
        }
        else if (rngtype == "linear") {
          if (reptype == "tick") {
            // Proportional to distance
            modelconf[jmod].stickinit[jstick] = RNGTYPE_LIN;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_LIN);
            try {
              // scale
              modelconf[jmod].stickparam[jstick][0] = state[j]["scale"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state linear scale: %s\n", e.what());
              return 1;
            }
            try {
              // offset
              modelconf[jmod].stickparam[jstick][1] = state[j]["offset"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  warning: state linear offset not defined,\n"
                  "           defaulting to 0\n");
              modelconf[jmod].stickparam[jstick][1] = 0.0;
            }
            ++jstick;
          }
          else {
            // Proportional to distance
            modelconf[jmod].stateinit[jstate] = RNGTYPE_LIN;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_LIN);
            try {
              // scale
              modelconf[jmod].stateparam[jstate][0] = state[j]["scale"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state linear scale: %s\n", e.what());
              return 1;
            }
            try {
              // offset
              modelconf[jmod].stateparam[jstate][1] = state[j]["offset"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  warning: state linear offset not defined,\n"
                  "           defaulting to 0\n");
              modelconf[jmod].stateparam[jstate][1] = 0.0;
            }
            ++jstate;
          }
        }
        else if (rngtype == "bounded linear") {
          if (reptype == "tick") {
            // Proportional to distance
            modelconf[jmod].stickinit[jstick] = RNGTYPE_BLIN;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_BLIN);
            try {
              // scale
              modelconf[jmod].stickparam[jstick][0] = state[j]["scale"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded linear scale: %s\n", e.what());
              return 1;
            }
            try {
              // offset
              modelconf[jmod].stickparam[jstick][1] = state[j]["offset"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  warning: state linear offset not defined,\n"
                  "           defaulting to 0\n");
              modelconf[jmod].stickparam[jstick][1] = 0.0;
            }
            try {
              // min value
              modelconf[jmod].stickparam[jstick][2] = state[j]["min"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded linear min: %s\n", e.what());
              return 1;
            }
            try {
              // max value
              modelconf[jmod].stickparam[jstick][3] = state[j]["max"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded linear max: %s\n", e.what());
              return 1;
            }
            ++jstick;
          }
          else {
            // Proportional to distance
            modelconf[jmod].stateinit[jstate] = RNGTYPE_BLIN;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_BLIN);
            try {
              // scale
              modelconf[jmod].stateparam[jstate][0] = state[j]["scale"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded linear scale: %s\n", e.what());
              return 1;
            }
            try {
              // offset
              modelconf[jmod].stateparam[jstate][1] = state[j]["offset"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  warning: state linear offset not defined,\n"
                  "           defaulting to 0\n");
              modelconf[jmod].stateparam[jstate][1] = 0.0;
            }
            try {
              // min value
              modelconf[jmod].stateparam[jstate][2] = state[j]["min"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded linear min: %s\n", e.what());
              return 1;
            }
            try {
              // max value
              modelconf[jmod].stateparam[jstate][3] = state[j]["max"].as<real_t>();
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state bounded linear max: %s\n", e.what());
              return 1;
            }
            ++jstate;
          }
        }
        else if (rngtype == "file") {
          // File type information
          std::string filetype;
          try {
            filetype = state[j]["filetype"].as<std::string>();
          } catch (YAML::RepresentationException& e) {
            // Default
            filetype = FT_DEFAULT;
          }
          if (filetype == "csv-sparse") {
            datatypes.push_back(FT_CSV_SPARSE);
          }
          else if (filetype == "csv-dense") {
            datatypes.push_back(FT_CSV_DENSE);
          }
          else {
            CkPrintf("  error: '%s' unsupported file type\n", filetype.c_str());
            return 1;
          }
          // Get the filename and index into the datafile list
          if (reptype == "tick") {
            // From file
            modelconf[jmod].stickinit[jstick] = RNGTYPE_FILE;
            modelconf[jmod].stickparam[jstick].resize(RNGPARAM_FILE);
            try {
              // filename
              // Add this filename to a list, and set the
              // model parameter to the filename's index in that list
              // TODO: check that the filename isn't already on the list
              //       (e.g. using the same matrix for connections as weights)
              modelconf[jmod].stickparam[jstick][0] = (real_t) datafiles.size();
              datafiles.push_back(state[j]["filename"].as<std::string>());
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state file name: %s\n", e.what());
              return 1;
            }
            ++jstick;
          }
          else {
            // From file
            modelconf[jmod].stateinit[jstate] = RNGTYPE_FILE;
            modelconf[jmod].stateparam[jstate].resize(RNGPARAM_FILE);
            try {
              // filename
              modelconf[jmod].stateparam[jstate][0] = (real_t) datafiles.size();
              datafiles.push_back(state[j]["filename"].as<std::string>());
            } catch (YAML::RepresentationException& e) {
              CkPrintf("  state file name: %s\n", e.what());
              return 1;
            }
            ++jstate;
          }
        }
        else {
          CkPrintf("  error: '%s' unknown state init\n", rngtype.c_str());
          return 1;
        }
      }
      // some sanity checking
      CkAssert(jstate == modelconf[jmod].nstate);
      CkAssert(jstick == modelconf[jmod].nstick);

      // Ports are their own 'node'
      YAML::Node port = modfile[i]["port"];
      modelconf[jmod].port.resize(port.size());
      for (std::size_t j = 0; j < port.size(); ++j) {
        try {
          modelconf[jmod].port[j] = port[j]["value"].as<std::string>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  port: %s\n", e.what());
          return 1;
        }
      }
    }

    // Polychronization (active models)
    // TODO: move non-standard simulation stuff to a different function?
    modelconf[jmod].grpactive = false;
    for (std::size_t j = 0; j < grpactives.size(); ++j) {
      if (modelconf[jmod].modname == grpactives[j]) {
        modelconf[jmod].grpactive = true;
        break;
      }
    }
    // Polychronization (mother vertices)
    modelconf[jmod].grpmother = false;
    for (std::size_t j = 0; j < grpmothers.size(); ++j) {
      if (modelconf[jmod].modname == grpmothers[j]) {
        modelconf[jmod].grpmother = true;
        break;
      }
    }
    // Polychronization (anchor edges)
    modelconf[jmod].grpanchor = false;
    for (std::size_t j = 0; j < grpanchors.size(); ++j) {
      if (modelconf[jmod].modname == grpanchors[j]) {
        modelconf[jmod].grpanchor = true;
        break;
      }
    }
    // update for next model
    ++jmod;
  }
  CkAssert(jmod == nmodels);

  // Fill in composite models and print out information
  for (std::size_t i = 0; i < modelconf.size(); ++i) {
    if (modelconf[i].modtype == 0) {
      // Print out model information
      std::string modelparts;
      for (std::size_t j = 0; j < modelconf[i].port.size(); ++j) {
        // Get first model name out of string
        std::string modname = modelconf[i].port[j].substr(0, modelconf[i].port[j].find(' '));
        modelconf[i].nstate += modelconf[modmap[modname]].nstate;
        modelconf[i].nstick += modelconf[modmap[modname]].nstick;
        std::ostringstream part;
        part << " " << modname;
        modelparts.append(part.str());
      }
      CkPrintf("  Model: %d   Name: %s   Type: %d   States: %u   Parts:%s\n",
          i+1, modelconf[i].modname.c_str(), modelconf[i].modtype, modelconf[i].nstate + modelconf[i].nstick,
          modelparts.c_str());
    }
    else {
      // Print out model information
      std::string modelports;
      // collect ports
      for (std::size_t j = 0; j < modelconf[i].port.size(); ++j) {
        std::ostringstream port;
        port << " " << modelconf[i].port[j];
        modelports.append(port.str());
      }
      // TODO: modtype to name for base model (may move into netdata for this?)
      CkPrintf("  Model: %d   Name: %s   Type: %d   States: %u   Params: %u   Ports:%s\n",
          i+1, modelconf[i].modname.c_str(), modelconf[i].modtype, modelconf[i].nstate + modelconf[i].nstick,
          modelconf[i].param.size(), (modelports == "") ? " None" : modelports.c_str());
    }
  }

  // Loop through modfile again to collect metamodels (e.g. records)
  idx_t jrec = 0;
  evtloglist.clear();
  recordlist.clear();
  for (std::size_t i = 0; i < modfile.size(); ++i) {
    // graph type
    std::string graphtype;
    try {
      // type
      graphtype = modfile[i]["type"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  type: %s\n", e.what());
      return 1;
    }
    if (graphtype == "record") {
      std::vector<std::string> evtlogstring;
      // Events to log are presented as a list of strings
      try {
        evtlogstring = modfile[i]["events"].as<std::vector<std::string>>();
      } catch (YAML::RepresentationException& e) {
        // Do nothing
      }
      if (evtlogstring.size()) {
        // check if already defined
        if (evtloglist.size()) {
          CkPrintf("  warning: event logging already defined\n");
        }
        // TODO: allow these to be recorded w.r.t. modname too
        // populate event log list
        evtloglist.resize(evtlogstring.size());
        for (std::size_t j = 0; j < evtlogstring.size(); ++j) {
          if (evtlogstring[j] == "spike") {
            evtloglist[j] = EVENT_SPIKE;
          }
          else if (evtlogstring[j] == "stim") {
            evtloglist[j] = EVENT_STIM;
          }
          else if (evtlogstring[j] == "current") {
            evtloglist[j] = EVENT_CURRENT;
          }
          else if (evtlogstring[j] == "count") {
            evtloglist[j] = EVENT_COUNT;
          }
          else if (evtlogstring[j] == "clamp") {
            evtloglist[j] = EVENT_CLAMP;
          }
          else if (evtlogstring[j] == "set rate") {
            evtloglist[j] = EVENT_RATE;
          }
          else {
            CkPrintf("  event type '%s' not supported for logging\n", evtlogstring[j].c_str());
          }
        }
      }
      // Record tracking (state, stick, coord) w.r.t. modname
      YAML::Node probe = modfile[i]["probes"];
      for (std::size_t j = 0; j < probe.size(); ++j) {
        recordlist.push_back(probe_t());
        // recording interval
        real_t treal;
        try {
          treal= probe[j]["tfreq"].as<real_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  recording interval (tfreq): %s\n", e.what());
          return 1;
        }
        recordlist[jrec].tfreq = (tick_t)(treal*TICKS_PER_MS);
        // modname
        std::string name;
        try {
          name = probe[j]["modname"].as<std::string>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  recording model name: %s\n", e.what());
          return 1;
        }
        // convert modname to modidx
        if (modmap.find(name) == modmap.end()) {
          CkPrintf("  error: model %s not defined\n", name.c_str());
          return 1;
        } else {
          recordlist[jrec].modidx = modmap[name];
        }
        // states to track
        try {
          recordlist[jrec].statename = probe[j]["state"].as<std::string>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  recording model state: %s\n", e.what());
          return 1;
        }
        ++jrec;
      }
    }
  }

  // Return success
  return 0;
}


// Read in graph information
//
int Main::ReadGraph() {
  // Load model file
  CkPrintf("Loading graph from %s/%s.graph\n", netwkdir.c_str(), filebase.c_str());
  YAML::Node graphfile;
  try {
    graphfile = YAML::LoadFile(netwkdir + "/" + filebase + ".graph");
  } catch (YAML::BadFile& e) {
    CkPrintf("  %s\n", e.what());
    return 1;
  }

  // Streams and vertices
  YAML::Node stream = graphfile["stream"];
  YAML::Node vertex = graphfile["vertex"];
  if (vertex.size() + stream.size() == 0) {
    CkPrintf("  error: graph has no vertices\n");
    return 1;
  }

  // TODO: Collect vertices and edges into a graph_t structure?
  // preallocate space
  idx_t jvtx = 0;
  idx_t nvtx = vertex.size() + stream.size();
  vertices.clear();
  vertices.resize(nvtx);
  
  // loop through the streams
  for (std::size_t i = 0; i < stream.size(); ++i) {
    std::string name;
    try {
      // modname
      name = stream[i]["modname"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  stream modname: %s\n", e.what());
      return 1;
    }
    if (modmap.find(name) == modmap.end()) {
      CkPrintf("  error: model %s not defined\n", name.c_str());
      return 1;
    }
    else {
      vertices[jvtx].modidx = modmap[name];
    }
    // order
    vertices[jvtx].order = 1;
    // shape
    vertices[jvtx].shape = VTXSHAPE_POINT;
    vertices[jvtx].param.resize(VTXPARAM_POINT);
    try {
      // coord
      vertices[jvtx].coord = stream[i]["coord"].as<std::vector<real_t>>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  warning: coord not defined, defaulting to [0.0, 0.0, 0.0]\n");
      vertices[jvtx].coord = {0.0, 0.0, 0.0};
    }
    if (vertices[jvtx].coord.size() != 3) {
      CkPrintf("  error: stream coord dimensions\n");
      return 1;
    }
    ++jvtx;
  }

  // loop through the vertices
  for (std::size_t i = 0; i < vertex.size(); ++i) {
    std::string name;
    try {
      // modname
      name = vertex[i]["modname"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  vertex modname: %s\n", e.what());
      return 1;
    }
    if (modmap.find(name) == modmap.end()) {
      CkPrintf("  error: model %s not defined\n", name.c_str());
      return 1;
    }
    else {
      vertices[jvtx].modidx = modmap[name];
    }
    try {
      // order
      vertices[jvtx].order = vertex[i]["order"].as<idx_t>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  vertex order: %s\n", e.what());
      return 1;
    }
    try {
      // shape
      name = vertex[i]["shape"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  vertex shape: %s\n", e.what());
      return 1;
    }
    if (name == "point") {
      vertices[jvtx].shape = VTXSHAPE_POINT;
      vertices[jvtx].param.resize(VTXPARAM_POINT);
    }
    else if (name == "circle") {
      vertices[jvtx].shape = VTXSHAPE_CIRCLE;
      vertices[jvtx].param.resize(VTXPARAM_CIRCLE);
      try {
        // radius
        vertices[jvtx].param[0] = vertex[i]["radius"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  vertex circle radius: %s\n", e.what());
        return 1;
      }
    }
    else if (name == "sphere") {
      vertices[jvtx].shape = VTXSHAPE_SPHERE;
      vertices[jvtx].param.resize(VTXPARAM_SPHERE);
      try {
        // radius
        vertices[jvtx].param[0] = vertex[i]["radius"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  vertex sphere radius: %s\n", e.what());
        return 1;
      }
    }
    else if (name == "sphere surface") {
      vertices[jvtx].shape = VTXSHAPE_SPHERE_SURFACE;
      vertices[jvtx].param.resize(VTXPARAM_SPHERE_SURFACE);
      try {
        // radius
        vertices[jvtx].param[0] = vertex[i]["radius"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  vertex sphere surface radius: %s\n", e.what());
        return 1;
      }
    }
    else if (name == "line") {
      vertices[jvtx].shape = VTXSHAPE_LINE;
      vertices[jvtx].param.resize(VTXPARAM_LINE);
      try {
        // length
        vertices[jvtx].param[0] = vertex[i]["length"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  vertex line length: %s\n", e.what());
        return 1;
      }
    }
    else if (name == "rectangle") {
      vertices[jvtx].shape = VTXSHAPE_RECT;
      vertices[jvtx].param.resize(VTXPARAM_RECT);
      try {
        // width
        vertices[jvtx].param[0] = vertex[i]["width"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  vertex rectangle width: %s\n", e.what());
        return 1;
      }
      try {
        // height
        vertices[jvtx].param[1] = vertex[i]["height"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  vertex rectangle height: %s\n", e.what());
        return 1;
      }
    }
    else {
      CkPrintf("  error: '%s' unknown shape\n", name.c_str());
      return 1;
    }
    try {
      // coord
      vertices[jvtx].coord = vertex[i]["coord"].as<std::vector<real_t>>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  warning: coord not defined, defaulting to [0.0, 0.0, 0.0]\n");
      vertices[jvtx].coord = {0.0, 0.0, 0.0};
    }
    if (vertices[jvtx].coord.size() != 3) {
      CkPrintf("  error: vertex coord dimensions\n");
      return 1;
    }
    ++jvtx;
  }
  CkAssert(jvtx == nvtx);
  
  // Edges
  YAML::Node edge = graphfile["edge"];
  if (edge.size() == 0) {
    CkPrintf("  error: graph has no edges\n");
    return 1;
  }

  // preallocate space
  edges.clear();
  edges.resize(edge.size());

  // loop through the edges
  for (std::size_t i = 0; i < edges.size(); ++i) {
    std::string name;
    std::vector<std::string> names;
    try {
      // source
      name = edge[i]["source"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  edge source: %s\n", e.what());
      return 1;
    }
    if (modmap.find(name) == modmap.end()) {
      CkPrintf("  error: model %s not defined\n", name.c_str());
      return 1;
    }
    else {
      edges[i].source = modmap[name];
    }
    edges[i].target.clear();
    try {
      // target(s)
      names = edge[i]["target"].as<std::vector<std::string>>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  edge targets: %s\n", e.what());
      return 1;
    }
    for (std::size_t j = 0; j < names.size(); ++j) {
      if (modmap.find(names[j]) == modmap.end()) {
        CkPrintf("  error: model %s not defined\n", name.c_str());
        return 1;
      }
      else {
        edges[i].target.push_back(modmap[names[j]]);
      }
    }
    try {
      // modname
      name = edge[i]["modname"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  edge modname: %s\n", e.what());
      return 1;
    }
    if (modmap.find(name) == modmap.end()) {
      CkPrintf("  error: model %s not defined\n", name.c_str());
      return 1;
    }
    else {
      edges[i].modidx = modmap[name];
    }
    try {
      // cutoff
      edges[i].cutoff = edge[i]["cutoff"].as<real_t>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  warning: cutoff not defined, defaulting to none\n");
      edges[i].cutoff = 0.0;
    }
    // Distance function
    std::string distype;
    try {
      // distance computation
      distype = edge[i]["dist"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      // Default to euclidean distance
      distype = std::string("euclidean");
    }
    if (distype == "sphere") {
      edges[i].distype = DISTYPE_SPHERE;
      edges[i].distparam.resize(DISTPARAM_SPHERE);
      // Additional information for computing distances on spheres
      // TODO: there's no checks for if this radius matches the one that
      //       was used to instantiate the neural population (yet)
      try {
        // TODO: assumes the sphere is centered at 0,0,0
        edges[i].distparam[0] = edge[i]["radius"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  connect distance sphere radius: %s\n", e.what());
        return 1;
      }
    }
    else if (distype == "periodic rectangle") {
      edges[i].distype = DISTYPE_PERIRECT;
      edges[i].distparam.resize(DISTPARAM_PERIRECT);
      // Additional information for computing on periodic rectangles
      try {
        // TODO: assumes the bottom left (x, y) corner is at 0,0,0
        //       also doesn't allow for (x, z) or (y, z) rectangles
        edges[i].distparam[0] = edge[i]["width"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  connect distance rect width: %s\n", e.what());
        return 1;
      }
      try {
        edges[i].distparam[1] = edge[i]["height"].as<real_t>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  connect distance rect height: %s\n", e.what());
        return 1;
      }
    }
    else {
      // No additional information needed for Euclidean
      edges[i].distype = DISTYPE_EUCLIDEAN;
      edges[i].distparam.resize(DISTPARAM_EUCLIDEAN);
    }

    // Connection types are their own 'node'
    YAML::Node conn = edge[i]["connect"];
    if (conn.size() == 0) {
      CkPrintf("  error: edge %s has no connections\n", name.c_str());
    }

    // preallocate space
    edges[i].conntype.resize(conn.size());
    edges[i].probparam.resize(conn.size());
    edges[i].maskparam.resize(conn.size());

    // loop through the connections
    for (std::size_t j = 0; j < conn.size(); ++j) {
      std::string conntype;
      try {
        // conntype
        conntype = conn[j]["type"].as<std::string>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  connect type: %s\n", e.what());
        return 1;
      }
      // based on connection type, get params
      if (conntype == "uniform") {
        // Randomly connect
        edges[i].conntype[j] = CONNTYPE_UNIF;
        edges[i].probparam[j].resize(PROBPARAM_UNIF);
        edges[i].maskparam[j].resize(MASKPARAM_UNIF);
        try {
          // probability threshold
          edges[i].probparam[j][0] = conn[j]["prob"].as<real_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect random prob: %s\n", e.what());
          return 1;
        }
      }
      else if (conntype == "sigmoid") {
        // Connect according to sigmoid property
        edges[i].conntype[j] = CONNTYPE_SIG;
        edges[i].probparam[j].resize(PROBPARAM_SIG);
        edges[i].maskparam[j].resize(MASKPARAM_SIG);
        try {
          // maximum probability
          edges[i].probparam[j][0] = conn[j]["maxprob"].as<real_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect sigmoid maxprob: %s\n", e.what());
          return 1;
        }
        try {
          // sigmoid midpoint
          edges[i].probparam[j][1] = conn[j]["midpoint"].as<real_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect sigmoid midpoint: %s\n", e.what());
          return 1;
        }
        try {
          // sigmoid slope
          edges[i].probparam[j][2] = conn[j]["slope"].as<real_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect sigmoid slope: %s\n", e.what());
          return 1;
        }
      }
      else if (conntype == "index") {
        // Connect by index of source and target
        edges[i].conntype[j] = CONNTYPE_IDX;
        edges[i].probparam[j].resize(PROBPARAM_IDX);
        edges[i].maskparam[j].resize(MASKPARAM_IDX);
        // By index requires single source and target
        if (edges[i].target.size() > 1) {
          CkPrintf("  connect index is single target only\n");
          return 1;
        }
        // source and target order
        for (std::size_t v = 0; v < vertices.size(); ++v) {
          if (edges[i].source == vertices[v].modidx) {
            edges[i].maskparam[j][0] = vertices[v].order;
          }
          if (edges[i].target[0] == vertices[v].modidx) {
            edges[i].maskparam[j][1] = vertices[v].order;
          }
        }
        try {
          // connection source multiplier
          edges[i].maskparam[j][2] = conn[j]["srcmul"].as<idx_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect index srcmul: %s\n", e.what());
          return 1;
        }
        try {
          // connection source offset
          edges[i].maskparam[j][3] = conn[j]["srcoff"].as<idx_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect index srcoff: %s\n", e.what());
          return 1;
        }
      }
      else if (conntype == "sample") {
        // Connect by sampling index from source population
        edges[i].conntype[j] = CONNTYPE_SMPL;
        edges[i].probparam[j].resize(PROBPARAM_SMPL);
        edges[i].maskparam[j].resize(MASKPARAM_SMPL);
        // By index requires single source and target
        if (edges[i].target.size() > 1) {
          CkPrintf("  connect sample is single target only\n");
          return 1;
        }
        // source order
        for (std::size_t v = 0; v < vertices.size(); ++v) {
          if (edges[i].source == vertices[v].modidx) {
            edges[i].maskparam[j][0] = vertices[v].order;
          }
        }
        try {
          // connection source sample number
          edges[i].maskparam[j][1] = conn[j]["order"].as<idx_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect sample order: %s\n", e.what());
          return 1;
        }
      }
      else if (conntype == "sample norm") {
        // Connect by sampling index from source population
        edges[i].conntype[j] = CONNTYPE_SMPL_NORM;
        edges[i].probparam[j].resize(PROBPARAM_SMPL_NORM);
        edges[i].maskparam[j].resize(MASKPARAM_SMPL_NORM);
        // By index requires single source and target
        if (edges[i].target.size() > 1) {
          CkPrintf("  connect sample is single target only\n");
          return 1;
        }
        try {
          // maximum probability
          edges[i].probparam[j][0] = conn[j]["variance"].as<real_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect sample norm variance: %s\n", e.what());
          return 1;
        }
        // source order
        for (std::size_t v = 0; v < vertices.size(); ++v) {
          if (edges[i].source == vertices[v].modidx) {
            edges[i].maskparam[j][0] = vertices[v].order;
          }
        }
        try {
          // connection source sample number
          edges[i].maskparam[j][1] = conn[j]["order"].as<idx_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect sample norm order: %s\n", e.what());
          return 1;
        }
      }
      else if (conntype == "sample anti-norm") {
        // Connect by sampling index from source population
        edges[i].conntype[j] = CONNTYPE_SMPL_ANORM;
        edges[i].probparam[j].resize(PROBPARAM_SMPL_ANORM);
        edges[i].maskparam[j].resize(MASKPARAM_SMPL_ANORM);
        // By index requires single source and target
        if (edges[i].target.size() > 1) {
          CkPrintf("  connect sample is single target only\n");
          return 1;
        }
        try {
          // maximum probability
          edges[i].probparam[j][0] = conn[j]["variance"].as<real_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect sample anti-norm variance: %s\n", e.what());
          return 1;
        }
        try {
          // maximum probability
          edges[i].probparam[j][1] = conn[j]["variance-y"].as<real_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect sample anti-norm variance-y: %s\n", e.what());
          return 1;
        }
        // source order
        for (std::size_t v = 0; v < vertices.size(); ++v) {
          if (edges[i].source == vertices[v].modidx) {
            edges[i].maskparam[j][0] = vertices[v].order;
          }
        }
        try {
          // connection source sample number
          edges[i].maskparam[j][1] = conn[j]["order"].as<idx_t>();
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect sample anti-norm order: %s\n", e.what());
          return 1;
        }
      }
      else if (conntype == "file") {
        // Connect by whatever is in a datafile
        edges[i].conntype[j] = CONNTYPE_FILE;
        edges[i].probparam[j].resize(PROBPARAM_FILE);
        edges[i].maskparam[j].resize(MASKPARAM_FILE);
        // By file requires single source and target
        if (edges[i].target.size() > 1) {
          CkPrintf("  connect file is single target only\n");
          return 1;
        }
        // source and target order
        for (std::size_t v = 0; v < vertices.size(); ++v) {
          if (edges[i].source == vertices[v].modidx) {
            edges[i].maskparam[j][0] = vertices[v].order;
          }
          if (edges[i].target[0] == vertices[v].modidx) {
            edges[i].maskparam[j][1] = vertices[v].order;
          }
        }
        try {
          // file name
          edges[i].probparam[j][0] = (real_t) datafiles.size();
          datafiles.push_back(conn[j]["filename"].as<std::string>());
        } catch (YAML::RepresentationException& e) {
          CkPrintf("  connect file name: %s\n", e.what());
          return 1;
        }
        // File type information
        std::string filetype;
        try {
          filetype = conn[j]["filetype"].as<std::string>();
        } catch (YAML::RepresentationException& e) {
          // Default
          filetype = FT_DEFAULT;
        }
        if (filetype == "csv-sparse") {
          datatypes.push_back(FT_CSV_SPARSE);
        }
        else if (filetype == "csv-dense") {
          datatypes.push_back(FT_CSV_DENSE);
        }
        else {
          CkPrintf("  error: '%s' unsupported file type\n", filetype.c_str());
          return 1;
        }
      }
      else {
        CkPrintf("  error: '%s' unknown connection type\n", conntype.c_str());
        return 1;
      }
    }
  }

  // Check that only one type of edge may exist between any two given vertices
  // TODO: Enable multiple edges between vertices one day
  std::vector<std::vector<idx_t>> connections;
  connections.resize(modelconf.size()+1);
  for (std::size_t i = 0; i < edges.size(); ++i) {
    for (std::size_t j = 0; j < edges[i].target.size(); ++j) {
      // add source target pairs
      connections[edges[i].source].push_back(edges[i].target[j]);
      // Sanity check that there are no 'none' sources or targets
      CkAssert(edges[i].source);
      CkAssert(edges[i].target[j]);
    }
  }
  // check for duplicates (connections[0] is 'none' model)
  for (std::size_t i = 1; i < connections.size(); ++i) {
    std::sort(connections[i].begin(), connections[i].end());
    for (std::size_t j = 1; j < connections[i].size(); ++j) {
      if (connections[i][j] == connections[i][j-1]) {
        CkPrintf("  error: multiple connection types between vertices\n"
                 "         %s to %s not allowed (yet)\n",
                 modelconf[i-1].modname.c_str(), modelconf[connections[i][j]-1].modname.c_str());
        return 1;
      }
    }
  }
  
  // Return success
  return 0;
}
