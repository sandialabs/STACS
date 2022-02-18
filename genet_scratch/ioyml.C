/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 *
 * ioyml.C
 * Handles YAML Ain't Markup Language format
 */

#include "genet.h"

// Using yaml-cpp (specification version 1.2)
#include "yaml-cpp/yaml.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ unsigned int randseed;
extern /*readonly*/ std::string netwkdir;
extern /*readonly*/ idx_t netparts;
extern /*readonly*/ int netfiles;
extern /*readonly*/ std::string filebase;
extern /*readonly*/ std::string filesave;


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
  // Random number seed
  try {
     randseed = config["randseed"].as<unsigned int>();
  } catch (YAML::RepresentationException& e) {
    std::random_device rd;
    randseed = rd();
    CkPrintf("  randseed not defined, seeding with: %u\n", randseed);
  }
  // Network data directory
  try {
    netwkdir = config["netwkdir"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  netwkdir: %s\n", e.what());
    return 1;
  }
  // Number of network parts
  try {
    netparts = config["netparts"].as<idx_t>();
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
  // Network data file
  try {
    filebase = config["filebase"].as<std::string>();
  } catch (YAML::RepresentationException& e) {
    CkPrintf("  filebase: %s\n", e.what());
    return 1;
  }
  if (mode == "build") {
    filesave = std::string("");
  }
  else {
    try {
      filesave = config["filesave"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  filesave not defined, defaulting to: \".o\"\n");
      filesave = std::string(".o");
    }
  }

  // Return success
  return 0;
}


/**************************************************************************
* Read in Models
**************************************************************************/

// Parse model file (multiple yaml docs in one file)
// TODO: make it possible to read from multiple files
//
int Main::ReadModel() {
  // Load model file
  CkPrintf("Loading models from %s/%s.model\n", netwkdir.c_str(), filebase.c_str());
  YAML::Node modfile;
  try {
    modfile = YAML::LoadAllFromFile(netwkdir + "/" + filebase + ".model");
  } catch (YAML::BadFile& e) {
    CkPrintf("  %s\n", e.what());
    return 1;
  }

  // Setup model data
  models.resize(modfile.size());
  // Data file names
  datafiles.clear();
  // Model map
  modmap.clear();
  modmap[std::string("none")] = 0;
  
  // Graph types
  graphtype.resize(GRAPHTYPE_NTYPE);
  graphtype[GRAPHTYPE_STR] = std::string("stream");
  graphtype[GRAPHTYPE_VTX] = std::string("vertex");
  graphtype[GRAPHTYPE_EDG] = std::string("edge");

  // Get model data
  for (std::size_t i = 0; i < modfile.size(); ++i) {
    std::string type;

    try {
      // type
      type = modfile[i]["type"].as<std::string>();
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  type: %s\n", e.what());
      return 1;
    }
    // Set type
    if (type == "stream") {
      models[i].type = GRAPHTYPE_STR;
    }
    else if (type == "vertex") {
      models[i].type = GRAPHTYPE_VTX;
    }
    else if (type == "edge") {
      models[i].type = GRAPHTYPE_EDG;
    }
    else {
      CkPrintf("  type: '%s' unknown type\n", type.c_str());
      return 1;
    }
    try {
      // modname
      models[i].modname = modfile[i]["modname"].as<std::string>();
      // modmap
      modmap[models[i].modname] = i+1;
    } catch (YAML::RepresentationException& e) {
      CkPrintf("  modname: %s\n", e.what());
      return 1;
    }

    // States are their own 'node'
    YAML::Node state = modfile[i]["state"];
    if (state.size() == 0) {
      CkPrintf("  warning: %s has no state\n", models[i].modname.c_str());
    }

    // Count sticks
    idx_t nstate = 0;
    idx_t nstick = 0;
    for (std::size_t j = 0; j < state.size(); ++j) {
      std::string reptype;
      try {
        // reptype
        reptype = state[j]["rep"].as<std::string>();
      } catch (YAML::RepresentationException& e) {
        reptype = std::string("real");
      }
      if (reptype == "tick") {
        ++nstick;
      }
      else {
        ++nstate;
      }
    }

    // preallocate space
    models[i].statetype.resize(nstate);
    models[i].stateparam.resize(nstate);
    models[i].sticktype.resize(nstick);
    models[i].stickparam.resize(nstick);
    idx_t jstate = 0;
    idx_t jstick = 0;

    // loop through the states
    for (std::size_t j = 0; j < state.size(); ++j) {
      std::string rngtype;
      std::string reptype;
      try {
        // rngtype
        rngtype = state[j]["type"].as<std::string>();
      } catch (YAML::RepresentationException& e) {
        CkPrintf("  state type: %s\n", e.what());
        return 1;
      }

      try {
        // reptype
        reptype = state[j]["rep"].as<std::string>();
      } catch (YAML::RepresentationException& e) {
        reptype = std::string("real");
      }
      if (reptype != "tick") {
        reptype = std::string("real");
      }

      // based on rng type, get params
      if (rngtype == "constant") {
        if (reptype == "tick") {
          // Constant value
          models[i].sticktype[jstick] = RNGTYPE_CONST;
          models[i].stickparam[jstick].resize(RNGPARAM_CONST);
          try {
            // value
            models[i].stickparam[jstick][0] = state[j]["value"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state constant value: %s\n", e.what());
            return 1;
          }
          ++jstick;
        }
        else {
          // Constant value
          models[i].statetype[jstate] = RNGTYPE_CONST;
          models[i].stateparam[jstate].resize(RNGPARAM_CONST);
          try {
            // value
            models[i].stateparam[jstate][0] = state[j]["value"].as<real_t>();
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
          models[i].sticktype[jstick] = RNGTYPE_UNIF;
          models[i].stickparam[jstick].resize(RNGPARAM_UNIF);
          try {
            // min value
            models[i].stickparam[jstick][0] = state[j]["min"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state uniform min: %s\n", e.what());
            return 1;
          }
          try {
            // max value
            models[i].stickparam[jstick][1] = state[j]["max"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state uniform max: %s\n", e.what());
            return 1;
          }
          ++jstick;
        }
        else {
          // Uniform distribution
          models[i].statetype[jstate] = RNGTYPE_UNIF;
          models[i].stateparam[jstate].resize(RNGPARAM_UNIF);
          try {
            // min value
            models[i].stateparam[jstate][0] = state[j]["min"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state uniform min: %s\n", e.what());
            return 1;
          }
          try {
            // max value
            models[i].stateparam[jstate][1] = state[j]["max"].as<real_t>();
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
          models[i].sticktype[jstick] = RNGTYPE_UNINT;
          models[i].stickparam[jstick].resize(RNGPARAM_UNINT);
          try {
            // min value
            models[i].stickparam[jstick][0] = state[j]["min"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state uniform min: %s\n", e.what());
            return 1;
          }
          try {
            // max value
            models[i].stickparam[jstick][1] = state[j]["max"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state uniform max: %s\n", e.what());
            return 1;
          }
          try {
            // int value
            models[i].stickparam[jstick][2] = state[j]["int"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state uniform int: %s\n", e.what());
            return 1;
          }
          ++jstick;
        }
        else {
          // Uniform distribution (intervalled)
          models[i].statetype[jstate] = RNGTYPE_UNINT;
          models[i].stateparam[jstate].resize(RNGPARAM_UNINT);
          try {
            // min value
            models[i].stateparam[jstate][0] = state[j]["min"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state uniform min: %s\n", e.what());
            return 1;
          }
          try {
            // max value
            models[i].stateparam[jstate][1] = state[j]["max"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state uniform max: %s\n", e.what());
            return 1;
          }
          try {
            // int value
            models[i].stateparam[jstate][2] = state[j]["int"].as<real_t>();
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
          models[i].sticktype[jstick] = RNGTYPE_NORM;
          models[i].stickparam[jstick].resize(RNGPARAM_NORM);
          try {
            // mean
            models[i].stickparam[jstick][0] = state[j]["mean"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state normal mean: %s\n", e.what());
            return 1;
          }
          try {
            // standard deviation
            models[i].stickparam[jstick][1] = state[j]["std"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state normal std: %s\n", e.what());
            return 1;
          }
          ++jstick;
        }
        else {
          // Normal distribution
          models[i].statetype[jstate] = RNGTYPE_NORM;
          models[i].stateparam[jstate].resize(RNGPARAM_NORM);
          try {
            // mean
            models[i].stateparam[jstate][0] = state[j]["mean"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state normal mean: %s\n", e.what());
            return 1;
          }
          try {
            // standard deviation
            models[i].stateparam[jstate][1] = state[j]["std"].as<real_t>();
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
          models[i].sticktype[jstick] = RNGTYPE_BNORM;
          models[i].stickparam[jstick].resize(RNGPARAM_BNORM);
          try {
            // mean
            models[i].stickparam[jstick][0] = state[j]["mean"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded normal mean: %s\n", e.what());
            return 1;
          }
          try {
            // standard deviation
            models[i].stickparam[jstick][1] = state[j]["std"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded normal std: %s\n", e.what());
            return 1;
          }
          try {
            // bounds
            models[i].stickparam[jstick][2] = state[j]["bound"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded normal bound: %s\n", e.what());
            return 1;
          }
          ++jstick;
        }
        else {
          // Normal distribution
          models[i].statetype[jstate] = RNGTYPE_BNORM;
          models[i].stateparam[jstate].resize(RNGPARAM_BNORM);
          try {
            // mean
            models[i].stateparam[jstate][0] = state[j]["mean"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded normal mean: %s\n", e.what());
            return 1;
          }
          try {
            // standard deviation
            models[i].stateparam[jstate][1] = state[j]["std"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded normal std: %s\n", e.what());
            return 1;
          }
          try {
            // bounds
            models[i].stateparam[jstate][2] = state[j]["bound"].as<real_t>();
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
          models[i].sticktype[jstick] = RNGTYPE_LBNORM;
          models[i].stickparam[jstick].resize(RNGPARAM_LBNORM);
          try {
            // mean
            models[i].stickparam[jstick][0] = state[j]["mean"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state lower bounded normal mean: %s\n", e.what());
            return 1;
          }
          try {
            // standard deviation
            models[i].stickparam[jstick][1] = state[j]["std"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state lower bounded normal std: %s\n", e.what());
            return 1;
          }
          try {
            // bounds
            models[i].stickparam[jstick][2] = state[j]["bound"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state lower bounded normal bound: %s\n", e.what());
            return 1;
          }
          ++jstick;
        }
        else {
          // Normal distribution
          models[i].statetype[jstate] = RNGTYPE_LBNORM;
          models[i].stateparam[jstate].resize(RNGPARAM_LBNORM);
          try {
            // mean
            models[i].stateparam[jstate][0] = state[j]["mean"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state lower bounded normal mean: %s\n", e.what());
            return 1;
          }
          try {
            // standard deviation
            models[i].stateparam[jstate][1] = state[j]["std"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state lower bounded normal std: %s\n", e.what());
            return 1;
          }
          try {
            // bounds
            models[i].stateparam[jstate][2] = state[j]["bound"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state lower bounded normal bound: %s\n", e.what());
            return 1;
          }
          ++jstate;
        }
      }
      else if (rngtype == "linear") {
        if (reptype == "tick") {
          // Proportional to distance
          models[i].sticktype[jstick] = RNGTYPE_LIN;
          models[i].stickparam[jstick].resize(RNGPARAM_LIN);
          try {
            // scale
            models[i].stickparam[jstick][0] = state[j]["scale"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state linear scale: %s\n", e.what());
            return 1;
          }
          try {
            // offset
            models[i].stickparam[jstick][1] = state[j]["offset"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  warning: state linear offset not defined,\n"
                "           defaulting to 0\n");
            models[i].stickparam[jstick][1] = 0.0;
          }
          ++jstick;
        }
        else {
          // Proportional to distance
          models[i].statetype[jstate] = RNGTYPE_LIN;
          models[i].stateparam[jstate].resize(RNGPARAM_LIN);
          try {
            // scale
            models[i].stateparam[jstate][0] = state[j]["scale"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state linear scale: %s\n", e.what());
            return 1;
          }
          try {
            // offset
            models[i].stateparam[jstate][1] = state[j]["offset"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  warning: state linear offset not defined,\n"
                "           defaulting to 0\n");
            models[i].stateparam[jstate][1] = 0.0;
          }
          ++jstate;
        }
      }
      else if (rngtype == "bounded linear") {
        if (reptype == "tick") {
          // Proportional to distance
          models[i].sticktype[jstick] = RNGTYPE_BLIN;
          models[i].stickparam[jstick].resize(RNGPARAM_BLIN);
          try {
            // scale
            models[i].stickparam[jstick][0] = state[j]["scale"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded linear scale: %s\n", e.what());
            return 1;
          }
          try {
            // offset
            models[i].stickparam[jstick][1] = state[j]["offset"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  warning: state linear offset not defined,\n"
                "           defaulting to 0\n");
            models[i].stickparam[jstick][1] = 0.0;
          }
          try {
            // min value
            models[i].stickparam[jstick][2] = state[j]["min"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded linear min: %s\n", e.what());
            return 1;
          }
          try {
            // max value
            models[i].stickparam[jstick][3] = state[j]["max"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded linear max: %s\n", e.what());
            return 1;
          }
          ++jstick;
        }
        else {
          // Proportional to distance
          models[i].statetype[jstate] = RNGTYPE_BLIN;
          models[i].stateparam[jstate].resize(RNGPARAM_BLIN);
          try {
            // scale
            models[i].stateparam[jstate][0] = state[j]["scale"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded linear scale: %s\n", e.what());
            return 1;
          }
          try {
            // offset
            models[i].stateparam[jstate][1] = state[j]["offset"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  warning: state linear offset not defined,\n"
                "           defaulting to 0\n");
            models[i].stateparam[jstate][1] = 0.0;
          }
          try {
            // min value
            models[i].stateparam[jstate][2] = state[j]["min"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded linear min: %s\n", e.what());
            return 1;
          }
          try {
            // max value
            models[i].stateparam[jstate][3] = state[j]["max"].as<real_t>();
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state bounded linear max: %s\n", e.what());
            return 1;
          }
          ++jstate;
        }
      }
      else if (rngtype == "file") {
        if (reptype == "tick") {
          // From file
          models[i].sticktype[jstick] = RNGTYPE_FILE;
          models[i].stickparam[jstick].resize(RNGPARAM_FILE);
          try {
            // filename
            // Add this filename to a list, and set the
            // model parameter to the filename's index in that list
            // TODO: check that the filename isn't already on the list
            //       (e.g. using the same matrix for connections as weights)
            models[i].stickparam[jstick][0] = (real_t) datafiles.size();
            datafiles.push_back(state[j]["filename"].as<std::string>());
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state file name: %s\n", e.what());
            return 1;
          }
          ++jstick;
        }
        else {
          // From file
          models[i].statetype[jstate] = RNGTYPE_FILE;
          models[i].stateparam[jstate].resize(RNGPARAM_FILE);
          try {
            // filename
            models[i].stateparam[jstate][0] = (real_t) datafiles.size();
            datafiles.push_back(state[j]["filename"].as<std::string>());
          } catch (YAML::RepresentationException& e) {
            CkPrintf("  state file name: %s\n", e.what());
            return 1;
          }
          ++jstate;
        }
      }
      else {
        CkPrintf("  error: '%s' unknown state type\n", rngtype.c_str());
        return 1;
      }
    }
    CkAssert(jstate == nstate);
    CkAssert(jstick == nstick);
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
        CkPrintf("  vertex circle radius: %s\n", e.what());
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
  connections.resize(models.size()+1);
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
                 models[i-1].modname.c_str(), models[connections[i][j]-1].modname.c_str());
        return 1;
      }
    }
  }
  
  // Return success
  return 0;
}
