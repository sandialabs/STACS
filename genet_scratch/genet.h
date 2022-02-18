/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_GENET_H__
#define __STACS_GENET_H__

#include <mpi.h>

#include <algorithm>
#include <list>
#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include <cctype>
#include <vector>
#include <unordered_map>

#include "typedefs.h"
#include "timing.h"

#include <mpi.h>
#include "mpi-interoperate.h"
#include "genet-mpi.h"

#include "pup_stl.h"
#include "genet.decl.h"


#define GRAPHTYPE_NTYPE 3

#define GRAPHTYPE_STR   0
#define GRAPHTYPE_VTX   1
#define GRAPHTYPE_EDG   2

#define RNGTYPE_NRNG    11
#define REPTYPE_REAL    0
#define REPTYPE_TICK    1

// TODO: reorder these numbers sometime
#define RNGTYPE_CONST   0
#define RNGTYPE_UNIF    1
#define RNGTYPE_UNINT   2
#define RNGTYPE_NORM    3
#define RNGTYPE_BNORM   4
#define RNGTYPE_LBNORM  5
#define RNGTYPE_LIN     6
#define RNGTYPE_LBLIN   7
#define RNGTYPE_UBLIN   8
#define RNGTYPE_BLIN    9
#define RNGTYPE_FILE    10

#define RNGPARAM_CONST  1
#define RNGPARAM_UNIF   2
#define RNGPARAM_UNINT  3
#define RNGPARAM_NORM   2
#define RNGPARAM_BNORM  3
#define RNGPARAM_LBNORM 3
#define RNGPARAM_LIN    2
#define RNGPARAM_LBLIN  3
#define RNGPARAM_UBLIN  3
#define RNGPARAM_BLIN   4
#define RNGPARAM_FILE   1

#define VTXSHAPE_POINT  0
#define VTXSHAPE_CIRCLE 1
#define VTXSHAPE_SPHERE 2
#define VTXSHAPE_RECT   3

#define VTXPARAM_POINT  0
#define VTXPARAM_CIRCLE 1
#define VTXPARAM_SPHERE 1
#define VTXPARAM_RECT   2

#define CONNTYPE_UNIF   0
#define PROBPARAM_UNIF  1
#define MASKPARAM_UNIF  0

#define CONNTYPE_SIG    1
#define PROBPARAM_SIG   3
#define MASKPARAM_SIG   0

#define CONNTYPE_IDX    2
#define PROBPARAM_IDX   0
#define MASKPARAM_IDX   4

#define CONNTYPE_FILE   3
#define PROBPARAM_FILE  1
#define MASKPARAM_FILE  2

#define EVENT_SPIKE     0


/**************************************************************************
* Charm++ Messages
**************************************************************************/

// Network distribution
//
void registerNetDist(void);
CkReductionMsg *netDist(int nMsg, CkReductionMsg **msgs);

// Model Information
//
#define MSG_Model 11
class mModel : public CMessage_mModel {
  public:
    idx_t *type;        // type of model (vertex/edge)
    idx_t *xmodname;    // model name prefix
    char *modname;      // model name
    idx_t *xstatetype;    // state generation prefix
    idx_t *xsticktype;    // stick generation prefix
    idx_t *statetype;     // state generation type
    idx_t *sticktype;     // stick representation type
    real_t *stateparam;   // state generation parameters
    real_t *stickparam;   // stick generation parameters
    idx_t *xdatafiles;  // prefix sum for filenames
    char *datafiles;    // filenames (concatenated)
    idx_t nmodel;
    idx_t nstateparam;
    idx_t nstickparam;
    idx_t ndatafiles;
};

#define MSG_Graph 17
class mGraph : public CMessage_mGraph {
  public:
    idx_t *vtxmodidx;     // Which vertex to build (modidx from modmap)
    idx_t *vtxorder;      // How many of each vertex to build
    idx_t *vtxshape;      // What shape to build vertices
    idx_t *xvtxparam;     // shape parameters prefix
    real_t *vtxparam;     // shape parameters
    real_t *vtxcoord;     // coordinate parameters
    idx_t *edgsource;     // source modidx of edge
    idx_t *xedgtarget;    // target prefix
    idx_t *edgtarget;     // target modidx of edge
    idx_t *edgmodidx;     // modidx of edge
    real_t *edgcutoff;    // cutoff distance of connection
    idx_t *xedgconntype;  // connection type prefix
    idx_t *edgconntype;   // connection type (computes probability threshold)
    idx_t *medgprobparam; // connection probability sizes
    real_t *edgprobparam; // connection probability parameters
    idx_t *medgmaskparam; // connection mask sizes
    idx_t *edgmaskparam;  // connection mask parameters
    idx_t nvtx;
    idx_t nvtxparam;
    idx_t nedg;
    idx_t nedgtarget;
    idx_t nedgconntype;
    idx_t nedgprobparam;
    idx_t nedgmaskparam;
};

#define MSG_Conn 6
class mConn : public CMessage_mConn {
  public:
    idx_t *vtxmodidx;   // vertex model
    idx_t *vtxordidx;   // vertex model order
    real_t *xyz;        // vertex coordinates
    idx_t *xadj;        // prefix for adjacency
    idx_t *adjcy;       // adjacent vertices
    idx_t *edgmodidx;   // edge models
    idx_t datidx;
    idx_t nvtx;
};

#define MSG_Metis 2
class mMetis : public CMessage_mMetis {
  public:
    idx_t *vtxdist; // number of vertices in data
    idx_t *edgdist; // number of edges in data
};
  
#define MSG_Part 14
class mPart : public CMessage_mPart {
  public:
    idx_t *vtxidx;
    idx_t *vtxmodidx;
    real_t *xyz;
    idx_t *xadj;
    idx_t *adjcy;
    idx_t *edgmodidx;
    real_t *state;
    tick_t *stick;
    idx_t *xevent;
    tick_t *diffuse;
    idx_t *type;
    idx_t *source;
    idx_t *index;
    real_t *data;
    idx_t datidx;
    idx_t prtidx;
    idx_t nvtx;
    idx_t nstate;
    idx_t nstick;
    idx_t nevent;
};
  
#define MSG_Order 2
class mOrder : public CMessage_mOrder {
  public:
    idx_t *vtxidxold;
    idx_t *vtxidxnew;
    idx_t datidx;
    idx_t nvtx;
};
  


/**************************************************************************
* Data Structures
**************************************************************************/

// Model
//
struct model_t {
  idx_t type;
  std::string modname;
  std::vector<idx_t> statetype;
  std::vector<std::vector<real_t>> stateparam;
  std::vector<idx_t> sticktype;
  std::vector<std::vector<real_t>> stickparam;
};

// Vertices
//
struct vertex_t {
  idx_t modidx;
  idx_t order;
  idx_t shape;
  std::vector<real_t> param;
  std::vector<real_t> coord;
};

// Edges
//
struct edge_t {
  idx_t source;
  std::vector<idx_t> target;
  idx_t modidx;
  real_t cutoff;
  std::vector<idx_t> conntype;
  std::vector<std::vector<real_t>> probparam;
  std::vector<std::vector<idx_t>> maskparam;
};

// Data files
//
struct datafile_t {
  std::string filename;
  // Sparse matrix (can also be used as a vector)
  std::vector<std::unordered_map<idx_t, real_t>> matrix;
};

// Size Distributions
//
struct dist_t {
  idx_t prtidx;
  idx_t nvtx;
  idx_t nedg;
  idx_t nstate;
  idx_t nstick;
  idx_t nevent;

  bool operator<(const dist_t& dist) const {
    return prtidx < dist.prtidx;
  }
};

// Events
//
struct event_t {
  tick_t diffuse;
  idx_t type;
  idx_t source;
  idx_t index;
  real_t data;
  
  bool operator<(const event_t& event) const {
    return diffuse < event.diffuse;
  }
};

// Vertex ordering
//
struct vtxorder_t {
  idx_t modidx;
  idx_t vtxidx;
  idx_t vtxidxloc; // local index of vertex
  bool operator < (const vtxorder_t& vtx) const {
    return (modidx < vtx.modidx);
  }
};

// Edge ordering
//
struct edgorder_t {
  idx_t edgidx;
  idx_t modidx;
  std::vector<real_t> state;
  std::vector<tick_t> stick;
  idx_t evtidx;
  bool operator < (const edgorder_t& edg) const {
    return (edgidx < edg.edgidx);
  }
};

/**************************************************************************
* Charm++ Mainchare
**************************************************************************/

// Main
//
class Main : public CBase_Main {
  public:
    /* Constructors and Destructors */
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);

    /* Control */
    void Control();
    void ReturnControl();
    void Halt(CkReductionMsg *msg);

    /* Persistence */
    int ParseConfig(std::string configfile);
    int ReadModel();
    int ReadGraph();
    int ReadMetis();
    int WriteDist();

    mModel* BuildModel();
    mGraph* BuildGraph();
    mMetis* BuildMetis();

  private:
    /* Chare Proxy */
    CProxy_GeNet genet;
    /* Models */
    std::vector<model_t> models;
    std::unordered_map<std::string, idx_t> modmap; // maps model name to object index
    std::vector<std::string> datafiles;
    /* Graph information */
    std::vector<vertex_t> vertices;
    std::vector<edge_t> edges;
    std::vector<std::string> graphtype;
    /* Persistence */
    std::vector<dist_t> netdist;
    std::vector<idx_t> vtxdist;
    std::vector<idx_t> edgdist;
    /* Bookkeeping */
    std::string mode;
    bool buildflag;
    bool partsflag;
    bool metisflag;
    bool orderflag;
    bool writeflag;
};


/**************************************************************************
* Charm++ Generate Network
**************************************************************************/

// Generate Network
//
class GeNet : public CBase_GeNet {
  public:
    /* Constructors and Destructors */
    GeNet(mModel *msg);
    GeNet(CkMigrateMessage* msg);
    ~GeNet();

    /* Build Network */
    void Build(mGraph *msg);
    void Connect(mConn *msg);

    void ConnRequest(idx_t reqidx);
    /* Reading Datafiles */
    int ReadDataCSV(datafile_t &datafile);

    /* Partition Network */
    void SetPartition();

    /* Reorder Network */
    void Read(mMetis *msg);
    void ScatterPart();
    void GatherPart(mPart *msg);
    void Order(mOrder *msg);
    void Reorder(mOrder *msg);
    mOrder* BuildOrder();

    /* Write Network */
    void Write(const CkCallback &cb);

    /* Connections */
    mConn* BuildPrevConn(idx_t reqidx);
    mConn* BuildCurrConn();
    mConn* BuildNextConn();
    idx_t MakeConnection(idx_t source, idx_t target, idx_t sourceidx, idx_t targetidx, real_t dist);
    std::vector<real_t> BuildEdgState(idx_t modidx, real_t dist, idx_t sourceidx, idx_t targetidx);
    std::vector<tick_t> BuildEdgStick(idx_t modidx, real_t dist, idx_t sourceidx, idx_t targetidx);

    /* Helper Functions */
    idx_t strtomodidx(const char* nptr, char** endptr) {
      const char *s;
      char c;
      int any;
      idx_t modidx = IDX_T_MAX;
      std::string name;
      // get rid of leading spaces
      s = nptr;
      do {
        c = *s++;
      } while (isspace(c));
      // read until space or end of line
      for (;; c = *s++) {
        // check for valid characters
        if (std::isalnum(c) || c == '_')
          name.append(&c,1);
        else
          break;
      }
      // If model exists, set it
      any = 0;
      if (modmap.find(name) != modmap.end()) {
        modidx = modmap[name];
        any = 1;
      }
      if (endptr != NULL)
        *endptr = (char *)(any ? s - 1 : nptr);
      return (modidx);
    }
    // Compute sigmoid
    real_t sigmoid(real_t x, real_t maxprob, real_t midpoint, real_t slope) {
      return maxprob * (1.0 - 1.0/(1.0 + std::exp( -slope * (x - midpoint) )));
    }
    // RNG State constant
    real_t rngconst(real_t *param) {
      return param[0];
    }
    // RNG State uniform
    real_t rngunif(real_t *param) {
      return param[0] + (param[1] - param[0])*((*unifdist)(rngine));
    }
    // RNG State uniform interval
    real_t rngunint(real_t *param) {
      return param[0] + param[2] * std::floor((((param[1] - param[0])/param[2])+1)*((*unifdist)(rngine)));
    }
    // RNG State normal
    real_t rngnorm(real_t *param) {
      return param[0] + (std::abs(param[1]))*((*normdist)(rngine));
    }
    // RNG State bounded normal
    real_t rngbnorm(real_t *param) {
      real_t state = (*normdist)(rngine);
      real_t bound = std::abs(param[2]);
      if (state > bound) { state = bound; }
      else if (state < -bound) { state = -bound; }
      return param[0] + (std::abs(param[1]))*state;
    }
    // RNG State lower bounded normal
    real_t rnglbnorm(real_t *param) {
      real_t state = (*normdist)(rngine);
      state = param[0] + (std::abs(param[1]))*state;
      if (state < param[2]) { state = param[2]; }
      return state;
    }
    // RNG State linear
    real_t rnglin(real_t *param, real_t dist) {
      return dist*param[0] + param[1];
    }
    // RNG State lower bounded linear
    real_t rnglblin(real_t *param, real_t dist) {
      real_t state = dist*param[0] + param[1];
      if (state < param[2]) {
        state = param[3];
      }
      return state;
    }
    // RNG State bounded linear
    real_t rngblin(real_t *param, real_t dist) {
      real_t state = dist*param[0] + param[1];
      if (state < param[2]) {
        state = param[2];
      }
      if (state > param[3]) {
        state = param[3];
      }
      return state;
    }
    // From datafile
    // (currently conforms to numpy savetxt format for csv)
    // Dimensions are stored: targetdim x sourcedim
    real_t rngfile(real_t *param, idx_t sourceidx, idx_t targetidx) {
      real_t state = 0.0;
      if (targetidx >= datafiles[(idx_t) (param[0])].matrix.size() ||
          datafiles[(idx_t) (param[0])].matrix[targetidx].find(sourceidx) ==
          datafiles[(idx_t) (param[0])].matrix[targetidx].end()) {
        // TODO: Throw an error if element doesn't exist
        CkPrintf("  error: datafile %s does not have element %" PRIidx ", %" PRIidx "\n",
                 datafiles[(idx_t) (param[0])].filename.c_str(), sourceidx, targetidx);
      } else {
        state = datafiles[(idx_t) (param[0])].matrix[targetidx][sourceidx];
      }
      return state;
    }

  private:
    /* Network Data */
    std::vector<idx_t> vtxdist;
    std::vector<real_t> xyz;
    std::vector<std::vector<idx_t>> adjcy;
    std::vector<std::vector<std::vector<real_t>>> state;
        // first level is the vertex, second level is the models, third is state data
    std::vector<std::vector<std::vector<tick_t>>> stick;
    std::vector<std::vector<event_t>> event;
    /* Models */
    std::vector<model_t> models;
    std::vector<std::string> modname;     // model names in order of object index
    std::unordered_map<std::string, idx_t> modmap; // maps model name to object index
    std::vector<std::string> rngtype;     // rng types in order of definitions
    std::vector<idx_t> vtxmodidx; // vertex model index into netmodel
    std::vector<idx_t> vtxordidx; // vertex index within model order
    std::vector<std::vector<idx_t>> edgmodidx; // edge model index into netmodel
    std::vector<datafile_t> datafiles;
    /* Connection information */
    std::vector<std::vector<std::vector<idx_t>>> adjcyconn;
        // first level is the data parts, second level are per vertex, third level is edges
    std::list<idx_t> adjcyreq; // parts that are requesting thier adjacency info
    std::vector<std::vector<std::vector<idx_t>>> edgmodidxconn; // edge model index into netmodel
        // first level is the data parts, second level are per vertex, third level is edges
    /* Graph information */
    std::vector<vertex_t> vertices; // vertex models and build information
    std::vector<edge_t> edges; // edge models and connection information
    /* Metis */
    std::vector<idx_t> vtxdistmetis; // distribution of vertices on data
    std::vector<idx_t> edgdistmetis; // distribution of edges on data
    std::vector<idx_t> partmetis; // which vertex goes to which part
    std::vector<std::vector<idx_t>> vtxidxpart; // vertex indices to go to a part
    std::vector<std::vector<idx_t>> vtxmodidxpart; // vertex models to go to a part
    std::vector<std::vector<idx_t>> xyzpart; // vertex models to go to a part
    std::vector<std::vector<std::vector<idx_t>>> adjcypart; // edge indices (by vtxidx) to go to a part
    std::vector<std::vector<std::vector<idx_t>>> edgmodidxpart; // edge models to go to a part
    std::vector<std::vector<std::vector<real_t>>> statepart;
        // first level is the part, second is the models, thrid is state data
    std::vector<std::vector<std::vector<tick_t>>> stickpart;
    std::vector<std::vector<std::vector<event_t>>> eventpart;
    /* Reordering */
    std::vector<std::vector<vtxorder_t>> vtxorder; // modidx and vtxidx for sorting
    std::vector<edgorder_t> edgorder; // edgidx and states for sorting
    std::vector<std::vector<real_t>> xyzorder; // coordinates by vertex
    std::vector<std::vector<std::vector<idx_t>>> adjcyorder; // adjacency by vertex
    std::vector<std::vector<std::vector<idx_t>>> adjcyreorder; // adjacency by vertex
    std::vector<std::vector<std::vector<idx_t>>> edgmodidxorder; // edge models by vertex
    std::vector<std::vector<std::vector<idx_t>>> edgmodidxreorder; // edge models by vertex
    std::vector<std::vector<std::vector<std::vector<real_t>>>> stateorder; // state by vertex
    std::vector<std::vector<std::vector<std::vector<real_t>>>> statereorder; // state by vertex
    std::vector<std::vector<std::vector<std::vector<tick_t>>>> stickorder; // stick by vertex
    std::vector<std::vector<std::vector<std::vector<tick_t>>>> stickreorder; // stick by vertex
    std::vector<std::vector<std::vector<event_t>>> eventorder; // event by vertex
    std::vector<std::vector<idx_t>> eventsourceorder; // source reordering
    std::vector<std::vector<idx_t>> eventindexorder; // index reordering
    std::list<mOrder *> ordering;
    /* Bookkeeping */
    int datidx;
    int cpdat;
    idx_t cpprt;
    idx_t nprt, xprt;
    idx_t norder; // total order
    idx_t norderdat; // order per data
    std::vector<idx_t> norderprt;  // order of vertices per network part
    std::vector<std::vector<idx_t>> nordervtx;  // order of vertex models 
    std::vector<std::vector<idx_t>> xordervtx;  // prefix of vertex models
    /* Random Number Generation */
    std::mt19937 rngine;
    std::uniform_real_distribution<real_t> *unifdist;
    std::normal_distribution<real_t> *normdist;
};


#endif //__STACS_GENET_H__
