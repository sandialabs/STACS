/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#ifndef __STACS_NETWORK_H__
#define __STACS_NETWORK_H__
#include <algorithm>
#include <array>
#include <deque>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "typedefs.h"
#include "timing.h"
#include "build.h"
#include "event.h"
#include "record.h"
#include "ckmulticast.h"
#include "network.decl.h"

#ifdef STACS_WITH_YARP
#include <yarp/os/all.h>
#endif


/**************************************************************************
* Data structures
**************************************************************************/

// Network Size Distributions
//
struct dist_t;

struct model_t;
struct vertex_t;
struct edge_t;

// Data files
//
struct datafile_t {
  std::string filename;
  idx_t filetype;
  // Sparse matrix (can also be used as a vector)
  std::vector<std::unordered_map<idx_t, real_t>> matrix;
};

// Vertex reordering
//
struct vtxreord_t {
  idx_t modidx;
  idx_t vtxidx;
  idx_t vtxidxloc; // local index of vertex
  bool operator < (const vtxreord_t& vtx) const {
    return (modidx < vtx.modidx);
  }
};

// Edge reordering
//
struct edgreord_t {
  idx_t edgidx;
  idx_t modidx;
  std::vector<real_t> state;
  std::vector<tick_t> stick;
  idx_t evtidx;
  bool operator < (const edgreord_t& edg) const {
    return (edgidx < edg.edgidx);
  }
};

// Auxiliary indices
//
struct auxidx_t {
  idx_t index;
  std::vector<idx_t> stateidx;
  std::vector<idx_t> stickidx;
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
    return (diffuse < event.diffuse ||
        (diffuse == event.diffuse && source < event.source));
  }
};

// Records
//
struct record_t {
  idx_t recidx;
  tick_t trec;
  std::vector<real_t> state;
  std::vector<tick_t> stick;
  std::vector<idx_t> index;
};

// Record list
//
struct track_t {
  tick_t trec;
  tick_t tfreq;
  idx_t rectype;
  idx_t recmodidx;
  idx_t recsttidx;
  // indices of where tracked state is
  std::vector<idx_t> recvtxidx;
  std::vector<std::vector<idx_t>> recedgidx;
};

// Spike-timing events (for Groups)
//
struct stamp_t {
  tick_t diffuse; // timestamp of spike
  idx_t source; // index of vertex that spiked

  bool operator<(const stamp_t & stamp) const {
    return (diffuse < stamp.diffuse ||
        (diffuse == stamp.diffuse && source < stamp.source));
  }
  bool operator==(const stamp_t & stamp) const {
    return (diffuse == stamp.diffuse &&
            source == stamp.source);
  }
};

// Spike-timing routes (for Groups)
//
struct route_t {
  tick_t diffuse; // timestamp of spike
  idx_t source; // index of vertex that spiked
  idx_t origin; // index of vertex that contributed to spike
  tick_t departure; // timestamp of vertex spike departure
  tick_t arrival; // timestamp of vertex spike arrival
  
  bool operator<(const route_t & route) const {
    return (diffuse < route.diffuse ||
        (diffuse == route.diffuse && source < route.source));
  }
};

// Spike-timing history (for Groups)
//
struct trace_t {
  idx_t origin; // index of vertex that (possibly) contributed to spike
  tick_t departure; // timestamp of vertex spike departure
  tick_t arrival; // timestamp of vertex spike arrival

  bool operator<(const trace_t& trace) const {
    return arrival < trace.arrival;
  }
};


/**************************************************************************
* Charm++ Init Nodes
**************************************************************************/

// Tick reduction
//
void registerMinTick(void);
CkReductionMsg *minTick(int nMsg, CkReductionMsg **msgs);

// Idx reduction
//
void registerMaxIdx(void);
CkReductionMsg *maxIdx(int nMsg, CkReductionMsg **msgs);

// Network Distribution Reduction
//
void registerNetDist(void);
CkReductionMsg *netDist(int nMsg, CkReductionMsg **msgs);

// Polychronous Neuronal Group Reduction
//
void registerNetGroup(void);
CkReductionMsg *netGroup(int nMsg, CkReductionMsg **msgs);

// Events List Reduction
//
void registerNetEvent(void);
CkReductionMsg *netEvent(int nMsg, CkReductionMsg **msgs);


/**************************************************************************
* Charm++ Messages
**************************************************************************/

// Graph adjacency distribution
//
#define MSG_Dist 5
class mDist : public CMessage_mDist {
  public:
    idx_t *vtxdist; // number of vertices in partitions
    idx_t *edgdist; // number of edges in partitions
    idx_t *statedist; // number of states in partitions
    idx_t *stickdist; // number of ticks in partitions
    idx_t *eventdist; // number of events in partitions
};

// Network model information
//
#define MSG_Model 34
class mModel : public CMessage_mModel {
  public:
    idx_t *modtype;     // model index identifier
    idx_t *graphtype;   // type of model (vertex/edge)
    idx_t *xmodname;    // model name prefix
    char *modname;      // model name
    idx_t *nstate;      // number of states  per model
    idx_t *nstick;      // number of sticks per model
    idx_t *nparam;      // number of params per model
    idx_t *xstatename;    // state names prefix
    idx_t *xstickname;    // stick names prefix
    char *statename;    // state specification names
    char *stickname;    // stick specification names
    idx_t *xstateinit;    // state generation prefix
    idx_t *xstickinit;    // stick generation prefix
    idx_t *stateinit;     // state generation type
    idx_t *stickinit;     // stick representation type
    real_t *stateparam;   // state generation parameters
    real_t *stickparam;   // stick generation parameters
    idx_t *xparamname;  // parameter names prefix
    char *paramname;  // parameter specification names
    idx_t *xparam;      // network model prefix
    real_t *param;      // network model parameters
    idx_t *xport;       // network port prefix
    char *port;        // network port name
    idx_t *xdatafiles;  // prefix sum for filenames
    char *datafiles;    // filenames (concatenated)
    idx_t *datatypes;   // data filetypes
    idx_t *evtloglist;  // event types for logging
    idx_t *recmodidx;    // model to record from
    tick_t *rectfreq;    // record interval
    idx_t *xrecstate;    // state names prefix
    char *recstate;      // state names of model to record
    // TODO: remove polychron from base sim/model and create separate header/files for it
    bool *grpactive;   // polychronization (active)
    bool *grpmother;   // polychronization (mother)
    bool *grpanchor;   // polychronization (anchor)
    idx_t nmodel;      // number of models
    idx_t nstateparam;   // number of state generation parameters (for prefix)
    idx_t nstickparam;   // number of stick generation parameters (for prefix)
    idx_t ndatafiles;    // number of datafiles (for prefix)
    idx_t nevtlog;      // number of events to log
    idx_t nrecord;      // number of records to track
    bool plastic;      // toggle for plasticity
    bool episodic;     // toggle for episodic simulation
    bool loadbal;      // toggle for periodic load balancing
};

// Network graph information
//
#define MSG_Graph 20
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
    idx_t *edgdistype;    // distance function type
    idx_t *medgdistparam; // distance function sizes
    real_t *edgdistparam; // distance function parameters
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
    idx_t nedgdistparam;
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
    idx_t *vtxidx;      // vertex global idx
    idx_t *xadj;        // prefix for adjacency
    idx_t *adjcy;       // adjacent vertices
    int datidx;
    idx_t nvtx;
};

#define MSG_Reorder 3
class mReorder : public CMessage_mReorder {
  public:
    idx_t *vtxdist;
    idx_t *vtxidxold;
    idx_t *vtxidxnew;
    int datidx;
    idx_t nvtx;
    idx_t nprt;
};


// Network partition data
//
#define MSG_Part 15
class mPart : public CMessage_mPart {
  public:
    /* Data */
    idx_t *vtxdist;
    idx_t *vtxidx;
    idx_t *vtxmodidx;
    real_t *xyz;
    idx_t *xadj;
    idx_t *adjcy;
    idx_t *edgmodidx;
    real_t *state;
    tick_t *stick;
    /* Events */
    idx_t *xevent;
    tick_t *diffuse;
    idx_t *type;
    idx_t *source;
    idx_t *index;
    real_t *data;
    /* Sizes */
    idx_t nvtx;
    idx_t nedg;
    idx_t nstate;
    idx_t nstick;
    idx_t nevent;
    /* Bookkeeping */
    int prtidx;
};

// Network record data
//
#define MSG_Record 9
class mRecord : public CMessage_mRecord {
  public:
    tick_t *diffuse;
    idx_t *type;
    idx_t *source;
    idx_t *index;
    real_t *data;
    tick_t *drift;
    idx_t *xdata;
    idx_t *xdiffuse;
    idx_t *xindex;
    /* Bookkeeping */
    int prtidx;
    idx_t iter;
    idx_t nevtlog;
    idx_t nrecord;
};

// Go-ahead
//
class mGo : public CkMcastBaseMsg, public CMessage_mGo {
  public:
    /* Constructors */
    mGo(idx_t i) : iter(i) { }
    /* Bookkeeping */
    idx_t iter;
};

// Network event data
//
#define MSG_Event 5
class mEvent : public CkMcastBaseMsg, public CMessage_mEvent {
  public:
    /* Data */
    tick_t *diffuse;
    idx_t *type;
    idx_t *source;
    idx_t *index;
    real_t *data;
    /* Bookkeeping */
    idx_t iter;
    idx_t nevent;
};

// RPC Command to Network
//
#define MSG_RPC 1
class mRPC : public CMessage_mRPC {
  public:
    real_t *rpcdata; // data sent with rpc message
    idx_t nrpcdata; // size of data
    int command;  // command
};


/**************************************************************************
* Network Models (abstract factory pattern)
**************************************************************************/

// Network Model
//
class NetModel {
  public:
    virtual ~NetModel() { }
    /* Getters */
    idx_t getModType() const { return modtype; }
    std::vector<std::string> getParamList() const { return paramlist; }
    std::vector<std::string> getStateList() const { return statelist; }
    std::vector<std::string> getStickList() const { return sticklist; }
    idx_t getNParam() const { return paramlist.size(); }
    idx_t getNState() const { return statelist.size(); }
    idx_t getNStick() const { return sticklist.size(); }
    idx_t getNAux() const { return auxstate.size() + auxstick.size(); }
    idx_t getNPort() const { return portlist.size(); }
    std::vector<real_t> getParam() const { return param; }
    std::vector<std::string> getAuxState() const { return auxstate; }
    std::vector<std::string> getAuxStick() const { return auxstick; }
    idx_t getParamIdx(const std::string& name) const {
      idx_t index = std::find(paramlist.begin(), paramlist.end(), name) - paramlist.begin();
      return (index < paramlist.size() ? index : -1);
    }
    idx_t getStateIdx(const std::string& name) const {
      idx_t index = std::find(statelist.begin(), statelist.end(), name) - statelist.begin();
      return (index < statelist.size() ? index : -1);
    }
    idx_t getStickIdx(const std::string& name) const {
      idx_t index = std::find(sticklist.begin(), sticklist.end(), name) - sticklist.begin();
      return (index < sticklist.size() ? index : -1);
    }
    // TODO: Move polychronization stuff out of model
    bool getPlastic() const { return plastic; }
    bool getActive() const { return active; }
    bool getMother() const { return mother; }
    bool getAnchor() const { return anchor; }
    /* Setters */
    void setRandom(std::uniform_real_distribution<real_t> *u, std::mt19937 *r) {
      unifdist = u;
      rngine = r;
    }
    void setParam(real_t *p) {
      param = std::vector<real_t>(p, p+paramlist.size());
    }
    void setPort(char *p) {
      portname.resize(portlist.size());
      for (std::size_t i = 0; i < portlist.size(); ++i) {
        portname[i] = std::string(p);
        p += portname[i].size() + 1;
      }
    }
    void setPlastic(bool plas) { plastic = plas; }
    void setActive(bool grpactive) { active = grpactive; }
    void setMother(bool grpmother) { mother = grpmother; }
    void setAnchor(bool grpanchor) { anchor = grpanchor; }
    /* Computational Abstract Functions */
    virtual tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) = 0;
    virtual void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) = 0;
    /* Periodic Computation */
    virtual void Leap(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) { }
    virtual void getLeap(std::vector<event_t>& events) { }
    /* Protocol Functions */
    virtual void OpenPorts() { }
    virtual void ClosePorts() { }
    /* Control Flow */
    virtual void Renew(std::vector<real_t>& state, std::vector<tick_t>& stick) { }
    virtual void Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) { }
    /* Default Parameters */
    virtual idx_t getDefaultStateType(idx_t paramidx) { return RNGTYPE_CONST; }
    virtual std::vector<real_t> getDefaultStateParam(idx_t paramidx) { std::vector<real_t> zero (1,0.0); return zero; }
    virtual idx_t getDefaultStickType(idx_t paramidx) { return RNGTYPE_CONST; }
    virtual std::vector<real_t> getDefaultStickParam(idx_t paramidx) { std::vector<real_t> zero (1,0.0); return zero; }
    virtual real_t getDefaultParam(idx_t paramidx) { return 0.0; }
  protected:
    /* Random Number Generation */
    std::mt19937 *rngine;
    std::uniform_real_distribution<real_t> *unifdist;
    /* Model Information */
    idx_t modtype;
    std::vector<std::string> paramlist;
    std::vector<std::string> statelist;
    std::vector<std::string> sticklist;
    std::vector<std::string> auxstate;
    std::vector<std::string> auxstick;
    std::vector<std::string> portlist;
    /* Model Data */
    std::vector<real_t> param;
    /* Protocol */
    std::vector<std::string> portname;
    /* Control Flow */
    bool plastic;
    /* Polychronization */
    bool active;
    bool mother;
    bool anchor;
};

// Network model template
//
template <idx_t TYPE, typename IMPL>
class ModelTmpl : public NetModel {
  public:
    static NetModel* Create() { return new IMPL(); }
    static const idx_t MODELTYPE; // for registration
    static void Enable() { volatile idx_t x = MODELTYPE; }
  protected:
    ModelTmpl() { modtype = MODELTYPE; }
  private:
    /* Bring type into scope */
    enum { _MODELTYPE = TYPE };
};

// Network model factory singleton
//
class ModelFactory {
  public:
    typedef NetModel* (*t_pfFactory)();

    static ModelFactory *newModel() {
      static ModelFactory factory;
      return &factory;
    }

    /* Register concrete factory */
    idx_t Register(idx_t modtype, t_pfFactory model) {
      //CkPrintf("Registering constructor for model %" PRIidx "\n", modtype);
      modellist[modtype] = model;
      return modtype;
    }

    /* Create concrete object */
    NetModel *Create(idx_t modtype) {
      return modellist[modtype]();
    }

    /* Model List*/
    std::unordered_map<idx_t, t_pfFactory> modellist;

  private:
    /* Empty */
    ModelFactory() { }
    ~ModelFactory() { }
    /* Prevent copies */
    ModelFactory(ModelFactory const&) { }
    ModelFactory& operator=(ModelFactory const&);
};

// Network model factory registration
//
template <idx_t TYPE, typename IMPL>
const idx_t ModelTmpl<TYPE, IMPL>::MODELTYPE = ModelFactory::newModel()->Register(
    ModelTmpl<TYPE, IMPL >::_MODELTYPE, &ModelTmpl<TYPE, IMPL >::Create);


// "none" model
//
class NoneModel : public ModelTmpl < 0, NoneModel > {
  public:
    NoneModel() { paramlist.resize(0); statelist.resize(0); sticklist.resize(0); auxstate.resize(0); auxstick.resize(0); portlist.resize(0); }
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) { return tdiff; }
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) { }
};


/**************************************************************************
* Charm++ Network Arrays
**************************************************************************/

// Network data
//
class Netdata : public CBase_Netdata {
  public:
    /* Constructors and Destructors */
    Netdata(mModel *msg);
    Netdata(CkMigrateMessage *msg);
    ~Netdata();

    /* Building */
    void Build(mGraph *msg);
    void BuildParts();
    void SaveBuild();
    void SaveCloseBuild();
    std::vector<real_t> BuildEdgState(idx_t modidx, real_t dist, idx_t sourceidx, idx_t targetidx);
    std::vector<tick_t> BuildEdgStick(idx_t modidx, real_t dist, idx_t sourceidx, idx_t targetidx);
    void ReBuildEdgState(idx_t modidx, real_t dist, std::vector<real_t>& rngstate);
    void ReBuildEdgStick(idx_t modidx, real_t dist, std::vector<tick_t>& rngstick);

    /* Connecting */
    idx_t MakeConnection(idx_t edg, idx_t sourceidx, idx_t targetidx, real_t dist);
    void ConnectVtx(mConn *msg);
    void RequestConnVtx(int reqidx);
    mConn* BuildConnVtx(int reqidx);
    void ConnectHandover();
    void ConnectEdg(mConn *msg);
    void RequestConnEdg(int reqidx);
    mConn* BuildConnEdg(int reqidx);

    /* Repartitioning */
    void LoadPart(mDist *msg);
    void ReadPart();
    void ScatterPart();
    void GatherPart(mPart *msg);
    void Reorder(mReorder *msg);
    void ReordEdg(mReorder *msg);
    mReorder* BuildReorder();
    void LoadRepart(mPart *msg);
    void BuildRepart();
    
    /* Reading Datafiles */
    int ReadDataCSV(datafile_t &datafile);
    int ReadDataCSVSparse(datafile_t &datafile);
    
    /* Loading */
    void LoadData(mDist *msg);
    void LoadNetwork(int prtidx, const CkCallback &cbnet);
    void ReadNetwork();

    /* Saving */
    void SaveNetwork(mPart *msg);
    void SaveCloseNetwork(mPart *msg);
    void CloseNetwork();
    void WriteNetwork(int checkflag=1);
    
    /* Recording */
    void SaveRecord(mRecord *msg);
    void SaveFinalRecord(mRecord *msg);
    void WriteRecord();
    void SaveEstimate(CkReductionMsg *msg);
    void SaveFinalEstimate(CkReductionMsg *msg);
    void WriteEstimate();
    
    /* Closing */
    void FinalizeNetwork();
    
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
    // RNG State lower bounded lognorm
    real_t rnglblognorm(real_t *param) {
      real_t mu = param[3]*param[0];
      real_t var = std::sqrt(std::log(param[1]/(mu*mu) + 1));
      real_t state = std::log(mu*mu/std::sqrt(mu*mu + param[1]));
      state = std::exp(state + var*(*normdist)(rngine));
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
    // Distance Functions
    real_t distfunc(idx_t distype, real_t *xyz1, real_t *xyz2, real_t *param) {
      if (distype == DISTYPE_EUCLIDEAN) {
        // Normal Euclidean (just compute inline)
        return sqrt((xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0])+
                    (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1])+
                    (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]));
      }
      else if (distype == DISTYPE_PERIRECT) {
        return dist_perirect(xyz1, xyz2, param);
      }
      else if (distype == DISTYPE_SPHERE) {
        return dist_spheresurf(xyz1, xyz2, param);
      }
      else {
        // This shoudn't happen
        CkPrintf("  error: unknown distance error type\n");
        return 0.0;
      }
    }
    // Periodic rectangle
    real_t dist_perirect(real_t *xyz1, real_t *xyz2, real_t *param) {
      // TODO: assumes bottom right corner (x,y) starts at 0,0,0
      // width/height subject to periodic boundary conditions
      real_t dx = std::abs(xyz1[0] - xyz2[0]);
      real_t dy = std::abs(xyz1[1] - xyz2[1]);
      if (dx > param[0]/2.0) { dx -= param[0]/2.0; }
      if (dy > param[1]/2.0) { dy -= param[1]/2.0; }
      // Distance currently doesn't depend on z
      return sqrt((dx*dx)+(dy*dy));
    }
    // Surface of sphere
    real_t dist_spheresurf(real_t *xyz1, real_t *xyz2, real_t *param) {
      // TODO: assumes the sphere is centered at 0,0,0
      // dist = r*theta = r*acos(dot(a,b)/r^2)
      real_t theta = std::acos((xyz1[0]*xyz2[0] + xyz1[1]*xyz2[1] + xyz1[2]*xyz2[2])/(param[0]*param[0]));
      return param[0] * theta;
    }

  private:
    /* Network Data */
    std::vector<mPart*> parts;
    std::vector<mRecord*> records;
    std::vector<event_t> grplog;
    /* Distributions */
    std::vector<dist_t> netdist; // over local partitions
    std::vector<idx_t> vtxdist; // over all partitions
    std::vector<idx_t> edgdist;
    std::vector<idx_t> statedist;
    std::vector<idx_t> stickdist;
    std::vector<idx_t> eventdist;
    CkCallback maindist;
    /* Network Model */
    std::vector<NetModel*> model;      // collection of model objects
    std::vector<model_t> modelconf;    // model information from config and implementation defaults
    std::vector<std::string> modname;     // model names in order of of object index
    std::unordered_map<std::string, idx_t> modmap; // maps model name to object index
    std::vector<std::string> rngtype; // rng types in order of definitions
    std::vector<datafile_t> datafiles;
    /* Network Build */
    std::vector<vertex_t> vertices; // vertex models and build information
    std::vector<edge_t> edges; // edge models and connection information
    idx_t norder; // total order
    idx_t norderdat; // order for this netdata chare
    std::vector<idx_t> xorderdat; // prefix order for netdata
    std::vector<idx_t> norderprt; // order per network partition
    std::vector<std::vector<idx_t>> xpopvtxidxprt;  // prefix of populations over partitions
    std::vector<std::vector<idx_t>> xglbvtxidxprt;  // prefix of vertex over partitions (population global)
    /* Network State */
    std::vector<idx_t> vtxmodidx; // vertex model index into netmodel
    std::vector<idx_t> vtxordidx; // vertex index within model order
    std::vector<real_t> xyz;
    std::vector<std::vector<idx_t>> adjcy;
    std::vector<std::vector<idx_t>> edgmodidx; // edge model index into netmodel
    std::vector<std::vector<std::vector<real_t>>> state;
        // first level is the vertex, second level is the models, third is state data
    std::vector<std::vector<std::vector<tick_t>>> stick;
    std::vector<std::vector<event_t>> event;
    /* Connecting */
    std::unordered_map<idx_t,idx_t> connmodmap; // (source*edges.size + target) to edge
    std::vector<std::set<idx_t>> connsampleset; // target to edges that are sample-based for it
    std::vector<std::set<idx_t>> adjcyset; // set of directed afferent edges
    std::vector<idx_t> nadjcysample; // counts of number of adjcy for index-based connections
    std::list<int> connvtxreq; // parts that are requesting their vertex into
    /* Repartitioning */
    std::vector<idx_t> vtxmetis; // distribution of vertices over files
    std::vector<idx_t> edgmetis; // distribution of edges over files
    std::vector<std::vector<idx_t>> vtxidxreprt; // vertex indices to go to a part
    std::vector<std::vector<idx_t>> vtxmodidxreprt; // vertex models to go to a part
    std::vector<std::vector<real_t>> xyzreprt; // coordinates to go to a part
    std::vector<std::vector<std::vector<idx_t>>> adjcyreprt; // edge indices (by vtxidx) to go to a part
    std::vector<std::vector<std::vector<idx_t>>> edgmodidxreprt; // edge models to go to a part
    std::vector<std::vector<std::vector<real_t>>> statereprt;
        // first level is the part, second is the models, thrid is state data
    std::vector<std::vector<std::vector<tick_t>>> stickreprt;
    std::vector<std::vector<std::vector<event_t>>> eventreprt;
    std::vector<std::vector<vtxreord_t>> vtxprted; // modidx and vtxidx for sorting
    std::vector<std::vector<real_t>> xyzprted; // coordinates by vertex
    std::vector<std::vector<std::vector<idx_t>>> adjcyprted; // adjacency by vertex
    std::vector<std::vector<std::vector<idx_t>>> edgmodidxprted; // edge models by vertex
    std::vector<std::vector<std::vector<std::vector<real_t>>>> stateprted; // state by vertex
    std::vector<std::vector<std::vector<std::vector<tick_t>>>> stickprted; // stick by vertex
    std::vector<std::vector<std::vector<event_t>>> eventprted; // event by vertex
    /* Reordering */
    std::list<mReorder *> reordlist;
    std::vector<edgreord_t> edgreord; // edgidx and states for sorting
    std::vector<std::vector<std::vector<idx_t>>> adjcyreord; // adjacency by vertex
    std::vector<std::vector<std::vector<idx_t>>> edgmodidxreord; // edge models by vertex
    std::vector<std::vector<std::vector<std::vector<real_t>>>> statereord; // state by vertex
    std::vector<std::vector<std::vector<std::vector<tick_t>>>> stickreord; // stick by vertex
    std::vector<std::vector<idx_t>> eventsourcereord; // source reordering
    std::vector<std::vector<idx_t>> eventindexreord; // index reordering
    /* Random Number Generation */
    std::mt19937 rngine;
    std::uniform_real_distribution<real_t> *unifdist;
    std::normal_distribution<real_t> *normdist;
    /* Bookkeeping */
    int datidx;
    int cprt, rprt;
    int nprt, xprt;
    int cpdat, cphnd, cpprt;
#ifdef STACS_WITH_YARP
    /* YARP */
    yarp::os::Network yarp;
#endif
};

// Network partitions
//
class Network : public CBase_Network {
  public:
    /* Constructors and Destructors */
    Network(mModel *msg);
    Network(CkMigrateMessage *msg);
    ~Network();

    /* Loading */
    mPart* BuildPart();
    void LoadNetwork(mPart *msg);
    void ReloadNetwork(mPart *msg);
    
    /* Simulation */
    void InitSim(CProxy_Netdata cpdata);
    void ContSim();
    void CycleSim();

    void InitSimGPU(CProxy_Netdata cpdata);
    void CycleSimGPU();
    
    /* Estimation */
    void InitEst(CProxy_Netdata cpdata);
    void CycleEst();

    /* Communication */
    void CreateComm();
    void GoAhead(mGo *msg);
    mEvent* BuildEvent();
    void CommEvent(mEvent *msg);
    void CommStamp(mEvent *msg);
    
    /* Computation */
    void SortEventCalendar();
    void LeapEvent();
    void HandleEvent(event_t &event, const idx_t i);
    void EstimateGroup(const idx_t i);
    
    /* Saving */
    mPart* BuildRepart();
    void RebalNetwork();
    void SaveNetwork();
    void SaveCloseNetwork();
    void CloseNetwork();
    void ResetNetwork();
    
    /* Recording */
    mRecord* BuildRecord();
    void AddRecord();
    void SaveRecord();
    void SaveFinalRecord();
    void SaveEstimate();
    void SaveFinalEstimate();

    /* Polychronization */
    void InitGroup(CProxy_Netdata cpdata);
    void FindGroup();
    void ComputeGroup(idx_t nseeds, int grpart);
    void ComputeGroup();
    mEvent* BuildGroupSeed(std::vector<event_t>& seed);
    void SeedGroup(mEvent *msg);
    void CycleGroup();
    void EvalGroup(CkReductionMsg *msg);
    void ReadGroup(idx_t groupidx);
    void WriteGroup(idx_t groupidx);

#ifdef STACS_WITH_YARP
    /* RPC Control */
    void CommRPC(mRPC *msg);
#else
    void CommRPC(mRPC *msg) { delete msg; }
#endif
    
  private:
    /* Multicast */
    CProxySection_Network netcomm;
    std::unordered_map<int, int> srcpart; // source partitions
    std::vector<int> trgpart; // target partitions
    /* Network Adjacency */
    std::vector<idx_t> vtxdist; // vertex distribution per part
    std::vector<idx_t> vtxidx; // local vertex index to global index
    std::unordered_map<idx_t, idx_t> vtxmap; // global vertex index to local index
    std::vector<idx_t> vtxmodidx; // vertex model index
    std::vector<std::vector<idx_t>> adjcy; // local vertex adjacency to global index
    std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>> adjmap; // mapping from global index to vector of target vertices and their adjcency
    /* Network Models */
    std::vector<NetModel*> model; // collection of model objects
    std::vector<std::vector<idx_t>> edgmodidx; // edge model index
    /* Network Data */
    std::vector<real_t> xyz; // spatial coordinates per vertex
    std::vector<std::vector<std::vector<real_t>>> state; // first level is the vertex, second level is the models, third is the state data
    std::vector<std::vector<std::vector<tick_t>>> stick; // first level is the vertex, second level is the models, third is the stick data
    std::vector<std::vector<auxidx_t>> vtxaux; // auxiliary indices per vertex (for vertex/edge cross-modification of state)
    std::vector<std::vector<std::vector<auxidx_t>>> edgaux; // auxiliary indices  per edge model
    /* Network Events */
    std::vector<std::vector<std::vector<event_t>>> evtcal; // event queue per vertex per iteration (similar to calendar queue)
    std::vector<std::vector<event_t>> evtcol; // collection of overflow/spillover event queue
    std::vector<event_t> events; // event buffer for generated events
    std::vector<event_t> evtext; // external events (and extra spillover)
    std::vector<event_t> evtrpc; // events generated from RPC
    /* Periodic Events */
    tick_t tleap;
    std::vector<event_t> leapevt; // set of periodic events
    std::vector<bool> leaplist; // models with periodic events
    std::vector<std::vector<std::array<idx_t, 2>>> leapidx; // indices into models
    /* Repartitioning */
    std::vector<idx_t> reprtidx; // partition that the vertex should move to
    /* Recording */
    std::vector<event_t> evtlog; // event logging
    std::vector<bool> evtloglist; // types of events to log
    std::vector<record_t> record; // record keeping
    std::vector<track_t> recordlist; // what to record
    std::vector<std::set<idx_t>> recordmodset;
    /* Polychronization */
    std::vector<std::vector<std::vector<stamp_t>>> grpstamps; // Groups per vertex (as mother)
    std::vector<std::vector<tick_t>> grpdur; // Group duration per vertex (as mother)
    std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>> grpmap; // mapping from global index to vector of target Groups
    std::vector<std::vector<std::deque<stamp_t>>> grpwindow; // Sliding window of stamps per group vertex (as mother)
    std::vector<event_t> grplog; // logging group activation
    std::vector<std::vector<event_t>> grpseeds; // Potential groups of the vertex (seed events)
    std::vector<std::deque<trace_t>> grptraces; // sliding window of contributing spike-timing events per vertex
    std::vector<route_t> grpleg; // Generated spike-timing routes during computation per partition
    std::vector<route_t> grproute; // Candidate group (evaluated on reduction)
    std::vector<std::vector<route_t>> grproutes; // Collection of groups for a given vertex
    /* Random Number Generation */
    std::mt19937 rngine;
    std::uniform_real_distribution<real_t> *unifdist;
    /* Bookkeeping */
    CProxy_Netdata netdata;
    CkCallback cyclepart;
    int prtidx, datidx;
    /* Configuration */
    bool plastic;
    bool episodic;
    bool loadbal;
    /* Timing */
    tick_t tsim;
    tick_t teps;
    idx_t epsidx;
    /* Computation */
    tick_t tcomp;
    idx_t ncomp, ccomp;
    idx_t compidx;
    int compart;
    /* Coordination */
    idx_t iter;
    idx_t commiter;
    idx_t dispiter;
    idx_t saveiter;
    idx_t reciter;
    idx_t baliter;
    idx_t nadjpart;
    idx_t cadjpart[2], partiter;
#ifdef STACS_WITH_YARP
    /* RPC Control */
    CkCallback netrpc;
    idx_t synciter;
    bool syncing;
#endif
};

#endif //__STACS_NETWORK_H__
