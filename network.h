/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_NETWORK_H__
#define __STACS_NETWORK_H__
#include <algorithm>
#include <array>
#include <deque>
#include <random>
#include <set>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "typedefs.h"
#include "timing.h"
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
  tick_t drift;
  idx_t type;
  std::vector<real_t> data;
  std::vector<tick_t> diffuse;
  std::vector<idx_t> index;
};

// Record list
//
struct track_t {
  tick_t trec;
  tick_t tfreq;
  std::vector<idx_t> type;
  std::vector<idx_t> index;
  std::vector<idx_t> model;
  std::vector<idx_t> value;
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
#define MSG_Dist 8
class mDist : public CMessage_mDist {
  public:
    idx_t *vtxdist; // number of vertices in partitions
    idx_t *edgdist; // number of edges in partitions
    idx_t *statedist; // number of states in partitions
    idx_t *stickdist; // number of ticks in partitions
    idx_t *eventdist; // number of events in partitions
    idx_t *modtype;
    idx_t *xmodname;
    char *modname;
    idx_t nmodel;
};

// Network model information
//
#define MSG_Model 10
class mModel : public CMessage_mModel {
  public:
    idx_t *modtype;
    idx_t *nstate;
    idx_t *nstick;
    idx_t *xparam;
    real_t *param;
    idx_t *xport;
    char *port;
    bool *grpactive;
    bool *grpmother;
    bool *grpanchor;
    idx_t nmodel;
    bool plastic;
    bool episodic;
};

// Network partition data
//
#define MSG_Part 14
class mPart : public CMessage_mPart {
  public:
    /* Data */
    idx_t *vtxdist;
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
    int partidx;
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
    int partidx;
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
    idx_t getNParam() const { return paramlist.size(); }
    idx_t getNState() const { return statelist.size(); }
    idx_t getNStick() const { return sticklist.size(); }
    idx_t getNAux() const { return auxstate.size() + auxstick.size(); }
    idx_t getNPort() const { return portlist.size(); }
    std::vector<real_t> getParam() const { return param; }
    std::vector<std::string> getAuxState() const { return auxstate; }
    std::vector<std::string> getAuxStick() const { return auxstick; }
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
    Netdata(mDist *msg);
    Netdata(CkMigrateMessage *msg);
    ~Netdata();

    /* Loading */
    void LoadNetwork(int partidx, const CkCallback &cbnet);
    void ReadNetwork();

    /* Saving */
    void SaveNetwork(mPart *msg);
    void SaveCloseNetwork(mPart *msg);
    void CloseNetwork();
    void WriteNetwork();
    
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

  private:
    /* Network Data */
    std::vector<mPart*> parts;
    std::vector<mRecord*> records;
    std::vector<event_t> grplog;
    /* Distributions */
    std::vector<dist_t> netdist;
    std::vector<idx_t> vtxdist;
    std::vector<idx_t> edgdist;
    std::vector<idx_t> statedist;
    std::vector<idx_t> stickdist;
    std::vector<idx_t> eventdist;
    CkCallback maindist;
    /* Network Model */
    std::vector<NetModel*> model;      // collection of model objects (empty)
    std::vector<std::string> modname;     // model names in order of of object index
    std::unordered_map<std::string, idx_t> modmap; // maps model name to object index
    /* Bookkeeping */
    int fileidx;
    int cpart, rpart;
    int npart, xpart;
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
    
    /* Simulation */
    void InitSim(CProxy_Netdata cpdata);
    void CycleSim();
    
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
    /* Recording */
    std::vector<event_t> evtlog; // event logging
    std::vector<bool> evtloglist; // types of events to log
    std::vector<record_t> record; // record keeping
    std::vector<track_t> recordlist; // what to record
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
    int partidx, fileidx;
    /* Configuration */
    bool plastic;
    bool episodic;
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
