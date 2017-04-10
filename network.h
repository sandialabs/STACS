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
    return diffuse < event.diffuse;
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

// Spike-timing events (for PNGs)
//
struct stamp_t {
  tick_t diffuse; // timestamp of spike
  idx_t source; // index of vertex that spiked

  bool operator<(const stamp_t & stamp) const {
    return diffuse < stamp.diffuse;
  }
  bool operator==(const stamp_t & stamp) const {
    return (diffuse == stamp.diffuse &&
            source == stamp.source);
  }
};

// Spike-timing routes (for PNGs)
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

// Spike-timing history (for PNGs)
//
struct trail_t {
  idx_t origin; // index of vertex that (possibly) contributed to spike
  tick_t departure; // timestamp of vertex spike departure
  tick_t arrival; // timestamp of vertex spike arrival

  bool operator<(const trail_t& trail) const {
    return arrival < trail.arrival;
  }
};


/**************************************************************************
* Charm++ Init Nodes
**************************************************************************/

// Tick reduction
//
void registerMinTick(void);
CkReductionMsg *minTick(int nMsg, CkReductionMsg **msgs);

// Network Distribution Reduction
//
void registerNetDist(void);
CkReductionMsg *netDist(int nMsg, CkReductionMsg **msgs);

// Polychronous Neuronal Group Reduction
//
void registerNetPNG(void);
CkReductionMsg *netPNG(int nMsg, CkReductionMsg **msgs);


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
    bool *pngactive;
    bool *pngmother;
    bool *pnganchor;
    idx_t nmodel;
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
    idx_t prtidx;
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
    idx_t prtidx;
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
    idx_t command;  // command
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
    bool getPNGActive() const { return pngactive; }
    bool getPNGMother() const { return pngmother; }
    bool getPNGAnchor() const { return pnganchor; }
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
    void setPNGActive(bool pactive) {
      pngactive = pactive;
    }
    void setPNGMother(bool pmother) {
      pngmother = pmother;
    }
    void setPNGAnchor(bool panchor) {
      pnganchor = panchor;
    }
    /* Protocol Functions */
    virtual void OpenPorts() { }
    virtual void ClosePorts() { }
    /* Abstract Functions */
    virtual void Reset(std::vector<real_t>& state, std::vector<tick_t>& stick) { }
    virtual tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) = 0;
    virtual void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) = 0;
    virtual void Hop(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) { }
    virtual void Leap(std::vector<event_t>& events) { }
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
    /* Polychronization */
    bool pngactive;
    bool pngmother;
    bool pnganchor;
};

// Network model template
//
template <idx_t TYPE, typename IMPL>
class NetModelTmpl : public NetModel {
  public:
    static NetModel* Create() { return new IMPL(); }
    static const idx_t MODELTYPE; // for registration
    static void Enable() { volatile idx_t x = MODELTYPE; }
  protected:
    NetModelTmpl() { modtype = MODELTYPE; }
  private:
    /* Bring type into scope */
    enum { _MODELTYPE = TYPE };
};

// Network model factory singleton
//
class NetModelFactory {
  public:
    typedef NetModel* (*t_pfFactory)();

    static NetModelFactory *getNetModel() {
      static NetModelFactory factory;
      return &factory;
    }

    /* Register concrete factory */
    idx_t Register(idx_t modtype, t_pfFactory netmodel) {
      //CkPrintf("Registering constructor for vtx id %" PRIidx "\n", modtype);
      modellist[modtype] = netmodel;
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
    NetModelFactory() { }
    ~NetModelFactory() { }
    /* Prevent copies */
    NetModelFactory(NetModelFactory const&) { }
    NetModelFactory& operator=(NetModelFactory const&);
};

// Network model factory registration
//
template <idx_t TYPE, typename IMPL>
const idx_t NetModelTmpl<TYPE, IMPL>::MODELTYPE = NetModelFactory::getNetModel()->Register(
    NetModelTmpl<TYPE, IMPL >::_MODELTYPE, &NetModelTmpl<TYPE, IMPL >::Create);


// "none" model
//
class NoneModel : public NetModelTmpl < 0, NoneModel > {
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
    void LoadNetwork(idx_t prtidx, const CkCallback &cbnet);
    void ReadNetwork();

    /* Saving */
    void SaveNetwork(mPart *msg);
    void SaveFinalNetwork(mPart *msg);
    void FinalizeNetwork();
    void WriteNetwork();

    /* Recording */
    void SaveRecord(mRecord *msg);
    void SaveFinalRecord(mRecord *msg);
    void WriteRecord();
    
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
    /* Distributions */
    std::vector<dist_t> netdist;
    std::vector<idx_t> vtxdist;
    std::vector<idx_t> edgdist;
    std::vector<idx_t> statedist;
    std::vector<idx_t> stickdist;
    std::vector<idx_t> eventdist;
    /* Network Model */
    std::vector<NetModel*> netmodel;      // collection of model objects (empty)
    std::vector<std::string> modname;     // model names in order of of object index
    std::unordered_map<std::string, idx_t> modmap; // maps model name to object index
    /* Bookkeeping */
    idx_t datidx;
    idx_t cprt, rprt;
    idx_t nprt, xprt;
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
    void InitSimPlastic(CProxy_Netdata cpdat);
    void InitSimStatic(CProxy_Netdata cpdat);
    void CycleSimPlastic();
    void CycleSimStatic();

    /* Communication */
    void CreateGroup();
    void GoAhead(mGo *msg);
    mEvent* BuildEvent();
    void CommEvent(mEvent *msg);
    void MarkEvent();
    
    /* Saving */
    void SaveNetwork();
    void SaveFinalNetwork();
    void FinalizeNetwork();
    void ResetNetwork();
    
    /* Recording */
    mRecord* BuildRecord();
    void StoreRecord();
    void SaveRecord();
    void SaveFinalRecord();

    /* Polychronization */
    void InitPNG(CProxy_Netdata cpdat);
    void FindPNG();
    void ComputePNG(idx_t nseeds, idx_t pngidx);
    void ComputePNG();
    mEvent* BuildPNGSeed(std::vector<event_t>& pngseed);
    void SeedPNG(mEvent *msg);
    void CyclePNG();
    void EvalPNG(CkReductionMsg *msg);
    void WritePNG(idx_t pngidx);

#ifdef STACS_WITH_YARP
    /* RPC Control */
    void CommRPC(mRPC *msg);
#else
    void CommRPC(mRPC *msg) { delete msg; }
#endif
    
  private:
    /* Multicast */
    CProxySection_Network netgroup;
    std::unordered_map<idx_t, idx_t> srcprt; // source partitions
    std::vector<idx_t> trgprt; // target partitions
    /* Network Adjacency */
    std::vector<idx_t> vtxdist; // vertex distribution per part
    std::vector<idx_t> vtxidx; // local vertex index to global index
    std::unordered_map<idx_t, idx_t> vtxmap; // global vertex index to local index
    std::vector<idx_t> vtxmodidx; // vertex model index into netmodel
    std::vector<std::vector<idx_t>> adjcy; // local vertex adjacency to global index
    std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>> adjmap; // mapping from global index to vector of target vertices and their adjcency
    /* Network Model */
    std::vector<NetModel*> netmodel; // collection of model objects
    std::vector<std::vector<idx_t>> edgmodidx; // edge model index into netmodel
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
    std::vector<event_t> evtleap; // set of periodic events
    std::vector<bool> leaplist; // models with periodic events
    std::vector<std::vector<std::array<idx_t, 2>>> leapidx; // indices into models
    /* Recording */
    std::vector<event_t> evtlog; // event logging
    std::vector<bool> evtloglist; // types of events to log
    std::vector<record_t> record; // record keeping
    std::vector<track_t> recordlist; // what to record
    /* Polychronization */
    std::vector<std::vector<std::vector<stamp_t>>> pngs; // PNGs per vertex (as mother)
    std::vector<std::vector<event_t>> pngseeds; // Potential PNGs of the vertex (seed events)
    std::vector<std::deque<trail_t>> pngtrail; // sliding window of contributing spike-timing events per vertex
    std::vector<route_t> pnglog; // Generated spike-timing routes during computation per partition
    std::vector<route_t> pngmap; // Candidate PNG (evaluated on reduction)
    std::vector<std::vector<route_t>> pngmaps; // Collection of PNG routes for a given vertex
    /* Random Number Generation */
    std::mt19937 rngine;
    std::uniform_real_distribution<real_t> *unifdist;
    /* Bookkeeping */
    CProxy_Netdata netdata;
    CkCallback cbcycleprt;
    idx_t prtidx, datidx;
    /* Timing */
    tick_t tsim;
    tick_t tdisp;
    tick_t tleap;
    /* Coordination */
    idx_t iter;
    idx_t nadjprt;
    idx_t cadjprt[2], prtiter;
#ifdef STACS_WITH_YARP
    /* RPC Control */
    CkCallback cbrpc;
    idx_t synciter;
#endif
    /* Checkpointing */
    idx_t checkiter;
    idx_t reciter;
    /* Computation */
    idx_t compidx;
    idx_t compendx;
    idx_t evalidx;
    idx_t ncomp, ccomp;
    tick_t tcomp;
};

#endif //__STACS_NETWORK_H__
