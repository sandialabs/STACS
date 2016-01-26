/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_NETWORK_H__
#define __STACS_NETWORK_H__
#include <algorithm>
#include <array>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "typedefs.h"
#include "timing.h"
#include "messages.h"
#include "ckmulticast.h"
#include "network.decl.h"


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


/**************************************************************************
* Charm++ Messages
**************************************************************************/

// Go-ahead
//
class mGo : public CkMcastBaseMsg, public CMessage_mGo {
  public:
    /* Constructors */
    mGo(idx_t i) : iter(i) {}
    /* Bookkeeping */
    idx_t iter;
};

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
#define MSG_Model 5
class mModel : public CMessage_mModel {
  public:
    idx_t *modtype;
    idx_t *nstate;
    idx_t *nstick;
    idx_t *xparam;
    real_t *param;
    idx_t nmodel;
};

// Network partition data
//
#define MSG_Part 13
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
    idx_t *target;
    idx_t *type;
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

// Network event data
//
#define MSG_Event 4
class mEvent : public CkMcastBaseMsg, public CMessage_mEvent {
  public:
    /* Data */
    tick_t *diffuse;
    idx_t *index;
    idx_t *type;
    real_t *data;
    /* Bookkeeping */
    idx_t iter;
    idx_t nevent;
};

// Network record data
//
#define MSG_Record 6
class mRecord : public CMessage_mRecord {
  public:
    tick_t *diffuse;
    idx_t *index;
    idx_t *type;
    real_t *data;
    tick_t *drift;
    idx_t *xdata;
    /* Bookkeeping */
    idx_t prtidx;
    idx_t iter;
    idx_t nrecevt;
    idx_t nrecord;
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
    idx_t getNParam() const { return nparam; }
    idx_t getNState() const { return nstate; }
    idx_t getNStick() const { return nstick; }
    std::vector<real_t> getParam() const { return param; }
    /* Setters */
    void setParam(real_t *p) {
      param = std::vector<real_t>(p, p+nparam);
    }
    /* Abstract Functions */
    virtual void Step(tick_t tdrift, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) = 0;
  protected:
    /* Bookkeeping */
    idx_t modtype;
    idx_t nparam;
    idx_t nstate;
    idx_t nstick;
    /* Model Data */
    std::vector<real_t> param;
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
    NetModelFactory() { };
    ~NetModelFactory() { };
    /* Prevent copies */
    NetModelFactory(NetModelFactory const&) { };
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
    NoneModel() { nparam = 0; nstate = 0; nstick = 0; }
    void Step(tick_t tdrift, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) { }
};


/**************************************************************************
* Charm++ Network Arrays
**************************************************************************/

// Network data
//
class NetData : public CBase_NetData {
  public:
    /* Constructors and Destructors */
    NetData(mDist *msg);
    NetData(CkMigrateMessage *msg);
    ~NetData();

    /* Persistence */
    void LoadNetwork(idx_t prtidx, const CkCallback &cb);
    void ReadCSR();
    void CheckNetwork(mPart *msg);
    void SaveNetwork(mPart *msg);
    void WriteCSR(bool check = false);

    /* Recording */
    void CheckRecord(mRecord *msg);
    void SaveRecord(mRecord *msg);
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
};

// Network partitions
//
class Network : public CBase_Network {
  public:
    /* Constructors and Destructors */
    Network(mModel *msg);
    Network(CkMigrateMessage *msg);
    ~Network();

    /* Multicast */
    void CreateGroup();

    /* Persistence */
    void LoadNetwork(CProxy_NetData cpdat);
    void LoadNetwork(mPart *msg);
    void SaveNetwork();

    /* Recording */
    void StoreRecord();
    void SaveRecord();

    /* Computation */
    void LoadNetModel();

    /* Communication */
    mEvent* BuildEvent();
    void CommEvent(mEvent *msg);
    void RedisEvent();

    /* Simulation */
    void GoAhead(mGo *msg);
    void Cycle();

#ifdef STACS_WITH_YARP
    /* RPC Control */
    void RPCMsg(mRPC *msg);
#else
    void RPCMsg(mRPC *msg) { delete msg; }
#endif
    
  private:
    /* Multicast */
    CProxySection_Network netgroup;
    std::unordered_map<idx_t, idx_t> srcprt;
    std::vector<idx_t> trgprt;
    /* Network Adjacency */
    std::vector<idx_t> vtxdist;             // vertex distribution per part
    std::vector<idx_t> vtxidx;              // local vertex index to global index
    std::vector<idx_t> vtxmodidx; // vertex model index into netmodel
    std::vector<std::vector<idx_t>> adjcy;  // local vertex adjacency to global index
    std::unordered_map<idx_t, std::vector<std::array<idx_t, 2>>> adjmap;
        // mapping from global index to vector of target vertices and their adjcency
    /* Network Model */
    std::vector<NetModel*> netmodel;               // collection of model objects
    std::vector<std::vector<idx_t>> edgmodidx; // edge model index into netmodel
    /* Network Data */
    std::vector<real_t> xyz;
    std::vector<std::vector<std::vector<real_t>>> state;
        // first level is the vertex, second level is the models, third is state data
    std::vector<std::vector<std::vector<tick_t>>> stick;
    /* Network Events */
    std::vector<std::vector<std::vector<event_t>>> event;
        // event queue per vertex per iteration (use as calendar queue)
    std::vector<std::vector<event_t>> evtaux; // spillover event queue
    std::vector<event_t> evtreaux;  // spillover spillover event queue
    std::vector<event_t> evtlog; // event buffer for generated events
    /* Recording */
    std::vector<record_t> record; // record keeping
    std::vector<recentry_t> recordlist; // what to record
    std::vector<event_t> recevt; // record keeping for events
    idx_t recevtlist; // types of events to record
    /* Timing */
    tick_t tsim;
    idx_t iter;
    /* Bookkeeping */
    CProxy_NetData netdata;
    idx_t prtidx, datidx;
    /* Coordination */
    idx_t nadjprt;
    idx_t cadjprt[2], prtiter;
    /* Checkpointing */
    bool cpflag;
    idx_t checkiter;
    bool recflag;
    idx_t reciter;
#ifdef STACS_WITH_YARP
    /* RPC Control */
    CkCallback cbrpc;
    idx_t synciter;
#endif
};

#endif //__STACS_NETWORK_H__
