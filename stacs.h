/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_STACS_H__
#define __STACS_STACS_H__

#include <string>
#include <sstream>
#include <unordered_map>
#include "typedefs.h"
#include "timing.h"
#include "ckmulticast.h"
#include "stacs.decl.h"

#ifdef STACS_WITH_YARP
#include <yarp/os/all.h>
#endif

#define CONFIG_DEFAULT "configs/config.yml"

#define RUNMODE_SIMULATE  "simulate"
#define RUNMODE_SIMGPU    "simgpu"
#define RUNMODE_BUILD     "build"
#define RUNMODE_FINDGROUP "findgroup"
#define RUNMODE_ESTIMATE  "estimate"
#define RUNMODE_EMPTY     ""
#define RUNMODE_DEFAULT   "simulate"

#define PLASTIC_DEFAULT  true // Plasticity on
#define EPISODIC_DEFAULT false // Episodes off
#define RPCPORT_DEFAULT "/stacs/rpc"
#define RPCPAUSE_DEFAULT true // Start paused

#define FILELOAD_DEFAULT ""
#define FILESAVE_DEFAULT ".out"
#define RECORDIR_DEFAULT "record"
#define GROUPDIR_DEFAULT "group"

#define GRPMINLEN_DEFAULT 7
#define GRPMAXDUR_DEFAULT 150.0

/**************************************************************************
* Data structures
**************************************************************************/

// Network Size Distributions
//
struct dist_t {
  int partidx;
  idx_t nvtx;
  idx_t nedg;
  idx_t nstate;
  idx_t nstick;
  idx_t nevent;

  bool operator<(const dist_t& dist) const {
    return partidx < dist.partidx;
  }
};

// Network Models
//
struct model_t {
  std::string modname;
  idx_t modtype;
  idx_t graphtype;
  idx_t nstate;
  idx_t nstick;
  idx_t nparam;
  std::vector<std::string> statename;
  std::vector<std::string> stickname;
  std::vector<idx_t> statetype;
  std::vector<idx_t> sticktype;
  std::vector<std::vector<real_t>> stateparam;
  std::vector<std::vector<real_t>> stickparam;
  std::vector<std::string> paramname;
  std::vector<real_t> param;
  std::vector<std::string> port;
  bool grpactive;
  bool grpmother;
  bool grpanchor;
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

/**************************************************************************
* Charm++ Messages
**************************************************************************/

// Graph adjacency distribution (vertex only)
//
#define MSG_Vtxs 2
class mVtxs : public CMessage_mVtxs {
  public:
    idx_t *vtxdist; // number of vertices in partitions
    char *rpcport; // yarp rpc port name
    int xrpcport;
};


/**************************************************************************
* Charm++ Mainchare
**************************************************************************/

// Main
//
class Main : public CBase_Main {
  public:
    /* Constructors and Destructors */
    Main(CkArgMsg *msg);
    Main(CkMigrateMessage *msg);

    /* Configuration */
    int ReadConfig(std::string configfile);
    int ReadDist();
    int ReadModel();
    int ReadGraph();

    /* Chare Messages */
    mDist* BuildDist();
    mModel* BuildModel();
    mGraph* BuildGraph();
#ifdef STACS_WITH_YARP
    mVtxs* BuildVtxs();
#endif

    /* Stacs */
    void Control();
    void Init();
    void Start();
    void Stop();
    void Halt();
    
    /* Network Distribution */
    void SaveDist(CkReductionMsg *msg);
    void SaveFinalDist(CkReductionMsg *msg);
    int WriteDist();

  private:
    /* Chare Arrays */
    CProxy_Netdata netdata;
    CProxy_Network network;
    CkCallback netcycle;
    /* Network */
    std::vector<dist_t> netdist;
    /* Model Information */
    // TODO: rename to modelconfig and combine edges/vertices into a graphconfig variable
    std::vector<model_t> models;
    std::unordered_map<std::string, std::size_t> modmap; // maps model name to object index
    std::vector<std::string> datafiles;
    /* Graph information */
    std::vector<vertex_t> vertices;
    std::vector<edge_t> edges;
    /* Configuration */
    std::string runmode;
    bool plastic;
    bool episodic;
    real_t teventq;
    real_t tdisplay;
    real_t trecord;
    real_t tsave;
    /* Groups */
    std::vector<std::string> grpactives;
    std::vector<std::string> grpmothers;
    std::vector<std::string> grpanchors;
    realidx_t grpvtxminreal, grpvtxmaxreal;
    /* Timing */
    std::chrono::system_clock::time_point wcstart, wcstop;
    /* Bookkeeping */
    int ninit, cinit;
    int nhalt, chalt;
    bool buildflag;
    bool writeflag;
#ifdef STACS_WITH_YARP
    /* YARP */
    CProxy_Stream stream;
    yarp::os::Network yarpnet;
    std::string rpcport;
    bool rpcpause;
#endif
};


/**************************************************************************
* Stream (YARP Remote Procedure Call)
**************************************************************************/

#ifndef STACS_WITH_YARP
// Empty classes if no YARP
class Stream : public CBase_Stream {
  public:
    /* Constructors and Destructors */
    Stream(mVtxs *msg) { delete msg; }
    Stream(CkMigrateMessage *msg) { delete msg; }
    ~Stream() { }
    /* Stream */
    void OpenRPC(CProxy_Network cpnet, const CkCallback &cbcycle, bool paused) { }
    void CloseRPC() { }
    /* Synchronization */
    void Sync(idx_t synciter) { };
    void Pause() { }
};
#else

// Port reader (client)
//
class RPCReader : public yarp::os::PortReader {
  public:
    /* Constructors and Destructors */
    RPCReader(std::vector<idx_t>& vtxs, CProxy_Network cpnet,
        const CkCallback &cbcycle, const CkCallback &cbstop,
        const CkCallback &cbpause, const CkCallback &cbsync) :
      vtxdist(vtxs), network(cpnet), netcycle(cbcycle),
      netstop(cbstop), netpause(cbpause), netsync(cbsync) { }

    /* Reader */
    virtual bool read(yarp::os::ConnectionReader& connection);

    /* RPC Messages */
    mRPC* BuildRPC(idx_t command, yarp::os::Bottle message);

    /* Coordination */
    void SetSyncFlag(idx_t sf) { syncflag = sf; }

  private:
    /* Network */
    std::vector<idx_t> vtxdist;
    CProxy_Network network;
    CkCallback netcycle;
    CkCallback netstop;
    CkCallback netpause;
    CkCallback netsync;
    /* Coordination */
    int syncflag;
};

// Port server (server)
//
class Stream : public yarp::os:: RpcServer, public CBase_Stream {
  public:
    /* Constructors and Destructors */
    Stream(mVtxs *msg);
    Stream(CkMigrateMessage *msg);
    ~Stream();

    /* Stream */
    void OpenRPC(CProxy_Network cpnet, const CkCallback &cbcycle, bool paused);
    void CloseRPC();

    /* Synchronization */
    void Sync(idx_t synciter);
    mRPC* BuildRPCSync(idx_t synciter);
    void Pause();

  private:
    /* Network */
    std::vector<idx_t> vtxdist;
    CProxy_Network network;
    CkCallback netcycle;
    CkCallback netstop;
    CkCallback netpause;
    /* Remote Procedure Call */
    std::string rpcport;
    RPCReader *rpcreader;
};

#endif //STACS_WITH_YARP

#endif //__STACS_STACS_H__
