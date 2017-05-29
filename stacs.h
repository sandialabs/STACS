/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_STACS_H__
#define __STACS_STACS_H__

#include <string>
#include <sstream>
#include "typedefs.h"
#include "timing.h"
#include "ckmulticast.h"
#include "stacs.decl.h"

#ifdef STACS_WITH_YARP
#include <yarp/os/all.h>
#endif

#define RUNMODE_SIM     "sim"
#define RUNMODE_PNG     "png"
#define RUNMODE_EST     "est"
#define RUNMODE_DEFAULT "sim"

#define RPCPORTNAME_DEFAULT "/stacs/rpc"
#define STARTPAUSED_DEFAULT true // Start paused
#define PLASTICITY_DEFAULT  true // Plasticity on

/**************************************************************************
* Data structures
**************************************************************************/

// Network Size Distributions
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

// Network Models
//
struct model_t {
  std::string modname;
  idx_t modtype;
  idx_t nstate;
  idx_t nstick;
  std::vector<real_t> param;
  std::vector<std::string> port;
  bool pngactive;
  bool pngmother;
  bool pnganchor;
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
    idx_t xrpcport;
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

    /* Chare Messages */
    mDist* BuildDist();
    mModel* BuildModel();
#ifdef STACS_WITH_YARP
    mVtxs* BuildVtxs();
#endif

    /* Stacs */
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
    CkCallback cbcycle;
    /* Configuration */
    std::vector<dist_t> netdist;
    std::vector<model_t> models;
    std::string runmode;
    bool plasticity;
    /* Polychronization */
    std::vector<std::string> pngactives;
    std::vector<std::string> pngmothers;
    std::vector<std::string> pnganchors;
    /* Timing */
    std::chrono::system_clock::time_point tstart, tstop;
    /* Bookkeeping */
    int ninit, cinit;
    int nhalt, chalt;
#ifdef STACS_WITH_YARP
    /* YARP */
    CProxy_Stream stream;
    yarp::os::Network yarp;
    std::string rpcport;
    bool startpaused;
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
    void OpenRPC(CProxy_Network cpnet, const CkCallback &cbcyc, bool paused) { }
    void CloseRPC() { }
    /* Synchronization */
    void Sync() { }
    void Pause() { }
};
#else

// Port reader (client)
//
class RPCReader : public yarp::os::PortReader {
  public:
    /* Constructors and Destructors */
    RPCReader(CProxy_Network cpn, std::vector<idx_t>& vtxs,
        const CkCallback &cbp, const CkCallback &cbs,
        const CkCallback &cbm, const CkCallback &cbc) :
      network(cpn), vtxdist(vtxs), cbpause(cbp), cbsync(cbs),
      cbmain(cbm), cbcycle(cbc) { }

    /* Reader */
    virtual bool read(yarp::os::ConnectionReader& connection);

    /* RPC Messages */
    mRPC* BuildRPC(idx_t command, yarp::os::Bottle message);

    /* Coordination */
    void SetSyncFlag(idx_t sf) { syncflag = sf; }

  private:
    /* Network */
    CProxy_Network network;
    std::vector<idx_t> vtxdist;
    CkCallback cbpause;
    CkCallback cbsync;
    CkCallback cbmain;
    CkCallback cbcycle;
    /* Coordination */
    idx_t syncflag;
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
    void OpenRPC(CProxy_Network cpnet, const CkCallback &cbcyc, bool paused);
    void CloseRPC();

    /* Synchronization */
    void Sync();
    void Pause();

  private:
    /* Network */
    CProxy_Network network;
    std::vector<idx_t> vtxdist;
    /* Callbacks */
    CkCallback cbmain;
    CkCallback cbcycle;
    /* Remote Procedure Call */
    std::string rpcport;
    RPCReader *rpcreader;
};

#endif //STACS_WITH_YARP

#endif //__STACS_STACS_H__
