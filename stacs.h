/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_STACS_H__
#define __STACS_STACS_H__

#include <string>
#include "typedefs.h"
#include "timing.h"
#include "stream.h"
#include "ckmulticast.h"
#include "stacs.decl.h"

#ifdef STACS_WITH_YARP
#include <yarp/os/all.h>
#endif

#define STARTPAUSED_DEFAULT true // Start paused
#define PLASTICITY_DEFAULT  true // Plasticity on


/**************************************************************************
* Data structures
**************************************************************************/

// Model
//
struct model_t {
  std::string modname;
  idx_t modtype;
  flag_t modflag;
  idx_t nstate;
  idx_t nstick;
  std::vector<real_t> param;
  std::vector<std::string> port;
};


/**************************************************************************
* Charm++ Messages
**************************************************************************/

// Graph adjacency distribution (vertex only)
//
#define MSG_VtxDist 1
class mVtxDist : public CMessage_mVtxDist {
  public:
    idx_t *vtxdist; // number of vertices in partitions
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
    int ParseConfig(std::string configfile);
    int ReadDist();
    int ReadModel();

    /* Chare Messages */
    mDist* BuildDist();
    mModel* BuildModel();
#ifdef STACS_WITH_YARP
    mVtxDist* BuildVtxDist();
#endif

    /* Stacs */
    void Init();
    void Start();
    void Stop();
    void Halt();
    
    /* Simulation */
    void CheckNetwork(CkReductionMsg *msg);
    void SaveNetwork(CkReductionMsg *msg);
    int WriteDist(bool check = false);

    /* Polychronization */
    void ResetPNG();

  private:
    /* Persistence */
    std::vector<dist_t> netdist;
    std::vector<model_t> models;
    /* Simulation */
    bool startpaused;
    bool plasticity;
    /* Polychronization */
    std::vector<std::string> actives;
    std::vector<std::string> pngmods;
    /* Chare Arrays */
    CProxy_NetData netdata;
    CProxy_Network network;
    /* Timing */
    std::chrono::system_clock::time_point tstart, tstop;
#ifdef STACS_WITH_YARP
    /* YARP */
    CProxy_StreamRPC streamrpc;
    yarp::os::Network yarp;
#endif
    /* Bookkeeping */
    int ninit, cinit;
    int nhalt, chalt;
};


/**************************************************************************
* Stream (YARP Remote Procedure Call)
**************************************************************************/

#ifndef STACS_WITH_YARP
// Empty classes if no YARP
class StreamRPC : public CBase_StreamRPC {
  public:
    /* Constructors and Destructors */
    StreamRPC(mVtxDist *msg) { delete msg; }
    StreamRPC(CkMigrateMessage *msg) { delete msg; }
    ~StreamRPC() { }
    /* Stream */
    void Open(CProxy_Network cpnet, bool startpaused) { }
    void Close() { }
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
    RPCReader(CProxy_Network n, std::vector<idx_t>& v,
        const CkCallback &cbp, const CkCallback &cbs, const CkCallback &cbm) :
      network(n), vtxdist(v), cbpause(cbp), cbsync(cbs), cbmain(cbm) { }

    /* Reader */
    virtual bool read(yarp::os::ConnectionReader& connection);

    /* RPC Messages */
    mRPC* BuildRPCMsg(idx_t command, yarp::os::Bottle message);

    /* Coordination */
    void SetSyncFlag(idx_t sf) { syncflag = sf; }

  private:
    /* Network */
    CProxy_Network network;
    std::vector<idx_t> vtxdist;
    CkCallback cbpause;
    CkCallback cbsync;
    CkCallback cbmain;
    /* Coordination */
    idx_t syncflag;
};

// Port server (server)
//
class StreamRPC : public yarp::os:: RpcServer, public CBase_StreamRPC {
  public:
    /* Constructors and Destructors */
    StreamRPC(mVtxDist *msg);
    StreamRPC(CkMigrateMessage *msg);
    ~StreamRPC();

    /* Stream */
    void Open(CProxy_Network cpnet, bool startpaused);
    void Close();

    /* Synchronization */
    void Sync();
    void Pause();

  private:
    /* Network */
    CProxy_Network network;
    std::vector<idx_t> vtxdist;
    /* Reader */
    RPCReader *rpcreader;
    /* Coordination */
    CkCallback *cbmain;
};

#endif //STACS_WITH_YARP

#endif //__STACS_STACS_H__
