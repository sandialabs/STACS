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
#include "messages.h"
#include "ckmulticast.h"
#include "stacs.decl.h"

#ifdef STACS_WITH_YARP
#include <yarp/os/all.h>
#endif


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
* Data structures
**************************************************************************/

// Model
//
struct model_t {
  std::string modname;
  idx_t modtype;
  idx_t nstate;
  idx_t nstick;
  std::vector<real_t> param;
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

    /* Persistence */
    int ReadDist();
    int ReadModel();
    int WriteDist(bool check = false);

    /* Chare Messages */
    mDist* BuildDist();
    mModel* BuildModel();
#ifdef STACS_WITH_YARP
    mVtxDist* BuildVtxDist();
#endif

    /* Simulation */
    void InitSim();
    void StartSim();
    void CheckSim(CkReductionMsg *msg);
    void SaveSim();
    void FiniSim(CkReductionMsg *msg);
    void Halt();

  private:
    /* Persistence */
    std::vector<dist_t> netdist;
    std::vector<model_t> models;
    /* Chare Arrays */
    CProxy_NetData netdata;
    CProxy_Network network;
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
* Remote Procedure Call
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
    void Open(CProxy_Network cpnet) { }
    void Close() { }
    void Callback() { }
    /* Synchronization */
    void RPCSync() { }
    void RPCPause() { }
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
    void Open(CProxy_Network cpnet);
    void Close();

    /* Synchronization */
    void RPCSync();
    void RPCPause();

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
