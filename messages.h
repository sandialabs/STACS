/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_MESSAGES_H__
#define __STACS_MESSAGES_H__

#define RPCPORT_DEFAULT "/stacs/rpc"

#define RPCCOMMAND_NONE    0
#define RPCCOMMAND_PAUSE   1
#define RPCCOMMAND_UNPAUSE 2
#define RPCCOMMAND_STOP    3
#define RPCCOMMAND_CHECK   4
#define RPCCOMMAND_STEP    5
#define RPCCOMMAND_STIM    6
#define RPCCOMMAND_PSTIM   7
#define RPCCOMMAND_OPEN    8
#define RPCCOMMAND_CLOSE   9

#define RPCSYNC_UNSYNCED   0
#define RPCSYNC_SYNCING    1
#define RPCSYNC_SYNCED     2

#define EVENT_EXTERNAL     1
#define EVENT_LOCALEDG     2
#define EVENT_LOCALVTX     4

#define EVTYPE_TOTAL       2
#define EVTYPE_RECORD      0
#define EVTYPE_SPIKE       1

#define RECORD_STATE       0
#define RECORD_STICK       1
#define RECORD_COORD       2

#endif //__STACS_MESSAGES_H__
