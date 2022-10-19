/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#ifndef __STACS_STREAM_H__
#define __STACS_STREAM_H__

#define RPCCOMMAND_NONE     0
#define RPCCOMMAND_PAUSE    1
#define RPCCOMMAND_PAUSED   2
#define RPCCOMMAND_UNPAUSE  3
#define RPCCOMMAND_STOP     4
#define RPCCOMMAND_CHECK    5
#define RPCCOMMAND_STEP     6
#define RPCCOMMAND_STIM     7
#define RPCCOMMAND_PSTIM    8
#define RPCCOMMAND_OPEN     9
#define RPCCOMMAND_CLOSE    10

#define RPCSYNC_UNSYNCED    0
#define RPCSYNC_SYNCING     1
#define RPCSYNC_SYNCED      2

#define RPCSTIM_FILE        0
#define RPCSTIM_POINT       1
#define RPCSTIM_CIRCLE      2
#define RPCSTIM_SPHERE      3

#endif //__STACS_STREAM_H__
