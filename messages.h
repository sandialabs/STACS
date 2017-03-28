/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_MESSAGES_H__
#define __STACS_MESSAGES_H__

#define RPCPORT_DEFAULT "/stacs/rpc"

#define RPCCOMMAND_NONE     0
#define RPCCOMMAND_PAUSE    1
#define RPCCOMMAND_UNPAUSE  2
#define RPCCOMMAND_STOP     3
#define RPCCOMMAND_CHECK    4
#define RPCCOMMAND_STEP     5
#define RPCCOMMAND_STIM     6
#define RPCCOMMAND_PSTIM    7
#define RPCCOMMAND_OPEN     8
#define RPCCOMMAND_CLOSE    9

#define RPCSYNC_UNSYNCED    0
#define RPCSYNC_SYNCING     1
#define RPCSYNC_SYNCED      2

#define RPCSTIM_FILE        0
#define RPCSTIM_POINT       1
#define RPCSTIM_CIRCLE      2
#define RPCSTIM_SPHERE      3

#define MODFLAG_DEFAULT     0x00
#define MODFLAG_ACTIVE      0x01
#define MODFLAG_PNGMOD      0x10

#define REMOTE_EDGES        0x00001
#define REMOTE_EDGE         0x00010
#define REMOTE_VERTEX       0x00100
#define LOCAL_EDGES         0x01000
#define LOCAL_VERTEX        0x10000

#define EVENT_TOTAL         4
#define EVENT_RECORD        0
#define EVENT_SPIKE         1
#define EVENT_EDGUP         2
#define EVENT_STIM          3

#define RECORD_STATE        0
#define RECORD_STICK        1
#define RECORD_COORD        2

#endif //__STACS_MESSAGES_H__
