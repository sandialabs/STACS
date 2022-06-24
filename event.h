/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_EVENT_H__
#define __STACS_EVENT_H__

#define EVENT_TOTAL         7
#define EVENT_SPIKE         0
#define EVENT_STIM          1
#define EVENT_SYNUP         2
#define EVENT_GROUP         3
#define EVENT_CURRENT       4
#define EVENT_COUNT         5
#define EVENT_CHGRATE       6

#define REMOTE_EDGES        0x00001
#define REMOTE_EDGE         0x00010
#define REMOTE_VERTEX       0x00100
#define LOCAL_EDGES         0x01000
#define LOCAL_VERTEX        0x10000

#endif //__STACS_EVENT_H__
