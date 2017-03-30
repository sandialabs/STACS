/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_EVENT_H__
#define __STACS_EVENT_H__

#define EVENT_TOTAL         4
#define EVENT_RECORD        0
#define EVENT_SPIKE         1
#define EVENT_EDGUP         2
#define EVENT_STIM          3

#define REMOTE_EDGES        0x00001
#define REMOTE_EDGE         0x00010
#define REMOTE_VERTEX       0x00100
#define LOCAL_EDGES         0x01000
#define LOCAL_VERTEX        0x10000

#define RECORD_STATE        0
#define RECORD_STICK        1
#define RECORD_COORD        2

#endif //__STACS_EVENT_H__
