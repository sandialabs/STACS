/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#ifndef __STACS_EVENT_H__
#define __STACS_EVENT_H__

#define EVENT_TOTAL         8
#define EVENT_SPIKE         0  // Regular spike
#define EVENT_STIM          1  // Data payload
#define EVENT_CURRENT       2  // Spike with data payload
#define EVENT_CLAMP         3  // Set to data payload
#define EVENT_SYNUP         4  // Synapse updates
#define EVENT_COUNT         5  // Multiple spikes
#define EVENT_RATE          6  // Spike rate
#define EVENT_GROUP         7  // Spatio-temporal group activation

#define REMOTE_EDGES        0x00001
#define REMOTE_EDGE         0x00010
#define REMOTE_VERTEX       0x00100
#define LOCAL_EDGES         0x01000
#define LOCAL_VERTEX        0x10000

#endif //__STACS_EVENT_H__
