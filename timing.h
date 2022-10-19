/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_TIMING_H__
#define __STACS_TIMING_H__
#include <chrono>
#include <thread>

#define TICKS_PER_MS     1000000 // 1 million
#define TSTEP_DEFAULT    1.0 // 1ms
#define TEVENTQ_DEFAULT  20.0 // 20ms
#define TDISPLAY_DEFAULT 100.0 // 100ms
#define TRECORD_DEFAULT  500.0 // 500ms
#define TBALANCE_DEFAULT 5000.0 // 5s
#define TSAVE_DEFAULT    2000.0 // 2s
#define TMAX_DEFAULT     10000.0 // 10s
#define EPISODES_DEFAULT 20 // 20 episodes
#define TEPISODE_DEFAULT 500.0 // 500ms

#endif //__STACS_TIMING_H__
