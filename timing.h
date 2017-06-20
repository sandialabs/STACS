/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_TIMING_H__
#define __STACS_TIMING_H__
#include <chrono>
#include <thread>

#define TICKS_PER_MS    1000000 // 1 million
#define TMAX_DEFAULT    10000.0 // 10s
#define TSTEP_DEFAULT       1.0 // 1ms
#define TQUEUE_DEFAULT     20.0 // 20ms
#define TCHECK_DEFAULT     50.0 // 50ms
#define TRECORD_DEFAULT    25.0 // 25ms
#define TDISPLAY_DEFAULT    1.0 // 1ms
#define TTRIAL_DEFAULT    500.0 // 500ms

#endif //__STACS_TIMING_H__
