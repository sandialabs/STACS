/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_TIMING_H__
#define __STACS_TIMING_H__
#include <chrono>
#include <thread>

#define TMAX_DEFAULT    10000.0 // 10s
#define TDISP_DEFAULT       1.0 // 1ms
#define TSTEP_DEFAULT       1.0 // 1ms
#define TCHECK_DEFAULT     50.0 // 50ms
#define TRECORD_DEFAULT    25.0 // 25ms
#define TICKS_PER_MS    1000000 // 1 million
#define EVTCAL_DEFAULT       20 // Days per year
#define SIMPAUSE_DEFAULT  false // Start unpaused

#define RUNMODE_DEFAULT       0
#define RUNMODE_SIM           0
#define RUNMODE_PNG           1

#endif //__STACS_TIMING_H__
