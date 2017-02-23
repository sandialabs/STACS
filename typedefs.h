/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_TYPEDEFS_H__
#define __STACS_TYPEDEFS_H__
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <vector>
//#include <math.h>

typedef float real4;
typedef double real8;
typedef long double real10;

typedef int64_t  idx_t;
typedef uint64_t uidx_t;
typedef uint64_t tick_t;
typedef uint16_t flag_t;
typedef real8 real_t;
typedef int16_t sample_t;

// Reading from file
#define strtoidx strtoll
#define strtouidx strtoull
#define strtotick strtoull
#define strtoreal strtod

// Printing to string
#define PRIidx PRId64
#define PRIuidx PRIu64
#define PRItick PRIu64
#define PRItickhex PRIx64
#define PRIreal ".2g"
#define PRIrealms ".1f"
#define PRIrealsec ".3f"
#define PRIrealfull " .7e"

// Some limits
#define IDX_T_MAX   0x7FFFFFFFFFFFFFFF
#define UIDX_T_MAX  0xFFFFFFFFFFFFFFFF
#define TICK_T_MAX  0xFFFFFFFFFFFFFFFF

/**************************************************************************
* Data Structures
**************************************************************************/

// Size Distributions
//
struct dist_t {
  idx_t prtidx;
  idx_t nvtx;
  idx_t nedg;
  idx_t nstate;
  idx_t nstick;
  idx_t nevent;

  bool operator<(const dist_t& dist) const {
    return prtidx < dist.prtidx;
  }
};

// Events
//
struct event_t {
  tick_t diffuse;
  idx_t type;
  idx_t source;
  idx_t index;
  real_t data;
  
  bool operator<(const event_t& event) const {
    return diffuse < event.diffuse;
  }
};

// Records
//
struct record_t {
  tick_t drift;
  std::vector<real_t> data;
};

// Record list
//
struct recentry_t {
  tick_t trec;
  tick_t tfreq;
  std::vector<idx_t> type;
  std::vector<idx_t> index;
  std::vector<idx_t> model;
  std::vector<idx_t> value;
};

// Auxiliary state
//
struct aux_t {
  idx_t index;
  std::vector<idx_t> stateidx;
  std::vector<idx_t> stickidx;
};

#endif //__STACS_TYPEDEFS_H__
