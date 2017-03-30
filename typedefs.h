/**
 * Copyright (C) 2017 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_TYPEDEFS_H__
#define __STACS_TYPEDEFS_H__
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

typedef float       real4;
typedef double      real8;
typedef long double real10;

typedef int64_t     idx_t;
typedef uint64_t    uidx_t;
typedef uint64_t    tick_t;
typedef uint16_t    flag_t;
typedef real8       real_t;
typedef int16_t     sample_t;

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

#endif //__STACS_TYPEDEFS_H__
