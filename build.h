/**
 * Copyright (C) 2022 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#ifndef __STACS_BUILD_H__
#define __STACS_BUILD_H__

#define GRAPHTYPE_NTYPE 3

#define GRAPHTYPE_STR   0
#define GRAPHTYPE_VTX   1
#define GRAPHTYPE_EDG   2

#define REPTYPE_REAL    0
#define REPTYPE_TICK    1

// TODO: reorder these numbers sometime
#define RNGTYPE_NRNG    11
#define RNGTYPE_CONST   0
#define RNGTYPE_UNIF    1
#define RNGTYPE_UNINT   2
#define RNGTYPE_NORM    3
#define RNGTYPE_BNORM   4
#define RNGTYPE_LBNORM  5
#define RNGTYPE_LIN     6
#define RNGTYPE_LBLIN   7
#define RNGTYPE_UBLIN   8
#define RNGTYPE_BLIN    9
#define RNGTYPE_FILE    10

#define RNGPARAM_CONST  1
#define RNGPARAM_UNIF   2
#define RNGPARAM_UNINT  3
#define RNGPARAM_NORM   2
#define RNGPARAM_BNORM  3
#define RNGPARAM_LBNORM 3
#define RNGPARAM_LIN    2
#define RNGPARAM_LBLIN  3
#define RNGPARAM_UBLIN  3
#define RNGPARAM_BLIN   4
#define RNGPARAM_FILE   1

#define VTXSHAPE_POINT  0
#define VTXSHAPE_CIRCLE 1
#define VTXSHAPE_SPHERE 2
#define VTXSHAPE_RECT   3

#define VTXPARAM_POINT  0
#define VTXPARAM_CIRCLE 1
#define VTXPARAM_SPHERE 1
#define VTXPARAM_RECT   2

#define CONNTYPE_UNIF   0
#define PROBPARAM_UNIF  1
#define MASKPARAM_UNIF  0

#define CONNTYPE_SIG    1
#define PROBPARAM_SIG   3
#define MASKPARAM_SIG   0

#define CONNTYPE_IDX    2
#define PROBPARAM_IDX   0
#define MASKPARAM_IDX   4

#define CONNTYPE_FILE   3
#define PROBPARAM_FILE  1
#define MASKPARAM_FILE  2

#define CONNTYPE_SMPL   4
#define PROBPARAM_SMPL  1
#define MASKPARAM_SMPL  2

#endif //__STACS_BUILD_H__
