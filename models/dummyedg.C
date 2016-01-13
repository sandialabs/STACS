/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DummyEdg : public NetModelTmpl < 2, DummyEdg > {
  public:
    /* Constructor */
    DummyEdg() { nparam = 2; nstate = 2; nstick = 1; }
    
    /* Simulation */
    void Step(tick_t tdrift);
};

/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
void DummyEdg::Step(tick_t tdrift) {
  CkPrintf("Stepping Edg\n");
}

