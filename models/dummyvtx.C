/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include <unistd.h>
#include "network.h"


/**************************************************************************
* Class declaration
**************************************************************************/
class DummyVtx : public NetModelTmpl < 1, DummyVtx > {
  public:
    /* Constructor */
    DummyVtx() { nparam = 1; nstate = 2; nstick = 0; }

    /* Simulation */
    void Step(tick_t tdrift);
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
void DummyVtx::Step(tick_t tdrift) {
  usleep(tdrift/1000); // tdrift ms
}
