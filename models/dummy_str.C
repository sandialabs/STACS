/**
 * Copyright (C) 2016 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Class declaration
**************************************************************************/
class DummyStr : public NetModelTmpl < 3, DummyStr > {
  public:
    /* Constructor */
    DummyStr() {
      // parameters
      paramlist.resize(1);
      paramlist[0] = "dp0";
      // states
      statelist.resize(0);
      // sticks
      sticklist.resize(0);
      // auxiliary states
      auxstate.resize(0);
      // auxiliary sticks
      auxstick.resize(0);
      // ports
      portlist.resize(1);
      portlist[0] = "dport";
#ifdef STACS_WITH_YARP
      port.resize(1);
#endif
    }

    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) { return tdiff; }
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) { }

    /* Protocol */
    void OpenPorts();
    void ClosePorts();

  private:
#ifdef STACS_WITH_YARP
    std::vector<yarp::os::Port> port;
#endif
};


/**************************************************************************
* Class methods
**************************************************************************/

// Open ports
//
void DummyStr::OpenPorts() {
#ifdef STACS_WITH_YARP
  for (std::size_t p = 0; p < portlist.size(); ++p) {
    port[p].open(portname[p].c_str());
  }
#endif
}

// Close ports
//
void DummyStr::ClosePorts() {
#ifdef STACS_WITH_YARP
  for (std::size_t p = 0; p < portlist.size(); ++p) {
    port[p].close();
  }
#endif
}
