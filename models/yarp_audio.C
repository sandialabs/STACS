/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

#ifdef STACS_WITH_YARP
#include <yarp/sig/Sound.h>

/**************************************************************************
* Ports
**************************************************************************/

class YarpAudioPort : public yarp::os::BufferedPort<yarp::sig::Sound> {
  public:
    YarpAudioPort();
    //using yarp::os::BufferedPort<yarp::sig::Sound>::onRead;
    virtual void onRead (yarp::sig::Sound& s) {
      // Receive sound
      CkPrintf("Received %d samples\n", s.getSamples());
    }
};

YarpAudioPort::YarpAudioPort() {
    CkPrintf("Receiving samples\n");
}
#endif

/**************************************************************************
* Class declaration
**************************************************************************/
class YarpAudio : public NetModelTmpl < 101, YarpAudio > {
  public:
    /* Constructor */
    YarpAudio() {
      // parameters
      paramlist.resize(0);
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
      portlist[0] = "in";
    }

    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) { return tdiff; }
    void Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) { }

    /* Protocol */
    void OpenPorts();
    void ClosePorts();

  private:
#ifdef STACS_WITH_YARP
    YarpAudioPort* port;
#endif
};


/**************************************************************************
* Class methods
**************************************************************************/

// Open ports
//
void YarpAudio::OpenPorts() {
#ifdef STACS_WITH_YARP
  port = new YarpAudioPort();
  port->useCallback();
  port->open(portname[0].c_str());
#endif
}

// Close ports
//
void YarpAudio::ClosePorts() {
#ifdef STACS_WITH_YARP
  port->close();
#endif
}

