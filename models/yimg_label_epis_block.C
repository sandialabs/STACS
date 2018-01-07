/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

#ifdef STACS_WITH_YARP
#include <yarp/sig/ImageFile.h>
#include <yarp/sig/Image.h>
#include <mutex>
#include <condition_variable>

/**************************************************************************
* Ports
**************************************************************************/

class YimgLabelEpisBlockPort : public yarp::os::BufferedPort<yarp::sig::ImageOf<yarp::sig::PixelInt> > {
  public:
    std::mutex mtx;
    std::condition_variable cvar;
    bool ready;
    std::vector<real_t> labels;
    YimgLabelEpisBlockPort(int np);
    //using yarp::os::BufferedPort<yarp::sig::Image>::onRead;
    virtual void onRead (yarp::sig::ImageOf<yarp::sig::PixelInt>& label) {
      // Receive image and store
      // TODO: Check that incoming image size is equal to number of pixels
      printf("Label received\n");
      labels.clear();
      // label.height should equal number of pixels (it's a vector)
      for (std::size_t i = 0; i < label.height(); ++i) {
        labels.push_back(label.getRawImage()[i]);
      }
      std::unique_lock<std::mutex> lck(mtx);
      ready = true;
      cvar.notify_all();
    }
  private:
    yarp::os::Property conf;
};

YimgLabelEpisBlockPort::YimgLabelEpisBlockPort(int np) {
  // Get an image write device
  CkPrintf("Streaming labels to stacs\n");
  ready = false;
}
#endif

/**************************************************************************
* Class declaration
**************************************************************************/
class YimgLabelEpisBlock : public ModelTmpl < 109, YimgLabelEpisBlock > {
  public:
    /* Constructor */
    YimgLabelEpisBlock() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "points";
      paramlist[1] = "ampl";
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

      // Local variables
      tsample = 0;
    }

    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
    void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) { }

    /* Protocol */
    void OpenPorts();
    void ClosePorts();
    
    void Rerun(std::vector<real_t>& state, std::vector<tick_t>& stick);

  private:
#ifdef STACS_WITH_YARP
    YimgLabelEpisBlockPort* port;
    yarp::os::BufferedPort<yarp::os::Bottle>* preq;
#endif
    tick_t tsample;
};


/**************************************************************************
* Class methods
**************************************************************************/

// Rerun model
//
void YimgLabelEpisBlock::Rerun(std::vector<real_t>& state, std::vector<tick_t>& stick) {
#ifdef STACS_WITH_YARP
  if (tsample == 0 && !(port->labels.size())) {
    // Request new image
    yarp::os::Bottle& b = preq->prepare();
    b.clear();
    b.addInt((int) port->labels.size());
    preq->write();
  }
#endif
  tsample = 0;
}

// Simulation step
//
tick_t YimgLabelEpisBlock::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
#ifdef STACS_WITH_YARP
  // Grab data from port if available
  if (tsample == 0) {
    std::unique_lock<std::mutex> lck(port->mtx);
    while (!(port->ready)) {
      port->cvar.wait(lck);
    }
    port->ready = false;
    tsample = TICKS_PER_MS;
    // generate events
    event_t event;
    event.type = EVENT_STIM;
    event.source = REMOTE_EDGE;
    for (idx_t i = 0; i < (idx_t) param[0]; ++i) {
      event.index = i;
      // label is raw millisecond delays
      int delay = port->labels[i];
      if (delay >= 0) {
        event.diffuse = tdrift + delay * TICKS_PER_MS;
        event.data = param[1];
        events.push_back(event);
        event.diffuse = tdrift + (delay+1) * TICKS_PER_MS;
        event.data = -param[1];
        events.push_back(event);
      }
    }
    // Request new label
    yarp::os::Bottle& b = preq->prepare();
    b.clear();
    b.addInt((int) port->labels.size());
    preq->write();
  }
#endif
  return tdiff;
}

// Open ports
//
void YimgLabelEpisBlock::OpenPorts() {
#ifdef STACS_WITH_YARP
  std::string recv = portname[0] + "/recv";
  std::string request = portname[0] + "/request";
  port = new YimgLabelEpisBlockPort((int) param[0]);
  port->setStrict(); // try not to drop input
  port->useCallback();
  port->open(recv.c_str());
  CkPrintf("Receiving labels on %s\n", recv.c_str());
  preq = new yarp::os::BufferedPort<yarp::os::Bottle>();
  preq->open(request.c_str());
  CkPrintf("Requesting labels from %s\n", request.c_str());
#endif
}

// Close ports
//
void YimgLabelEpisBlock::ClosePorts() {
#ifdef STACS_WITH_YARP
  port->close();
  yarp::os::Bottle& b = preq->prepare();
  b.clear();
  b.addInt(-1);
  preq->write();
  preq->close();
#endif
}

