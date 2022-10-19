/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

#ifdef STACS_WITH_YARP
#include <yarp/sig/ImageFile.h>
#include <yarp/sig/Image.h>

/**************************************************************************
* Ports
**************************************************************************/

class YimgUspsEpisPort : public yarp::os::BufferedPort<yarp::sig::ImageOf<yarp::sig::PixelMono> > {
  public:
    std::vector<real_t> imgs;
    YimgUspsEpisPort(int np);
    //using yarp::os::BufferedPort<yarp::sig::Image>::onRead;
    virtual void onRead (yarp::sig::ImageOf<yarp::sig::PixelMono>& img) {
      // Receive image and store
      // TODO: Check that incoming image size is equal to number of pixels
      printf("Image received\n");
      imgs.clear();
      for (std::size_t i = 0; i < img.height(); ++i) {
        for (std::size_t j = 0; j < img.width(); ++j) {
          //printf("%3d ", img.getRawImage()[i*img.height() + j]);
          imgs.push_back(img.getRawImage()[i*img.height() + j]);
        }
        //printf("\n");
      }
    }
  private:
    yarp::os::Property conf;
};

YimgUspsEpisPort::YimgUspsEpisPort(int np) {
  // Get an image write device
  CkPrintf("Streaming to stacs\n");
}
#endif

/**************************************************************************
* Class declaration
**************************************************************************/
class YimgUspsEpis : public ModelTmpl < 202, YimgUspsEpis > {
  public:
    /* Constructor */
    YimgUspsEpis() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "npxl";
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
    
    /* Control Flow */
    void Renew(std::vector<real_t>& state, std::vector<tick_t>& stick);

  private:
#ifdef STACS_WITH_YARP
    YimgUspsEpisPort* port;
    yarp::os::BufferedPort<yarp::os::Bottle>* preq;
#endif
    tick_t tsample;
};


/**************************************************************************
* Class methods
**************************************************************************/

// Renew model
//
void YimgUspsEpis::Renew(std::vector<real_t>& state, std::vector<tick_t>& stick) {
#ifdef STACS_WITH_YARP
  if (tsample == 0 && !(port->imgs.size())) {
    // Request new image
    yarp::os::Bottle& b = preq->prepare();
    b.clear();
    b.addInt((int) port->imgs.size());
    preq->write();
  }
#endif
  tsample = 0;
}

// Simulation step
//
tick_t YimgUspsEpis::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
#ifdef STACS_WITH_YARP
  // Grab data from port if available
  if (tsample == 0 && port->imgs.size()) {
    tsample = TICKS_PER_MS;
    // generate events
    event_t event;
    event.type = EVENT_STIM;
    event.source = REMOTE_EDGE;
    for (idx_t i = 0; i < (idx_t) param[0]; ++i) {
      event.index = i;
      // Multiplier goes from 0 to 255
      int delay = (int)((256.0 - port->imgs[i]) * 24.0 / 256.0);
      if (delay < 20) {
        event.diffuse = tdrift + delay * TICKS_PER_MS;
        event.data = param[1];
        events.push_back(event);
        event.diffuse = tdrift + (delay+1) * TICKS_PER_MS;
        event.data = -param[1];
        events.push_back(event);
      }
    }
    // Request new phone
    yarp::os::Bottle& b = preq->prepare();
    b.clear();
    b.addInt((int) port->imgs.size());
    preq->write();
  }
#endif
  return tdiff;
}

// Open ports
//
void YimgUspsEpis::OpenPorts() {
#ifdef STACS_WITH_YARP
  std::string recv = portname[0] + "/recv";
  std::string request = portname[0] + "/request";
  port = new YimgUspsEpisPort((int) param[0]);
  port->setStrict(); // try not to drop input
  port->useCallback();
  port->open(recv.c_str());
  CkPrintf("Receiving USPS digits on %s\n", recv.c_str());
  preq = new yarp::os::BufferedPort<yarp::os::Bottle>();
  preq->open(request.c_str());
  CkPrintf("Requesting USPS digits from %s\n", request.c_str());
#endif
}

// Close ports
//
void YimgUspsEpis::ClosePorts() {
#ifdef STACS_WITH_YARP
  port->close();
  yarp::os::Bottle& b = preq->prepare();
  b.clear();
  b.addInt(-1);
  preq->write();
  preq->close();
#endif
}

