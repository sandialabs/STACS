/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

#ifdef STACS_WITH_YARP
#include <yarp/sig/Sound.h>
#include <fftw3.h>

/**************************************************************************
* Ports
**************************************************************************/

class YaudWinMFCCEpisPort : public yarp::os::BufferedPort<yarp::sig::Sound> {
  private:
    // Sound samples
    real_t samplefreq;
    std::vector<real_t> samplebuffer;
    // Power spectral density
    int ninput;
    int noverlap;
    int noutput;
    std::vector<real_t> window;
    real_t windownormsquare;
    std::vector<real_t> psd;
    // Mel freqencies
    int nmel;
    std::vector<real_t> mel;
    std::vector<std::vector<real_t>> filters;
    //FILE *pMEL;
    //char melfile[100];
  public:
    std::deque<std::vector<real_t>> mels;
    YaudWinMFCCEpisPort(int nm);
    //using yarp::os::BufferedPort<yarp::sig::Sound>::onRead;
    virtual void onRead (yarp::sig::Sound& ssig) {
      // Receive sound
      CkPrintf("Received %d samples\n", ssig.getSamples());
      std::size_t nsample = ssig.getSamples();
      samplebuffer.resize(nsample);
      for (std::size_t s = 0; s < nsample; ++s) {
        samplebuffer[s] = ((real_t)((sample_t)ssig.get(s)))/32768.0;
      }
      //CkPrintf("Sample samples: %f, %f, %f, %f, %f\n", samplebuffer[50], samplebuffer[100], samplebuffer[150], samplebuffer[200], samplebuffer[250]);
      // Split samples into windows
      int nwindow = (((int)samplebuffer.size()) - noverlap)/(ninput - noverlap);
      for (int w = 0; w < nwindow; ++w) {
        // Find power spectral density of window
        double* inputbuffer = static_cast<double*>(fftw_malloc(ninput * sizeof(double)));
        fftw_complex* outputbuffer = static_cast<fftw_complex*>(fftw_malloc(noutput * sizeof(fftw_complex)));
        for (int n = 0; n < ninput; ++n) {
          inputbuffer[n] = samplebuffer[w * (ninput - noverlap) + n] * window[n];
        }
        fftw_plan plan = fftw_plan_dft_r2c_1d(ninput, inputbuffer, outputbuffer, FFTW_ESTIMATE);
        fftw_execute(plan);
        for (int n = 0; n < noutput; ++n) {
          psd[n] = (outputbuffer[n][0] * outputbuffer[n][0] + outputbuffer[n][1] * outputbuffer[n][1]) / (samplefreq * windownormsquare);
        }
        // Filtering to mel scale
        for (int m = 0; m < nmel; ++m) {
          mel[m] = 0.0;
          for (int n = 0; n < noutput; ++n) {
            mel[m] += filters[m][n] * psd[n];
          }
          mel[m] = (mel[m] == 0.0) ? -16.0 : log10(mel[m]);
        }
        mels.push_back(mel);
        // Free allocated things
        fftw_free(inputbuffer);
        fftw_free(outputbuffer);
        fftw_destroy_plan(plan);
        // Print to file for testing
        //for (int m = 0; m < nmelfreq; ++m) {
        //  fprintf(pMEL, "%f,", mel[m]);
        //}
        //fprintf(pMEL, "\n");
        //mels.pop_front();
      }
      // Delete to the end of the processed sound samples
      samplebuffer.clear();//erase(samplebuffer.begin(), samplebuffer.begin() + nwindow * (ninput - noverlap));
    }
};

YaudWinMFCCEpisPort::YaudWinMFCCEpisPort(int nm) {
  // Set up parameters
  // Sound samples
  samplefreq = 16000.0; // Hz
  samplebuffer.clear();
  // Power spectral density
  ninput = 512; // 32ms
  noverlap = 256; // 50%
  noutput = ninput/2+1;
  // Hamming window
  window.resize(ninput);
  windownormsquare = 0.0;
  for (int n = 0; n < ninput; ++n) {
    window[n] = 0.54 - 0.46 * cos(6.28318530718 * n / (ninput - 1));
    windownormsquare += window[n] * window[n];
  }
  psd.resize(noutput);
  // Mel frequencies
  nmel = nm;
  mel.resize(nmel);
  mels.clear();
  // Filterbank
  filters.clear();
  real_t mellow = 1125.0 * log(1.0 + 300.0/700.0);
  real_t melhigh = 1125.0 * log(1.0 + 8000.0/700.0);
  std::vector<int> freqcenters;
  freqcenters.clear();
  freqcenters.push_back(0);
  for (int m = 0; m < nmel; ++m) {
    real_t melcenter = mellow + m * ((melhigh - mellow)/(nmel));
    real_t freqcenter = 700.0 * (exp(melcenter/1125) -1);
    freqcenters.push_back(std::floor((ninput+1)*freqcenter/samplefreq));
  }
  freqcenters.push_back(noutput-1);
  for (int m = 1; m <= nmel; ++m) {
    std::vector<real_t> filter;
    filter.resize(noutput);
    for (int k = 0; k < freqcenters[m-1]; ++k) {
      filter[k] = 0.0;
    }
    for (int k = freqcenters[m-1]; k < freqcenters[m]; ++k) {
      filter[k] = ((real_t)(k - freqcenters[m-1]))/(freqcenters[m] - freqcenters[m-1]);
    }
    for (int k = freqcenters[m]; k <= freqcenters[m+1]; ++k) {
      filter[k] = ((real_t)(freqcenters[m+1] - k))/(freqcenters[m+1] - freqcenters[m]);
    }
    for (int k = freqcenters[m+1]+1; k < noutput; ++k) {
      filter[k] = 0.0;
    }
    filters.push_back(filter);
  }
  // Open File
  //sprintf(psdfile, "melout.csv");
  //pMEL = fopen(melfile,"w");
  //if (pMEL == NULL) {
  //  CkPrintf("Error opening file for writing\n");
  //  CkExit();
  //}
}
#endif

/**************************************************************************
* Class declaration
**************************************************************************/
class YaudWinMFCCEpis : public ModelTmpl < 102, YaudWinMFCCEpis > {
  public:
    /* Constructor */
    YaudWinMFCCEpis() {
      // parameters
      paramlist.resize(1);
      paramlist[0] = "nmel";
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
      tsample = 16 * TICKS_PER_MS;
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
    YaudWinMFCCEpisPort* port;
    yarp::os::BufferedPort<yarp::os::Bottle>* preq;
#endif
    tick_t tsample;
};


/**************************************************************************
* Class methods
**************************************************************************/

// Renew model
//
void YaudWinMFCCEpis::Renew(std::vector<real_t>& state, std::vector<tick_t>& stick) {
#ifdef STACS_WITH_YARP
  // Request new phone
  yarp::os::Bottle& b = preq->prepare();
  b.clear();
  b.addInt((int) port->mels.size());
  preq->write();
#endif
}

// Simulation step
//
tick_t YaudWinMFCCEpis::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
#ifdef STACS_WITH_YARP
  // Check every 16 ms
  // TODO: move to parameters
  // Grab data from port if available
  if ((tsample >= 16 * TICKS_PER_MS) && port->mels.size()) {
    tsample = 0;
    std::vector<real_t> mel = port->mels.front();
    port->mels.pop_front();
    //CkPrintf("Popping front of deque.\n");
    // generate events
    event_t event;
    event.type = EVENT_STIM;
    event.source = REMOTE_EDGE;
    // overlap of 16ms
    // total time of 32ms
    for (idx_t i = 0; i < (idx_t) param[0]; ++i) {
      event.index = i;
      event.diffuse = tdrift + 4 * TICKS_PER_MS;
      event.data = (mel[i] + 16.0)/3;
      events.push_back(event);
      event.diffuse = tdrift + 8 * TICKS_PER_MS;
      event.data = (mel[i] + 16.0)/3;
      events.push_back(event);
      event.diffuse = tdrift + 12 * TICKS_PER_MS;
      event.data = (mel[i] + 16.0)/3;
      events.push_back(event);
      event.diffuse = tdrift + 20 * TICKS_PER_MS;
      event.data = -(mel[i] + 16.0)/3;
      events.push_back(event);
      event.diffuse = tdrift + 24 * TICKS_PER_MS;
      event.data = -(mel[i] + 16.0)/3;
      events.push_back(event);
      event.diffuse = tdrift + 28 * TICKS_PER_MS;
      event.data = -(mel[i] + 16.0)/3;
      events.push_back(event);
    }
  }
  tsample += tdiff;
#endif
  return tdiff;
}

// Open ports
//
void YaudWinMFCCEpis::OpenPorts() {
#ifdef STACS_WITH_YARP
  std::string recv = portname[0] + "/recv";
  std::string request = portname[0] + "/request";
  port = new YaudWinMFCCEpisPort((int) param[0]);
  port->setStrict(); // try not to drop input
  port->useCallback();
  port->open(recv.c_str());
  CkPrintf("Receiving audio samples on %s\n", recv.c_str());
  preq = new yarp::os::BufferedPort<yarp::os::Bottle>();
  preq->open(request.c_str());
  CkPrintf("Requesting audio samples from %s\n", request.c_str());
#endif
}

// Close ports
//
void YaudWinMFCCEpis::ClosePorts() {
#ifdef STACS_WITH_YARP
  port->close();
  yarp::os::Bottle& b = preq->prepare();
  b.clear();
  b.addInt(-1);
  preq->write();
  preq->close();
#endif
}

