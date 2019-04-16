/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

#ifdef STACS_WITH_YARP
#include <yarp/sig/Sound.h>
#include <fftw3.h>
#include <mutex>
#include <condition_variable>

/**************************************************************************
* Ports
**************************************************************************/

class YaudPhnMFCCImgTunedEpisBlockPort : public yarp::os::BufferedPort<yarp::sig::Sound> {
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
    std::mutex mtx;
    std::condition_variable cvar;
    bool ready;
    std::vector<std::vector<real_t>> mels;
    YaudPhnMFCCImgTunedEpisBlockPort(int nm);
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
      mels.clear();
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
      std::unique_lock<std::mutex> lck(mtx);
      ready = true;
      cvar.notify_all();
    }
};

YaudPhnMFCCImgTunedEpisBlockPort::YaudPhnMFCCImgTunedEpisBlockPort(int nm) {
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
  ready = false;
}
#endif

/**************************************************************************
* Class declaration
**************************************************************************/
class YaudPhnMFCCImgTunedEpisBlock : public ModelTmpl < 111, YaudPhnMFCCImgTunedEpisBlock > {
  public:
    /* Constructor */
    YaudPhnMFCCImgTunedEpisBlock() {
      // parameters
      paramlist.resize(2);
      paramlist[0] = "nmel";
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

      // Assuming nmel = 30
      mfccticks.push_back(std::vector<real_t>{-9.418,-8.7922,-7.9317,-7.4624,-7.1495,-6.9148,-6.7583,-6.6019,-6.3672,-6.1326});
      mfccticks.push_back(std::vector<real_t>{-9.785,-9.1635,-8.6308,-8.0094,-7.5655,-7.2104,-6.8553,-6.5889,-6.3226,-6.0563});
      mfccticks.push_back(std::vector<real_t>{-9.5814,-9.0663,-8.4482,-7.9331,-7.5211,-7.006,-6.697,-6.3879,-6.0789,-5.6668});
      mfccticks.push_back(std::vector<real_t>{-9.5549,-9.0515,-8.5481,-8.1454,-7.642,-7.1386,-6.7359,-6.4338,-6.0311,-5.5277});
      mfccticks.push_back(std::vector<real_t>{-9.6563,-9.1657,-8.7732,-8.2826,-7.8901,-7.4976,-7.007,-6.6145,-6.222,-5.6333});
      mfccticks.push_back(std::vector<real_t>{-9.676,-9.2817,-8.8874,-8.493,-8.0987,-7.7044,-7.31,-6.9157,-6.5214,-5.9299});
      mfccticks.push_back(std::vector<real_t>{-9.6867,-9.2234,-8.9455,-8.5749,-8.2043,-7.8337,-7.4631,-7.0925,-6.6293,-6.0734});
      mfccticks.push_back(std::vector<real_t>{-9.7805,-9.3414,-8.9902,-8.639,-8.3756,-8.0243,-7.6731,-7.3219,-6.8829,-6.2682});
      mfccticks.push_back(std::vector<real_t>{-9.8999,-9.5223,-9.1446,-8.8613,-8.5781,-8.2948,-7.9172,-7.5395,-7.0674,-6.4065});
      mfccticks.push_back(std::vector<real_t>{-10.007,-9.6468,-9.2868,-9.0168,-8.7467,-8.3867,-8.0267,-7.6667,-7.1266,-6.4966});
      mfccticks.push_back(std::vector<real_t>{-10.037,-9.6725,-9.3078,-9.0343,-8.7607,-8.4872,-8.1225,-7.6667,-7.2108,-6.6638});
      mfccticks.push_back(std::vector<real_t>{-10.14,-9.7001,-9.3482,-9.0843,-8.8204,-8.4686,-8.1167,-7.7649,-7.325,-6.7093});
      mfccticks.push_back(std::vector<real_t>{-10.158,-9.664,-9.3349,-9.088,-8.7588,-8.4297,-8.1005,-7.7714,-7.3599,-6.7839});
      mfccticks.push_back(std::vector<real_t>{-10.215,-9.7154,-9.3822,-9.1322,-8.799,-8.5491,-8.2158,-7.8826,-7.4661,-6.9662});
      mfccticks.push_back(std::vector<real_t>{-10.273,-9.7991,-9.4832,-9.1673,-8.8514,-8.6145,-8.2986,-7.9827,-7.6667,-7.1139});
      mfccticks.push_back(std::vector<real_t>{-10.188,-9.7161,-9.4013,-9.1653,-8.9292,-8.6145,-8.2997,-8.0636,-7.6702,-7.1981});
      mfccticks.push_back(std::vector<real_t>{-10.267,-9.8002,-9.4115,-9.1006,-8.8674,-8.6341,-8.3232,-8.0123,-7.6236,-7.1572});
      mfccticks.push_back(std::vector<real_t>{-10.391,-9.8527,-9.4682,-9.1606,-8.9299,-8.6222,-8.3915,-8.0839,-7.6994,-7.238});
      mfccticks.push_back(std::vector<real_t>{-10.543,-9.9421,-9.5664,-9.1907,-8.9652,-8.6647,-8.4393,-8.1387,-7.8381,-7.3873});
      mfccticks.push_back(std::vector<real_t>{-10.557,-10.017,-9.5534,-9.1673,-8.8584,-8.6267,-8.3951,-8.0862,-7.7773,-7.314});
      mfccticks.push_back(std::vector<real_t>{-10.583,-9.9594,-9.5695,-9.1796,-8.8676,-8.6337,-8.3218,-8.0098,-7.6979,-7.23});
      mfccticks.push_back(std::vector<real_t>{-10.587,-10.043,-9.576,-9.2648,-8.9536,-8.7202,-8.409,-8.0978,-7.7866,-7.3198});
      mfccticks.push_back(std::vector<real_t>{-10.741,-10.26,-9.8604,-9.5404,-9.2203,-8.9003,-8.5803,-8.2602,-7.8602,-7.3801});
      mfccticks.push_back(std::vector<real_t>{-10.971,-10.498,-10.104,-9.7882,-9.4727,-9.1573,-8.8419,-8.4476,-8.0533,-7.4225});
      mfccticks.push_back(std::vector<real_t>{-11.1,-10.684,-10.352,-10.019,-9.6866,-9.354,-9.0214,-8.6056,-8.1067,-7.5247});
      mfccticks.push_back(std::vector<real_t>{-11.239,-10.901,-10.563,-10.226,-9.9723,-9.6345,-9.2124,-8.8747,-8.2837,-7.6083});
      mfccticks.push_back(std::vector<real_t>{-11.27,-10.945,-10.701,-10.376,-10.052,-9.7267,-9.4019,-8.9959,-8.5087,-7.7779});
      mfccticks.push_back(std::vector<real_t>{-11.25,-10.934,-10.617,-10.3,-10.062,-9.7457,-9.4289,-9.0329,-8.5578,-7.8451});
      mfccticks.push_back(std::vector<real_t>{-11.22,-10.91,-10.601,-10.291,-10.059,-9.7489,-9.4392,-9.1295,-8.5874,-7.8906});
      mfccticks.push_back(std::vector<real_t>{-11.185,-10.865,-10.544,-10.304,-9.9832,-9.7427,-9.4222,-9.1016,-8.6207,-7.8994});
      CkAssert(mfccticks.size() == 30);
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
    YaudPhnMFCCImgTunedEpisBlockPort* port;
    yarp::os::BufferedPort<yarp::os::Bottle>* preq;
#endif
    tick_t tsample;
    std::vector<std::vector<real_t>> mfccticks;
};


/**************************************************************************
* Class methods
**************************************************************************/

// Renew model
//
void YaudPhnMFCCImgTunedEpisBlock::Renew(std::vector<real_t>& state, std::vector<tick_t>& stick) {
#ifdef STACS_WITH_YARP
  if (tsample == 0 && !(port->mels.size())) {
    // Request new phone
    yarp::os::Bottle& b = preq->prepare();
    b.clear();
    b.addInt((int) port->mels.size());
    preq->write();
  }
#endif
  tsample = 0;
}

// Simulation step
//
tick_t YaudPhnMFCCImgTunedEpisBlock::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {
#ifdef STACS_WITH_YARP
  // Grab data from port if available
  if (tsample == 0) {
    std::unique_lock<std::mutex> lck(port->mtx);
    while (!(port->ready)) {
      port->cvar.wait(lck);
    }
    port->ready = false;
    tsample = TICKS_PER_MS;
    std::vector<real_t> mel = port->mels.front();
    for (std::size_t ms = 1; ms < port->mels.size(); ++ms) {
      for (std::size_t m = 0; m < mel.size(); ++m) {
        mel[m] = (mel[m]*ms + port->mels[ms][m])/(ms+1);
      }
    }
    //CkPrintf("Popping front of deque.\n");
    // generate events
    event_t event;
    event.type = EVENT_STIM;
    event.source = REMOTE_EDGE;
    // total time of 20ms
    for (idx_t i = 0; i < (idx_t) param[0]; ++i) {
      event.index = i;
      if (mel[i] < mfccticks[i][0]) {
        continue;
      }
      for (std::size_t mt = 1; mt < mfccticks[i].size(); ++mt) {
        if (mel[i] > mfccticks[i][mfccticks[i].size() - mt]) {
          event.diffuse = tdrift + (mt-1) * TICKS_PER_MS;
          event.data = param[1];
          events.push_back(event);
          event.diffuse = tdrift + (mt) * TICKS_PER_MS;
          event.data = -param[1];
          events.push_back(event);
          break;
        }
      }
    }
    // Request new phone
    yarp::os::Bottle& b = preq->prepare();
    b.clear();
    b.addInt((int) port->mels.size());
    preq->write();
  }
#endif
  return tdiff;
}

// Open ports
//
void YaudPhnMFCCImgTunedEpisBlock::OpenPorts() {
#ifdef STACS_WITH_YARP
  std::string recv = portname[0] + "/recv";
  std::string request = portname[0] + "/request";
  port = new YaudPhnMFCCImgTunedEpisBlockPort((int) param[0]);
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
void YaudPhnMFCCImgTunedEpisBlock::ClosePorts() {
#ifdef STACS_WITH_YARP
  port->close();
  yarp::os::Bottle& b = preq->prepare();
  b.clear();
  b.addInt(-1);
  preq->write();
  preq->close();
#endif
}

