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

class YarpAudioPort : public yarp::os::BufferedPort<yarp::sig::Sound> {
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
    YarpAudioPort(int nm);
    //using yarp::os::BufferedPort<yarp::sig::Sound>::onRead;
    virtual void onRead (yarp::sig::Sound& ssig) {
      // Receive sound
      CkPrintf("Received %d samples\n", ssig.getSamples());
      std::size_t noldsample = samplebuffer.size();
      std::size_t nnewsample = ssig.getSamples();
      samplebuffer.resize(noldsample + nnewsample);
      for (std::size_t s = 0; s < nnewsample; ++s) {
        samplebuffer[noldsample+s] = ((real_t)((sample_t)ssig.get(s)))/32768.0;
      }
        CkPrintf("Sample samples: %f, %f, %f, %f, %f\n", samplebuffer[50], samplebuffer[100], samplebuffer[150], samplebuffer[200], samplebuffer[250]);
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
      samplebuffer.erase(samplebuffer.begin(), samplebuffer.begin() + nwindow * (ninput - noverlap));
    }
};

YarpAudioPort::YarpAudioPort(int nm) {
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
class YarpAudio : public NetModelTmpl < 101, YarpAudio > {
  public:
    /* Constructor */
    YarpAudio() {
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
      tsample = 32 * TICKS_PER_MS;
    }

    /* Simulation */
    tick_t Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog);
    void Jump(const event_t& evt, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<aux_t>& aux) { }

    /* Protocol */
    void OpenPorts();
    void ClosePorts();

  private:
#ifdef STACS_WITH_YARP
    YarpAudioPort* port;
#endif
    tick_t tsample;
};


/**************************************************************************
* Class methods
**************************************************************************/

// Simulation step
//
tick_t YarpAudio::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& evtlog) {
  // Check every 32 ms
  // TODO: move to parameters
  // Grab data from port if available
  if ((tsample >= 32 * TICKS_PER_MS) && port->mels.size()) {
    tsample = 0;
    std::vector<real_t> mel = port->mels.front();
    port->mels.pop_front();
    //CkPrintf("Popping front of deque.\n");
    // generate events
    event_t evtpre;
    evtpre.type = EVENT_STIM;
    evtpre.source = REMOTE_EDGE;
    for (idx_t i = 0; i < (idx_t) param[0]; ++i) {
      evtpre.index = i;
      evtpre.diffuse = tdrift;
      evtpre.data = (mel[i] + 16.0)/2;
      evtlog.push_back(evtpre);
      evtpre.diffuse = tdrift + 32 * TICKS_PER_MS;
      evtpre.data = -(mel[i] + 16.0)/2;
      evtlog.push_back(evtpre);
    }
  }
  tsample += tdiff;
  return tdiff;
}

// Open ports
//
void YarpAudio::OpenPorts() {
#ifdef STACS_WITH_YARP
  port = new YarpAudioPort((int) param[0]);
  port->useCallback();
  port->open(portname[0].c_str());
  CkPrintf("Receiving audio samples on %s\n", portname[0].c_str());
#endif
}

// Close ports
//
void YarpAudio::ClosePorts() {
#ifdef STACS_WITH_YARP
  port->close();
#endif
}

