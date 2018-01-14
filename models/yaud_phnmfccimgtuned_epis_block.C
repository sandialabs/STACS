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
      mfccticks.push_back(std::vector<real_t>{-9.6526,-9.1833,-8.714,-8.1664,-7.697,-7.3841,-7.2277,-7.0712,-6.9148,-6.7583,-6.6801,-6.5237,-6.3672,-6.2108,-5.9761});
      mfccticks.push_back(std::vector<real_t>{-9.9625,-9.5186,-9.1635,-8.7196,-8.3645,-7.9206,-7.5655,-7.2992,-7.1216,-6.9441,-6.7665,-6.5889,-6.4114,-6.1451,-5.8787});
      mfccticks.push_back(std::vector<real_t>{-9.8904,-9.3753,-8.9633,-8.6542,-8.2422,-7.9331,-7.5211,-7.2121,-7.006,-6.697,-6.491,-6.2849,-6.0789,-5.8729,-5.5638});
      mfccticks.push_back(std::vector<real_t>{-9.7563,-9.3536,-8.9509,-8.6488,-8.3468,-8.0447,-7.7427,-7.34,-7.0379,-6.8366,-6.5345,-6.3332,-6.0311,-5.7291,-5.427});
      mfccticks.push_back(std::vector<real_t>{-9.8526,-9.4601,-9.1657,-8.8713,-8.577,-8.2826,-7.9882,-7.6939,-7.3995,-7.1051,-6.8108,-6.5164,-6.222,-5.9277,-5.4371});
      mfccticks.push_back(std::vector<real_t>{-9.9718,-9.5775,-9.2817,-8.986,-8.6902,-8.3945,-8.1973,-7.9015,-7.6058,-7.4086,-7.1129,-6.8171,-6.5214,-6.1271,-5.6341});
      mfccticks.push_back(std::vector<real_t>{-9.9646,-9.5014,-9.2234,-8.9455,-8.7602,-8.4822,-8.2969,-8.019,-7.741,-7.5558,-7.2778,-6.9999,-6.7219,-6.3513,-5.8881});
      mfccticks.push_back(std::vector<real_t>{-9.9561,-9.6049,-9.3414,-9.078,-8.8146,-8.639,-8.3756,-8.2,-7.9365,-7.7609,-7.4975,-7.2341,-6.8829,-6.5316,-6.0048});
      mfccticks.push_back(std::vector<real_t>{-10.089,-9.7111,-9.5223,-9.239,-9.0502,-8.8613,-8.5781,-8.3893,-8.2004,-7.9172,-7.7283,-7.4451,-7.0674,-6.6897,-6.2177});
      mfccticks.push_back(std::vector<real_t>{-10.187,-9.8268,-9.5568,-9.3768,-9.1968,-8.9268,-8.7467,-8.5667,-8.2967,-8.1167,-7.8467,-7.4866,-7.2166,-6.7666,-6.3166});
      mfccticks.push_back(std::vector<real_t>{-10.22,-9.8548,-9.6725,-9.399,-9.2166,-9.0343,-8.8519,-8.5784,-8.3961,-8.1225,-7.849,-7.5755,-7.302,-6.8461,-6.3903});
      mfccticks.push_back(std::vector<real_t>{-10.316,-9.964,-9.7001,-9.4362,-9.2603,-9.0843,-8.8204,-8.6445,-8.3806,-8.2047,-7.9408,-7.6769,-7.325,-6.9732,-6.5333});
      mfccticks.push_back(std::vector<real_t>{-10.405,-9.9932,-9.664,-9.4172,-9.2526,-9.0057,-8.8411,-8.5943,-8.4297,-8.1828,-7.936,-7.6891,-7.3599,-7.0308,-6.6193});
      mfccticks.push_back(std::vector<real_t>{-10.465,-10.049,-9.7154,-9.4655,-9.2988,-9.0489,-8.8823,-8.6324,-8.4658,-8.2158,-8.0492,-7.7993,-7.4661,-7.1328,-6.7163});
      mfccticks.push_back(std::vector<real_t>{-10.431,-10.036,-9.7202,-9.4832,-9.3253,-9.0883,-8.9304,-8.7724,-8.5355,-8.3775,-8.1406,-7.9037,-7.6667,-7.3508,-6.9559});
      mfccticks.push_back(std::vector<real_t>{-10.424,-10.031,-9.7161,-9.48,-9.3226,-9.0866,-8.9292,-8.7718,-8.5358,-8.3784,-8.1423,-7.985,-7.6702,-7.4342,-7.0407});
      mfccticks.push_back(std::vector<real_t>{-10.5,-10.033,-9.7225,-9.4892,-9.256,-9.1006,-8.9451,-8.7119,-8.5564,-8.4009,-8.1677,-7.9345,-7.7013,-7.3904,-7.0017});
      mfccticks.push_back(std::vector<real_t>{-10.622,-10.16,-9.8527,-9.5451,-9.3144,-9.1606,-8.9299,-8.7761,-8.6222,-8.3915,-8.2377,-8.007,-7.7763,-7.4687,-7.0842});
      mfccticks.push_back(std::vector<real_t>{-10.769,-10.243,-9.9421,-9.6415,-9.4161,-9.1907,-8.9652,-8.815,-8.6647,-8.4393,-8.289,-8.0636,-7.8381,-7.5376,-7.1619});
      mfccticks.push_back(std::vector<real_t>{-10.866,-10.326,-9.9395,-9.6306,-9.3989,-9.1673,-8.9356,-8.7812,-8.5495,-8.3951,-8.2406,-8.009,-7.7773,-7.4684,-7.0823});
      mfccticks.push_back(std::vector<real_t>{-10.895,-10.349,-9.9594,-9.6475,-9.4135,-9.1796,-8.9456,-8.7117,-8.5557,-8.3997,-8.1658,-7.9318,-7.6979,-7.3859,-7.074});
      mfccticks.push_back(std::vector<real_t>{-10.899,-10.354,-10.043,-9.7317,-9.4204,-9.2648,-9.0314,-8.798,-8.6424,-8.409,-8.2534,-8.02,-7.7866,-7.4754,-7.0864});
      mfccticks.push_back(std::vector<real_t>{-10.981,-10.5,-10.18,-9.9404,-9.7004,-9.4604,-9.2203,-9.0603,-8.8203,-8.6603,-8.4203,-8.1802,-7.9402,-7.5402,-7.1401});
      mfccticks.push_back(std::vector<real_t>{-11.129,-10.813,-10.498,-10.182,-9.9459,-9.7882,-9.5516,-9.315,-9.0784,-8.8419,-8.6053,-8.3687,-8.0533,-7.659,-7.1859});
      mfccticks.push_back(std::vector<real_t>{-11.266,-10.934,-10.684,-10.435,-10.185,-9.936,-9.7697,-9.5203,-9.2708,-9.0214,-8.7719,-8.5225,-8.1899,-7.7741,-7.2752});
      mfccticks.push_back(std::vector<real_t>{-11.408,-11.07,-10.817,-10.648,-10.394,-10.226,-9.9723,-9.8034,-9.5501,-9.2968,-9.0435,-8.7058,-8.3681,-7.946,-7.355});
      mfccticks.push_back(std::vector<real_t>{-11.432,-11.188,-10.945,-10.782,-10.539,-10.376,-10.133,-9.8891,-9.7267,-9.4831,-9.1583,-8.9147,-8.5087,-8.1027,-7.4531});
      mfccticks.push_back(std::vector<real_t>{-11.409,-11.092,-10.854,-10.696,-10.458,-10.3,-10.062,-9.9041,-9.6665,-9.5081,-9.1913,-8.9538,-8.5578,-8.0826,-7.5283});
      mfccticks.push_back(std::vector<real_t>{-11.375,-11.065,-10.833,-10.678,-10.446,-10.291,-10.059,-9.9038,-9.6715,-9.5166,-9.2843,-8.9746,-8.6649,-8.2003,-7.5808});
      mfccticks.push_back(std::vector<real_t>{-11.346,-11.025,-10.785,-10.624,-10.464,-10.224,-10.063,-9.903,-9.6626,-9.5023,-9.2619,-9.0214,-8.6207,-8.22,-7.5788});
      CkAssert(mfccticks.size() == 30);
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
    YaudPhnMFCCImgTunedEpisBlockPort* port;
    yarp::os::BufferedPort<yarp::os::Bottle>* preq;
#endif
    tick_t tsample;
    std::vector<std::vector<real_t>> mfccticks;
};


/**************************************************************************
* Class methods
**************************************************************************/

// Rerun model
//
void YaudPhnMFCCImgTunedEpisBlock::Rerun(std::vector<real_t>& state, std::vector<tick_t>& stick) {
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

