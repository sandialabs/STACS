/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/


/**************************************************************************
* Network Recording (local)
**************************************************************************/

// Store records
//
void Network::AddRecord() {
  // Periodic Records
  for (std::size_t r = 0; r < recordlist.size(); ++r) {
    if (recordlist[r].trec <= tsim) {
      // Set next recording
      recordlist[r].trec = tsim + recordlist[r].tfreq;
      // record time
      record.push_back(record_t());
      record.back().recidx = r;
      record.back().trec = tsim;
      record.back().state.clear();
      record.back().stick.clear();
      record.back().index.clear();
      // add data accordingly
      if (recordlist[r].rectype == RECORD_STATE) {
        for (std::size_t i = 0; i < recordlist[r].recvtxidx.size(); ++i) {
          for (std::size_t j = 0; j < recordlist[r].recedgidx[recordlist[r].recvtxidx[i]].size(); ++j) {
            record.back().state.push_back(state[recordlist[r].recvtxidx[i]][recordlist[r].recedgidx[recordlist[r].recvtxidx[i]][j]][recordlist[r].recsttidx]);
          }
        }
      }
      else if (recordlist[r].rectype == RECORD_STICK) {
        for (std::size_t i = 0; i < recordlist[r].recvtxidx.size(); ++i) {
          for (std::size_t j = 0; j < recordlist[r].recedgidx[recordlist[r].recvtxidx[i]].size(); ++j) {
            record.back().stick.push_back(stick[recordlist[r].recvtxidx[i]][recordlist[r].recedgidx[recordlist[r].recvtxidx[i]][j]][recordlist[r].recsttidx]);
          }
        }
      }
      /*
      else if (recordlist[r].rectype == RECORD_COORD) {
        record.back().data.push_back(xyz[recordlist[r].index[i]*3+0]);
        record.back().data.push_back(xyz[recordlist[r].index[i]*3+1]);
        record.back().data.push_back(xyz[recordlist[r].index[i]*3+2]);
      }
      */
    }
  }
}


/**************************************************************************
* Network Recording (to file)
**************************************************************************/

// Send Records for writing
//
void Network::SaveRecord() {
  // Build record message for saving
  mRecord* mrecord = BuildRecord();
  netdata(datidx).SaveRecord(mrecord);
  
  // Start a new cycle (checked data sent)
  cyclepart.send();
}

// Send Records for writing (final)
//
void Network::SaveFinalRecord() {
  // Build record message for saving
  mRecord* mrecord = BuildRecord();
  netdata(datidx).SaveFinalRecord(mrecord);
}

// Write records to file
//
void Netdata::SaveRecord(mRecord *msg) {
  // Stash record
  records[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++rprt == nprt) {
    rprt = 0;
    
    // Write data
    //CkPrintf("  Writing records %" PRIidx "\n", datidx);
    WriteRecord();

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete records[i];
    }
  }
}

// Write records to file (final)
//
void Netdata::SaveFinalRecord(mRecord *msg) {
  // Stash record
  records[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++rprt == nprt) {
    rprt = 0;
    
    // Write data
    //CkPrintf("  Writing records %" PRIidx "\n", datidx);
    WriteRecord();

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete records[i];
    }
    
    // Return control to main (halting)
    contribute(0, NULL, CkReduction::nop);
  }
}


/**************************************************************************
* Build Messages
**************************************************************************/

// Build record message for writing
//
mRecord* Network::BuildRecord() {
  /* Bookkeeping */
  idx_t ndata = 0;
  idx_t ndiffuse = 0;
  idx_t nindex = 0;

  // Count data points
  for (std::size_t i = 0; i < record.size(); ++i) {
    ndata += record[i].state.size();
    ndiffuse += record[i].stick.size();
    nindex += record[i].index.size();
  }

  // Initialize distribution message
  int msgSize[MSG_Record];
  msgSize[0] = evtlog.size()+ndiffuse;      // diffuse
  msgSize[1] = evtlog.size()+record.size(); // type
  msgSize[2] = evtlog.size();               // source
  msgSize[3] = evtlog.size()+nindex;        // index
  msgSize[4] = evtlog.size()+ndata;         // data
  msgSize[5] = record.size();               // drift
  msgSize[6] = record.size()+1;             // xdata
  msgSize[7] = record.size()+1;             // xdiffuse
  msgSize[8] = record.size()+1;             // xindex
  mRecord *mrecord = new(msgSize, 0) mRecord;
  mrecord->nevtlog = evtlog.size();
  mrecord->nrecord = record.size();
  mrecord->iter = iter;
  mrecord->prtidx = prtidx;

  // Pack event information
  for (std::size_t i = 0; i < evtlog.size(); ++i) {
    mrecord->diffuse[i] = evtlog[i].diffuse;
    mrecord->type[i] = evtlog[i].type;
    mrecord->source[i] = evtlog[i].source;
    mrecord->index[i] = evtlog[i].index;
    mrecord->data[i] = evtlog[i].data;
  }
  
  // Counters
  std::size_t jdata = evtlog.size();
  std::size_t jdiffuse = evtlog.size();
  std::size_t jindex = evtlog.size();
  
  // Prefix starts at end of event data
  mrecord->xdata[0] = evtlog.size();
  mrecord->xdiffuse[0] = evtlog.size();
  mrecord->xindex[0] = evtlog.size();

  // Pack record information
  for (std::size_t i = 0; i < record.size(); ++i) {
    mrecord->drift[i] = record[i].trec;
    mrecord->type[evtlog.size()+i] = record[i].recidx;
    // data
    mrecord->xdata[i+1] = mrecord->xdata[i] + record[i].state.size();
    for (std::size_t j = 0; j < record[i].state.size(); ++j) {
      mrecord->data[jdata++] = record[i].state[j];
    }
    // diffuse
    mrecord->xdiffuse[i+1] = mrecord->xdiffuse[i] + record[i].stick.size();
    for (std::size_t j = 0; j < record[i].stick.size(); ++j) {
      mrecord->diffuse[jdiffuse++] = record[i].stick[j];
    }
    // index
    mrecord->xindex[i+1] = mrecord->xindex[i] + record[i].index.size();
    for (std::size_t j = 0; j < record[i].index.size(); ++j) {
      mrecord->index[jindex++] = record[i].index[j];
    }
  }
  CkAssert(jdata == evtlog.size() + ndata);
  CkAssert(jdiffuse == evtlog.size() + ndiffuse);
  CkAssert(jindex == evtlog.size() + nindex);
  
  // Clear records
  evtlog.clear();
  record.clear();

  return mrecord;
}
