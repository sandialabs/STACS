/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 */

#include "network.h"

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ std::string filebase;


/**************************************************************************
* Network Recording
**************************************************************************/

// Store records
//
void Network::StoreRecord() {
  // Periodic Records
  for (std::size_t r = 0; r < recordlist.size(); ++r) {
    if (recordlist[r].trec <= tsim) {
      // Set next recording
      recordlist[r].trec = tsim + recordlist[r].tfreq;
      // record time
      record.push_back(record_t());
      record.back().drift = tsim;
      // add data
      for (std::size_t i = 0; i < recordlist[r].index.size(); ++i) {
        if (recordlist[r].type[i] == RECORD_STATE) {
          record.back().data.push_back(state[recordlist[r].index[i]][recordlist[r].model[i]][recordlist[r].value[i]]);
        }
        else if (recordlist[r].type[i] == RECORD_STICK) {
          record.back().data.push_back((real_t)(stick[recordlist[r].index[i]][recordlist[r].model[i]][recordlist[r].value[i]]));
        }
        else if (recordlist[r].type[i] == RECORD_COORD) {
          record.back().data.push_back(xyz[recordlist[r].index[i]*3+0]);
          record.back().data.push_back(xyz[recordlist[r].index[i]*3+1]);
          record.back().data.push_back(xyz[recordlist[r].index[i]*3+2]);
        }
      }
    }
  }
  // Event records are recorded in event handling
}

// Build record message for writing
//
mRecord* Network::BuildRecord() {
  /* Bookkeeping */
  idx_t ndata = 0;

  // Count data points
  for (std::size_t i = 0; i < record.size(); ++i) {
    ndata += record[i].data.size();
  }

  // Initialize distribution message
  int msgSize[MSG_Record];
  msgSize[0] = recevt.size();       // diffuse
  msgSize[1] = recevt.size();       // index
  msgSize[2] = recevt.size();       // type
  msgSize[3] = recevt.size()+ndata; // data
  msgSize[4] = record.size();       // drift
  msgSize[5] = record.size()+1;     // xdata
  mRecord *mrecord = new(msgSize, 0) mRecord;
  mrecord->nrecevt = recevt.size();
  mrecord->nrecord = record.size();
  mrecord->iter = iter;
  mrecord->prtidx = prtidx;

  // Pack event information
  for (std::size_t i = 0; i < recevt.size(); ++i) {
    mrecord->diffuse[i] = recevt[i].diffuse;
    mrecord->index[i] = recevt[i].index;
    mrecord->type[i] = recevt[i].type;
    mrecord->data[i] = recevt[i].data;
  }
  
  // Counters
  idx_t jdata = recevt.size();
  
  // Prefix starts at the end of event data
  mrecord->xdata[0] = recevt.size();

  // Pack record information
  for (std::size_t i = 0; i < record.size(); ++i) {
    mrecord->drift[i] = record[i].drift;
    // data
    mrecord->xdata[i+1] = mrecord->xdata[i] + record[i].data.size();
    for (std::size_t j = 0; j < record[i].data.size(); ++j) {
      mrecord->data[jdata++] = record[i].data[j];
    }
  }
  CkAssert(jdata == recevt.size() + ndata);
  
  // Clear records
  record.clear();
  recevt.clear();

  return mrecord;
}

// Send Records for writing (checking)
//
void Network::CheckRecord() {
  // Build record message for saving
  mRecord* mrecord = BuildRecord();
  netdata(datidx).CheckRecord(mrecord);
  
  // Start a new cycle (checked data sent)
  thisProxy(prtidx).Cycle();
}

// Send Records for writing
//
void Network::SaveRecord() {
  // Build record message for saving
  mRecord* mrecord = BuildRecord();
  netdata(datidx).SaveRecord(mrecord);
}


/**************************************************************************
* Netdata recording to file
**************************************************************************/

// Write periodic records to file
//
void NetData::CheckRecord(mRecord *msg) {
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

// Write periodic records to file (final)
//
void NetData::SaveRecord(mRecord *msg) {
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
* Writing Record Information
**************************************************************************/

// Writing Records to file
//
void NetData::WriteRecord() {
  /* File operations */
  FILE *pRecord;
  char recfile[100];

  // Open File
  sprintf(recfile, "%s.record.%" PRIidx ".%" PRIidx "", filebase.c_str(), datidx, records[0]->iter);
  pRecord = fopen(recfile,"w");
  if (pRecord == NULL) {
    CkPrintf("Error opening files for recording %" PRIidx "\n", datidx);
    CkExit();
  }

  // Loop through parts
  for (idx_t k = 0; k < nprt; ++k) {
    // Loop through events
    for (idx_t e = 0; e < records[k]->nrecevt; ++e) {
      // types lacking data
      if (records[k]->type[e] == EVTYPE_SPIKE) {
        fprintf(pRecord, "%" PRIidx " %" PRItickhex " %" PRIidx "\n",
            records[k]->type[e], records[k]->diffuse[e], records[k]->index[e]);
      }
      // events with data
      else {
        fprintf(pRecord, "%" PRIidx " %" PRItickhex " %" PRIidx " %" PRIrealfull "\n",
            records[k]->type[e], records[k]->diffuse[e], records[k]->index[e], records[k]->data[e]);
      }
    }
    // Loop through records
    for (idx_t r = 0; r < records[k]->nrecord; ++r) {
      // 'event type' 0 followed by amount of data
      fprintf(pRecord, "0 %" PRItickhex " %" PRIidx "", records[k]->drift[r], (records[k]->xdata[r+1] - records[k]->xdata[r]));
      // data
      for (idx_t d = records[k]->xdata[r]; d < records[k]->xdata[r+1]; ++d) {
        fprintf(pRecord, " %" PRIrealfull "", records[k]->data[d]);
      }
      // one line per record
      fprintf(pRecord, "\n");
    }
  }
  
  // Cleanup
  fclose(pRecord);
}
