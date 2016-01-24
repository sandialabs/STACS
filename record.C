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

// Store periodic records
//
void Network::StoreRecord() {
  // record time
  record.push_back(record_t());
  record.back().drift = tsim;

  // Add data
  for (std::size_t i = 0; i < recordlist.size(); ++i) {
    if (recordlist[i].type == RECORD_STATE) {
      record.back().data.push_back(state[recordlist[i].vertex][recordlist[i].model][recordlist[i].index]);
    }
    else if (recordlist[i].type == RECORD_STICK) {
      record.back().data.push_back((real_t)(stick[recordlist[i].vertex][recordlist[i].model][recordlist[i].index]));
    }
    else if (recordlist[i].type == RECORD_COORD) {
      record.back().data.push_back(xyz[((recordlist[i].vertex)*3)+0]);
      record.back().data.push_back(xyz[((recordlist[i].vertex)*3)+1]);
      record.back().data.push_back(xyz[((recordlist[i].vertex)*3)+2]);
    }
  }
}

// Build message for recording (events)
//
mEvent* Network::BuildRecevt() {
  // Initialize distribution message
  int msgSize[MSG_Event];
  msgSize[0] = recevt.size();     // diffuse
  msgSize[1] = recevt.size();     // index
  msgSize[2] = recevt.size();     // type
  msgSize[3] = recevt.size();     // data
  mEvent *mrecevt = new(msgSize, 0) mEvent;
  mrecevt->nevent = recevt.size();
  mrecevt->iter = iter;
  mrecevt->prtidx = prtidx;
  
  // Pack record information
  for (std::size_t i = 0; i < recevt.size(); ++i) {
    mrecevt->diffuse[i] = recevt[i].diffuse;
    mrecevt->index[i] = recevt[i].index;
    mrecevt->type[i] = recevt[i].type;
    mrecevt->data[i] = recevt[i].data;
  }

  return mrecevt;
}

// Build message for recording (periodic)
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
  msgSize[0] = record.size();   // drift
  msgSize[1] = record.size()+1; // xdata
  msgSize[2] = ndata;           // data
  mRecord *mrecord = new(msgSize, 0) mRecord;
  mrecord->nrecord = record.size();
  mrecord->iter = iter;
  mrecord->prtidx = prtidx;

  // Counters
  idx_t jdata = 0;
  
  // Prefixes start at 0
  mrecord->xdata[0] = 0;

  // Pack record information
  for (std::size_t i = 0; i < record.size(); ++i) {
    mrecord->drift[i] = record[i].drift;
    // data
    mrecord->xdata[i+1] = mrecord->xdata[i] + record[i].data.size();
    for (std::size_t j = 0; j < record[i].data.size(); ++j) {
      mrecord->data[jdata++] = record[i].data[j];
    }
  }
  CkAssert(jdata == ndata);

  return mrecord;
}


/**************************************************************************
* Netdata recording to file
**************************************************************************/

// Write irregular events to file
//
void NetData::CheckRecevt(mEvent *msg) {
  // Stash record
  recevts[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++eprt == nprt) {
    eprt = 0;
    
    // Write data
    CkPrintf("  Writing events %" PRIidx "\n", datidx);
    WriteRecevt();

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete recevts[i];
    }
  }
}

// Write irregular events to file (final)
//
void NetData::SaveRecevt(mEvent *msg) {
  // Stash record
  recevts[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++eprt == nprt) {
    eprt = 0;
    
    // Write data
    CkPrintf("  Writing events %" PRIidx "\n", datidx);
    WriteRecevt();

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete recevts[i];
    }
    
    // Return control to main
    contribute(0, NULL, CkReduction::nop);
  }
}

// Write periodic records to file
//
void NetData::CheckRecord(mRecord *msg) {
  // Stash record
  records[msg->prtidx - xprt] = msg;
  
  // Wait for all parts
  if (++rprt == nprt) {
    rprt = 0;
    
    // Write data
    CkPrintf("  Writing records %" PRIidx "\n", datidx);
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
    CkPrintf("  Writing records %" PRIidx "\n", datidx);
    WriteRecord();

    // Cleanup stash
    for (idx_t i = 0; i < nprt; ++i) {
      delete records[i];
    }
    
    // Return control to main
    contribute(0, NULL, CkReduction::nop);
  }
}


/**************************************************************************
* Writing Record Information
**************************************************************************/

// Writing Events to file
//
void NetData::WriteRecevt() {
  /* File operations */
  FILE *pRecevt;
  char recfile[100];

  // Open File
  sprintf(recfile, "%s.recevt.%" PRIidx ".%" PRIidx "", filebase.c_str(), datidx, recevts[0]->iter);
  pRecevt = fopen(recfile,"w");
  if (pRecevt == NULL) {
    CkPrintf("Error opening files for recording %" PRIidx "\n", datidx);
    CkExit();
  }

  // Loop through parts
  for (idx_t k = 0; k < nprt; ++k) {
    // Loop through records
    for (idx_t e = 0; e < recevts[k]->nevent; ++e) {
      if (recevts[k]->type[e] == EVENT_SPIKE) {
        fprintf(pRecevt, "%" PRItickhex " %" PRIidx " %" PRIidx "\n",
            recevts[k]->diffuse[e], recevts[k]->index[e], recevts[k]->type[e]);
      }
      else {
        fprintf(pRecevt, "%" PRItickhex " %" PRIidx " %" PRIidx " %" PRIrealfull "\n",
            recevts[k]->diffuse[e], recevts[k]->index[e], recevts[k]->type[e], recevts[k]->data[e]);
      }
    }
  }
  
  // Cleanup
  fclose(pRecevt);
}

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
    // Loop through records
    for (idx_t r = 0; r < records[k]->nrecord; ++r) {
      fprintf(pRecord, "%" PRItickhex " %" PRIidx "", records[k]->drift[r], (records[k]->xdata[r+1] - records[k]->xdata[r]));
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
