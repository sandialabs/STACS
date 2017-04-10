/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchrnous Cortical Streams (stacs)
 *
 * iocsr.C
 * Handles condensed sparse row format
 */

#include "stacs.h"
#include "network.h"

// Maximum size of input line (bytes)
#define MAXLINE 1280000

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ idx_t npnet;
extern /*readonly*/ std::string filebase;
extern /*readonly*/ std::string filemod;
extern /*readonly*/ std::string filedir;
extern /*readonly*/ std::string recdir;
extern /*readonly*/ std::string pngdir;


/**************************************************************************
* Network Distribution
**************************************************************************/

// Read graph distribution
//
int Main::ReadDist() {
  /* File operations */
  FILE *pDist;
  char csrfile[100];
  char *line;
  char *oldstr, *newstr;

  // Prepare buffer
  line = new char[MAXLINE];
  
  // Open file for reading
  CkPrintf("Reading network distribution\n");//from %s/%s.dist\n", filedir.c_str(), filebase.c_str());
  sprintf(csrfile, "%s/%s.dist", filedir.c_str(), filebase.c_str());
  pDist = fopen(csrfile,"r");
  if (pDist == NULL || line == NULL) {
    return 1;
  }

  netdist.resize(npnet+1);

  // Get distribution info
  for (idx_t i = 0; i < npnet+1; ++i) {
    while(fgets(line, MAXLINE, pDist) && line[0] == '%');
    oldstr = line;
    newstr = NULL;
    // vtxdist
    netdist[i].nvtx = strtoidx(oldstr, &newstr, 10);
    oldstr = newstr;
    // edgdist
    netdist[i].nedg = strtoidx(oldstr, &newstr, 10);
    oldstr = newstr;
    // statedist
    netdist[i].nstate = strtoidx(oldstr, &newstr, 10);
    oldstr = newstr;
    // stickdist
    netdist[i].nstick = strtoidx(oldstr, &newstr, 10);
    oldstr = newstr;
    // eventdist
    netdist[i].nevent = strtoidx(oldstr, &newstr, 10);
    oldstr = newstr;
    // prtidx
    netdist[i].prtidx = 0;
  }

  // Cleanup
  fclose(pDist);
  delete[] line;

  return 0;
}

// Write graph adjacency distribution
//
int Main::WriteDist() {
  /* File operations */
  FILE *pDist;
  char csrfile[100];
  /* Bookkeeping */
  idx_t nvtx;
  idx_t nedg;
  idx_t nstate;
  idx_t nstick;
  idx_t nevent;

  // Open File
  CkPrintf("Writing network distribution\n");
  sprintf(csrfile, "%s/%s%s.dist", filedir.c_str(), filebase.c_str(), filemod.c_str());//(check ? ".check" : filemod.c_str()));
  pDist = fopen(csrfile,"w");
  if (pDist == NULL) {
    printf("Error opening file for writing\n");
    return 1;
  }

  // Sort distribution
  std::sort(netdist.begin(), netdist.end());

  // Write to file
  nvtx = nedg = nstate = nstick = nevent = 0;
  fprintf(pDist, "%" PRIidx " %" PRIidx " %" PRIidx " %" PRIidx " %" PRIidx "\n", nvtx, nedg, nstate, nstick, nevent);
  for (std::size_t i = 0; i < netdist.size(); ++i) {
    nvtx += netdist[i].nvtx;
    nedg += netdist[i].nedg;
    nstate += netdist[i].nstate;
    nstick += netdist[i].nstick;
    nevent += netdist[i].nevent;
    fprintf(pDist, "%" PRIidx " %" PRIidx " %" PRIidx " %" PRIidx " %" PRIidx "\n", nvtx, nedg, nstate, nstick, nevent);
  }

  // Cleanup
  fclose(pDist);

  return 0;
}


/**************************************************************************
* Network data files
**************************************************************************/

// Read in network files
//
void Netdata::ReadNetwork() {
  /* File operations */
  FILE *pCoord;
  FILE *pAdjcy;
  FILE *pState;
  FILE *pEvent;
  char csrfile[100];
  char *line;
  char *oldstr, *newstr;
  /* Bookkeeping */
  idx_t jadjcy; // edges
  idx_t jedgmodidx; // edgmodidx
  idx_t jstate; // state
  idx_t jstick; // time state
  idx_t jevent; // event
  
  // Prepare buffer
  line = new char[MAXLINE];
  
  // Open files for reading
  sprintf(csrfile, "%s/%s.coord.%" PRIidx "", filedir.c_str(), filebase.c_str(), datidx);
  pCoord = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s.adjcy.%" PRIidx "", filedir.c_str(), filebase.c_str(), datidx);
  pAdjcy = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s.state.%" PRIidx "", filedir.c_str(), filebase.c_str(), datidx);
  pState = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s.event.%" PRIidx "", filedir.c_str(), filebase.c_str(), datidx);
  pEvent = fopen(csrfile,"r");
  if (pCoord == NULL || pAdjcy == NULL || pState == NULL ||
      pEvent == NULL || line == NULL) {
    CkPrintf("Error opening graph files on %" PRIidx "\n", datidx);
    CkExit();
  }

  // Read in Parts
  for (idx_t k = 0; k < nprt; ++k) {
    idx_t nvtx = vtxdist[xprt+k+1] - vtxdist[xprt+k];
    idx_t nedg = edgdist[xprt+k+1] - edgdist[xprt+k];
    idx_t nstate = statedist[xprt+k+1] - statedist[xprt+k];
    idx_t nstick = stickdist[xprt+k+1] - stickdist[xprt+k];
    idx_t nevent = eventdist[xprt+k+1] - eventdist[xprt+k];

    // Initialize partition data message
    int msgSize[MSG_Part];
    msgSize[0] = npnet+1;       // vtxdist
    msgSize[1] = nvtx;          // vtxmodidx
    msgSize[2] = nvtx*3;        // xyz
    msgSize[3] = nvtx+1;        // xadj
    msgSize[4] = nedg;          // adjcy
    msgSize[5] = nedg;          // edgmodidx
    msgSize[6] = nstate;        // state
    msgSize[7] = nstick;        // stick
    msgSize[8] = nvtx+1;        // xevent
    msgSize[9] = nevent;        // diffuse
    msgSize[10] = nevent;       // type
    msgSize[11] = nevent;       // source
    msgSize[12] = nevent;       // index
    msgSize[13] = nevent;       // data
    parts[k] = new(msgSize, 0) mPart;
    
    // Data sizes
    parts[k]->nvtx = nvtx;
    parts[k]->nedg = nedg;
    parts[k]->nstate = nstate;
    parts[k]->nstick = nstick;
    parts[k]->nevent = nevent;
    parts[k]->prtidx = xprt+k;

    // vtxdist
    for (idx_t i = 0; i < npnet+1; ++i) {
      parts[k]->vtxdist[i] = vtxdist[i];
    }

    // first prefix entry is zero
    parts[k]->xadj[0] = 0;
    parts[k]->xevent[0] = 0;
    
    // initialize counters
    jadjcy = 0;
    jedgmodidx = 0;
    jstate = 0;
    jstick = 0;
    jevent = 0;

    // Extract Graph Information
    for (idx_t i = 0; i < nvtx; ++i) {
      // Read in line (coordinates)
      while(fgets(line, MAXLINE, pCoord) && line[0] == '%');
      oldstr = line;
      newstr = NULL;
      // xyz
      for (idx_t j = 0; j < 3; ++j) {
        real_t coord = strtoreal(oldstr, &newstr);
        oldstr = newstr;
        parts[k]->xyz[i*3+j] = coord;
      }

      // Read line (vertex adjacency)
      while(fgets(line, MAXLINE, pAdjcy) && line[0] == '%');
      oldstr = line;
      newstr = NULL;
      for(;;) {
        idx_t edg = strtoidx(oldstr, &newstr, 10);
        // check for end of line
        if (edg == 0 && oldstr != line)
          break;
        oldstr = newstr;
        // adjcy
        parts[k]->adjcy[jadjcy++] = edg;
      }
      // xadj
      parts[k]->xadj[i+1] = jadjcy;
      
      // Read line (models and states)
      while(fgets(line, MAXLINE, pState) && line[0] == '%');
      oldstr = line;
      newstr = NULL;
      // extract model index from name
      idx_t modidx = strtomodidx(oldstr, &newstr);
      oldstr = newstr;

      // first one is the vertex
      CkAssert(modidx != IDX_T_MAX);
      // modidx
      parts[k]->vtxmodidx[i] = modidx;
      // TODO: update models
      for(idx_t s = 0; s < netmodel[modidx]->getNState(); ++s) {
        real_t state = strtoreal(oldstr, &newstr);
        oldstr = newstr;
        // state
        parts[k]->state[jstate++] = state;
      }
      for(idx_t s = 0; s < netmodel[modidx]->getNStick(); ++s) {
        tick_t stick = strtotick(oldstr, &newstr, 16);
        oldstr = newstr;
        // state
        parts[k]->stick[jstick++] = stick;
      }
      
      // then the edges
      modidx = strtomodidx(oldstr, &newstr);
      oldstr = newstr;
      //CkPrintf(" Modidx: %s\n", oldstr);
      while (modidx != IDX_T_MAX) {
        // modidx
        parts[k]->edgmodidx[jedgmodidx++] = modidx;
        for(idx_t s = 0; s < netmodel[modidx]->getNState(); ++s) {
          real_t state = strtoreal(oldstr, &newstr);
          oldstr = newstr;
          // state
          parts[k]->state[jstate++] = state;
        }
        for(idx_t s = 0; s < netmodel[modidx]->getNStick(); ++s) {
          tick_t stick = strtotick(oldstr, &newstr, 16);
          oldstr = newstr;
          // state
          parts[k]->stick[jstick++] = stick;
        }
        modidx = strtomodidx(oldstr, &newstr);
        oldstr = newstr;
      }

      // Read line (events)
      while(fgets(line, MAXLINE, pEvent) && line[0] == '%');
      oldstr = line;
      newstr = NULL;
      // xevent
      jevent += strtoidx(oldstr, &newstr, 10);
      parts[k]->xevent[i+1] = jevent;
      oldstr = newstr;
      for(idx_t e = parts[k]->xevent[i]; e < parts[k]->xevent[i+1]; ++e) {
        // diffuse
        parts[k]->diffuse[e] = strtotick(oldstr, &newstr, 16);
        oldstr = newstr;
        // type
        idx_t type = strtoidx(oldstr, &newstr, 10);
        oldstr = newstr;
        parts[k]->type[e] = type;
        // source
        parts[k]->source[e] = strtoidx(oldstr, &newstr, 10);
        oldstr = newstr;
        // index
        parts[k]->index[e] = strtoidx(oldstr, &newstr, 10);
        oldstr = newstr;
        // event types lacking data
        if (type == EVENT_SPIKE) {
          // data
          parts[k]->data[e] = 0.0;
        }
        // event types with data
        else {
          // data
          parts[k]->data[e] = strtoreal(oldstr, &newstr);
          oldstr = newstr;
        }
      }
    }
    // Sanity check
    CkAssert(jadjcy == nedg);
    CkAssert(jedgmodidx == nedg);
    CkAssert(jstate == nstate);
    CkAssert(jstick == nstick);
    CkAssert(jevent == nevent);
  }

  // Cleanup
  fclose(pCoord);
  fclose(pAdjcy);
  fclose(pState);
  fclose(pEvent);
  delete[] line;
}
    
// Write out network files
//
void Netdata::WriteNetwork() {
  /* File operations */
  FILE *pCoord;
  FILE *pAdjcy;
  FILE *pState;
  FILE *pEvent;
  char csrfile[100];

  // Open files for writing
  CkPrintf("Writing network data files %" PRIidx "\n", datidx);
  sprintf(csrfile, "%s/%s%s.coord.%" PRIidx "", filedir.c_str(), filebase.c_str(), filemod.c_str(), datidx);
  pCoord = fopen(csrfile,"w");
  sprintf(csrfile, "%s/%s%s.adjcy.%" PRIidx "", filedir.c_str(), filebase.c_str(), filemod.c_str(), datidx);
  pAdjcy = fopen(csrfile,"w");
  sprintf(csrfile, "%s/%s%s.state.%" PRIidx "", filedir.c_str(), filebase.c_str(), filemod.c_str(), datidx);
  pState = fopen(csrfile,"w");
  sprintf(csrfile, "%s/%s%s.event.%" PRIidx "", filedir.c_str(), filebase.c_str(), filemod.c_str(), datidx);
  pEvent = fopen(csrfile,"w");
  if (pCoord == NULL || pAdjcy == NULL || pState == NULL || pEvent == NULL) {
    CkPrintf("Error opening files for writing %" PRIidx "\n", datidx);
    CkExit();
  }

  // Initialize distributions
  netdist.resize(nprt);

  // Loop through parts
  for (idx_t k = 0; k < nprt; ++k) {
    // Set up netdist
    netdist[k].prtidx = xprt + k;
    netdist[k].nvtx = parts[k]->nvtx;
    netdist[k].nedg = parts[k]->nedg;
    netdist[k].nstate = parts[k]->nstate;
    netdist[k].nstick = parts[k]->nstick;
    netdist[k].nevent = parts[k]->nevent;

    // initialize counters
    idx_t jstate = 0;
    idx_t jstick = 0;

    // Graph adjacency information
    for (idx_t i = 0; i < parts[k]->nvtx; ++i) {
      // xyz
      // vertex coordinates
      fprintf(pCoord, " %" PRIrealfull " %" PRIrealfull " %" PRIrealfull "\n",
          parts[k]->xyz[i*3+0], parts[k]->xyz[i*3+1], parts[k]->xyz[i*3+2]);

      // vertex state
      CkAssert(parts[k]->vtxmodidx[i] > 0);
      fprintf(pState, " %s", modname[parts[k]->vtxmodidx[i]].c_str());
      for(idx_t s = 0; s < netmodel[parts[k]->vtxmodidx[i]]->getNState(); ++s) {
        fprintf(pState, " %" PRIrealfull "", parts[k]->state[jstate++]);
      }
      for(idx_t s = 0; s < netmodel[parts[k]->vtxmodidx[i]]->getNStick(); ++s) {
        fprintf(pState, " %" PRItickhex "", parts[k]->stick[jstick++]);
      }

      // Edges
      for (idx_t j = parts[k]->xadj[i]; j < parts[k]->xadj[i+1]; ++j) {
        // adjcy
        fprintf(pAdjcy, " %" PRIidx "", parts[k]->adjcy[j]);
        // edge state
        fprintf(pState, " %s", modname[parts[k]->edgmodidx[j]].c_str());
        for(idx_t s = 0; s < netmodel[parts[k]->edgmodidx[j]]->getNState(); ++s) {
          fprintf(pState, " %" PRIrealfull "", parts[k]->state[jstate++]);
        }
        for(idx_t s = 0; s < netmodel[parts[k]->edgmodidx[j]]->getNStick(); ++s) {
          fprintf(pState, " %" PRItickhex "", parts[k]->stick[jstick++]);
        }
      }

      // Events
      fprintf(pEvent, " %" PRIidx "", parts[k]->xevent[i+1] - parts[k]->xevent[i]);
      for (idx_t e = parts[k]->xevent[i]; e < parts[k]->xevent[i+1]; ++e) {
        // event
        if (parts[k]->type[e] == EVENT_SPIKE) {
          fprintf(pEvent, " %" PRItickhex " %" PRIidx " %" PRIidx " %" PRIidx "",
              parts[k]->diffuse[e], parts[k]->type[e], parts[k]->source[e], parts[k]->index[e]);
        }
        else {
          fprintf(pEvent, " %" PRItickhex " %" PRIidx " %" PRIidx " %" PRIidx " %" PRIrealfull "",
              parts[k]->diffuse[e], parts[k]->type[e], parts[k]->source[e], parts[k]->index[e], parts[k]->data[e]);
        }
      }

      // xadj
      fprintf(pAdjcy, "\n");
      fprintf(pState, "\n");
      fprintf(pEvent, "\n");
    }
    // Sanity check
    CkAssert(jstate == parts[k]->nstate);
    CkAssert(jstick == parts[k]->nstick);
  }

  // Cleanup
  fclose(pCoord);
  fclose(pAdjcy);
  fclose(pState);
  fclose(pEvent);
}


/**************************************************************************
* Record Information
**************************************************************************/

// Writing Records to file
//
void Netdata::WriteRecord() {
  /* File operations */
  FILE *pEvtlog;
  FILE *pRecord;
  char recfile[100];

  // Only save when data exists
  idx_t nevtlog = 0;
  idx_t nrecord = 0;
  for (idx_t k = 0; k < nprt; ++k) {
    nevtlog += records[k]->nevtlog;
    nrecord += records[k]->nrecord;
  }

  // Open File
  if (nevtlog) {
    sprintf(recfile, "%s/%s/%s.log.%" PRIidx ".%" PRIidx "", filedir.c_str(), recdir.c_str(), filebase.c_str(), datidx, records[0]->iter);
    pEvtlog = fopen(recfile,"w");
    if (pEvtlog == NULL) {
      CkPrintf("Error opening files for recording %" PRIidx "\n", datidx);
      CkExit();
    }
  }
  if (nrecord) {
    sprintf(recfile, "%s/%s/%s.record.%" PRIidx ".%" PRIidx "", filedir.c_str(), recdir.c_str(), filebase.c_str(), datidx, records[0]->iter);
    pRecord = fopen(recfile,"w");
    if (pRecord == NULL) {
      CkPrintf("Error opening files for recording %" PRIidx "\n", datidx);
      CkExit();
    }
  }

  // TODO: Store records indexed by time and then by type
  // Loop through parts
  for (idx_t k = 0; k < nprt; ++k) {
    // Loop through events
    for (idx_t e = 0; e < records[k]->nevtlog; ++e) {
      // event types lacking data
      if (records[k]->type[e] == EVENT_SPIKE) {
        fprintf(pEvtlog, "%" PRIidx " %" PRItickhex " %" PRIidx "\n",
            records[k]->type[e], records[k]->diffuse[e], records[k]->source[e]);
      }
      // event types with data
      else {
        fprintf(pEvtlog, "%" PRIidx " %" PRItickhex " %" PRIidx " %" PRIidx " %" PRIrealfull "\n",
            records[k]->type[e], records[k]->diffuse[e], records[k]->source[e], records[k]->index[e], records[k]->data[e]);
      }
    }
    // Loop through other records
    // TODO: add record type and modify accordingly
    for (idx_t r = 0; r < records[k]->nrecord; ++r) {
      // timestamp followed by number of data entries
      fprintf(pRecord, "%" PRIidx " %" PRItickhex " %" PRIidx " %" PRIidx " %" PRIidx "",
          records[k]->type[records[k]->nevtlog+r], records[k]->drift[r],
          (records[k]->xdata[r+1] - records[k]->xdata[r]),
          (records[k]->xdiffuse[r+1] - records[k]->xdiffuse[r]),
          (records[k]->xindex[r+1] - records[k]->xindex[r]));
      // real data
      for (idx_t s = records[k]->xdata[r]; s < records[k]->xdata[r+1]; ++s) {
        fprintf(pRecord, " %" PRIrealfull "", records[k]->data[s]);
      }
      // tick data
      for (idx_t s = records[k]->xdiffuse[r]; s < records[k]->xdiffuse[r+1]; ++s) {
        fprintf(pRecord, " %" PRItickhex "", records[k]->diffuse[s]);
      }
      // idx data
      for (idx_t s = records[k]->xindex[r]; s < records[k]->xindex[r+1]; ++s) {
        fprintf(pRecord, " %" PRIidx "", records[k]->index[s]);
      }
      // one line per record
      fprintf(pRecord, "\n");
    }
  }
  
  // Cleanup
  if (nevtlog) { fclose(pEvtlog); }
  if (nrecord) { fclose(pRecord); }
}


/**************************************************************************
* PNG Information
**************************************************************************/

// Writing PNGs to file
//
// TODO: Move this to Netdata?
void Network::WritePNG(idx_t pngidx) {
  /* File operations */
  FILE *pPNG;
  FILE *pMap;
  char pngfile[100];

  // Open File
  sprintf(pngfile, "%s/%s/%s.png.%" PRIidx "", filedir.c_str(), pngdir.c_str(), filebase.c_str(), vtxidx[pngidx]);
  pPNG = fopen(pngfile,"w");
  sprintf(pngfile, "%s/%s/%s.map.%" PRIidx "", filedir.c_str(), pngdir.c_str(), filebase.c_str(), vtxidx[pngidx]);
  pMap = fopen(pngfile,"w");
  if (pPNG == NULL || pMap == NULL) {
    CkPrintf("Error opening files for PNG output %" PRIidx "\n", vtxidx[pngidx]);
    CkExit();
  }

  // Loop through pngs
  CkAssert(pngmaps.size() == pngs[pngidx].size());
  for (std::size_t p = 0; p < pngmaps.size(); ++p) {
    // Loop through maps
    for (std::size_t s = 0; s < pngmaps[p].size(); ++s) {
      fprintf(pMap, "%" PRItickhex " %" PRIidx " %" PRIidx " %" PRItickhex " %" PRItickhex "\n",
          pngmaps[p][s].diffuse, pngmaps[p][s].source, pngmaps[p][s].origin, pngmaps[p][s].departure, pngmaps[p][s].arrival);
    }
    // empty line between maps
    fprintf(pMap, "\n");
    // Loop through stamps
    for (std::size_t s = 0; s < pngs[pngidx][p].size(); ++s) {
      fprintf(pPNG, " %" PRItickhex " %" PRIidx "",
          pngs[pngidx][p][s].diffuse, pngs[pngidx][p][s].source);
    }
    // newline between stamps
    fprintf(pPNG, "\n");
  }

  // Cleanup
  fclose(pPNG);
  fclose(pMap);
}


