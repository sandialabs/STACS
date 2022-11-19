/**
 * Copyright (C) 2015 Felix Wang
 *
 * Simulation Tool for Asynchronous Cortical Streams (stacs)
 *
 * iocsr.C
 * Handles condensed sparse row format
 */

#include "stacs.h"
#include "network.h"

// Maximum size of input line (bytes)
#define MAXLINE 2560000

/**************************************************************************
* Charm++ Read-Only Variables
**************************************************************************/
extern /*readonly*/ std::string netwkdir;
extern /*readonly*/ int netparts;
extern /*readonly*/ int netfiles;
extern /*readonly*/ std::string filebase;
extern /*readonly*/ std::string fileload;
extern /*readonly*/ std::string filesave;
extern /*readonly*/ std::string recordir;
extern /*readonly*/ std::string groupdir;


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
  CkPrintf("Reading network distribution\n");//from %s/%s.dist\n", netwkdir.c_str(), filebase.c_str());
  sprintf(csrfile, "%s/%s.dist", netwkdir.c_str(), filebase.c_str());
  pDist = fopen(csrfile,"r");
  if (pDist == NULL || line == NULL) {
    return 1;
  }

  netdist.resize(netparts+1);

  // Get distribution info
  for (int i = 0; i < netparts+1; ++i) {
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
int Main::WriteDist(int checkflag) {
  /* File operations */
  FILE *pDist;
  FILE *pMetis;
  char csrfile[100];
  /* Bookkeeping */
  idx_t nvtx;
  idx_t nedg;
  idx_t nstate;
  idx_t nstick;
  idx_t nevent;

  std::string filecheck = (checkflag ? filesave : "");
  // Open File
  CkPrintf("Writing network distribution\n");
  sprintf(csrfile, "%s/%s%s.dist", netwkdir.c_str(), filebase.c_str(), filecheck.c_str());
  pDist = fopen(csrfile,"w");
  sprintf(csrfile, "%s/%s%s.metis", netwkdir.c_str(), filebase.c_str(), filecheck.c_str());
  pMetis = fopen(csrfile,"w");
  if (pDist == NULL || pMetis == NULL) {
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
  
  // For Metis
  nvtx = nedg = 0;
  fprintf(pMetis, "%" PRIidx " %" PRIidx "\n", nvtx, nedg);
  for (int datidx = 0; datidx < netfiles; ++datidx) {
    idx_t ndiv = netparts/netfiles;
    idx_t nrem = netparts%netfiles;
    idx_t nprt = ndiv + (datidx < nrem);
    idx_t xprt = datidx*ndiv + (datidx < nrem ? datidx : nrem);
    for (idx_t jprt = 0; jprt < nprt; ++jprt) {
      nvtx += netdist[xprt+jprt].nvtx;
      nedg += netdist[xprt+jprt].nedg;
    }
    fprintf(pMetis, "%" PRIidx " %" PRIidx "\n", nvtx, nedg);
  }

  // Cleanup
  fclose(pDist);
  fclose(pMetis);

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
  sprintf(csrfile, "%s/%s%s.coord.%d", netwkdir.c_str(), filebase.c_str(), fileload.c_str(), datidx);
  pCoord = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s%s.adjcy.%d", netwkdir.c_str(), filebase.c_str(), fileload.c_str(), datidx);
  pAdjcy = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s%s.state.%d", netwkdir.c_str(), filebase.c_str(), fileload.c_str(), datidx);
  pState = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s%s.event.%d", netwkdir.c_str(), filebase.c_str(), fileload.c_str(), datidx);
  pEvent = fopen(csrfile,"r");
  if (pCoord == NULL || pAdjcy == NULL || pState == NULL ||
      pEvent == NULL || line == NULL) {
    CkPrintf("Error opening graph files on %d\n", datidx);
    CkExit();
  }

  // Read in Parts
  for (int k = 0; k < nprt; ++k) {
    idx_t nvtx = vtxdist[xprt+k+1] - vtxdist[xprt+k];
    idx_t nedg = edgdist[xprt+k+1] - edgdist[xprt+k];
    idx_t nstate = statedist[xprt+k+1] - statedist[xprt+k];
    idx_t nstick = stickdist[xprt+k+1] - stickdist[xprt+k];
    idx_t nevent = eventdist[xprt+k+1] - eventdist[xprt+k];

    // Initialize partition data message
    int msgSize[MSG_Part];
    msgSize[0] = netparts+1;    // vtxdist
    msgSize[1] = 0;             // vtxidx (implicit)
    msgSize[2] = nvtx;          // vtxmodidx
    msgSize[3] = nvtx*3;        // xyz
    msgSize[4] = nvtx+1;        // xadj
    msgSize[5] = nedg;          // adjcy
    msgSize[6] = nedg;          // edgmodidx
    msgSize[7] = nstate;        // state
    msgSize[8] = nstick;        // stick
    msgSize[9] = nvtx+1;        // xevent
    msgSize[10] = nevent;       // diffuse
    msgSize[11] = nevent;       // type
    msgSize[12] = nevent;       // source
    msgSize[13] = nevent;       // index
    msgSize[14] = nevent;       // data
    parts[k] = new(msgSize, 0) mPart;
    
    // Data sizes
    parts[k]->nvtx = nvtx;
    parts[k]->nedg = nedg;
    parts[k]->nstate = nstate;
    parts[k]->nstick = nstick;
    parts[k]->nevent = nevent;
    parts[k]->prtidx = xprt+k;

    // vtxdist
    for (int i = 0; i < netparts+1; ++i) {
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
      for(idx_t s = 0; s < model[modidx]->getNState(); ++s) {
        real_t state = strtoreal(oldstr, &newstr);
        oldstr = newstr;
        // state
        parts[k]->state[jstate++] = state;
      }
      for(idx_t s = 0; s < model[modidx]->getNStick(); ++s) {
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
        for(idx_t s = 0; s < model[modidx]->getNState(); ++s) {
          real_t state = strtoreal(oldstr, &newstr);
          oldstr = newstr;
          // state
          parts[k]->state[jstate++] = state;
        }
        for(idx_t s = 0; s < model[modidx]->getNStick(); ++s) {
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

// Read network parts (reordering)
//
void Netdata::ReadPart() {
  /* Bookkeeping */
  idx_t nsizedat;
  idx_t nstatedat;
  idx_t nstickdat;
  idx_t neventdat;
  /* File operations */
  FILE *pPart;
  FILE *pCoord;
  FILE *pAdjcy;
  FILE *pState;
  FILE *pEvent;
  char csrfile[100];
  char *line;
  char *oldstr, *newstr;

  // Prepare buffer
  line = new char[MAXLINE];

  // Open files for reading
  sprintf(csrfile, "%s/%s.part.%d", netwkdir.c_str(), filebase.c_str(), datidx);
  pPart = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s.coord.%d", netwkdir.c_str(), filebase.c_str(), datidx);
  pCoord = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s.adjcy.%d", netwkdir.c_str(), filebase.c_str(), datidx);
  pAdjcy = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s.state.%d", netwkdir.c_str(), filebase.c_str(), datidx);
  pState = fopen(csrfile,"r");
  sprintf(csrfile, "%s/%s.event.%d", netwkdir.c_str(), filebase.c_str(), datidx);
  pEvent = fopen(csrfile,"r");
  if (pPart == NULL || pCoord == NULL || pAdjcy == NULL ||
      pState == NULL || pEvent == NULL) {
    CkPrintf("Error opening files for reading\n");
    CkExit();
  }
  if (line == NULL) {
    CkPrintf("Could not allocate memory for lines\n");
    CkExit();
  }

  // Initialize sizes
  std::vector<idx_t> reprt(vtxmetis[datidx+1] - vtxmetis[datidx]);
  nsizedat = 0;
  nstatedat = 0;
  nstickdat = 0;
  neventdat = 0;

  // Read in graph information
  for (std::size_t i = 0; i < reprt.size(); ++i) {
    while(fgets(line, MAXLINE, pPart) && line[0] == '%');
    oldstr = line;
    newstr = NULL;
    // reprt
    reprt[i] = strtoidx(oldstr, &newstr, 10);
    CkAssert(reprt[i] < netparts);
    vtxidxreprt[reprt[i]].push_back(vtxmetis[datidx]+i);
    adjcyreprt[reprt[i]].push_back(std::vector<idx_t>());

    // Read in line (coordinates)
    while(fgets(line, MAXLINE, pCoord) && line[0] == '%');
    oldstr = line;
    newstr = NULL;
    // xyz
    for (idx_t j = 0; j < 3; ++j) {
      real_t coord = strtoreal(oldstr, &newstr);
      oldstr = newstr;
      xyzreprt[reprt[i]].push_back(coord);
    }

    // Read line (per vertex)
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
      adjcyreprt[reprt[i]].back().push_back(edg);
      ++nsizedat;
    }

    // Extract State Information
    // Read line (vertex)
    while(fgets(line, MAXLINE, pState) && line[0] == '%');
    oldstr = line;
    newstr = NULL;
    // extract model index from name
    idx_t modidx = strtomodidx(oldstr, &newstr);
    CkAssert(modidx != IDX_T_MAX);
    oldstr = newstr;
    // vtxmodidx
    vtxmodidxreprt[reprt[i]].push_back(modidx);
    statereprt[reprt[i]].push_back(std::vector<real_t>());
    stickreprt[reprt[i]].push_back(std::vector<tick_t>());
    CkAssert(modidx > 0);
    for(std::size_t s = 0; s < modelconf[modidx-1].stateinit.size(); ++s) {
      real_t stt = strtoreal(oldstr, &newstr);
      oldstr = newstr;
      // state
      statereprt[reprt[i]].back().push_back(stt);
      ++nstatedat;
    }
    for(std::size_t s = 0; s < modelconf[modidx-1].stickinit.size(); ++s) {
      tick_t stt = strtotick(oldstr, &newstr, 16);
      oldstr = newstr;
      // state
      stickreprt[reprt[i]].back().push_back(stt);
      ++nstickdat;
    }

    // edgmodidx
    modidx = strtomodidx(oldstr, &newstr);
    oldstr = newstr;
    edgmodidxreprt[reprt[i]].push_back(std::vector<idx_t>());
    while (modidx != IDX_T_MAX) {
      // modidx
      edgmodidxreprt[reprt[i]].back().push_back(modidx);
      statereprt[reprt[i]].push_back(std::vector<real_t>());
      stickreprt[reprt[i]].push_back(std::vector<tick_t>());
      // only push edge state if model and not 'none'
      if (modidx > 0) {
        for(std::size_t s = 0; s < modelconf[modidx-1].stateinit.size(); ++s) {
          real_t stt = strtoreal(oldstr, &newstr);
          oldstr = newstr;
          // state
          statereprt[reprt[i]].back().push_back(stt);
          ++nstatedat;
        }
        for(std::size_t s = 0; s < modelconf[modidx-1].stickinit.size(); ++s) {
          tick_t stt = strtotick(oldstr, &newstr, 16);
          oldstr = newstr;
          // state
          stickreprt[reprt[i]].back().push_back(stt);
          ++nstickdat;
        }
      }
      modidx = strtomodidx(oldstr, &newstr);
      oldstr = newstr;
    }

    // Extract event information
    // Read line (per vertex)
    while(fgets(line, MAXLINE, pEvent) && line[0] == '%');
    oldstr = line;
    newstr = NULL;
    // number of events
    idx_t jevent = strtoidx(oldstr, &newstr, 10);
    oldstr = newstr;
    eventreprt[reprt[i]].push_back(std::vector<event_t>());
    event_t eventpre;
    for (idx_t j = 0; j < jevent; ++j) {
      // diffuse
      eventpre.diffuse = strtotick(oldstr, &newstr, 16);
      oldstr = newstr;
      // type
      idx_t type = strtoidx(oldstr, &newstr, 10);
      oldstr = newstr;
      eventpre.type = type;
      // source
      eventpre.source = strtoidx(oldstr, &newstr, 10);
      oldstr = newstr;
      // index
      eventpre.index = strtoidx(oldstr, &newstr, 10);
      oldstr = newstr;
      // data
      if (type == EVENT_SPIKE) {
        eventpre.data = 0.0;
      }
      else {
        eventpre.data = strtoreal(oldstr, &newstr);
        oldstr = newstr;
      }
      eventreprt[reprt[i]].back().push_back(eventpre);
    }
    neventdat += jevent;
    CkAssert(eventreprt[reprt[i]].back().size() == (std::size_t) jevent);
  }

  // Cleanup
  fclose(pPart);
  fclose(pCoord);
  fclose(pAdjcy);
  fclose(pState);
  fclose(pEvent);
  delete[] line;

  // Print out some information
  CkPrintf("  File: %d   Vertices: %zu   Edges: %" PRIidx "   States: %" PRIidx "   Sticks: %" PRIidx"   Events: %" PRIidx"\n",
      datidx, reprt.size(), nsizedat, nstatedat, nstickdat, neventdat);

  // Prepare for reording partitions
  cpdat = 0;
  cpprt = 0;
  norderdat = 0;
  vtxprted.resize(nprt);
  xyzprted.resize(nprt);
  adjcyprted.resize(nprt);
  edgmodidxprted.resize(nprt);
  stateprted.resize(nprt);
  stickprted.resize(nprt);
  eventprted.resize(nprt);
  adjcyreord.resize(nprt);
  edgmodidxreord.resize(nprt);
  statereord.resize(nprt);
  stickreord.resize(nprt);
  norderprt.resize(nprt);
  for (idx_t i = 0; i < nprt; ++i) {
    norderprt[i] = 0;
  }
  vtxdist.resize(netfiles+1);
  vtxdist[0] = 0;
}
    
// Write out network files
//
void Netdata::WriteNetwork(int checkflag) {
  /* File operations */
  FILE *pCoord;
  FILE *pAdjcy;
  FILE *pState;
  FILE *pEvent;
  char csrfile[100];

  std::string filecheck = (checkflag ? filesave : "");
  // Open files for writing
  CkPrintf("Writing network data files %d\n", datidx);
  sprintf(csrfile, "%s/%s%s.coord.%d", netwkdir.c_str(), filebase.c_str(), filecheck.c_str(), datidx);
  pCoord = fopen(csrfile,"w");
  sprintf(csrfile, "%s/%s%s.adjcy.%d", netwkdir.c_str(), filebase.c_str(), filecheck.c_str(), datidx);
  pAdjcy = fopen(csrfile,"w");
  sprintf(csrfile, "%s/%s%s.state.%d", netwkdir.c_str(), filebase.c_str(), filecheck.c_str(), datidx);
  pState = fopen(csrfile,"w");
  sprintf(csrfile, "%s/%s%s.event.%d", netwkdir.c_str(), filebase.c_str(), filecheck.c_str(), datidx);
  pEvent = fopen(csrfile,"w");
  if (pCoord == NULL || pAdjcy == NULL || pState == NULL || pEvent == NULL) {
    CkPrintf("Error opening files for writing %d\n", datidx);
    CkExit();
  }

  // Initialize distributions
  netdist.resize(nprt);

  // Loop through parts
  for (int k = 0; k < nprt; ++k) {
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
      for(idx_t s = 0; s < model[parts[k]->vtxmodidx[i]]->getNState(); ++s) {
        fprintf(pState, " %" PRIrealfull "", parts[k]->state[jstate++]);
      }
      for(idx_t s = 0; s < model[parts[k]->vtxmodidx[i]]->getNStick(); ++s) {
        fprintf(pState, " %" PRItickhex "", parts[k]->stick[jstick++]);
      }

      // Edges
      for (idx_t j = parts[k]->xadj[i]; j < parts[k]->xadj[i+1]; ++j) {
        // adjcy
        fprintf(pAdjcy, " %" PRIidx "", parts[k]->adjcy[j]);
        // edge state
        fprintf(pState, " %s", modname[parts[k]->edgmodidx[j]].c_str());
        for(idx_t s = 0; s < model[parts[k]->edgmodidx[j]]->getNState(); ++s) {
          fprintf(pState, " %" PRIrealfull "", parts[k]->state[jstate++]);
        }
        for(idx_t s = 0; s < model[parts[k]->edgmodidx[j]]->getNStick(); ++s) {
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
  for (int k = 0; k < nprt; ++k) {
    nevtlog += records[k]->nevtlog;
    nrecord += records[k]->nrecord;
  }

  // Open File
  if (nevtlog) {
    sprintf(recfile, "%s/%s/%s.evtlog.%" PRIidx ".%d", netwkdir.c_str(), recordir.c_str(), filebase.c_str(), records[0]->iter, datidx);
    pEvtlog = fopen(recfile,"w");
    if (pEvtlog == NULL) {
      CkPrintf("Error opening files for recording %d\n", datidx);
      CkExit();
    }
  }
  if (nrecord) {
    sprintf(recfile, "%s/%s/%s.record.%" PRIidx ".%d", netwkdir.c_str(), recordir.c_str(), filebase.c_str(), records[0]->iter, datidx);
    pRecord = fopen(recfile,"w");
    if (pRecord == NULL) {
      CkPrintf("Error opening files for recording %d\n", datidx);
      CkExit();
    }
  }

  // TODO: Store records indexed by time and then by type
  // Loop through parts
  for (int k = 0; k < nprt; ++k) {
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
      // real data (states)
      for (idx_t s = records[k]->xdata[r]; s < records[k]->xdata[r+1]; ++s) {
        fprintf(pRecord, " %" PRIrealfull "", records[k]->data[s]);
      }
      // tick data (sticks)
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
* Group Information
**************************************************************************/

// Writing Groups to file
//
// TODO: Move this to Netdata?
void Network::WriteGroup(idx_t groupidx) {
  /* File operations */
  FILE *pStamp;
  FILE *pRoute;
  char grpfile[100];

  // Open File
  sprintf(grpfile, "%s/%s/%s.stamp.%" PRIidx "", netwkdir.c_str(), groupdir.c_str(), filebase.c_str(), vtxidx[groupidx]);
  pStamp = fopen(grpfile,"w");
  sprintf(grpfile, "%s/%s/%s.route.%" PRIidx "", netwkdir.c_str(), groupdir.c_str(), filebase.c_str(), vtxidx[groupidx]);
  pRoute = fopen(grpfile,"w");
  if (pStamp == NULL || pRoute == NULL) {
    CkPrintf("Error opening files for group output %" PRIidx "\n", vtxidx[groupidx]);
    CkExit();
  }

  // Loop through groups
  CkAssert(grproutes.size() == grpstamps[groupidx].size());
  for (std::size_t p = 0; p < grproutes.size(); ++p) {
    // Loop through maps
    for (std::size_t s = 0; s < grproutes[p].size(); ++s) {
      fprintf(pRoute, "%" PRItickhex " %" PRIidx " %" PRIidx " %" PRItickhex " %" PRItickhex "\n",
          grproutes[p][s].diffuse, grproutes[p][s].source, grproutes[p][s].origin, grproutes[p][s].departure, grproutes[p][s].arrival);
    }
    // empty line between maps
    fprintf(pRoute, "\n");
    // Loop through stamps
    for (std::size_t s = 0; s < grpstamps[groupidx][p].size(); ++s) {
      fprintf(pStamp, " %" PRItickhex " %" PRIidx "",
          grpstamps[groupidx][p][s].diffuse, grpstamps[groupidx][p][s].source);
    }
    // newline between stamps
    fprintf(pStamp, "\n");
  }

  // Cleanup
  fclose(pStamp);
  fclose(pRoute);
}

// Reading Groups from file
//
void Network::ReadGroup(idx_t groupidx) {
  /* File operations */
  FILE *pStamp;
  char grpfile[100];
  char *line;
  char *oldstr, *newstr;

  // Prepare buffer
  line = new char[MAXLINE];
  
  sprintf(grpfile, "%s/%s/%s.stamp.%" PRIidx "", netwkdir.c_str(), groupdir.c_str(), filebase.c_str(), vtxidx[groupidx]);
  pStamp = fopen(grpfile,"r");
  if (line == NULL || pStamp == NULL) {
    //CkPrintf("Warning: Group file does not exist %" PRIidx "\n", vtxidx[groupidx]);
    return;
  }

  // Each line is a Group
  for (idx_t p = 0;; ++p) {
    // Read line (stamps)
    while(fgets(line, MAXLINE, pStamp) && line[0] == '%');
    if (feof(pStamp)) { break; }
    oldstr = line;
    newstr = NULL;
    std::vector<stamp_t> stamps;
    std::set<idx_t> grpsource;
    for(;;) {
      stamp_t stamp;
      // diffuse
      stamp.diffuse = strtotick(oldstr, &newstr, 16);
      // check for end of line
      if (stamp.diffuse == 0 && oldstr != line)
        break;
      oldstr = newstr;
      // source
      stamp.source = strtoidx(oldstr, &newstr, 10);
      oldstr = newstr;
      // Add to group
      stamps.push_back(stamp);
      // Add to source set
      grpsource.insert(stamp.source);
    }
    // Add to groups
    grpstamps[groupidx].push_back(stamps);
    // Add to groups duration
    grpdur[groupidx].push_back(stamps.back().diffuse + 10*TICKS_PER_MS);
    // Add to map
    for (std::set<idx_t>::iterator source = grpsource.begin(); source != grpsource.end(); ++source) {
      grpmap[(*source)].push_back(std::array<idx_t, 2>{{groupidx, p}});
    }
  }
  grpwindow[groupidx].resize(grpstamps[groupidx].size());

  // Cleanup
  fclose(pStamp);
  delete[] line;
}


/**************************************************************************
* Estimation information
**************************************************************************/

// Writing Records to file
//
void Netdata::WriteEstimate() {
  /* File operations */
  FILE *pGrplog;
  char recfile[100];

  // Only save when data exists
  if (grplog.size() > 1) {
    // Open File
    sprintf(recfile, "%s/%s/%s.grplog.%" PRIidx "", netwkdir.c_str(), recordir.c_str(), filebase.c_str(), grplog[0].index);
    pGrplog = fopen(recfile,"w");
    if (pGrplog == NULL) {
      CkPrintf("Error opening files for recording\n");
      CkExit();
    }

    // Loop through events
    for (std::size_t e = 1; e < grplog.size(); ++e) {
      fprintf(pGrplog, "%" PRIidx " %" PRItickhex " %" PRIidx " %" PRIidx " %" PRIrealsec "\n",
          grplog[e].type, grplog[e].diffuse, grplog[e].source, grplog[e].index, grplog[e].data);
    }

    // Cleanup
    fclose(pGrplog);
  }
}

/**************************************************************************
* GeNet (network data files) (dense)
**************************************************************************/


// Read data file (csv)
//
int Netdata::ReadDataCSV(datafile_t &datafile) {
  FILE *pData;
  char csvfile[100];
  char *line;
  char *oldstr, *newstr;

  // Prepare buffer
  line = new char[MAXLINE];

  // Open files for reading
  // TODO: single-node file reads instead of per-process
  //       integrate this with MPI-IO?
  sprintf(csvfile, "%s/%s", netwkdir.c_str(), datafile.filename.c_str());
  pData = fopen(csvfile,"r");
  if (pData == NULL || line == NULL) {
    CkPrintf("Error opening file for reading\n");
    return 1;
  }

  // Initialize matrix
  datafile.matrix.clear();

  // Read csv into matrix
  // Dimensions are stored: targetdim x sourcedim
  // TODO: transpose the input file when reading?
  //       storage in csr-target-major order makes a
  //       single-threaded read-distribute more practical
  for (idx_t j = 0;; ++j) {
    // read in row
    while(fgets(line, MAXLINE, pData) && line[0] == '%');
    if (feof(pData)) { break; }
    oldstr = line;
    newstr = NULL;
    std::unordered_map<idx_t, real_t> row;
    // read in columns (comma delimited)
    idx_t i = 0;
    for (;;) {
      // check for empty element at beginning of file
      // TODO: is this robust enough?
      while (isspace(oldstr[0])) { ++oldstr; }
      while (oldstr[0] == ',') { ++oldstr; ++i; }
      // check for end of line (added by fgets)
      if (oldstr[0] == '\0') { break; }
      // element
      real_t element;
      element = strtoreal(oldstr, &newstr);
      oldstr = newstr;
      // Add element to row
      row.emplace(i, element);
      //CkPrintf("  %" PRIidx ", %" PRIidx ": %" PRIreal "\n", i, j, element);
      // check for empty element (again)
      while (isspace(oldstr[0])) { ++oldstr; }
      while (oldstr[0] == ',') { ++oldstr; ++i; }
      // check for end of line (added by fgets)
      if (oldstr[0] == '\0') { break; }
    }
    // Add to matrix
    datafile.matrix.push_back(row);
  }

  // Cleanup
  fclose(pData);
  delete[] line;

  return 0;
}

/**************************************************************************
* GeNet (network data files) (sparse)
**************************************************************************/


// Read data file (csv)
//
int Netdata::ReadDataCSVSparse(datafile_t &datafile) {
  FILE *pData;
  char csvfile[100];
  char *line;
  char *oldstr, *newstr;

  // Prepare buffer
  line = new char[MAXLINE];

  // Open files for reading
  // TODO: single-node file reads instead of per-process
  //       integrate this with MPI-IO?
  sprintf(csvfile, "%s/%s", netwkdir.c_str(), datafile.filename.c_str());
  pData = fopen(csvfile,"r");
  if (pData == NULL || line == NULL) {
    CkPrintf("Error opening file for reading\n");
    return 1;
  }

  // Initialize matrix
  datafile.matrix.clear();

  // Read csv into matrix
  // Dimensions are stored: targetdim x condensed sparse rows of source
  //                        as sourceidx:datavalue (data is optional)
  // TODO: transpose the input file when reading?
  //       storage in csr-target-major order makes a
  //       single-threaded read-distribute more practical
  for (idx_t j = 0;; ++j) {
    // read in row
    while(fgets(line, MAXLINE, pData) && line[0] == '%');
    if (feof(pData)) { break; }
    oldstr = line;
    newstr = NULL;
    std::unordered_map<idx_t, real_t> row;
    for (;;) {
      // check for empty element at beginning of file
      // TODO: is this robust enough?
      while (isspace(oldstr[0])) { ++oldstr; }
      while (oldstr[0] == ',') { ++oldstr; }
      // check for end of line (added by fgets)
      if (oldstr[0] == '\0') { break; }
      // source index
      idx_t sourceidx;
      sourceidx = strtoidx(oldstr, &newstr, 10);
      oldstr = newstr;
      // Skip the colon
      while (isspace(oldstr[0])) { ++oldstr; }
      while (oldstr[0] == ':') { ++oldstr; }
      while (isspace(oldstr[0])) { ++oldstr; }
      // element
      real_t element;
      // also handle no element case
      if (oldstr[0] == ',' || oldstr[0] == '\0') {
        element = 0.0;
      }
      else {
        element = strtoreal(oldstr, &newstr);
        oldstr = newstr;
      }
      // Add element to row
      row.emplace(sourceidx, element);
      //CkPrintf("  %" PRIidx ", %" PRIidx ": %" PRIreal "\n", sourceidx, j, element);
      // check for empty element (again)
      while (isspace(oldstr[0])) { ++oldstr; }
      while (oldstr[0] == ',') { ++oldstr; }
      // check for end of line (added by fgets)
      if (oldstr[0] == '\0') { break; }
    }
    // Add to matrix
    datafile.matrix.push_back(row);
  }

  // Cleanup
  fclose(pData);
  delete[] line;

  return 0;
}
