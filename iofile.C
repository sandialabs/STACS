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


/**************************************************************************
* File Data
**************************************************************************/

// Load data from files into partitions
//
void Netdata::LoadFile() {
  // TODO: perhaps split which netdata are reading which files?
  // Read in data files
  for (std::size_t i = 0; i < datafiles.size(); ++i) {
    // Read in files
    if (datafiles[i].filetype == FT_CSV_SPARSE) {
      // read in data (as csr if sparse flag set)
      if (ReadFileCSVSparse(datafiles[i])) {
        CkPrintf("Error reading data file %s...\n", datafiles[i].filename.c_str());
        CkExit();
      }
    }
    else if (datafiles[i].filetype == FT_CSV_DENSE) {
      // read in data (as matrix)
      if (ReadFileCSV(datafiles[i])) {
        CkPrintf("Error reading data file %s...\n", datafiles[i].filename.c_str());
        CkExit();
      }
    }
  }
  
  // Return control to main
  contribute(0, NULL, CkReduction::nop);
}

// Coordination with NetData chare array
//
void Network::LoadFile() {
  // No need to load if no files
  nfile = cfile = 0;
  if (datafiles.size() == 0) {
    // Return control to main
    contribute(0, NULL, CkReduction::nop);
  }
  else {
    std::vector<idx_t> datamodidx;
    datamodidx.resize(datafiles.size());
    // Figure out which files and rows to request
    // This is done by walking through the different models and determining the target vertex model
    // For vertices, this is simply through the model index in modelconf where the init is through a file
    // For edges, we need to additionally search through the edges list to get the target model
    // TODO: is there a more optimal way to do this?
    for (std::size_t i = 0; i < modelconf.size(); ++i) {
      // Vertex state
      if (modelconf[i].graphtype == GRAPHTYPE_VTX) {
        for (std::size_t s = 0; s < modelconf[i].stateinit.size(); ++s) {
          if (modelconf[i].stateinit[s] == RNGTYPE_FILE) {
            for (std::size_t v = 0; v < vertices.size(); ++v) {
              if (vertices[v].modidx == modmap[modelconf[i].modname]) {
                datamodidx[(idx_t) modelconf[i].stateparam[s][0]] = v;
              }
            }
          }
        }
        for (std::size_t s = 0; s < modelconf[i].stickinit.size(); ++s) {
          if (modelconf[i].stickinit[s] == RNGTYPE_FILE) {
            for (std::size_t v = 0; v < vertices.size(); ++v) {
              if (vertices[v].modidx == modmap[modelconf[i].modname]) {
                datamodidx[(idx_t) modelconf[i].stickparam[s][0]] = v;
              }
            }
          }
        }
      }
      // Edge state
      else if (modelconf[i].graphtype == GRAPHTYPE_EDG) {
        for (std::size_t s = 0; s < modelconf[i].stateinit.size(); ++s) {
          if (modelconf[i].stateinit[s] == RNGTYPE_FILE) {
            for (std::size_t e = 0; e < edges.size(); ++e) {
              if (edges[e].modidx == modmap[modelconf[i].modname]) {
                for (std::size_t v = 0; v < vertices.size(); ++v) {
                  if (vertices[v].modidx == edges[e].target[0]) {
                    datamodidx[(idx_t) modelconf[i].stateparam[s][0]] = v;
                  }
                }
              }
            }
          }
        }
        for (std::size_t s = 0; s < modelconf[i].stickinit.size(); ++s) {
          if (modelconf[i].stickinit[s] == RNGTYPE_FILE) {
            for (std::size_t e = 0; e < edges.size(); ++e) {
              if (edges[e].modidx == modmap[modelconf[i].modname]) {
                for (std::size_t v = 0; v < vertices.size(); ++v) {
                  if (vertices[v].modidx == edges[e].target[0]) {
                    datamodidx[(idx_t) modelconf[i].stickparam[s][0]] = v;
                  }
                }
              }
            }
          }
        }
      }
    }
    // Edge connection
    for (std::size_t e = 0; e < edges.size(); ++e) {
      for (std::size_t k = 0; k < edges[e].conntype.size(); ++k) {
        if (edges[e].conntype[k] == CONNTYPE_FILE) {
          for (std::size_t v = 0; v < vertices.size(); ++v) {
            if (vertices[v].modidx == edges[e].target[0]) {
              datamodidx[(idx_t) edges[e].probparam[k][0]] = v;
            }
          }
        }
      }
    }

    // Request rows from files
    nfile = datafiles.size();
    for (std::size_t dfidx = 0; dfidx < datafiles.size(); ++dfidx) {
      //CkPrintf("  Datafile %zu (%s)is for target model %" PRIidx " (%s)\n",
      //    dfidx, datafiles[dfidx].filename.c_str(), datamodidx[dfidx],
      //    modelconf[datamodidx[dfidx]].modname.c_str());
      int xrow = xordervtx[datamodidx[dfidx]];
      int nrow = nordervtx[datamodidx[dfidx]];

      // Request matrix rows from datafiles
      netdata(datidx).LoadMatrix(dfidx, xrow, nrow,
          CkCallback(CkIndex_Network::LoadMatrix(NULL), thisProxy(prtidx)));
    }
  }
}

void Netdata::LoadMatrix(int dfidx, int xrow, int nrow, const CkCallback &cbpart) {
  // Compute number of columns
  int ncol = 0;
  for (int i = xrow; i < xrow+nrow; ++i) {
    ncol += datafiles[dfidx].matrix[i].size();
  }
  //CkPrintf("  Loading datafile %d (%d - %d) (%d) on %d\n",
  //    dfidx, xrow, xrow+nrow, ncol, datidx);
  
  // Initialize matrix message
  int msgSize[MSG_Matrix];
  msgSize[0] = nrow+1; // rows
  msgSize[1] = ncol; // cols
  msgSize[2] = ncol; // values
  mMatrix *mmatrix = new(msgSize, 0) mMatrix;
  mmatrix->dfidx = dfidx;
  mmatrix->xrow = xrow;
  mmatrix->nrow = nrow;
  
  // Build a compressed sparse row of selected rows
  mmatrix->rows[0] = 0;
  idx_t jcol = 0;
  for (int i = 0; i < nrow; ++i) {
    mmatrix->rows[i+1] = mmatrix->rows[i] + datafiles[dfidx].matrix[xrow+i].size();
    std::unordered_map<idx_t, real_t>::iterator icol;
    for (icol = datafiles[dfidx].matrix[xrow+i].begin();
        icol != datafiles[dfidx].matrix[xrow+i].end(); ++icol) {
      mmatrix->cols[jcol] = icol->first;
      mmatrix->values[jcol++] = icol->second;
    }
  }

  // Send part to network
  cbpart.send(mmatrix);
}

void Network::LoadMatrix(mMatrix *msg) {
  // Copy over rows and columns
  datafiles[msg->dfidx].xrow = msg->xrow;
  datafiles[msg->dfidx].matrix.clear();
  for (idx_t i = 0; i < msg->nrow; ++i) {
    std::unordered_map<idx_t, real_t> row;
    for (idx_t j = msg->rows[i]; j < msg->rows[i+1]; ++j) {
      row.emplace(msg->cols[j], msg->values[j]);
    }
    datafiles[msg->dfidx].matrix.push_back(row);
  }
  //CkPrintf("  Loaded datafile %d (%d - %d) on %d\n",
  //    msg->dfidx, msg->xrow, msg->xrow+datafiles[msg->dfidx].matrix.size(), prtidx);

  delete msg;

  // Return control to main when all files loaded
  if (++cfile == nfile) {
    contribute(0, NULL, CkReduction::nop);
  }
}


/**************************************************************************
* Network data files (dense)
**************************************************************************/

// Read data file (csv)
//
int Netdata::ReadFileCSV(datafile_t &datafile) {
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
* Network data files (sparse)
**************************************************************************/

// Read data file (csv)
//
int Netdata::ReadFileCSVSparse(datafile_t &datafile) {
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
