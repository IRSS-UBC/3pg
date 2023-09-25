// 3PG main program.  ANL 6/3/2000.  

/*
All source code remains the property and copyright of CSIRO. 

CSIRO accepts no responsibility for the use of 3PG(S) or of the model 3-PG in
the form supplied or as subsequently modified by third parties. CSIRO disclaims
liability for all losses, damages and costs incurred by any person as a result
of relying on this software. 
Use of this software assumes agreement to this condition of use
*/

static char rcsid[] = "$Id: 3pg.cpp,v 1.10 2001/08/02 06:34:10 lou026 Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gdal.h"
#include "gdal_priv.h"
#include "GDALRasterImage.hpp"
#include "Data_io.hpp"
#include "The_3PG_Model.hpp"

// Need to provide getopt on MSVC. 
#ifdef WIN32
extern "C"
{
  extern int getopt(int argc, char **argv, char *opts);
  extern char *optarg;
}
#endif

// Maximum file path length. 
#define MAXFILE 1000

FILE *logfp;
char usage[] = 
"-d <default parameter file> -s <site parameter file>\n";
char program[] = "3pg";

//----------------------------------------------------------------------------------------

void copyright(FILE *fp)
{ 
  char copymessage[]=
    "This version of 3-PG(S) has been developed by:\n"
    //"Nicholas Coops [Nicholas.Coops@csiro.au],\n"
    //"Anders Siggins [Anders.Siggins@csiro.au],\n"
    //"and Andrew Loughhead.\n"
    "CSIRO Forestry and Forest Products,\n"
    "Private Bag 10, Clayton South 3169, Victoria, Australia.\n"
    "For Aracruz.\n"
    "Contact (Aracruz): Auro Campi de Almeida [aca@aracruz.com.br]\n"
    "Contact (CSIRO)  : Anders.Siggins@csiro.au\n"
    "                 : Nicholas.Coops@csiro.au\n"
    "\n"
    "Please read Release Notes for additional infomation:\n"
    "\n"
    "Revision: 2.53 \n"
    "Date: 13 Dec 2002 \n"
    "Programmer:  Anders Siggins\n"
    "\n"
    "\"DISCLAIMER\" \n"
    "------------------------------------------\n"
    "CSIRO accepts no responsibility for the use of 3PG(S) or of the model 3-PG in\n"
    "the form supplied or as subsequently modified by third parties. CSIRO disclaims\n"
    "liability for all losses, damages and costs incurred by any person as a result\n"
    "of relying on this software. Use of this software assumes agreement to this\n"
    "condition of use.\n"
    "\n"
    "Removal of this statement violates the spirit in which 3PG(S) was released by\n"
    "CSIRO Forestry and Forest Products.\n";
  fprintf(fp, copymessage);
  return;
}

//----------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  std::string optarg;
  extern int optind;
  int c;
  int result;

  GDALRasterImage *refGrid; // Pointer variable refGrid pointing to GDALRasterImage 
  bool spatial=0;
  long nrows, ncols;
  MYDate spMinMY, spMaxMY; 
  std::string defParamFile = ""; 
  std::string siteParamFile = "";

  // Copyright message for stdout
  //copyright(stdout);
  
  /* Command line options */

  
  while (((c = getopt(argc, argv, "d:s:")) != EOF))
  {
    switch (c) {
    case 'd':
      defParamFile = optarg;
      break;
      
    case 's':
      siteParamFile = optarg;
      break;
      
    default:
      fprintf(stderr, "Usage: %s %s", program, usage);
      exit(0);
    }
  }

  // Usage
  if (defParamFile[0] == 0 || siteParamFile[0] == 0) {
    fprintf(stderr, "Usage: %s %s", program, usage);
    exit(0);
  }

  // Open the log file. 
  logfp = openLogFile(siteParamFile); 

  // Copyright message for log file
  //copyright(logfp);

  // result = _setmaxstdio(2048);
  if (result == -1)
  {
    fprintf(logfp, "Maximum files number could not be increased\n");
    fprintf(stdout, "Maximum files number could not be increased\n");
  }

  // Load the parameters and output variables. 
  InitInputParams();
  readParamFile(defParamFile);
  readParamFile(siteParamFile);
  if (!haveAllParams())
    exit(1);

  // Check for a spatial run, if so open input grids and define refGrid. 
  refGrid = openInputGrids();
  spatial = ( refGrid != NULL );

  // In point mode require the point mode output filename. 
  if (!spatial)
    if (!havePointOpFile())
      exit(1);

  if (spatial) {
    // Grid dimensions
    nrows = refGrid->nRows;
    ncols = refGrid->nCols;
  }
  else {
    nrows = 0;
    ncols = 0;
  }

  // Find the over all start year and end year. 
  findRunPeriod( refGrid, spMinMY, spMaxMY ); 
  fprintf( logfp, "first run mon/year = %2d/%4d, last run mon/year = %2d/%4d\n", spMinMY.mon, 
	  spMinMY.year, spMaxMY.mon, spMaxMY.year );

  fprintf( stdout, "Expected run period of simulation: %d - %d\n", spMinMY.year, spMaxMY.year );

  // NOTE: don't this this is necessary for GDAL stuff... but I guess we'll see
  // ResetGrids(); //findrun period has moved all grids to the end - move them back to the start...

  if (spatial) {
    // Open output grids. 
    openOutputGrids( refGrid );
    openRegularOutputGrids( refGrid, spMinMY, spMaxMY );

    // Open and read sample point file
    readSampleFile( refGrid ); 
    fprintf(logfp, "Processing %u cells...\n", nrows * ncols);
  }

 
//if (spatial) {
//    int cellsDone = 0; 
//    int cellsTotal = (nrows - 2) * (ncols - 2); 
//    int lastProgress = -1; 
//    fprintf(stdout, "                              \r");

// Run the model. 
  if (spatial) {

    int cellsDone = 0;
    int cellsTotal = (nrows) * (ncols); 
    int lastProgress = -1; 
    fprintf(stdout, "                              \r");
    for (int j = 0; j < cellsTotal ; j++)
    {
      //PrintGrids();
      if (j == 16080)
        bool test = true;
      int progress = ((100 * cellsDone) / cellsTotal);
      if (progress > lastProgress)
        fprintf(stdout, "Completed %2u%%\r", progress);
      runTreeModel( spMinMY, spMaxMY, spatial, j);
      cellsDone++;
      lastProgress = progress;
      
    //   fprintf(logfp, "   row %u\n", j); fflush(logfp);
    }

    fprintf(stdout, "Completed %2u%%\r", 100);

    // Write output grids. 
//    writeOutputGrids();
//  }

    // Write output grids. 
    writeOutputGrids();
  }
  else {
    runTreeModel(spMinMY, spMaxMY, spatial, 0);
  }

  CloseGrids(); //Close all files currently open...

  return 0;
}
