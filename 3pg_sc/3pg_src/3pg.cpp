// 3PG main program.  ANL 6/3/2000.  

/*
All source code remains the property and copyright of CSIRO. 

CSIRO accepts no responsibility for the use of 3PG(S) or of the model 3-PG in
the form supplied or as subsequently modified by third parties. CSIRO disclaims
liability for all losses, damages and costs incurred by any person as a result
of relying on this software. 
Use of this software assumes agreement to this condition of use
*/

// static char rcsid[] = "$Id: 3pg.cpp,v 1.10 2001/08/02 06:34:10 lou026 Exp $";

#include <cstdlib>
#include <cstring>
#include <iostream>
#include "gdal.h"
#include "gdal_priv.h"
#include "GDALRasterImage.hpp"
#include "Data_io.hpp"
#include "The_3PG_Model.hpp"
#include <boost/program_options.hpp>

// Need to provide getopt on MSVC. 
//#ifdef WIN32
//extern "C"
//{
//  extern int getopt(int argc, char **argv, char *opts);
//  extern char *optarg;
//}
//#endif

// Maximum file path length. 
#define MAXFILE 1000

// FILE *logfp;
// char usage[] = 
// "-d <default parameter file> -s <site parameter file>\n";
// char program[] = "3pg";

//----------------------------------------------------------------------------------------

void copyright()
{ 
    std::string copymessage =
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
  std::cout << copymessage << std::endl;
  // fprintf(fp, copymessage);
  return;
}

//function to parse -d and -s options from command line using boost::program_options
void parseCommandLine(int argc, char *argv[], std::string& defParamFile, std::string& siteParamFile)
{
  namespace po = boost::program_options;
  std::vector<std::string> inputFiles[2];
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("d", po::value<std::string>(&defParamFile), "default parameter file")
    ("s", po::value<std::string>(&siteParamFile), "site parameter file")
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(EXIT_SUCCESS);
  }
  if (vm.count("d") && vm.count("s")) {
    std::cout << "Default parameter file: " << defParamFile << std::endl;
    std::cout << "Site parameter file: " << siteParamFile << std::endl;
  } else {
    std::cout << "Failed to parse parameter and site parameter file" << std::endl;
    exit(EXIT_FAILURE);
  }
}
//----------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  std::string optarg;
  //extern int optind;
  //int c;
  //int result;

  GDALRasterImage *refGrid; // Pointer variable refGrid pointing to GDALRasterImage 
  bool spatial=0;
  long nrows, ncols;
  MYDate spMinMY, spMaxMY; 
  std::string defParamFile = "Y:\\Francois\\_Vaughan\\species_lpp.txt"; 
  std::string siteParamFile = "Y:\\Francois\\_Vaughan\\site_subarea.txt";

  // Copyright message for stdout
  //copyright();
  
  /* Command line options */
  //parseCommandLine(argc, argv, defParamFile, siteParamFile);
  // while (((c = getopt(argc, argv, "d:s:")) != EOF))
  // {
  //   switch (c) {
  //   case 'd':
  //     defParamFile = optarg;
  //     break;
      
  //   case 's':
  //     siteParamFile = optarg;
  //     break;
      
  //   default:
  //     std::cout << "Usage: " << program << " " << usage << std::endl;
  //     // fprintf(stderr, "Usage: %s %s", program, usage);
  //     exit(EXIT_FAILURE);
  //   }
  // }

  // Check for params files using isspace.
  if (defParamFile.empty() || siteParamFile.empty()) {

      std::cout << "Usage: 3pg <default parameter file> -s <site param%eter file>" << std::endl;
      // fprintf(stderr, "Usage: %s %s", program, usage);
      exit(EXIT_FAILURE);
  }
  // Open the log file. 
  // logfp = openLogFile(siteParamFile); 

  // Copyright message for log file
  //copyright(logfp);

  // Increase the number of files that can be open at once (Windows only)
  // result = _setmaxstdio(2048);
  // if (result == -1)
  // {
  //   fprintf(logfp, "Maximum files number could not be increased\n");
  //   fprintf(stdout, "Maximum files number could not be increased\n");
  // }

  // Load the parameters and output variables. 
  InitInputParams();
  readParamFile(defParamFile);
  readParamFile(siteParamFile);
  if (!haveAllParams())
    exit(EXIT_FAILURE);

  // Check for a spatial run, if so open input grids and define refGrid. 
  refGrid = openInputGrids();
  spatial = ( refGrid != NULL );

  // In point mode require the point mode output filename. 
  if (!spatial)
    if (!havePointOpFile())
      exit(EXIT_FAILURE);

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
  // TODO: findRunPeriod reads the entire input grid, which is unnecessary. Find some modern way to do this.
  std::cout << "Finding run period..." << std::endl;
  findRunPeriod( refGrid, spMinMY, spMaxMY ); 
  // fprintf( logfp, "first run mon/year = %2d/%4d, last run mon/year = %2d/%4d\n", spMinMY.mon, 
	//   spMinMY.year, spMaxMY.mon, spMaxMY.year );
  // fprintf( stdout, "Expected run period of simulation: %d - %d\n", spMinMY.year, spMaxMY.year );
  std::cout << "first run mon/year = " << spMinMY.mon << "/" << spMinMY.year << ", last run mon/year = " << spMaxMY.mon << "/" << spMaxMY.year << std::endl;
  std::cout << "Expected run period of simulation: " << spMinMY.year << " - " << spMaxMY.year << std::endl;

  // NOTE: don't think ResetGrids is necessary for GDAL stuff... but I guess we'll see
  // ResetGrids(); 

  if (spatial) {
    // Open output grids. 
    std::cout << "Opening output grids..." << std::endl;
    openOutputGrids( refGrid );
    std::cout << "Output grids opened." << std::endl;
    std::cout << "Opening regular output grids..." << std::endl;
    openRegularOutputGrids( refGrid, spMinMY, spMaxMY );

    std::cout << "Regular output grids opened." << std::endl;
    std::cout << "Reading points from sample file..." << std::endl;
    // Open and read sample point file
    readSampleFile( refGrid ); 
    std::cout << "Points read from sample file." << std::endl;
    // fprintf(logfp, "Processing %u cells...\n", nrows * ncols);
  }

 
  if (spatial) {
      int cellsDone = 0;
      int cellsTotal = (nrows - 2) * (ncols - 2);
      int lastProgress = -1;
      //fprintf(stdout, "                              \r");
  }
// Run the model. 
  if (spatial) {
     std::cout << "Running model..." << std::endl;
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

//     fprintf(stdout, "Completed %2u%%\r", 100);

//     // Write output grids. 
    //writeOutputGrids();
  }

//     // Write output grids. 
//     writeOutputGrids();
//   }
//   else {
//     runTreeModel(spMinMY, spMaxMY, spatial, 0);
//   }

  CloseGrids(); //Close all files currently open...

  return EXIT_SUCCESS;
}
