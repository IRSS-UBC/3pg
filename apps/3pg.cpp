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
#include "util.hpp"
#include <boost/program_options.hpp>
#include "DataOutput.hpp"

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
std::string VERSION = "0.1";
std::string COPYMSG = "This version of 3-PG has been revised by:\n"
                        //"Nicholas Coops [Nicholas.Coops@csiro.au],\n"
                        //"Anders Siggins [Anders.Siggins@csiro.au],\n"
                        //"and Andrew Loughhead.\n"
                        "Sarah (Vaughan) and Joe\n"
                        "Version: " + VERSION + "\n"
                        "Revisions based on Siggins' 2.53 version\n\n"
                        "Better message TBD. Enjoy!\n"
                        "--------------------------------------\n";

Logger logger("logfile.txt");

class InputParser {
public:
    InputParser(int& argc, char** argv) {
        for (int i = 1; i < argc; ++i)
            this->tokens.push_back(std::string(argv[i]));
    }
    const std::string& getCmdOption(const std::string& option) const {
        std::vector<std::string>::const_iterator itr;
        itr = std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
            return *itr;
        }
        static const std::string empty_string("");
        return empty_string;
    }
    bool cmdOptionExists(const std::string& option) const {
        return std::find(this->tokens.begin(), this->tokens.end(), option)
            != this->tokens.end();
    }
private:
    std::vector <std::string> tokens;
};

//----------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    std::string optarg;
    //extern int optind;
    //int c;
    //int result;

    GDALRasterImage* refGrid; // Pointer variable refGrid pointing to GDALRasterImage 
    DataOutput* dataOutput; //thread safe data output class
    bool spatial = 0;
    long nrows, ncols;
    MYDate spMinMY, spMaxMY;
    std::string defParamFile;
    std::string siteParamFile;

    /* Parse command line args */
    InputParser input(argc, argv);
    if (!input.cmdOptionExists("-d")) {
        std::cout << "Missing species definition file. Pass path with -d flag." << std::endl;
    }
    defParamFile = input.getCmdOption("-d");
    if (defParamFile.empty()) {
        std::cout << "Path to species definition file is empty. Exiting... " << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!input.cmdOptionExists("-s")) {
        std::cout << "Missing site parameter file. Pass path with -s flag." << std::endl;
    }
    siteParamFile = input.getCmdOption("-s");
    if (siteParamFile.empty()) {
        std::cout << "Path to site parameter file is empty. Exiting... " << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string outPath = getOutPathTMP(siteParamFile);
    logger.StartLog(outPath);

    /* Copyright */
    std::cout << COPYMSG << std::endl;
    logger.Log(COPYMSG);


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

    //initialize dataOutput class
    initDataOutput(refGrid);
  }
  else {
    nrows = 0;
    ncols = 0;
  }

  // Find the over all start year and end year. 
  // TODO: findRunPeriod reads the entire input grid, which is unnecessary. Find some modern way to do this.
  std::cout << "Finding run period..." << std::endl;
  findRunPeriod(spMinMY, spMaxMY); 

  // NOTE: don't think ResetGrids is necessary for GDAL stuff... but I guess we'll see
  // ResetGrids(); 

  if (spatial) {
    // Open output grids. 
    //openOutputGrids( refGrid );
    //std::cout << "Opening regular output grids..." << std::endl;
    //openRegularOutputGrids( refGrid, spMinMY, spMaxMY );
    //std::cout << "Regular output grids opened." << std::endl;
    //std::cout << "Reading points from sample file..." << std::endl;
    // Open and read sample point file
    readSampleFile( refGrid ); 
    std::cout << "Points read from sample file." << std::endl;
    // fprintf(logfp, "Processing %u cells...\n", nrows * ncols);
  }

 
// Run the model. 
  if (spatial) {
     int cellsDone = 0;
     int cellsTotal = (nrows) * (ncols); 
     int lastProgress = -1; 
     std::cout << "Processing..." << cellsTotal << " cells... " << std::endl;
     logger.Log("Processing..." + to_string(cellsTotal) + " cells... ");

     for (int i = 0; i < nrows; i++) {
         for (int j = 0; j < ncols; j++) {

             //calculate/print progress
             int progress = (100 * cellsDone / cellsTotal);
             if (progress > lastProgress) {
                 fprintf(stdout, "Completed %2u%%\r", progress);
             }
            
             int cellIndex = i * ncols + j;
             runTreeModel(spMinMY, spMaxMY, spatial, cellIndex);

             //increment progress
             cellsDone++;
             lastProgress = progress;
         }
     }

     deleteDataOutput();
  }

  return EXIT_SUCCESS;
}
