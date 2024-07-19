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
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <thread>
#include "DataOutput.hpp"
#include "DataInput.hpp"
#include "ParamStructs.hpp"

// Maximum file path length. 
#define MAXFILE 1000

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

extern bool modelMode3PGS;

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

class Progress {
private:
    int progress;
    int lastProgress;
    int rowsTotal;
    int rowsDone;
    std::mutex progressMutex;
public:
    Progress(int rowsTotal) {
        this->rowsTotal = rowsTotal;
        this->rowsDone = 0;
        this->progress = 0;
        this->lastProgress = -1;
    }
    void rowCompleted() {
        //lock mutex
        this->progressMutex.lock();

        //increment rows
        this->rowsDone++;

        //calculate updated percentage
        this->progress = (100 * this->rowsDone / this->rowsTotal);

        //if percentage has incremented, display
        if (this->progress > this->lastProgress) {
            fprintf(stdout, "Completed %2u%%\r", this->progress);
            this->lastProgress = this->progress;
        }

        //release mutex
        this->progressMutex.unlock();
    }
};

//----------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    GDALRasterImage* refGrid; // Pointer variable refGrid pointing to GDALRasterImage 
    DataOutput* dataOutput; //thread safe data output class
    bool spatial = 0;
    long nrows, ncols;
    MYDate spMinMY, spMaxMY;
    std::string defParamFile;
    std::string siteParamFile;
    std::unordered_map<std::string, PPPG_OP_VAR> opVars;

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
    DataInput *dataInput = new DataInput();

    /* Copyright */
    std::cout << COPYMSG << std::endl;
    logger.Log(COPYMSG);

    // Load the parameters
    readSpeciesParamFile(defParamFile, dataInput);
    opVars = readSiteParamFile(siteParamFile, dataInput);

    //check that we have all the correct parameters
    if (!dataInput->inputFinished(modelMode3PGS)) {
        exit(EXIT_FAILURE);
    }

    // Check for a spatial run, if so open input grids and define refGrid. 
    openInputGrids();
    refGrid = dataInput->getRefGrid();

    nrows = refGrid->nRows;
    ncols = refGrid->nCols;

    //initialize dataOutput class
    initDataOutput(refGrid);

    // Find the over all start year and end year. 
    // TODO: findRunPeriod reads the entire input grid, which is unnecessary. Find some modern way to do this.
    std::cout << "Finding run period..." << std::endl;
    dataInput->findRunPeriod(spMinMY, spMaxMY);

    readSampleFile(opVars, refGrid); 
    std::cout << "Points read from sample file." << std::endl;
 
    // Run the model. 
    long cellsTotal = nrows * ncols;
    std::cout << "Processing..." << cellsTotal << " cells... " << std::endl;
    logger.Log("Processing..." + to_string(cellsTotal) + " cells... ");

    //unsigned int numThreads = std::thread::hardware_concurrency();
    int nthreads = 4;
    boost::asio::thread_pool pool(nthreads);

    Progress progress(refGrid->nRows);

    for (int i = 0; i < nrows; i++) {
        boost::asio::post(pool, [opVars, spMinMY, spMaxMY, i, ncols, dataInput, &progress] {
            int cellIndexStart = i * ncols;
            for (int j = 0; j < ncols; j++) {
                int cellIndex = cellIndexStart + j;
                runTreeModel(opVars, spMinMY, spMaxMY, cellIndex, dataInput);
            }

            writeRowDataOutput(i);
            progress.rowCompleted();
        });
    }

    pool.join();

    deleteDataOutput();
    delete dataInput;

    return EXIT_SUCCESS;
}