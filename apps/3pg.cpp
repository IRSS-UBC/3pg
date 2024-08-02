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
//
// 
// target Windows 7 and above
#define _WIN32_WINNT 0x0601

#include <cstdlib>
#include <cstring>
#include <iostream>
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

//----------------------------------------------------------------------------------------
std::string VERSION = "4.0";
std::string COPYMSG = "3-PG: A Forest Growth Model for Natural Climate Solutions\n"
"Developed for use by BP\n\n"
"For more information, please see:\n"
"Website:\n"
"Repository:\n"
"----------------------------------------------------------------------------\n"
"This version of 3-PG has been developed by:\n"
"Integrated Remote Sensing Studio (IRSS),\n"
"Faculty of Forestry, University of British Columbia\n"
"2424 Main Mall, Vancouver, BC, Canada, V6T 1Z4\n"
"Developers : Nicholas Coops, Fracois du Toit, Sarah Zwiep, and Joe Meyer\n"
"Contact : nicholas.coops@ubc.ca\n\n"

"Revision : " + VERSION + "\n"
"Date : AUGUST 2024\n\n"
"3-PG is based on work by Ander Siggins and CSIRO Australia\n"
"----------------------------------------------------------------------------\n"
"\"DISCLAIMER\"\n"
"IRSS accepts no responsibility for the use of 3PG(S) or of the model 3-PG in\n"
"the form supplied or as subsequently modified by third parties. IRSS disclaims\n"
"liability for all losses, damages and costs incurred by any person as a result\n"
"of relying on this software. Use of this software assumes agreement to this\n"
"condition of use.\n\n"

"Removal of this statement violates the spirit in which 3-PG was released by\n"
"IRSS and CSIRO Australia\n\n"

"Licensed under CC BY 4.0\n"
"----------------------------------------------------------------------------\n\n";

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
    std::string prefix;
    const int outputWidth;

    void printProgress() {
        std::cout << '\r' << std::string(outputWidth, ' ') << '\r';
        std::cout << prefix << ' ' << std::setw(3) << progress << "% complete" << std::flush;
    }

public:
    Progress(int rowsTotal, const std::string& prefix = "Running 3PG model...")
        : rowsTotal(rowsTotal), rowsDone(0), progress(0), lastProgress(0),
        prefix(prefix), outputWidth(prefix.length() + 25)
    {  
        printProgress();
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
            printProgress();
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
    bool spatial = 0;
    long nrows, ncols;
    MYDate spMinMY, spMaxMY;
    std::string defParamFile;
    std::string siteParamFile;
    std::unordered_map<std::string, PPPG_OP_VAR> opVars;

    /* Copyright */
    std::cout << COPYMSG << std::endl;

    /* Parse command line args */
    InputParser input(argc, argv);
    if (!input.cmdOptionExists("-d")) {
        std::cout << "Missing species definition file. Pass path with -d flag" << std::endl;
        exit(EXIT_FAILURE);
    }
    defParamFile = input.getCmdOption("-d");
    if (defParamFile.empty()) {
        std::cout << "Path to species definition file is empty. Exiting... " << std::endl;
        exit(EXIT_FAILURE);
    }
    if (!input.cmdOptionExists("-s")) {
        std::cout << "Missing site parameter file. Pass path with -s flag" << std::endl;
        exit(EXIT_FAILURE);
    }
    siteParamFile = input.getCmdOption("-s");
    if (siteParamFile.empty()) {
        std::cout << "Path to site parameter file is empty. Exiting... " << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string outPath = getOutPathTMP(siteParamFile);
    logger.StartLog(outPath);
    logger.Log(COPYMSG);
    DataInput *dataInput = new DataInput();

    // Load the parameters
    std::cout << "Loading and validating parameters...";
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
    std::cout << "  Complete" << std::endl;
    dataInput->findRunPeriod(spMinMY, spMaxMY);

 
    // Run the model. 
    long cellsTotal = nrows * ncols;

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