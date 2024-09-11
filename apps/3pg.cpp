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
#include "util.hpp"

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
        prefix(prefix), outputWidth(static_cast<int>(prefix.length()) + 25)
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

class Logger {
private:
    std::string logName;
    std::string logLoc;
    std::ofstream log;
    bool logging = false;
public:
    Logger(const std::string& filename) {
        logName = filename;
    }

    ~Logger() {
        if (log.is_open())
        {
            log.close();
        }
    }

    void StartLog(const std::string& outPath) {
        this->logging = true;
        logLoc = outPath;
        log.open(logLoc + logName, std::ios::trunc);
        std::string currDate = GetCurrentDate();
        std::string currTime = GetCurrentTime();
        log << "-------------------\n";
        log << "OS date: " << std::setw(20) << currDate << "\n";
        log << "OS time: " << std::setw(20) << currTime << "\n";
        log << "-------------------\n";
        log.close();
    }

    std::string GetCurrentDate() {
        if (!this->logging) {
            return "ERROR logger not started";
        }
        auto now = std::chrono::system_clock::now();
        auto local_time = std::chrono::current_zone()->to_local(now);
        return std::format("{:%d-%m-%Y}", local_time);
    }

    std::string GetCurrentTime() {
        if (!this->logging) {
            return "ERROR logger not started";
        }
        auto now = std::chrono::system_clock::now();
        auto local_time = std::chrono::current_zone()->to_local(now);
        return std::format("{:%H:%M:%S}", local_time);
    }

    void Log(const std::string& logMsg){
        if (!this->logging) {
            return;
        }
        log.open(logLoc + logName, std::ios::app);
        log << logMsg + "\n";
        log.close();
    }

};

//----------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    //lambda function of logger, to be passed as parameter to other parts of program that log
    //I'm doing this because every other way I tried resulted in some encredibly bizzare build errors.
    //
    //A potential long term goal would be to have the logging class become a private member of the dataInput
    //class, although this seems to work for now.
    std::string filename = "logfile.txt";
    Logger logger(filename);
    std::function<void(std::string)> log = [&logger](std::string message) {
        logger.Log(message);
    };

    bool spatial = 0;
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

    setLogFunc(log);
    std::string outPath = getOutPathTMP(siteParamFile);
    logger.StartLog(outPath);
    DataInput dataInput(log);

    /* Copyright */
    std::cout << COPYMSG << std::endl;
    logger.Log(COPYMSG);

    // Load the parameters
    std::cout << "Loading and validating parameters...";
    readSpeciesParamFile(defParamFile, dataInput);
    readSiteParamFile(siteParamFile, dataInput);

    //check that we have all the correct parameters
    if (!dataInput.inputFinished(modelMode3PGS)) {
        exit(EXIT_FAILURE);
    }

    // Check for a spatial run, if so open input grids and define refGrid. 
    RefGridProperties refGrid = dataInput.getRefGrid();
    DataOutput dataOutput(refGrid, outPath, dataInput.getOpVars());

    // Find the over all start year and end year. 
    std::cout << "  Complete" << std::endl;
    dataInput.findRunPeriod(spMinMY, spMaxMY);
 
    // Run the model. 
    long cellsTotal = refGrid.nRows * refGrid.nCols;

    //unsigned int numThreads = std::thread::hardware_concurrency();
    int nthreads = 4;
    boost::asio::thread_pool pool(nthreads);

    Progress progress(refGrid.nRows);

    for (int i = 0; i < refGrid.nRows; i++) {
        boost::asio::post(pool, [spMinMY, spMaxMY, i, refGrid, &dataInput, &dataOutput, &progress] {
            int cellIndexStart = i * refGrid.nCols;
            for (int j = 0; j < refGrid.nCols; j++) {
                int cellIndex = cellIndexStart + j;
                runTreeModel(spMinMY, spMaxMY, cellIndex, dataInput, dataOutput);
            }

            dataOutput.writeRow(i);
            progress.rowCompleted();
        });
    }

    pool.join();

    return EXIT_SUCCESS;
}