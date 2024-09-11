// Data input routines for 3PG.  

/*
All source code remains the property and copyright of CSIRO. 

CSIRO accepts no responsibility for the use of 3PG(S) or of the model 3-PG in
the form supplied or as subsequently modified by third parties. CSIRO disclaims
liability for all losses, damages and costs incurred by any person as a result
of relying on this software. 
Use of this software assumes agreement to this condition of use
*/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdint>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include <filesystem>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include "GDALRasterImage.hpp"
#include "util.hpp"
#include "Data_io.hpp"
#include "ParamStructs.hpp"

// namespace fs = std::filesystem;
static std::string rcsid = "$Id: Data_io.cpp,v 1.10 2001/08/02 06:41:01 lou026 Exp $";

// Case insensitive string comparison is strncasecmp on Solaris and strnicmp on WIN32. 
// Use strncasecmp in the code. 
#ifdef WIN32
#define strncasecmp strnicmp
#endif

#define Pi 3.1415927
#ifndef MAXLINE
#define MAXLINE 1000
#endif

#define PPPG_MAX_SERIES_YEARS 150
#define PPPG_MAX_SERIES_LENGTH PPPG_MAX_SERIES_YEARS*12

//logger
//extern Logger logger; 

// Controls and counters
extern double DaysInMonth[13];                  // array for days in months
extern bool modelMode3PGS;

std::string outPath = "./";

//----------------------------------------------------------------------------------

std::function<void(std::string)> logMessage;
void setLogFunc(std::function<void(std::string)>& log) {
    logMessage = log;
}

//----------------------------------------------------------------------------------

string getOutPathTMP(const std::string& siteParamFile)
{
    /**
     * Get the output path for logging.
     *
     * This function exists ONLY because we need know where to write the log file
     * and the path is in the parameter file, but logging needs to begin before the
     * current parameter reading procedure begins. This function is just
     * frankenstein-ed code from 'readParamFile' and 'readOtherParams'.
     *
     * Should be made redundant when we refactor the parameter reading.
     */
    std::string line;
    std::string cp;
    std::ifstream inFile(siteParamFile);
    while (std::getline(inFile, line))
    {
        if (line.empty())
            continue;
        if (line[0] == '/' && line[1] == '/')
            continue;
        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of(","), boost::token_compress_on);
        for (int i = 0; i < tokens.size(); i++)
        {
            boost::trim(tokens[i]);
        }
        std::string pName = tokens.front();
        boost::trim_if(pName, boost::is_any_of("\""));
        std::vector<std::string> pValues;
        boost::split(pValues, tokens.at(1), boost::is_any_of(" \t"), boost::token_compress_on);

        if (pName == "Output directory")
        {
            if (pValues.empty())
            {
                std::cout << "No output directory specified." << std::endl;
                logMessage("No ouput directory specified.");
                exit(EXIT_FAILURE);
            }
            else if (pValues.size() > 1) {
                std::cout << "More than one value element detected in output directory specification." << std::endl;
                logMessage("More than one value element detected in output directory specification.");
                exit(EXIT_FAILURE);
            }
            else {
                cp = pValues.front();
                if (pValues.empty()) {
                    outPath = ".";
                }
                else {
                    if (cp.back() != '\\')
                        cp += '\\';
                    if (std::filesystem::exists(cp)) {
                        outPath = cp;
                    }
                    else {
                        std::cout << "Output directory " << cp << " does not exist." << std::endl;
                        logMessage("Output directory " + cp + " does not exist.");
                        exit(EXIT_FAILURE);
                    }
                }
            }
            break;
        }
    }
    return outPath;
}
//----------------------------------------------------------------------------------

bool readOtherParam(const std::string& pName, std::vector<std::string> pValue)
{
  // Set various miscellaneous parameters. 
  std::string cp;

  // Look for Output directory. 
  // Allow no directory to be specified, in which case force the current directory. 
  // Check that the directory exists.
  if (namesMatch("Output directory", pName)) {
    if (pValue.empty()) 
    {
        std::cout << "No output directory specified." << std::endl; 
        logMessage("No ouput directory specified.");
        exit(EXIT_FAILURE);
    }
    else if (pValue.size() > 1) {
      std::cout << "More than one value element detected in output directory specification." << std::endl;
      logMessage("More than one value element detected in output directory specification.");
      exit(EXIT_FAILURE);
    }
    else {
      cp = pValue.front();
      if (pValue.empty()) {
        outPath = ".";
        return true;
      }
      else {
        // Check that the directory exists. 
        // Make sure of the trailing /
        if (cp.back() != '\\')
          cp += '\\';
        if (std::filesystem::exists(cp)) {
          outPath = cp;
          return true;
        }
        else {
          std::cout << "Output directory " << cp << " does not exist." << std::endl;
          logMessage("Output directory " + cp + " does not exist.");
          exit(EXIT_FAILURE);
        }
      }
      std::cout << "   output path: " << outPath << std::endl; // "   output path: %s\n
      logMessage("   output path: " + outPath);
      return true;
    }
    
  }
  // Model mode (Standard 3PG or 3PGS)
  else if (namesMatch("Model mode", pName)) {
    if (pValue.empty()) {
      std::cout << "No model mode specified." << std::endl;
      logMessage("No model mode specified.");
      return false;
    }
    if (pValue.size() > 1) {
      std::cout << "More than one value element detected in model mode specification." << std::endl;
      logMessage("More than one value element detected in model mode specification.");
      exit(EXIT_FAILURE);
    }
    if ("3PGS" == pValue.front())
      modelMode3PGS = true;
    else if ("3PG" == pValue.front())
      modelMode3PGS = false;
    else {
      std::cout << "Invalid value for parameter 'Model mode': " << pValue.front() << std::endl;
      logMessage("Invalid value for parameter 'Model mode': " + pValue.front());
      exit(EXIT_FAILURE);
    }
    return true;
  }
  else
    return false;
}

//----------------------------------------------------------------------------------

void readSpeciesParamFile(const std::string& speciesFile, DataInput& dataInput) {
    std::string line, pName;
    std::string pValue;
    std::string cp;
    int lineNo = 0;

    auto isDoubleQuote = [](char c) { return c == '\"'; };
    std::ifstream inFile(speciesFile);
    logMessage("Reading  species parameter from file '" + speciesFile + "'...");
    while (std::getline(inFile, line)) {
        lineNo++;
        if (line.empty()) { continue; }
        if (line[0] == '/' && line[1] == '/') { continue; }
        // Parse the parameter name and parameter value from Species file
        // While still implemented at .txt, format must be:
        //          "paramName", paramValue
        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of(","), boost::token_compress_on);
        for (int i = 0; i < tokens.size(); i++) {
            boost::trim(tokens[i]);
        }
        std::string pName = tokens.front();
        boost::trim_if(pName, boost::is_any_of("\""));
        std::vector<std::string> pValues;
        boost::split(pValues, tokens.at(1), boost::is_any_of(" \t"), boost::token_compress_on);
        if (dataInput.tryAddInputParam(pName, pValues)) {
            continue; 
        }
        else {
            std::cout << "Invalid site parameter: " << pName << std::endl;
            logMessage("Invalid site parameter: " + pName);
            exit(EXIT_FAILURE);
        }
    }
}

void readSiteParamFile(const std::string& paramFile, DataInput& dataInput)
{
  // Read a text file containing 3PG parameters.  Comments are allowed
  // and must begin with C++ style '//'.  Comments can begin at any
  // position on a line.  Parameter lines can begin at any position,
  // and must start with the the string 'id' defined in the array
  // 'params' above, in double quotes. Extra data on a parameter line
  // will be ignored.  Missing data on a parameter line will cause a
  // warning.  Lines which are not comments and cannot be identified
  // as parameters will cause the program to exit.  
  std::string line, pName;
  std::string pValue;
  std::string cp;
  int lineNo=0;
  //std::unordered_map<std::string, PPPG_OP_VAR> opVars;

  auto isDoubleQuote = [](char c) { return c == '\"'; };
  std::ifstream inFile(paramFile);
  logMessage("Reading input parameters from file '" + paramFile + "'...");
  while (std::getline(inFile, line)) {
      lineNo++;
    // Skip blank lines
    if (line.empty())
      continue;
    // Skip comments
    if (line[0] == '/' && line[1] == '/')
      continue;
    // Tokenize the line
    std::vector<std::string> tokens;
    /*for (std::string& i : tokens)
        std::cout << i << ' ' << std::endl;*/
    boost::split(tokens, line, boost::is_any_of(","), boost::token_compress_on);
    // trim leading whitespace from each token using boost::trim
    for (int i = 0; i < tokens.size(); i++)
      boost::trim(tokens[i]);
    // First token is the parameter name
    std::string pName = tokens.front();
    // Trim double quotations from the name
    boost::trim_if(pName, boost::is_any_of("\""));

    // Second and subsequent tokens are the parameter values, put them all into a vector
    std::vector<std::string> pValues;

    //in the case where series parameters are given by year, there will only be one token on the line
    //and boost::split() will throw an error. Check to make sure the token size is safe to call this
    //function on.
    if (tokens.size() > 1) {
        boost::split(pValues, tokens.at(1), boost::is_any_of(" \t"), boost::token_compress_on);
    }

    if (dataInput.tryAddInputParam(pName, pValues)) { continue; }
    if (dataInput.tryAddOutputParam(pName, pValues, lineNo)) { continue; }
    else if (readOtherParam(pName, pValues)) { continue; }
    else if (dataInput.tryAddSeriesParam(pName, pValues, inFile, lineNo)) { continue; }
    else if (dataInput.tryAddManagementParam(pName, inFile, lineNo)) { continue; }
    else {
        std::cout << "Cannot read parameter in file " << paramFile << ", line: " << lineNo << ": " << pName << std::endl;
        logMessage("Cannot read parameter in file " + paramFile + ", line: " + to_string(lineNo) + ": " + pName);
        exit(EXIT_FAILURE);
    }
  }
}