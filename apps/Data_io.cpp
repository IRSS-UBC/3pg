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
extern Logger logger; 

// Controls and counters
extern double DaysInMonth[13];                  // array for days in months
extern bool modelMode3PGS;

// Time variant management factors
extern int nFertility;                          // size of site fertility array
//extern MANAGE_TABLE Fertility[1000];            // time-variant site fertility
extern int nMinAvailSW;                         // size of MinAvailSW array
//extern MANAGE_TABLE MinAvailSW[1000];           // time-variant MinAvailSW (mm)
extern int nIrrigation;                         // size of irrigation array
//extern MANAGE_TABLE Irrigation[1000];           // time-variant irrigation (ML/y)
extern double Irrig;                            // current annual irrigation (ML/y)

std::string outPath = "./";

//----------------------------------------------------------------------------------

// Initialisation of output variable array. This lists all possible output variables 
// and sets up the mapping of the output variable to its name, which is used in 
// parsing the parameter file.  

std::unordered_set<std::string> output_var_names {
    "StemNo",
    "WF",
    "WR",
    "WS",
    "TotalW",
    "LAI",
    "cLAI",
    "MAI",
    "avDBH",
    "BasArea",
    "StandVol",
    "GPP",
    "cGPP",
    "NPP",
    "cNPP",
    "delWAG",
    "cumWabv",
    "Transp",
    "cTransp",
    "ASW",
    "fSW",
    "fVPD",
    "fT",
    "fNutr",
    "fFrost",
    "APAR",
    "APARu",
    "EvapTransp",
    "cEvapTransp",
    "LAIx",
    "ageLAIx",
    "MAIx",
    "ageMAIx",
    "FR",
    "PhysMod",
    "alphaC",
    "fAge",
    "fracBB",
    "WUE",
    "cWUE",
    "CVI",
    "cCVI",
    "TotalLitter",
    "cLitter"
};

// 3PG management table parameters. At most one value per year. 
PPPG_MT_PARAM FertMT[PPPG_MAX_SERIES_YEARS+1];
PPPG_MT_PARAM IrrigMT[PPPG_MAX_SERIES_YEARS+1];
PPPG_MT_PARAM MinAswMT[PPPG_MAX_SERIES_YEARS+1];

//----------------------------------------------------------------------------------

bool openGrid(PPPG_VVAL& vval)
{
    if (vval.spType == pTif) {
        string openString = "   opening raster from " + vval.gridName + "...";
        string succesReadString = "read raster";
        string failReadString = "failed";

        try {
            vval.g = new GDALRasterImage(vval.gridName);
            logger.Log(openString + succesReadString);
        }
        catch (const std::exception&) {
            logger.Log(openString + failReadString);
            exit(EXIT_FAILURE);
        }

    }
    else {
        return false;
    }
    return true;
}

//----------------------------------------------------------------------------------

bool getVVal(double &val, PPPG_VVAL vval, int k)
{
  GDALRasterImage *fg;
  float result;
  
  if (vval.spType == pScalar)
    val = vval.sval;
  else if (vval.spType == pTif) {
    fg = vval.g;
    result = fg->GetVal(k);
    if (fg->IsNoData(result)) {
      return false; 
    }
    else
    {
       val = result;
    }
  }
  else { // pNull
    return false; 
  }
  return true; 
}

//-----------------------------------------------------------------------------

double lookupManageTable( int year, int table, double def, int k )
{
  // Lookup the value of a cell in a management table, for a year. 
  // Entries in management tables apply up to but not including the 
  // next year listed, for consistency with the VB version.  Ie, the 
  // management tables describe periods, not events.  def is the 'default' 
  // value, which will correspond to the general parameter value for the 
  // attribute (ie there can be both an "FR" parameter and an FR management 
  // table).  
  PPPG_MT_PARAM *mt;  
  GDALRasterImage *fg;
  int i;
  double val; 

  if ( table == MT_FERTILITY ) 
    mt = FertMT; 
  else if (table == MT_MINASW )
    mt = MinAswMT; 
  else if (table == MT_IRRIGATION )
    mt = IrrigMT; 
  else {
    std::cout << "Program error: called lookupManageTable with invalid table" << std::endl;
    logger.Log("Program error: called lookupManageTable with invalid table");
    exit(EXIT_FAILURE);
  }

  // Load the table entries.  Unfortunately we have to reload every 
  // entries value to allow values from earlier entries to persist through 
  // NODATA cells in later entries.  If we hit NODATA, we must look at the next 
  // table entry as, if it exists, its value will apply. 
  val = def; 
  bool hitdata = false; 
  for (i = 0; mt[i].year > 0; i++) {
    // Read earlier table entries. 
    if ( mt[i].data.spType == pScalar ) {
      val = mt[i].data.sval; 
    }
    else if ( mt[i].data.spType == pTif ) {
      fg = (GDALRasterImage *)mt[i].data.g;
      if (fg->IsNoData(fg->GetVal(k)))
        continue; 
      else
        val = fg->GetVal(k);
    }
    if ( mt[i].year >= year )
      break;
  }

  return val; 
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

PPPG_OP_VAR readOutputParam(const std::string& pName, const std::vector<std::string>& pValue, int lineNo)
{
  // For a parameter name pName and a parameter value, pValue, both as strings, read 
  // the value into an appropriate variable. The parameter name can be either the 
  // same as the variable name, or it can be a long descriptive name.  
  std::string cp;
  std::string outstr;
  bool yearlyOutput = false;

  PPPG_OP_VAR opVar;
  if (pValue.empty()) {
    outstr = "No grid name for param " + pName +  " on line: " + to_string(lineNo);
    std::cout << outstr << std::endl;
    logger.Log(outstr);
    exit(EXIT_FAILURE);
  }
  if (pValue.size() > 5) {
      outstr = "More than 5 value elements detected for param " + pName + " on line: " + to_string(lineNo);
      std::cout << outstr << std::endl;
      logger.Log(outstr);
      exit(EXIT_FAILURE);
  }
  // First token in the pValue is the output grid filename, outPath and filename are concatenated for the full path
  opVar.gridName = pValue.front();
  const std::filesystem::path filePath = opVar.gridName;
  if (filePath.extension() == ".tif") // Heed the dot.
  {
      opVar.spType = pTif;
  }
  else
  {
      outstr = filePath.filename().generic_string() + " is an invalid filename. Found " + filePath.extension().generic_string() + " but must be '.tif'";
      std::cout << outstr << std::endl;
      logger.Log(outstr);
      exit(EXIT_FAILURE);
  }
  // Check for optional second, third, fourth and fifth tokens; these are used to specify recurring output pattern.
  // The following parsing rules apply:
  //    - If second token exists, then a third and fourth token must also exist. A fifth token is optional.
  //    - The third token must be an integer, representing the start year of the recurrence pattern.
  //    - The fourth token must be 'monthly' or 'month'
  //    - If fourth token is 'monthly', then fifth token must not exist (assumed to be 1).
  //    - If fourth token is 'month', then fifth token must be an integer between 1 and 12.
  try {
    // NOTE: using .at() accessor here to ensure an exception is thrown if the index is out of range
    // There is likely a better way to do this whole flow of logic, but for now we are sticking as close
    // as possible to the original code.
    cp = pValue.at(1);
    yearlyOutput = true;
  }
  catch (const std::out_of_range&) {
      // No recurring year output
      yearlyOutput = false;
  }
  if (yearlyOutput == true) {
    // Look for start year
    try {
      opVar.recurStart = std::stoi(cp);
    }
    catch (std::invalid_argument const&) {
      outstr = "Expected an integer start year in recuring output specification on line " +  to_string(lineNo);
      std::cout << outstr << std::endl;
      logger.Log(outstr);
      exit(EXIT_FAILURE);
    }
    // Look for interval
    try {
      cp = pValue.at(2);
      try {
        const int interval = std::stoi(cp);
        if (interval == 0)
        {
            outstr = "Found interval of zero years in recuring output specification on line " + to_string(lineNo) + ". Expected non-zero";
            std::cout << outstr << std::endl;
            logger.Log(outstr);
            exit(EXIT_FAILURE);
        }
        opVar.recurYear = interval;
      }
      catch (std::invalid_argument const&) {
          outstr = "Expected an integer interval in recuring output specification on line " + to_string(lineNo);
          std::cout << outstr << std::endl;
          logger.Log(outstr);
          exit(EXIT_FAILURE);
      }
    }
    catch (const std::out_of_range&) {
        outstr = "Found start year but no interval in recuring output specification on line " + to_string(lineNo);
        std::cout << outstr << std::endl;
        logger.Log(outstr);
        exit(EXIT_FAILURE);
    }
    // Look for 'monthly' or 'month' keywords. 
    try {
      cp = pValue.at(3);
      if (cp == "monthly")
        opVar.recurMonthly = true;
      else if (cp == "month") {
        opVar.recurMonthly = false;
        // If 'month', look for the month interger
        try {
          cp = pValue.at(4);
          try {
            opVar.recurMonth = std::stoi(cp);
            if (opVar.recurMonth == 0)
            {
                outstr = "Found month of zero in recuring output specification on line " + to_string(lineNo) + ". Expected non-zero";
                std::cout << outstr << std::endl;
                logger.Log(outstr);
                exit(EXIT_FAILURE);
            }
          }
          catch (std::invalid_argument const&) {
              outstr = "Expected an integer month in recuring output specification on line " + to_string(lineNo);
              std::cout << outstr << std::endl;
              logger.Log(outstr);
              exit(EXIT_FAILURE);
          }
        }
        catch (const std::out_of_range&) {
            outstr = "Found 'month' keyword but no month in recuring output specification on line " + to_string(lineNo);
            std::cout << outstr << std::endl;
            logger.Log(outstr);
            exit(EXIT_FAILURE);
        }
      }
      else {
          outstr = "Unrecognised keyword '" + cp + "' on line " + to_string(lineNo);
          std::cout << outstr << std::endl;
          logger.Log(outstr);
          exit(EXIT_FAILURE);
      }
    }
    catch (const std::out_of_range&) {
        
        outstr = "Found start year and interval but no keyword in recuring output specification on line " + to_string(lineNo);
        std::cout << outstr << std::endl;
        logger.Log(outstr);
        exit(EXIT_FAILURE);
    }
  }
  // Mark the variable for later writing
    opVar.write = true;
    logger.Log("   variable: " + opVar.id + "   grid: " + opVar.gridName);
    if (opVar.recurStart)
    {
        string outputGridString = "      starting in " + to_string(opVar.recurStart) + ", writing every " + to_string(opVar.recurYear) + " years";
        if (opVar.recurMonthly)
        {
            outputGridString = outputGridString + ", with monthly values";
        }
        else if (opVar.recurMonth != 0)
        {
            outputGridString = outputGridString + ", on the " + to_string(opVar.recurMonth) + " month";

        }
        logger.Log(outputGridString);
    }
  return opVar;
}

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
                logger.Log("No ouput directory specified.");
                exit(EXIT_FAILURE);
            }
            else if (pValues.size() > 1) {
                std::cout << "More than one value element detected in output directory specification." << std::endl;
                logger.Log("More than one value element detected in output directory specification.");
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
                        logger.Log("Output directory " + cp + " does not exist.");
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
        logger.Log("No ouput directory specified.");
        exit(EXIT_FAILURE);
    }
    else if (pValue.size() > 1) {
      std::cout << "More than one value element detected in output directory specification." << std::endl;
      logger.Log("More than one value element detected in output directory specification.");
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
          logger.Log("Output directory " + cp + " does not exist.");
          exit(EXIT_FAILURE);
        }
      }
      std::cout << "   output path: " << outPath << std::endl; // "   output path: %s\n
      logger.Log("   output path: " + outPath);
      return true;
    }
    
  }
  // Model mode (Standard 3PG or 3PGS)
  else if (namesMatch("Model mode", pName)) {
    if (pValue.empty()) {
      std::cout << "No model mode specified." << std::endl;
      logger.Log("No model mode specified.");
      return false;
    }
    if (pValue.size() > 1) {
      std::cout << "More than one value element detected in model mode specification." << std::endl;
      logger.Log("More than one value element detected in model mode specification.");
      exit(EXIT_FAILURE);
    }
    if ("3PGS" == pValue.front())
      modelMode3PGS = true;
    else if ("3PG" == pValue.front())
      modelMode3PGS = false;
    else {
      std::cout << "Invalid value for parameter 'Model mode': " << pValue.front() << std::endl;
      logger.Log("Invalid value for parameter 'Model mode': " + pValue.front());
      exit(EXIT_FAILURE);
    }
    return true;
  }
  else
    return false;
}

//----------------------------------------------------------------------------------

bool readParam( PPPG_VVAL &vval, std::string pValue )
{
  std::string cp;

   cp = pValue;
  try {
    vval.sval = std::stod(cp);
    vval.spType = pScalar;
    return true;
  }
  catch (std::invalid_argument const&) {
    // If not, then try a grid name
    cp = pValue;
    const std::filesystem::path filePath = cp;
    // Check that the file is a TIF file
    if (filePath.extension() == ".tif") // Heed the dot.
    {
        vval.gridName = cp;
        vval.spType = pTif;
        return true;
    }
    else
    {
        std::cout << filePath.filename() << " is an invalid filetype (" << filePath.extension() << ")" << std::endl; 
        logger.Log(filePath.filename().generic_string() + " is an invalid filetype (" + filePath.extension().generic_string() + ")");
        exit(EXIT_FAILURE);
    }
    return false;
  }
}

//----------------------------------------------------------------------------------

bool readInputManageParam(const std::string pName, std::ifstream& inFile, int &lineNo)
{
  // Read management table input parameters.
  // A table must begin on the line following the keyword identifying it.
  // The table has one entry per line, each entry consists of a year and a value, seperated by whitespace.
  // A blank line terminates the table.  Each value can be either a constant or a grid name.  
  std::string line;
  std::string tok, cp;
  PPPG_MT_PARAM *tab;
  int i, *nRead; 
  std::string tabName;

  // Are we reading a managment table?
  if ( namesMatch( "Management: fertility", pName ) ) {
    tab = FertMT;
    tabName = "Fertility MT";
    nRead = &nFertility; 
  }
  else if ( namesMatch( "Management: irrigation", pName ) ) {
    tab = IrrigMT; 
    tabName = "Irrigation MT";
    nRead = &nIrrigation; 
  }
  else if ( namesMatch( "Management: MinASW", pName ) ) {
    tab = MinAswMT;
    tabName = "Min ASW MT";
    nRead = &nMinAvailSW;
  }
  else
    return false; 

  // Read the table using ifstreams.
  // The table has one entry per line, each entry consists of a year and a value, seperated by whitespace.
  // A blank line terminates the table.  Each value can be either a constant or a grid name.
  i = 0; 
  while (std::getline(inFile, line)) {
    lineNo++; // Passed-by-reference, will alter the value of lineNo in the calling function.
    if (line.empty())
      // Blank line terminates the table.
      tab[i].year = -1; 
      break;
    // Tokenize the line
    std::vector<std::string> tTokens;
    boost::split(tTokens, line, boost::is_any_of(", \n\t"));
    if (tTokens.size() != 2) {
      std::cout << "Could not read management table at line " << lineNo << std::endl;
      logger.Log("Could not read management table at line " + to_string(lineNo));
      exit(EXIT_FAILURE);
    }
    // trim leading whitespace from each token using boost::trim
    for (int i = 0; i < tTokens.size(); i++)
      boost::trim(tTokens[i]);

    // Read the year
    try {
      tab[i].year = std::stoi(tTokens.front());
    }
    catch (std::invalid_argument const&) {
      std::cout << "Expected an integer year in management table at line " << lineNo << std::endl;
      logger.Log("Expected an interger year in management table at line " + to_string(lineNo));
      exit(EXIT_FAILURE);
    }
    // Read the second token, which is either a constant or a grid name.
    if( !readParam( tab[i].data, tTokens.back() )) {
      std::cout << "Could not read management table value at line " << lineNo << std::endl;
      logger.Log("Could not read management table value at line " + to_string(lineNo));
      exit(EXIT_FAILURE);
    }
    else {
      tab[i].got = 1;
    }
    if (tab[i].data.spType == pScalar) {
      std::cout << "   " << tabName << " year: " << tab[i].year << "   value: " << tab[i].data.sval << std::endl;
      logger.Log("   " + tabName + " year: " + to_string(tab[i].year) + "   value:" + to_string(tab[i].data.sval));
    }
    else {
      std::cout << "   " << tabName << " year: " << tab[i].year << "   grid: " << tab[i].data.gridName << std::endl;
      logger.Log("   " + tabName + " year: " + to_string(tab[i].year) + "   grid: " + tab[i].data.gridName);
    } 
    i++; 
  }
  *nRead = i; 
  return true;  
}

//----------------------------------------------------------------------------------

void readSpeciesParamFile(const std::string& speciesFile, DataInput& dataInput) {
    std::string line, pName;
    std::string pValue;
    std::string cp;
    int lineNo = 0;

    auto isDoubleQuote = [](char c) { return c == '\"'; };
    std::ifstream inFile(speciesFile);
    logger.Log("Reading  species parameter from file '" + speciesFile + "'...");
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
            logger.Log("Invalid site parameter: " + pName);
            exit(EXIT_FAILURE);
        }
    }
}

std::unordered_map<std::string, PPPG_OP_VAR> readSiteParamFile(const std::string& paramFile, DataInput& dataInput)
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
  std::unordered_map<std::string, PPPG_OP_VAR> opVars;

  auto isDoubleQuote = [](char c) { return c == '\"'; };
  std::ifstream inFile(paramFile);
  logger.Log("Reading input parameters from file '" + paramFile + "'...");
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
    if (output_var_names.find(pName) != output_var_names.end()) {
        opVars.emplace(pName, readOutputParam(pName, pValues, lineNo));
    }
    else if (readOtherParam(pName, pValues)) { continue; }
    else if (dataInput.tryAddSeriesParam(pName, pValues, inFile, lineNo)) { continue; }
    else if (readInputManageParam(pName, inFile, lineNo)) { continue; }
    else {
        std::cout << "Cannot read parameter in file " << paramFile << ", line: " << lineNo << ": " << pName << std::endl;
        logger.Log("Cannot read parameter in file " + paramFile + ", line: " + to_string(lineNo) + ": " + pName);
        exit(EXIT_FAILURE);
    }
  }
  return opVars;
}

//----------------------------------------------------------------------------------

GDALRasterImage* openInputGrids( )
{
  // Open any grids in the params array, the climate and NDVI series arrays, and 
  // the management tables. 
  // Copy the grid parameters of the first grid opened to refGrid. 
  int j; 
  bool spatial = false, first = true;
  GDALRasterImage *refGrid;

  logger.Log("Opening input rasters...");

  // Open all management table grids. 
  PPPG_MT_PARAM *tab;
  PPPG_MT_PARAM *tablist[] = { FertMT, IrrigMT, MinAswMT, NULL }; 
  for (j = 0 ; tablist[j] != NULL; j++) {
    tab = tablist[j]; 
    for (int i = 0; tab[i].year > 0; i++) {
      if ( openGrid( tab[i].data ) ) {
        spatial = true; 
        if ( first ) {
          refGrid = (GDALRasterImage *)tab[i].data.g;
          first = false; 
        }
        else if ( ( fabs( refGrid->xMin - tab[i].data.g->xMin ) > 0.0001 ) 
          || ( fabs( refGrid->yMin - tab[i].data.g->yMin ) > 0.0001 )
          || ( fabs( refGrid->xMax - tab[i].data.g->xMax ) > 0.0001 )
          || ( fabs( refGrid->yMax - tab[i].data.g->yMax ) > 0.0001 ) 
          || ( refGrid->nRows != tab[i].data.g->nRows ) 
          || ( refGrid->nCols != tab[i].data.g->nCols ) ) {
            std::cout << "Grid dimensions must match, raster " << tab[i].data.gridName << " differs from first raster." << std::endl;
            logger.Log("Grid dimensions must match, raster " + tab[i].data.gridName + " differs from first raster.");
            exit(EXIT_FAILURE);
          // sprintf(outstr, "Grid dimensions must match, grid %s differs from first grid.\n", 
          //   tab[i].data.gridName ); 
          // logAndExit(logfp, outstr); 
        }
      }
    }
  }

  //if (!spatial) {
  //  std::cout << "None" << std::endl;
  //  // fprintf(logfp, "none\n");
  //  refGrid = NULL; 
  //}

  return nullptr;
}

//----------------------------------------------------------------------------------

int writeOutputGrids(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, DataOutput& dataOutput, long cellIndex) {
    //for each possible output variable
    for (auto& [pN, opV] : opVars) {
        //if it has been marked to be written
        if (opV.write) {
            //determine value, name, and tell dataOutput class to write
            float val = (float)(opV.v);
            std::string name = opV.gridName;
            name = name.substr(0, name.find_last_of("."));
            dataOutput.setVal(-1, -1, name, cellIndex, val);
        }
    }
    return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------------

void writeMonthlyOutputGrids(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, DataOutput& dataOutput, int calYear, int calMonth, MYDate minMY, MYDate maxMY, long cellIndex) {
    //for each possible output variable
    for (auto& [pN, opV] : opVars) {

        //skip output variable if it is not marked for recurring output
        if (opV.recurYear == -1) {
            continue;
        }

        //skip output variable if it is not marked for recurring output
        if (!opV.recurStart) {
            continue;
        }

        // skip output variable if it is not at the recur interval
        if (((calYear - opV.recurStart) % opV.recurYear) != 0) {
            continue;
        }

        //mx is no longer used to index an array, but is useful (for now) for checking
        //whether we've gone above or below the max or min allowed year/month combo.
        int mx = (calYear - minMY.year) * 12 + (calMonth - 1);
        int maxInd = (maxMY.year - minMY.year + 1) * 12 + 12;

        //ensure the year and month are not below the min
        if (mx < 0) {
            continue;
        }

        //ensure the year and month are not above the max
        if (mx > maxInd) {
            string outStr = "Program error, mx=" + to_string(mx) + " too high in writeMonthlyOutputGrids at month/year " + to_string(calMonth) + "/" + to_string(calYear);
            std::cout << outStr << std::endl;
            logger.Log(outStr);
            exit(EXIT_FAILURE);
        }

        //skip output variable if we're not meant to be printing every month AND we're not on the month we're meant to be printing
        if (!opV.recurMonthly && opV.recurMonth != calMonth) {
            continue;
        }

        //ensure output type is tif
        if (opV.spType != pTif) {
            std::string outStr = "output type must be tif";
            std::cout << outStr << std::endl;
            logger.Log(outStr);
            exit(EXIT_FAILURE);
        }

        //determine value, name, and tell dataOutput class to write
        float val = (float)(opV.v);
        std::string name = opV.gridName;
        name = name.substr(0, name.find_last_of("."));
        dataOutput.setVal(calYear, calMonth, name, cellIndex, val);
    }
}