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

#define GRID_NAME_LENGTH 300
#define PPPG_MAX_SERIES_YEARS 150
#define PPPG_MAX_SERIES_LENGTH PPPG_MAX_SERIES_YEARS*12

extern Logger logger;

//----------------------------------------------------------------------------------
// Global variables.  These are used to push values into the runTreeModel routine. 
// Where the VB version uses an integer variable the equivalent C declaration 
// has been left in but commented.  Apart from that this section is a direct 
// translation from the The_3PG_Model module in the VB version.  Note that the 
// arrays used directly by the models have been left '1' based, so there is an 
// unused array member at index 0.  

// Controls and counters
extern double DaysInMonth[13];                  // array for days in months
extern bool modelMode3PGS;

// Site characteristics, site specific parameters
extern std::string siteName;                      // name of site
extern double FR;                              // current site fertility rating
extern double SWconst, SWpower;                 // soil parameters for soil class

// Time variant management factors
extern int nFertility;                          // size of site fertility array
//extern MANAGE_TABLE Fertility[1000];            // time-variant site fertility
extern int nMinAvailSW;                         // size of MinAvailSW array
//extern MANAGE_TABLE MinAvailSW[1000];           // time-variant MinAvailSW (mm)
extern int nIrrigation;                         // size of irrigation array
//extern MANAGE_TABLE Irrigation[1000];           // time-variant irrigation (ML/y)
extern double Irrig;                            // current annual irrigation (ML/y)

// Mean monthly weather data
//int mYears;                                   // years of met data available
// ANL changed this from int to double
extern double mYears;                           // years of met data available
extern double mDayLength[13];                   // day length
//int mFrostDays[13];                           // frost days/month
// ANL changed this from int to double
extern double mFrostDays[13];                   // frost days/month
extern double mSolarRad[13];                    // solar radiation (MJ/m2/day)
extern double mTx[13];                          // maximum temperature
extern double mTn[13];                          // minimum temperature
extern double mTav[13];                         // mean daily temperature
extern double mVPD[13];                         // mean daily VPD
extern double mRain[13];                        // total monthly rain + irrigation
extern double mNDVI[13];                        // ANL monthly NDVI for 3PGS mode
extern double mNetRad[13];                      // ANL can use net instead of short wave

// Stand data
// extern char SpeciesName[100];                // name of species
// int StandAge;                                // stand age
// ANL changed StandAge from int to double
extern double StandAge;                         // stand age
extern double ASW;                        // available soil water
extern double StemNo;                  // stem numbers
extern double WF;                          // foliage biomass
extern double WR;                          // root biomass
extern double WS;                          // stem biomass
extern double LAIi, LAI;                        // canopy leaf area index
extern double MAIi, MAI;                        // mean annual volume increment
extern double avDBHi, avDBH;                    // average stem DBH
extern double TotalW;                           // total biomass
extern double BasArea;                          // basal area
extern double StandVol;                         // stem volume
extern double LAIx, ageLAIx;                    // peak LAI and age at peak LAI
extern double MAIx, ageMAIx;                    // peak MAI and age at peak MAI
extern double cumTransp;                        // annual stand transporation
extern double cumIrrig;                         // annual irrig. to maintain MinASW

// Stand factors that are specifically age dependent
extern double SLA;
extern double Littfall;
extern double fracBB;
extern double CanCover;

// Parameter values
// int MaxAge;
// ANL changed MaxAge from int to double
extern double Interception;
extern double Density;
extern double pfsConst, pfsPower;                     // derived from pFS2, pFS20

// Intermediate monthly results
extern double m, alphaC;
extern double RAD, PAR;
extern double lightIntcptn;
extern double fAge, fT, fFrost;
extern double fVPD, fSW, fNutr;
extern double CanCond;
extern double Transp, EvapTransp;
extern double AvStemMass;
extern double APAR, APARu;
extern double GPPmolc, GPPdm, NPP;
extern double pR, pS, pF, pFS;
extern double delWF, delWR, delWS;
extern double delFloss, delRloss;
extern double monthlyIrrig;

// Annual results
extern double cLAI, cGPP, cNPP, cCVI, cRainInt, cEvapTransp, cTransp, cWUE;
extern double cumGPP, cumWabv;
extern double abvgrndEpsilon, totalEpsilon;
extern double StemGrthRate;
extern double cLitter;
extern double CumdelWF, CumdelWR, CumdelWS;
extern double CumAPARU, cumARAD;
extern double CumStemLoss;
extern double CutStemMass1, CutStemMass2, CutStemMass3;

//----------------------------------------------------------------------------------

// 3PGS variables
extern double delWAG;

// ANL - other globals.
bool yearlyOutput, monthlyOutput; 
bool samplePointsYearly = false, samplePointsMonthly = false;
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

 
//----------------------------------------------------------------------------------

// 3PG series parameters, climate and NDVI. These either have 12 values, for 
// nominal values such as Esoclim data, or some multiple of 12 values, for time 
// series data.  
PPPG_SERIES_PARAM Tmax_vals;
PPPG_SERIES_PARAM Tmin_vals;
PPPG_SERIES_PARAM Tavg_vals;
PPPG_SERIES_PARAM Rain_vals;
PPPG_SERIES_PARAM SolarRad_vals;
PPPG_SERIES_PARAM FrostDays_vals;
PPPG_SERIES_PARAM NdviAvh_vals;
PPPG_SERIES_PARAM NetRad_vals;
PPPG_SERIES_PARAM Vpd_vals;

// 3PG management table parameters. At most one value per year. 
PPPG_MT_PARAM FertMT[PPPG_MAX_SERIES_YEARS+1];
PPPG_MT_PARAM IrrigMT[PPPG_MAX_SERIES_YEARS+1];
PPPG_MT_PARAM MinAswMT[PPPG_MAX_SERIES_YEARS+1];

//----------------------------------------------------------------------------------

// Sample points. 
#define MAX_SAMPLE_POINTS 50
std::string sampleIpFile;

FILE *sampleIpFp; 
struct {
  std::string id = "-1";
  FILE *fp; 
  double lat; 
  double lon; 
  int cellIndex;
} samplePoints[MAX_SAMPLE_POINTS + 1];

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//Define and provide initialization/deletion functions for data output:
// 
//I imagine this is not the best way of going about this. This will likely need to be
//changes as we refactor the way Data_io works.
DataOutput* dataOutput;
void initDataOutput(GDALRasterImage* refGrid) {
    dataOutput = new DataOutput(refGrid, outPath);
}
void deleteDataOutput() {
    delete dataOutput;
    dataOutput = nullptr;
}

bool openGrid(PPPG_VVAL& vval)
{
    if (vval.spType == pTif) {
        string openString = "   opening raster from " + vval.gridName + "...";
        string succesReadString = "read raster";
        string failReadString = "failed";
        std::cout << openString;

        try {
            vval.g = new GDALRasterImage(vval.gridName);
            std::cout << succesReadString << std::endl;
            logger.Log(openString + succesReadString);
        }
        catch (const std::exception& e) {
            std::cout << failReadString << std::endl;
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
  float val; 

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

bool getSeriesVal(double& val, int ser, int calMonth, int calYear, int k)
{
    PPPG_SERIES_PARAM* series;
    std::string errstr;

    // Which series. 
    switch (ser) {
    case SS_TMAX:      series = &Tmax_vals;      break;
    case SS_TMIN:      series = &Tmin_vals;      break;
    case SS_TAVG:      series = &Tavg_vals;      break;
    case SS_RAIN:      series = &Rain_vals;      break;
    case SS_SOLARRAD:  series = &SolarRad_vals;  break;
    case SS_FROSTDAYS: series = &FrostDays_vals; break;
    case SS_NDVI_AVH:  series = &NdviAvh_vals;   break;
    case SS_NETRAD:    series = &NetRad_vals;    break;
    case SS_VPD:       series = &Vpd_vals;       break;
    default: break;
    }

    // Long run or 1 year data.
    int i;
    if (series->oneYear) {
        // calYear is irrelevant, calMonth will range from 1 to 12, series elements will 
        // be indexed from 0 to 11. 
        i = calMonth - 1;
    }
    else {
        // Long run data. 
        i = (calYear - series->start) * 12 + calMonth - 1;
        if (i > series->vlen * 12 - 1) {
            errstr = "Attempted lookup of series element " + to_string(i) + " in series " + to_string(ser) + ", only " + to_string((series->vlen * 12 - 1)) + " entries in series.\nCheck Start age and End years.";
            std::cout << errstr << std::endl;
            logger.Log(errstr);
            exit(EXIT_FAILURE);
        }
        if (i < 0)
        {
            errstr = "Attempted lookup of year before start year";
            std::cout << errstr << std::endl;
            logger.Log(errstr);
            exit(EXIT_FAILURE);
        }
    }
    if (getVVal(val, series->data[i], k))
        if (series->data[i].g->IsNoData(val)) {
            return false;
        }
        else {
            return true;
        }  
    return false; 

}
//----------------------------------------------------------------------------------

void readSampleFile(std::unordered_map<std::string, PPPG_OP_VAR> &opVars, GDALRasterImage *refGrid)
{
  // Read a text file of sample points, one per line, in the format idstring, 
  // xcoord, ycoord; find the index number of the cell the points fall in. 
  char *line, *fname;
  char *id, *xstr, *ystr, *cp;
  int ind=0;
  int opn;
  double lat, lon;
  int cellIndex;
  std::string errstr;

  if (!samplePointsMonthly && !samplePointsYearly)
    return;

    if ((sampleIpFp = fopen(sampleIpFile.c_str(), "r")) == NULL) {
        std::string ipStr = sampleIpFile;
        errstr = "Could not open sample point file " + ipStr;
        std::cout << errstr << std::endl;
        logger.Log(errstr);
        exit(EXIT_FAILURE);
    }
  while (fgets(line, MAXLINE-1, sampleIpFp) != NULL) {
    if (sscanf(line, "%s %s %s", id, xstr, ystr) != 3)
      return;
    // Change 'D' to 'e'.  Arcinfo ungenerate writes exponents with D, 
    // atof and scanf only read e or E. 
    // TODO: Need to enforce sample points use e/E 
    // for (cp = xstr; *cp != '\0'; cp++)
    //   if (*cp == 'D')
    //     *cp = 'e';
    // for (cp = ystr; *cp != '\0'; cp++)
    //   if (*cp == 'D')
    //     *cp = 'e';
    
    // Read coordinates.
    // TODO: convert id to string
    samplePoints[ind].id = id;
    lat = atof(xstr);
    lon = atof(ystr);
    // TODO: Need to create a function that gets cell index from x/y with GDAL
    cellIndex = refGrid->IndexFrom(lat, lon);
    samplePoints[ind].lat = lat;
    samplePoints[ind].lon = lon;
    samplePoints[ind].cellIndex = cellIndex;

    // Open output file for this point. Make sure its in outPath. 
    // sprintf(fname, "%s3pg.%s.txt", outPath, id);
    if ((samplePoints[ind].fp = fopen(fname, "w")) == NULL)
    {
        string nameStr = fname;
        errstr = "Error opening output sample file " + nameStr;
        std::cout << errstr << std::endl;
        logger.Log(errstr);
        exit(EXIT_FAILURE);
    }

    // Write header line for each sample file. 
    // For each output sample file
    fprintf(samplePoints[ind].fp, "year, month, id, ");
    for (auto& [key, oV] : opVars) {
      fprintf(samplePoints[ind].fp, "%s, ", oV.id.c_str());
    };
    ind++;
  }
  // Make sure end of sample is marked. 
  samplePoints[ind].id[0] = 0;
}

//----------------------------------------------------------------------------------

PPPG_OP_VAR readOutputParam(const std::string& pName, const std::vector<std::string>& pValue, int lineNo)
{
  // For a parameter name pName and a parameter value, pValue, both as strings, read 
  // the value into an appropriate variable. The parameter name can be either the 
  // same as the variable name, or it can be a long descriptive name.  
   int pInd, pInd1, pInd2;
  std::string cp;
  std::string outstr;

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
  catch (const std::out_of_range& oor) {
    std::cout << "No recurring year output detected." << std::endl;
  }
  if (yearlyOutput == true) {
    // Look for start year
    try {
      opVar.recurStart = std::stoi(cp);
    }
    catch (std::invalid_argument const& e) {
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
      catch (std::invalid_argument const& e) {
          outstr = "Expected an integer interval in recuring output specification on line " + to_string(lineNo);
          std::cout << outstr << std::endl;
          logger.Log(outstr);
          exit(EXIT_FAILURE);
      }
    }
    catch (const std::out_of_range& oor) {
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
          catch (std::invalid_argument const& e) {
              outstr = "Expected an integer month in recuring output specification on line " + to_string(lineNo);
              std::cout << outstr << std::endl;
              logger.Log(outstr);
              exit(EXIT_FAILURE);
          }
        }
        catch (const std::out_of_range& oor) {
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
    catch (const std::out_of_range& oor) {
        
        outstr = "Found start year and interval but no keyword in recuring output specification on line " + to_string(lineNo);
        std::cout << outstr << std::endl;
        logger.Log(outstr);
        exit(EXIT_FAILURE);
    }
  }
  // Mark the variable for later writing
    opVar.write = true;
    std::cout << "   variable: " << opVar.id << "   grid: " << opVar.gridName << std::endl;
    logger.Log("   variable: " + opVar.id + "   grid: " + opVar.gridName);
    if (opVar.recurStart)
    {
        string outputGridString = "      starting in " + to_string(opVar.recurStart) + ", writing every " + to_string(opVar.recurYear) + " years";
        std::cout << outputGridString << " years";
        if (opVar.recurMonthly)
        {
            std::cout << ", with monthly values";
            outputGridString = outputGridString + ", with monthly values";
        }
        else if (opVar.recurMonth != 0)
        {
            std::cout << ", on the " << opVar.recurMonth << " month";
            outputGridString = outputGridString + ", on the " + to_string(opVar.recurMonth) + " month";

        }
        std::cout << std::endl;
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
  // Look for sample points file.
  // if pName matched but pValue is empty, bail.
  // If the parameter is present it must have a value. 
  // Allow optional "monthly" keyword following file name. Delimeted by space.
  else if (namesMatch("sample points file", pName)) {
    if (pValue.empty()) {
      std::cout << "No sample points file specified." << std::endl;
      logger.Log("No sample points file specified.");
      exit(EXIT_FAILURE);
    }
    else if (pValue.size() > 2) {
      std::cout << "More than two value elements detected in sample points file specification." << std::endl;
      logger.Log("More than two value elements detected in sample points file specification.");
      exit(EXIT_FAILURE);
    }
    
    // Set first token to sampleIpFile
    sampleIpFile = pValue.front();
    // If there is a second token, check if it is 'monthly' and set samplePointsMonthly to true
    // otherwise set samplePointsYearly to true
    if (pValue.size() == 2) {
      cp = pValue.at(1);
      if (cp == "monthly") {
        samplePointsMonthly = true;
      }
      else if (cp == "yearly") {
        samplePointsYearly = true;
      }
      else {
        std::cout << "Unrecognised keyword \"" << cp << "\" in sample points file specification." << std::endl;
        logger.Log("Unrecognised keyword '" + cp + "' in sample points file specification.");
        exit(EXIT_FAILURE);
      }
    }
    else {
      samplePointsYearly = true;
    }
    return true;
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
  double dv; 

   cp = pValue;
  try {
    vval.sval = std::stod(cp);
    vval.spType = pScalar;
    return true;
  }
  catch (std::invalid_argument const& inv) {
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
    catch (std::invalid_argument const& inv) {
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

bool readInputSeriesParam(std::string pName, std::vector<std::string> pValue, std::ifstream& paramFp, int &lineNo)
{
  // Read 'series' input parameters, ie climate and NDVI. 
  // Two styles of input are permitted. 
  // Firstly, the parameter name can be followed on the same line by 12 values, in this case the 12 values will 
  // be reused for each run year.
  // Secondly, the parameter name can be the only thing on the line (other than a comment),
  // on each following line must be a year followed by 12 values, until the sequence is terminated by a blank 
  // line.  The year values must be in ascending order. 
  PPPG_SERIES_PARAM *series; 
  int ser; 
  int series_yr, prev_yr; 
  std::string tok, line;

  // Which series. 
  if (      namesMatch ("Tmax",        pName ) ) ser = SS_TMAX; 
  else if ( namesMatch ("Tmin",        pName ) ) ser = SS_TMIN;
  else if ( namesMatch ("Tavg",        pName ) ) ser = SS_TAVG;
  else if ( namesMatch ("Rain",        pName ) ) ser = SS_RAIN; 
  else if ( namesMatch ("Solar Radtn", pName ) ) ser = SS_SOLARRAD; 
  else if ( namesMatch ("Frost days",  pName ) ) ser = SS_FROSTDAYS; 
  else if ( namesMatch ("NDVI_AVH",    pName ) ) ser = SS_NDVI_AVH; 
  else if ( namesMatch ("Net radtn",   pName ) ) ser = SS_NETRAD; 
  else if ( namesMatch ("VPD",         pName ) ) ser = SS_VPD; 
  else return false; 

  switch (ser) {
  case SS_TMAX:      series = &Tmax_vals; break; 
  case SS_TMIN:      series = &Tmin_vals; break; 
  case SS_TAVG:      series = &Tavg_vals; break;
  case SS_RAIN:      series = &Rain_vals; break; 
  case SS_SOLARRAD:  series = &SolarRad_vals; break; 
  case SS_FROSTDAYS: series = &FrostDays_vals; break; 
  case SS_NDVI_AVH:  series = &NdviAvh_vals; break; 
  case SS_NETRAD:    series = &NetRad_vals; break; 
  case SS_VPD:       series = &Vpd_vals; break; 
  default: break;
  }

  // Check that latitude is already set. 
  /*if (Lat > 400) {
    sprintf(outstr, "Must set Latitude parameter before series parameters\n");
    logAndExit(logfp, outstr);
  }
  */
  // See which style of input the series has, and where to start
  // populating the series arrays. In the northern hemisphere the
  // model will run from january to december, in the southern
  // hemisphere it will run from july to june. The series parameters
  // are indexed from january of the starting year in both case.  When
  // we only have one years data and are running in the southern
  // hemisphere put the Jan-June values after the July-Dec values.
  if (pValue.empty()) 
    // First style described above. 
    series->oneYear = false;
  else 
    // Second style above.  
    series->oneYear = true; 
  
  // Read values for one year style. 
  if ( series->oneYear ) {
    series->vlen = 1; 
    series->data = new PPPG_VVAL[series->vlen * 12]; 
    int i;

    for ( i = 0; i < 12; i++ ) {
      try{
        std::string mValue = pValue.at(i);
        if ( readParam( series->data[i], pValue[i] ) ) {
          if ( series->data[i].spType == pScalar ) {
              string monthlyConstantStr = "   " + pName + "           month " + to_string(i + 1) + " constant: " + to_string(series->data[i].sval);
              std::cout << monthlyConstantStr << std::endl;
              logger.Log(monthlyConstantStr);
          }
          else if (series->data[i].spType == pTif) {
              string monthlyGridStr = "   " + pName + "          month " + to_string(i + 1) + " grid: " + series->data[i].gridName;
              std::cout << monthlyGridStr << std::endl;
              logger.Log(monthlyGridStr);
          }
          }
          else {
            string monthlyFailStr = "Could not read parameter " + pName + " at month " + to_string(i + 1);
            std::cout << monthlyFailStr << std::endl;
            logger.Log(monthlyFailStr);
            exit(EXIT_FAILURE);
          }
      } catch (const std::out_of_range& oor) {
          std::cout << "No value for " << pName << " at month "  << i+1 << std::endl;
          logger.Log("No value for " + pName + " at month " + to_string(i + 1));
          exit(EXIT_FAILURE);
          // sprintf(outstr, "Incomplete series on line %d.\n", lineNo); 
          // logAndExit(logfp, outstr); 
      }
    }
    series->got = true;

    if (i < 12) {
        std::cout << "Incomplete series on line " << lineNo << std::endl;
        logger.Log("Incomplete series on line " + to_string(lineNo));
        exit(EXIT_FAILURE);
    }
  }
  // Time series style
  else {
    // Find out how many years in the series. 
    int place = paramFp.tellg(); 
    int ss_lineNo = lineNo; 
    prev_yr = -1; 
    while ( std::getline( paramFp, line )) {
      lineNo++;
      std::vector<std::string> sTokens;
      sTokens = boost::split(sTokens, line, boost::is_any_of(", \n\t"));
      if ( sTokens.empty() )
        break;
      try {
        series_yr = std::stoi(sTokens.front());
      }
      catch (const std::out_of_range& oor) {
        std::cout << "Could not read year in series data at line " << lineNo << std::endl;
        logger.Log("Couls not read year in series data at line " + to_string(lineNo));
        exit(EXIT_FAILURE);
        // logAndExit(logfp, outstr); 
      }
      if ( prev_yr < 0 ) {
        prev_yr = series_yr - 1;
        series->start = series_yr;
      }
      if (series_yr - 1 != prev_yr) {
        std::cout << "Series year on line " << lineNo << " is not consecutive." << std::endl;
        logger.Log("Series year on line " + to_string(lineNo) + " is not consecutive.");
        exit(EXIT_FAILURE);
        // sprintf(outstr, "series year on line %d is not consecutive\n", lineNo);
        // logAndExit(logfp, outstr);
      }
      prev_yr = series_yr; 
    }
    series->vlen = series_yr - series->start + 1; 
    
    // Allocate the space and read the series, have already checked the years. 
    paramFp.seekg(place);
    lineNo = ss_lineNo; 
    series->data = new PPPG_VVAL[series->vlen * 12]; 
    for ( int ss = 0; ss < series->vlen; ss++) {
      std::getline( paramFp, line );
      lineNo++; 
      std::vector<std::string> sTokens;
      sTokens = boost::split(sTokens, line, boost::is_any_of(", \n\t")); 
      if (sTokens.size() != 13) {
        std::cout << "Could not parse series format (year and months) on " << lineNo << std::endl;
        logger.Log("Could not parse series format (year and months) on " + to_string(lineNo));
        exit(EXIT_FAILURE);
      }
      try {
        series_yr = std::stoi(sTokens.front());
      }
      catch (std::invalid_argument const& inv) {
          std::cout << "Could not read series year on line " << lineNo << std::endl;
          logger.Log("Could not read series year on line " + to_string(lineNo));
          exit(EXIT_FAILURE);
      }
      // Read the monthly values. 
      int si; 
      for (int mn = 0; mn < 12; mn++) {
        si = ss * 12 + mn;
        tok = pValue.at(mn);
        if ( readParam(series->data[si], tok ) ) {
            if (series->data[si].spType == pScalar)
            {
                string monthVStr = "   " + pName + " year " + to_string(series->start + ss) + " month " + to_string(mn + 1) + " constant: " + to_string(series->data[si].sval);
                std::cout << monthVStr << std::endl;
                logger.Log(monthVStr);
            }
                
            else if (series->data[si].spType == pTif)
            {
                string monthTIFStr = "   " + pName + " year " + to_string(series->start + ss) + " month " + to_string(mn + 1) + " grid: " + series->data[si].gridName;
                std::cout << monthTIFStr << std::endl;
                logger.Log(monthTIFStr);
            }
            
            // fprintf(logfp, "   %-34s  %4d/%02d grid: %s\n", pName, series->start + ss, mn+1, series->data[si].gridName );
        }
        else {
          std::cout << "Could not read parameter " << pName << " at line " << lineNo << std::endl;
          logger.Log("Could not read parameter " + pName + " at line " + to_string(lineNo));
          exit(EXIT_FAILURE);

        }
      }
    }
    series->got = true; 
  }
  return true;
}

//----------------------------------------------------------------------------------

void readSpeciesParamFile(const std::string& speciesFile, DataInput *dataInput) {
    FILE* paramFp;
    std::string line, pName;
    std::string pValue;
    std::string cp;
    int lineNo = 0;
    int len;

    auto isDoubleQuote = [](char c) { return c == '\"'; };
    std::ifstream inFile(speciesFile);
    std::cout << "Reading species parameter from file '" << speciesFile << "'..." << std::endl;
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
        if (dataInput->tryAddParam(pName, pValues)) {
            continue; 
        }
        else {
            std::cout << "Invalid site parameter: " << pName << std::endl;
            logger.Log("Invalid site parameter: " + pName);
            exit(EXIT_FAILURE);
        }
    }
}

std::unordered_map<std::string, PPPG_OP_VAR> readSiteParamFile(const std::string& paramFile, DataInput *dataInput)
{
  // Read a text file containing 3PG parameters.  Comments are allowed
  // and must begin with C++ style '//'.  Comments can begin at any
  // position on a line.  Parameter lines can begin at any position,
  // and must start with the the string 'id' defined in the array
  // 'params' above, in double quotes. Extra data on a parameter line
  // will be ignored.  Missing data on a parameter line will cause a
  // warning.  Lines which are not comments and cannot be identified
  // as parameters will cause the program to exit.  
  FILE *paramFp;
  std::string line, pName;
  std::string pValue;
  std::string cp;
  int lineNo=0;
  int len;
  std::unordered_map<std::string, PPPG_OP_VAR> opVars;

  auto isDoubleQuote = [](char c) { return c == '\"'; };
  std::ifstream inFile(paramFile);
  std::cout << "Reading input parameters from file '" << paramFile << "'..." << std::endl;
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
    boost::split(pValues, tokens.at(1), boost::is_any_of(" \t"), boost::token_compress_on);
    if (dataInput->tryAddParam(pName, pValues)) { continue; }
    if (output_var_names.find(pName) != output_var_names.end()) {
        opVars.emplace(pName, readOutputParam(pName, pValues, lineNo));
    }
    else if (readOtherParam(pName, pValues)) { continue; }
    else if (readInputSeriesParam(pName, pValues, inFile, lineNo)) { continue; }
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

bool haveAllParams()
{ 
  // Check that we have read a value for all parameters. 
  int i, pInd;
  int indSeed, indWRi, indWFi, indWSi;  //indexes for values to check if required variable available
  bool missing = false;

  std::cout << "Checking all required parameters have been set.." << std::endl;
  logger.Log("Checking all required parameters have been set..");

  // Temperature series
  if ((!Tmax_vals.got) && (!userTavgSeries())){
    std::cout << "No Tmax or Tavg data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Tmax or Tavg data");
  }
  if ((!Tmin_vals.got) && (!userTavgSeries())) {
    std::cout << "No Tmin or Tavg data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Tmin or Tavg data");
  }
  if (userTavgSeries()) {
    std::cout << "Using Tavg series" << std::endl;
    if (!userVpdSeries())
      std::cout << "No VPD series but VPD series is required if Tavg series is used" << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, "A VPD series is required if an average temperature series is used\n"); 
  }
  if (!Rain_vals.got) {
    std::cout << "No Rain data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Rain data");
  }
  if (!SolarRad_vals.got && !NetRad_vals.got) {
    std::cout << "No Solar or Net Radiation data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Solar or Net Radiation data");
  }
  if (!FrostDays_vals.got) {
    std::cout << "No Frost data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Frost data");
  }

  // Check various optional parameters
  //if ( !modelMode3PGS && NdviAvh_vals->got) {
  if ( !modelMode3PGS && NdviAvh_vals.got) {
    std::cout << "NdviAvh_vals not used in 3PGS mode." << std::endl;
    logger.Log("NdviAvh_vals not used in 3PGS mode.");
  }
  return !missing;
}

//----------------------------------------------------------------------------------

GDALRasterImage* openInputGrids( )
{
  // Open any grids in the params array, the climate and NDVI series arrays, and 
  // the management tables. 
  // Copy the grid parameters of the first grid opened to refGrid. 
  int pn, j; 
  bool spatial = false, first = true;
  PPPG_SERIES_PARAM *ser; 
  GDALRasterImage *refGrid;

  std::cout << "Opening input rasters..." << std::endl;
  logger.Log("Opening input rasters...");

  // Open all series grids. 
  // For each series. 
  PPPG_SERIES_PARAM *serlist[] = { &Tmax_vals, &Tmin_vals, &Rain_vals, &SolarRad_vals, 
    &FrostDays_vals, &NdviAvh_vals, &NetRad_vals, &Vpd_vals, &Tavg_vals, NULL }; 
  for (j = 0; serlist[j] != NULL; j++) {
    ser = serlist[j];
    for (int i = 0; i < ser->vlen * 12; i++) {
      if ( openGrid( ser->data[i] ) ) {
        spatial = true; 
        if ( first ) {
          refGrid = (GDALRasterImage *)(ser->data[i].g);
          first = false;
        }
        else if ( ( fabs( refGrid->xMin - ser->data[i].g->xMin ) > 0.0001 ) 
          || ( fabs( refGrid->yMin - ser->data[i].g->yMin ) > 0.0001 )
          || ( fabs( refGrid->xMax - ser->data[i].g->xMax ) > 0.0001 )
          || ( fabs( refGrid->yMax - ser->data[i].g->yMax ) > 0.0001 ) 
          || ( refGrid->nRows != ser->data[i].g->nRows ) 
          || ( refGrid->nCols != ser->data[i].g->nCols ) ) {
            std::cout << "Grid dimensions must match, raster " << ser->data[i].gridName << " differs from first raster." << std::endl;
            logger.Log("Grid dimensions must match, raster " + ser->data[i].gridName + " differs from first raster.");
            exit(EXIT_FAILURE);
          // sprintf(outstr, "Grid dimensions must match, grid %s differs from first grid.\n", 
          //   ser->data[i].gridName ); 
          // logAndExit(logfp, outstr); 
        }
      }
    }
  }

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

  if (!spatial) {
    std::cout << "None" << std::endl;
    // fprintf(logfp, "none\n");
    refGrid = NULL; 
  }


  return refGrid;
}

//----------------------------------------------------------------------------------

int writeOutputGrids(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, bool hitNODATA, long cellIndex) {
    //for each possible output variable
    for (auto& [pN, opV] : opVars) {
        //if it has been marked to be written
        if (opV.write) {
            //determine value, name, and tell dataOutput class to write
            float val = (float)(opV.v);
            std::string name = opV.gridName;
            name = name.substr(0, name.find_last_of("."));
            dataOutput->write(-1, -1, name, cellIndex, val, hitNODATA);
        }
    }
    return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------------

void writeMonthlyOutputGrids(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, int calYear, int calMonth, bool hitNODATA, MYDate minMY, MYDate maxMY, long cellIndex) {
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
        dataOutput->write(calYear, calMonth, name, cellIndex, val, hitNODATA);
    }
}

//----------------------------------------------------------------------------------

void writeSampleFiles(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, int cellIndex, int month, long calYear)
{
    int opn, sInd, i;
    static bool firstTime = true;

    // Do we want to sample this point? 
    sInd = -1;
    for (i = 0; samplePoints[i].id[0] != 0; i++)
        if (samplePoints[i].cellIndex == cellIndex) {
            sInd = i;
            break;
        }
    if (sInd == -1)
        return;

    fprintf(samplePoints[sInd].fp, "%d, %d, %s, ", calYear,
        month, samplePoints[sInd].id.c_str());
    // For each variable
    for (auto& [pN, opV] : opVars) {
        fprintf(samplePoints[sInd].fp, "%f, ", (opV.v));
    }
    fprintf(samplePoints[sInd].fp, "\n");
}


//----------------------------------------------------------------------------------

void writeSampleFiles(std::unordered_map<std::string, PPPG_OP_VAR> opVars, int cellIndex, int month, long calYear)
{
  int opn, sInd, i;
  static bool firstTime = true;

  // Do we want to sample this point? 
  sInd = -1;
  for (i = 0; samplePoints[i].id[0] != 0; i++) 
    if (samplePoints[i].cellIndex == cellIndex) {
      sInd = i;
      break;
    }
  if (sInd == -1)
    return;

  fprintf(samplePoints[sInd].fp, "%d, %d, %s, ", calYear, 
    month, samplePoints[sInd].id.c_str());
  // For each variable
  for (auto& [pN, opV] : opVars) {
      fprintf(samplePoints[sInd].fp, "%f, ", (opV.v));
  }
  fprintf(samplePoints[sInd].fp, "\n");
}

//----------------------------------------------------------------------------------

bool userVpdSeries(void)
{
  return (Vpd_vals.data != NULL); 
}

//----------------------------------------------------------------------------------

bool userNetRadSeries(void)
{
  return (NetRad_vals.data != NULL);
}

//----------------------------------------------------------------------------------

bool userTavgSeries(void)
{
  return (Tavg_vals.data != NULL); 
}