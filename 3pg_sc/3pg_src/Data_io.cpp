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
// #include <format>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>

#include "GDALRasterImage.hpp"
#include "util.hpp"
#include "Data_io.hpp"

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

extern FILE *logfp;
char * outstr;

#define GRID_NAME_LENGTH 300
#define PPPG_MAX_SERIES_YEARS 150
#define PPPG_MAX_SERIES_LENGTH PPPG_MAX_SERIES_YEARS*12

// Possible types of parameter - null or not yet set, scalar (constant), 
// ByteGrid for a BIL, and FloatGrid for a floating point file.  The last 
// two are handled by the JG Grid class.  
typedef enum {pNull, pScalar, pTif} ParamSpatial;

typedef struct {
  ParamSpatial spType;                   // Scalar, grid, or null
  double sval;                           // Scalar value. 
  std::string gridName;                        // ptr to grid file name
  GDALRasterImage *g;                               // ptr to grid value
} PPPG_VVAL; 

// 3PG 'parameters'. These are all stored as double.  The 'id' string field 
// is set to the name of the variable in the initialisation below.  Within 
// the model itself we don't reference the parameters via this type.  Its 
// used to help identify parameter lines, and to help get grid values into 
// the model.  
typedef struct {
  std::string id = "-1";                        // String version of the variable name. 
  double *adr;                     // The address of the model variable. 
  bool got;                        // Has the parameter been set? 
  PPPG_VVAL data;                  // Variant value
} PPPG_PARAM; 

// 3PG output variables. In spatial mode output variables may be written 
// repeatedly, on a time step defined by recurStart, recurYear, and recurMonthly. 
typedef struct {
  std::string id;                       // String version of the variable name. 
  double *adr;                    // The address of the model variable. 
  ParamSpatial spType;            // If its a spatial parameter and what grid type. 
  std::string gridName;  // The gridname, in spatial mode. 
  GDALRasterImage *g;                        // The final output grid, in spatial mode. 
  bool write;                     // Whether the variable is wanted. 
  int recurStart;                 // First year to write regular output. 
  int recurYear;                  // Interval on which to write regular output. 
  int recurMonth;                 // Single month number we want output in. 
  bool recurMonthly;              // Whether to write every month in an output year. 
  FILE **RO;                      // The output files for regular output. 
} PPPG_OP_VAR;

// 3PG 'series' parameters.  This is any parameter with a time series for value, 
// in particular the climate parameters, and NDVI.  
typedef struct {
  int start;                             // Calendar year of first entry. 
  PPPG_VVAL *data;                       // Array of variant values.
  int vlen;                              // Number of entries in array. 
  bool oneYear;                          // Array is of a single 'average' year (eg esoclim). 
  bool got;                              // Have read the series. 
} PPPG_SERIES_PARAM; 

// 3PG 'management table' parameters.  Only one value per year is allowed.  
typedef struct {
  int year;                                   // Calendar year
  bool got;
  PPPG_VVAL data; 
} PPPG_MT_PARAM; 

//----------------------------------------------------------------------------------
// Global variables.  These are used to push values into the runTreeModel routine. 
// Where the VB version uses an integer variable the equivalent C declaration 
// has been left in but commented.  Apart from that this section is a direct 
// translation from the The_3PG_Model module in the VB version.  Note that the 
// arrays used directly by the models have been left '1' based, so there is an 
// unused array member at index 0.  

// Controls and counters
//int StartAge, EndAge;                         // age of trees at start/end of run
//int StartMonth;                               // month of year to start run
//int yearPlanted;                              // year trees planted
// ANL changed these three from int to double
extern double StartAge, EndAge;                 // age of trees at start/end of run
extern double StartMonth;                       // month of year to start run
extern double yearPlanted;                      // year trees planted
// int DaysInMonth[13];                         // array for days in months
extern double DaysInMonth[13];                  // array for days in months
//extern bool showDetailedResults;              // TRUE ==> show monthly results
//extern bool showStandSummary;                 // TRUE ==> show stand summary
extern bool modelMode3PGS;

// Site characteristics, site specific parameters
extern std::string siteName[100];                      // name of site
extern double Lat;                              // site latitude
extern double MaxASW, MinASWp;                  // maximum & minimum available soil water
extern double FRp, FR;                              // current site fertility rating
extern double FRstart, FRend, FRdec;            // Start, end and decrement % for fertility decrease with time
extern double soilIndex;                        // soil class index
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
extern double ASW, ASWi;                        // available soil water
extern double MinASWTG;                         // soil water modifier corrector
extern double StemNoi, StemNo;                  // stem numbers
extern double SeedlingMass;
extern double WFi, WF;                          // foliage biomass
extern double WRi, WR;                          // root biomass
extern double WSi, WS;                          // stem biomass
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
extern double MaxAge;
extern double gammaFx, gammaF0, tgammaF;
extern double Rttover;
extern double SLA0, SLA1, tSLA;
extern double fullCanAge;
extern double k;
extern double pFS2, pFS20;
extern double StemConst, StemPower;
extern double SWconst0, SWpower0;
extern double Interception;
extern double BLcond;
extern double MaxCond, CoeffCond;
extern double y;
extern double growthTmax, growthTmin, growthTopt;
extern double thinPower;           //Added 29-07-02
extern double mF, mR, mS;          //Added 29-07-02
extern double wSx1000;
extern double m0, fN0, fNn;
extern double alpha, alphaC;  //alphaC added 11/07/02
extern double pRx, pRn;
extern double nAge, rAge;
extern double kF;
extern double fracBB0, fracBB1, tBB;
extern double fracBB; //fracBB added 11/07/02
extern double Density;
extern double pfsConst, pfsPower;                     // derived from pFS2, pFS20
extern double rhoMin, rhoMax, tRho;             // Standage varying density 3-06-02 
extern double PhysMod;
extern double WUE;                              //Added 16/07/02
extern double CVI;                              //Added 16/07/02
extern double TotalLitter;                      //Added 16/07/02

//Conversion factors
extern double Qa, Qb; 
extern double gDM_mol; 
extern double molPAR_MJ; 

//Additional factors (conductance)
extern double LAIgcx;
extern double MaxIntcptn;
extern double LAImaxIntcptn;

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
extern double NDVI_FPAR_intercept, NDVI_FPAR_constant; 
extern double delWAG;

// ANL - other globals.
bool yearlyOutput, monthlyOutput; 
bool samplePointsYearly = false, samplePointsMonthly = false;
FILE *pointModeFp;   // Output file for point mode only.  
std::string outPath = "./";

//----------------------------------------------------------------------------------

// Initialisation of parameter array. This sets up the mapping between the variable 
// and its name, which is used in parsing the parameter files. 
PPPG_PARAM params[] = 
{
  {"paramError", NULL},
  {"pFS2",         &pFS2},
  {"pFS20",        &pFS20},
  {"StemConst",    &StemConst},
  {"StemPower",    &StemPower},
  {"pRx",          &pRx},
  {"pRn",          &pRn},

  // Temperature modifier (fT) | cardinal temperatures
  // ANL - these have been renamed from just Tmax etc, to avoid confusion with the 
  // climate variables. 
  {"growthTmin",   &growthTmin},
  {"growthTopt",   &growthTopt},
  {"growthTmax",   &growthTmax},
    
  // Frost modifier
  {"kF",           &kF},

  // Litterfall & root turnover
  {"gammaFx",      &gammaFx},
  {"gammaF0",      &gammaF0},
  {"tgammaF",      &tgammaF},
  {"Rttover",      &Rttover},

  // conductances
  {"MaxCond",      &MaxCond},
  {"CoeffCond",    &CoeffCond},
  {"BLcond",       &BLcond},

  // fertility effects
  {"m0",           &m0},
  {"fN0",          &fN0},
  {"fNn",          &fNn},

  //Thinning effects
  {"thinPower",    &thinPower},
  {"mF",           &mF},
  {"mR",           &mR},
  {"mS",           &mS},
  
  // Soil water modifier (fSW) | soil characteristics
  {"SWconst0",     &SWconst0},
  {"SWpower0",     &SWpower0},

  // stem numbers
  {"wSx1000",      &wSx1000},

  // Age modifier (fAge)
  {"MaxAge",       &MaxAge},
  {"nAge",         &nAge},
  {"rAge",         &rAge},

  // Canopy structure and processes | specific leaf area
  {"SLA0",         &SLA0},
  {"SLA1",         &SLA1},
  {"tSLA",         &tSLA},
  {"k",            &k},
  {"fullCanAge",   &fullCanAge},
  {"alpha",        &alpha},
  {"fracBB0",      &fracBB0},
  {"fracBB1",      &fracBB1},
  {"tBB",          &tBB},

  // various
  {"y",            &y},
  {"rhoMin",       &rhoMin},
  {"rhoMax",       &rhoMax},
  {"tRho",         &tRho},             // Standage varying density 3-06-02 
  
  //Conversions
  {"Qa",           &Qa},
  {"Qb",           &Qb},
  {"gDM_mol",      &gDM_mol},
  {"molPAR_MJ",    &molPAR_MJ},

  //Additional conversion factors 
  {"LAIgcx",           &LAIgcx},
  {"MaxIntcptn",       &MaxIntcptn},
  {"LAImaxIntcptn",    &LAImaxIntcptn},

  // 3PG site parameters. 
  {"Lat",          &Lat},
  {"FRp",          &FRp},
  {"FRstart",      &FRstart},  //These three variables relate to fertility decrease with age
  {"FRend",        &FRend},
  {"FRdec",        &FRdec},
  {"soilIndex",    &soilIndex}, 
  {"MaxASW",       &MaxASW},
  {"MinASWp",      &MinASWp},

  // Initial conditions. 
  {"StartAge",     &StartAge},
  {"EndAge",       &EndAge},
  {"StartMonth",   &StartMonth},
  {"yearPlanted",  &yearPlanted},  /* CHECK! do we still use this?*/  
  {"SeedlingMass", &SeedlingMass},
  {"WFi",          &WFi},
  {"WRi",          &WRi},
  {"WSi",          &WSi},
  {"StemNoi",      &StemNoi},
  {"ASWi",         &ASWi},
  {"MinASWTG",     &MinASWTG},
//  {"yearPlanted",  &yearPlanted},  /* This has been moved as strange errors were occuring with grids here*/

  // ANL - extras for 3PGS mode
  {"NDVI_FPAR_intercept", &NDVI_FPAR_intercept}, 
  {"NDVI_FPAR_constant",  &NDVI_FPAR_constant}, 

  {"", NULL}  // NULL entries used to mark array ends. 
};

//----------------------------------------------------------------------------------

// Initialisation of output variable array. This lists all possible output variables 
// and sets up the mapping of the output variable to its name, which is used in 
// parsing the parameter file.  
PPPG_OP_VAR opVars[] = {
  {"opVarError", NULL}, 
  {"StemNo",     &StemNo},
  {"WF",         &WF}, 
  {"WR",         &WR}, 
  {"WS",         &WS}, 
  {"TotalW",     &TotalW}, 
  {"LAI",        &LAI},
  {"cLAI",       &cLAI},
  {"MAI",        &MAI}, 
  {"avDBH",      &avDBH}, 
  {"BasArea",    &BasArea}, 
  {"StandVol",   &StandVol},
  {"GPP",        &GPPdm},
  {"cGPP",       &cGPP},
  {"NPP",        &NPP}, 
  {"cNPP",       &cNPP},
  {"delWAG",     &delWAG}, 
  {"cumWabv",    &cumWabv}, 
  {"Transp",     &Transp},
  {"cTransp",    &cTransp},
  {"ASW",        &ASW},
  {"fSW",        &fSW},
  {"fVPD",       &fVPD},
  {"fT",         &fT}, 
  {"fNutr",      &fNutr}, 
  {"fFrost",     &fAge}, 
  {"APAR",       &APAR}, 
  {"APARu",      &APARu}, 
  {"EvapTransp", &EvapTransp},
  {"cEvapTransp",&cEvapTransp}, //Added 08/11/02
  {"LAIx",       &LAIx},
  {"ageLAIx",    &ageLAIx},
  {"MAIx",       &MAIx},    //Added 29/07/2002
  {"ageMAIx",    &ageMAIx}, //Added 29/07/2002
  {"FR",         &FR},     //Added 11/07/2002
  {"PhysMod",    &PhysMod}, //Added 11/07/2002
  {"alphaC",     &alphaC},  //Added 11/07/2002
  {"fAge",       &fAge},    //Added 11/07/2002
  {"fracBB",     &fracBB},
  {"WUE",        &WUE},     //Added 16/07/02
  {"cWUE",       &cWUE},    //Added 08/11/02
  {"CVI",        &CVI},     //Added 16/07/02
  {"cCVI",       &cCVI},    //Added 08/11/02
  {"TotalLitter", &TotalLitter}, //Added 16/07/02
  {"cLitter",    &cLitter},
  {NULL,         NULL}
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
  std::string id;
  FILE *fp; 
  double x; 
  double y; 
  std::tuple<int, int> cellIndex;
} samplePoints[MAX_SAMPLE_POINTS + 1];

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

int pNameToInd(const std::string& id)
{
  // For a parameter name return its index in the parameter array. 
  // Return 0 on non-existant parameter name. Position zero in 
  // the parameter name array is occupied by an error marker. 
  //
  // We have to take a paranoid approach to comparing the strings, 
  // because we have parameter names which are substrings of other 
  // parameter names.  For example, 'Rain' is a substring of 
  // 'Rainfall Interception Factor (%)'.  strcmp seems to just compare 
  // to the length of the shorter string. 
  int pn;
  int w1, w2;

  w1 = id.length();
  for (pn=0; params[pn].id != "-1"; pn++) {
    w2 = params[pn].id.length();
    if (id.compare(params[pn].id) == 0)
      if (w1 == w2)
        return pn;
  }
  //fprintf(stderr, "Warning, lookup of non-existent parameter name: %s\n", 
  //        id);
  return 0;
}

//----------------------------------------------------------------------------------

int opNameToInd(const std::string& id)
{
  // Exactly the same as pNameToInd, but looking at the opVars array. 
  std::size_t pn;
  std::size_t w1, w2;

  w1 = id.length();
  for (pn=0; opVars[pn].id != "-1"; pn++) {
    w2 = opVars[pn].id.length();
    if (id.compare(opVars[pn].id) == 0)
      if (w1 == w2)
        return pn;
  }
  //fprintf(stderr, "Warning, lookup of non-existent output variable name: %s\n", 
  //        id);
  return 0;
}

//----------------------------------------------------------------------------------

bool getVVal(double &val, PPPG_VVAL vval, std::tuple<int,int> k)
{
  GDALRasterImage *fg;
  float result;
  
  if (vval.spType == pScalar)
    val = vval.sval;
  else if (vval.spType == pTif) {
    fg = (GDALRasterImage *)vval.g;
    result = fg->GetVal(std::get<0>(k), std::get<1>(k));
    if (result == fg->noData)
      return false; 
    val = result;
  }
  else { // pNull
    return false; 
  }
  return true; 
}

//-----------------------------------------------------------------------------

double lookupManageTable( int year, int table, double def, std::tuple<int,int> k )
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
    sprintf(outstr, "Program error: called lookupManageTable with invalid table\n");
    logAndExit(logfp, outstr);
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
      if ( fg->GetVal(std::get<0>(k), std::get<1>(k)) == fg->noData ) 
        continue; 
      else
        val = fg->GetVal(std::get<0>(k), std::get<1>(k));; 
    }
    if ( mt[i].year >= year )
      break;
  }

  return val; 
}

//----------------------------------------------------------------------------------

bool getSeriesVal(double &val, int ser, int calMonth, int calYear, std::tuple<int,int> k)
{
  PPPG_SERIES_PARAM *series;

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
    i = calMonth -1; 
  }
  else {
    // Long run data. 
    i = (calYear - series->start) * 12 + calMonth -1; 
    if (i > series->vlen * 12 - 1) {
      // Should not happen as we will sanity check series before running. 
      // sprintf(outstr.c_str(), "Attempted lookup of series element %d in series %d, only %d entries in series.\nCheck Start/End Ages.",
      //   i, ser, (series->vlen*12-1)); 
      logAndExit(logfp, outstr);
    }
    if (i < 0)
    {
      // sprintf(outstr, "Attempted lookup of year before start year\n");
      logAndExit(logfp, outstr);
    }
  }
  if ( getVVal(val, series->data[i], k) )
    return true;
  else 
    return false; 

}
//----------------------------------------------------------------------------------

void readSampleFile(GDALRasterImage *refGrid)
{
  // Read a text file of sample points, one per line, in the format idstring, 
  // xcoord, ycoord; find the index number of the cell the points fall in. 
  char *line, *fname;
  char *id, *xstr, *ystr, *cp;
  int ind=0;
  int opn;
  double lat, lon;
  std::tuple<int, int> cellIndex;

  if (!samplePointsMonthly && !samplePointsYearly)
    return;

    if ((sampleIpFp = fopen(sampleIpFile.c_str(), "r")) == NULL) {
      // std::cout << outstr << "\n";
      // sprintf(outstr.c_str(), "Could not open sample point file %s\n", sampleIpFile);
      fprintf(logfp, outstr);
      fprintf(stderr, outstr);
      exit(1);
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
    cellIndex = refGrid->XYfrom(lat, lon);
    samplePoints[ind].x = lat;
    samplePoints[ind].y = lon;
    samplePoints[ind].cellIndex = cellIndex;

    // Open output file for this point. Make sure its in outPath. 
    // sprintf(fname, "%s3pg.%s.txt", outPath, id);
    if ((samplePoints[ind].fp = fopen(fname, "w")) == NULL) {
      sprintf(outstr, "Error opening output sample file %s\n", fname);
      logAndExit(logfp, outstr);
    }

    // Write header line for each sample file. 
    // For each output sample file
    fprintf(samplePoints[ind].fp, "year, month, id, ");
    for (opn = 1; opVars[opn].adr != NULL; opn++)
    //for (opn = 1; opVars[opn].id != "-1"; opn++)
      fprintf(samplePoints[ind].fp, "%s, ", opVars[opn].id);
    fprintf(samplePoints[ind].fp, "\n");

    ind++;
  }
  // Make sure end of sample is marked. 
  samplePoints[ind].id[0] = 0;
}

//----------------------------------------------------------------------------------

bool readInputParam(const std::string& pName, std::string pValue)
{
  // For a parameter name pName and a parameter value, pValue, both as strings, read 
  // the value into an appropriate variable. The parameter name can be either the 
  // same as the variable name, or it can be a long descriptive name, which matches 
  // the description of the parameter given in the VB version. In some cases other 
  // synonymns are allowed also. If pName does not match with any of the defined 
  // parameter names, return false. 
  int pInd=0;
  const std::string cp;

  // Find the index within parameters array of n.  
  // Allometric relationships & partitioning
  if (namesMatch("pFS2", pName) ||
           namesMatch("Foliage:stem partitioning ratio @ D=2 cm", pName)) pInd = pNameToInd("pFS2");
  else if (namesMatch("pFS20", pName) ||
           namesMatch("Foliage:stem partitioning ratio @ D=20 cm", pName)) pInd = pNameToInd("pFS20");
  else if (namesMatch("StemConst", pName) ||
           namesMatch("Constant in the stem mass v. diam. relationship", pName)) 
    pInd = pNameToInd("StemConst");
  else if (namesMatch("StemPower", pName) ||
           namesMatch("Power in the stem mass v. diam. relationship", pName)) pInd = pNameToInd("StemPower");
  else if (namesMatch("pRx", pName) ||
           namesMatch("Maximum fraction of NPP to roots", pName)) pInd = pNameToInd("pRx");
  else if (namesMatch("pRn", pName) ||
           namesMatch("Minimum fraction of NPP to roots", pName)) pInd = pNameToInd("pRn");

  // Temperature modifier (fT) | cardinal temperatures
  else if (namesMatch("growthTmin", pName) ||
           namesMatch("Minimum temperature for growth", pName)) pInd = pNameToInd("growthTmin");
  else if (namesMatch("growthTopt", pName) ||
           namesMatch("Optimum temperature for growth", pName)) pInd = pNameToInd("growthTopt");
  else if (namesMatch("growthTmax", pName) ||
           namesMatch("Maximum temperature for growth", pName)) pInd = pNameToInd("growthTmax");
    
  // Frost modifier
  else if (namesMatch("kF", pName) ||
           namesMatch("Days production lost per frost day", pName)) pInd = pNameToInd("kF");

  // Litterfall & root turnover
  else if (namesMatch("gammaFx", pName) ||
           namesMatch("Maximum litterfall rate", pName)) pInd = pNameToInd("gammaFx");
  else if (namesMatch("gammaF0", pName) ||
           namesMatch("Litterfall rate at t = 0", pName)) pInd = pNameToInd("gammaF0");
  else if (namesMatch("tgammaF", pName) ||
           namesMatch("Age at which litterfall rate has median value", pName)) 
    pInd = pNameToInd("tgammaF");
  else if (namesMatch("Rttover", pName) ||
           namesMatch("Average monthly root turnover rate", pName)) pInd = pNameToInd("Rttover");

  // conductances
  else if (namesMatch("MaxCond", pName) ||
           namesMatch("Maximum canopy conductance", pName)) pInd = pNameToInd("MaxCond");
  else if (namesMatch("CoeffCond", pName) ||
           namesMatch("Defines stomatal response to VPD", pName)) pInd = pNameToInd("CoeffCond");
  else if (namesMatch("BLcond", pName) ||
           namesMatch("Canopy boundary layer conductance", pName)) pInd = pNameToInd("BLcond");

  // fertility effects
  else if (namesMatch("m0", pName) ||
           namesMatch("Value of 'm' when FR = 0", pName)) pInd = pNameToInd("m0");
  else if (namesMatch("fN0", pName) ||
           namesMatch("Value of 'fNutr' when FR = 0", pName)) pInd = pNameToInd("fN0");
  else if (namesMatch("fNn", pName) ||
           namesMatch("Power of (1-FR) in 'fNutr'", pName)) pInd = pNameToInd("fNn");  //added 22-07-02
  
  // Soil water modifier (fSW) | soil characteristics
  else if (namesMatch("SWconst0", pName) ||
           namesMatch("Moisture ratio deficit for fq = 0.5", pName)) pInd = pNameToInd("SWconst0");
  else if (namesMatch("SWpower0", pName) ||
           namesMatch("Power of moisture ratio deficit", pName)) pInd = pNameToInd("SWpower0");

  // stem numbers
  else if (namesMatch("wSx1000", pName) ||
           namesMatch("Max. stem mass per tree @ 1000 trees/hectare", pName)) 
    pInd = pNameToInd("wSx1000");

  //Thinning Parameters 29-07-02
  else if (namesMatch("thinPower", pName) ||
           namesMatch("Power in self-thinning rule", pName)) pInd = pNameToInd("thinPower");
  else if (namesMatch("mF", pName) ||
           namesMatch("Fraction mean single-tree foliage biomass lost per dead tree", pName)) pInd = pNameToInd("mF");
  else if (namesMatch("mR", pName) ||
           namesMatch("Fraction mean single-tree root biomass lost per dead tree", pName)) pInd = pNameToInd("mR");
  else if (namesMatch("mS", pName) ||
           namesMatch("Fraction mean single-tree stem biomass lost per dead tree", pName)) pInd = pNameToInd("mS");
           
  // Age modifier (fAge)
  else if (namesMatch("MaxAge", pName) ||
           namesMatch("Maximum stand age used in age modifier", pName)) pInd = pNameToInd("MaxAge");
  else if (namesMatch("nAge", pName) ||
           namesMatch("Power of relative age in function for fAge", pName)) pInd = pNameToInd("nAge");
  else if (namesMatch("rAge", pName) ||
           namesMatch("Relative age to give fAge = 0.5", pName)) pInd = pNameToInd("rAge");

  // Canopy structure and processes | specific leaf area
  else if (namesMatch("SLA0", pName) ||
           namesMatch("Specific leaf area at age 0", pName)) pInd = pNameToInd("SLA0");
  else if (namesMatch("SLA1", pName) ||
           namesMatch("Specific leaf area for mature leaves", pName)) pInd = pNameToInd("SLA1");
  else if (namesMatch("tSLA", pName) ||
           namesMatch("Age at which specific leaf area = (SLA0+SLA1)/2", pName)) pInd = pNameToInd("tSLA");
  else if (namesMatch("k", pName) ||
           namesMatch("Extinction coefficient for absorption of PAR by canopy", pName)) pInd = pNameToInd("k");
  else if (namesMatch("fullCanAge", pName) ||
           namesMatch("Age at canopy cover", pName)) pInd = pNameToInd("fullCanAge");
  else if (namesMatch("alpha", pName) ||
           namesMatch("Canopy quantum efficiency", pName)) pInd = pNameToInd("alpha");

  // Branch and bark fraction (fracBB)
  else if (namesMatch("fracBB0", pName) ||
           namesMatch("Branch and bark fraction at age 0", pName)) pInd = pNameToInd("fracBB0");
  else if (namesMatch("fracBB1", pName) ||
           namesMatch("Branch and bark fraction for mature stands", pName)) pInd = pNameToInd("fracBB1");
  else if (namesMatch("tBB", pName) ||
           namesMatch("Age at which fracBB = (fracBB0+fracBB1)/2", pName)) pInd = pNameToInd("tBB");

  // various
  else if (namesMatch("y", pName) ||
           namesMatch("Ratio NPP/GPP", pName)) pInd = pNameToInd("y");
  else if (namesMatch("Density", pName) ||
           namesMatch("Basic density", pName)) pInd = pNameToInd("Density");

  //conversion factors - Addition 29th November Anders Siggins
  else if (namesMatch("Qa", pName) ||
           namesMatch("Intercept of net v. solar radiation relationship", pName)) pInd = pNameToInd("Qa");
  else if (namesMatch("Qb", pName) ||
           namesMatch("Slope of net v. solar radiation relationship", pName)) pInd = pNameToInd("Qb");
  else if (namesMatch("gDM_mol", pName) ||
           namesMatch("Molecular weight of dry matter", pName)) pInd = pNameToInd("gDM_mol");
  else if (namesMatch("molPAR_MJ", pName) ||
           namesMatch("Conversion of solar radiation to PAR", pName)) pInd = pNameToInd("molPAR_MJ");

  //Additional conversion factors
  else if (namesMatch("LAIgcx", pName) ||
           namesMatch("LAI for maximum canopy conductance", pName)) pInd = pNameToInd("LAIgcx");
  else if (namesMatch("MaxIntcptn", pName) ||
           namesMatch("Maximum proportion of rainfall evaporated from canopy", pName)) pInd = pNameToInd("MaxIntcptn");
  else if (namesMatch("LAImaxIntcptn", pName) ||
           namesMatch("LAI for maximum rainfall interception", pName)) pInd = pNameToInd("LAImaxIntcptn");

  // 3PG site parameters. 
  else if (namesMatch("Lat", pName) ||
           namesMatch("Latitude", pName)) pInd = pNameToInd("Lat");
  else if (namesMatch("FR", pName) ||
           namesMatch("Fertility rating", pName)) pInd = pNameToInd("FRp");
  else if (namesMatch("soilIndex", pName) ||
           namesMatch("Soil Index", pName) ||
           namesMatch("soil class", pName)) pInd = pNameToInd("soilIndex");
  else if (namesMatch("MaxASW", pName) ||
           namesMatch("Maximum ASW", pName)) pInd = pNameToInd("MaxASW");
  else if (namesMatch("MinASW", pName) ||
           namesMatch("Minimum ASW", pName)) pInd = pNameToInd("MinASWp");

  // Initial conditions. 
  else if (namesMatch("StartAge", pName) || namesMatch("Initial age", pName) || 
           namesMatch("Start age", pName)) pInd = pNameToInd("StartAge");
  else if (namesMatch("EndAge", pName) ||
           namesMatch("End age", pName)) pInd = pNameToInd("EndAge");
  else if (namesMatch("StartMonth", pName) || namesMatch("Start Month", pName) || 
           namesMatch("Start month", pName)) pInd = pNameToInd("StartMonth");
  else if (namesMatch("SeedlingMass", pName) || namesMatch("Seedling Mass", pName) || 
           namesMatch("Seedling mass", pName)) pInd = pNameToInd("SeedlingMass");
  else if (namesMatch("FRstart", pName)) pInd = pNameToInd("FRstart");        //Fertility modifiers that depend on age
  else if (namesMatch("FRend", pName)) pInd = pNameToInd("FRend");
  else if (namesMatch("FRdec", pName)) pInd = pNameToInd("FRdec");
  else if (namesMatch("WFi", pName) ||
           namesMatch("W foliage", pName)) pInd = pNameToInd("WFi");
  else if (namesMatch("WRi", pName) ||
           namesMatch("W root", pName)) pInd = pNameToInd("WRi");
  else if (namesMatch("WSi", pName) ||
           namesMatch("W stem", pName)) pInd = pNameToInd("WSi");
  else if (namesMatch("StemNoi", pName) ||
           namesMatch("Stem no", pName)) pInd = pNameToInd("StemNoi");
  else if (namesMatch("ASWi", pName) ||
           namesMatch("Initial soil water", pName)) pInd = pNameToInd("ASWi");
  else if (namesMatch("MinASWTG", pName)) pInd = pNameToInd("MinASWTG");
  else if (namesMatch("rhoMin", pName) ||
           namesMatch("Minimum basic density - for young trees", pName)) pInd = pNameToInd("rhoMin");  //Standage varying density 15/07/2002
  else if (namesMatch("rhoMax", pName) ||   //Standage varying density 15/07/2002
           namesMatch("Maximum basic density - for older trees", pName)) pInd = pNameToInd("rhoMax");
  else if (namesMatch("tRho", pName)   ||
           namesMatch("Age at which rho = (rhoMin+rhoMax)/2", pName)) pInd = pNameToInd("tRho");      //Standage varying density 15/07/2002 
  else if (namesMatch("yearPlanted", pName) ||
           namesMatch("Year Planted", pName)) pInd = pNameToInd("yearPlanted");

  // 3PGS mode
  else if (namesMatch("NDVI_FPAR_intercept", pName)) pInd = pNameToInd("NDVI_FPAR_intercept");
  else if (namesMatch("NDVI_FPAR_constant", pName)) pInd = pNameToInd("NDVI_FPAR_constant");
    
  // If no index was found its not a basic input parameter. 
  if (pInd == 0)
    return false;

  // Soil Index/Class is a special case, if its specified with a character code, rewrite 
  // it as an integer.  Must be carefull not to match a grid name.  
  if ( pInd == pNameToInd( "soilIndex" ) ) {
      if (namesMatch("S", pValue))
          pValue = "1";
      else if (namesMatch("SL", pValue))
          pValue = "2";
      else if (namesMatch("CL", pValue))
          pValue = "3";
      else if (namesMatch("C", pValue))
          pValue = "4";
  }
  
  // Is the parameter a constant value (a float).
  if (sscanf(pValue.c_str(), "%lf", params[pInd].adr) == 1) {
    fprintf(logfp, "   %-40s constant:  % 9.3f\n", params[pInd].id, 
      *(params[pInd].adr)); 
    params[pInd].data.spType = pScalar;
    params[pInd].got = 1;
  }
  else {
    // Is the parameter a grid name (a string). 
    // REFERENCES: https://stackoverflow.com/questions/43114174/convert-a-string-to-std-filesystem-path
    // and https://stackoverflow.com/questions/51949/
    params[pInd].data.gridName = pValue;
    const std::filesystem::path filePath = params[pInd].data.gridName;
    if (filePath.extension() == ".tif") // Heed the dot.
    {
        std::cout << filePath.stem() << " is a valid typ  e.";
        params[pInd].data.spType = pTif;
    }
    else
    {
        std::cout << filePath.filename() << " is an invalid type. File extension must be '.tif'";
        logAndExit(logfp, outstr);
        // Output: "myFile.cfg is an invalid type"
    }
    fprintf(logfp, "   %-40s grid: %s\n", params[pInd].id, 
      params[pInd].data.gridName);
    params[pInd].got = 1;
  }
  return true;
}

//----------------------------------------------------------------------------------

bool readOutputParam(const std::string& pName, const std::string& pValue, int lineNo)
{
  // For a parameter name pName and a parameter value, pValue, both as strings, read 
  // the value into an appropriate variable. The parameter name can be either the 
  // same as the variable name, or it can be a long descriptive name.  
  int pInd, pInd1, pInd2;
  std::string cp;

  pInd = pInd1 = pInd2 = 0;

  // Find the index for the param in the opVars array. 
  // 3PG only output variables. 
  if (namesMatch("StemNo", pName) ||
    namesMatch("Stocking density", pName)) pInd1 = opNameToInd("StemNo");
  else if (namesMatch("WF", pName) ||
    namesMatch("Weight of foliage", pName)) pInd1 = opNameToInd("WF");
  else if (namesMatch("WR", pName) ||
    namesMatch("Weight of roots", pName)) pInd1 = opNameToInd("WR");
  else if (namesMatch("WS", pName) ||
    namesMatch("Weight of stems", pName)) pInd1 = opNameToInd("WS");
  else if (namesMatch("TotalW", pName) ||
      namesMatch("Total weight", pName)) pInd1 = opNameToInd("TotalW");
  else if (namesMatch("MAI", pName) ||
    namesMatch("Mean Annual Increment", pName)) pInd1 = opNameToInd("MAI");
  else if (namesMatch("avDBH", pName) ||
    namesMatch("Average DBH", pName)) pInd1 = opNameToInd("avDBH");
  else if (namesMatch("BasArea", pName) ||
    namesMatch("Basal Area", pName)) pInd1 = opNameToInd("BasArea");
  else if (namesMatch("StandVol", pName) ||
    namesMatch("Stand volume", pName)) pInd1 = opNameToInd("StandVol");
  else if (namesMatch("Transp", pName) ||
           namesMatch("Transpiration", pName)) pInd1 = opNameToInd("Transp");
  else if (namesMatch("cTransp", pName))  pInd1 = opNameToInd("cTransp");
  else if (namesMatch("ASW", pName) ||
           namesMatch("Available Soil Water", pName)) pInd1 = opNameToInd("ASW");
  else if (namesMatch("EvapTransp", pName)) pInd1 = opNameToInd("EvapTransp");
  else if (namesMatch("cEvapTransp", pName)) pInd1 = opNameToInd("cEvapTransp");
  else if (namesMatch("LAIx", pName)) pInd1 = opNameToInd("LAIx");
  else if (namesMatch("ageLAIx", pName)) pInd1 = opNameToInd("ageLAIx"); 
  else if (namesMatch("GPP", pName) ||
    namesMatch("Gross Primary Production (tDM/ha)", pName)) pInd1 = opNameToInd("GPP");
  else if (namesMatch("cGPP", pName))  pInd1 = opNameToInd("cGPP");
  else if (namesMatch("TotalLitter", pName)) pInd1 = opNameToInd("TotalLitter");
  else if (namesMatch("cLitter", pName)) pInd1 = opNameToInd("cLitter");
  
  
  // 3PGS only output variables
  if (namesMatch("delWAG", pName) ||
    namesMatch("change in aboveground biomass (tDM/ha)", pName)) pInd2 = opNameToInd("delWAG");
  else if (namesMatch("cumWabv", pName) ||
    namesMatch("accumulated aboveground biomass (tDM/ha)", pName)) pInd2 = opNameToInd("cumWabv");

  // Check we only matched from the appropriate set. 
  if (modelMode3PGS) {
    if (pInd1 != 0) {
      sprintf(outstr, "The output variable %s is not supported in 3PGS mode\n", pName);
      logAndExit(logfp, outstr);
    }
  }
  else {
    if (pInd2 != 0) {
      sprintf(outstr, "The output variable %s is not supported in 3PG mode\n", pName);
      logAndExit(logfp, outstr);
    }
  }
  if (modelMode3PGS)
    pInd = pInd2; 
  else
    pInd = pInd1; 

  // Output variables common to both modes. 
  if (namesMatch("NPP", pName) ||
    namesMatch("Net Primary Production (tDM/ha)", pName)) pInd = opNameToInd("NPP"); //Modified 26/07/02
  else if (namesMatch("cNPP", pName))  pInd = opNameToInd("cNPP");
  else if (namesMatch("LAI", pName) ||
    namesMatch("Leaf Area Index", pName)) pInd = opNameToInd("LAI");
  else if (namesMatch("cLAI", pName)) pInd = opNameToInd("cLAI"); //Added 7 November 2002
  else if (namesMatch("FRout", pName))    pInd = opNameToInd("FR"); //Added 11/07/2002
  else if (namesMatch("PhysMod", pName))  pInd = opNameToInd("PhysMod");//Added 11/07/2002
  else if (namesMatch("alphaC", pName))  pInd = opNameToInd("alphaC");//Added 11/07/2002
  else if (namesMatch("fAge", pName))  pInd = opNameToInd("fAge"); //Added 11/07/2002
  else if (namesMatch("fracBB", pName))  pInd = opNameToInd("fracBB"); //Added 11/07/2002
  else if (namesMatch("WUE", pName))  pInd = opNameToInd("WUE"); //Added 16/07/2002
  else if (namesMatch("cWUE", pName))  pInd = opNameToInd("cWUE");
  else if (namesMatch("CVI", pName))  pInd = opNameToInd("CVI"); //Added 16/07/2002
  else if (namesMatch("cCVI", pName))  pInd = opNameToInd("cCVI");
  else if (namesMatch("fSW", pName))      pInd = opNameToInd("fSW");
  else if (namesMatch("fVPD", pName))     pInd = opNameToInd("fVPD");
  else if (namesMatch("fT", pName))       pInd = opNameToInd("fT"); 
  else if (namesMatch("fNutr", pName))    pInd = opNameToInd("fNutr"); 
  else if (namesMatch("fFrost", pName))   pInd = opNameToInd("fFrost"); 
  else if (namesMatch("APAR", pName))     pInd = opNameToInd("APAR"); 
  else if (namesMatch("APARu", pName))    pInd = opNameToInd("APARu");
  else if (namesMatch("MAIx", pName))    pInd = opNameToInd("MAIx");
  else if (namesMatch("ageMAIx", pName))    pInd = opNameToInd("ageMAIx");

  // Did we match a name?  
  if (pInd == 0)
    return false;

  // First token in the pValue is the output grid name. 
  if (pValue.empty()) {
    sprintf(outstr, "Error: can't read grid name on line %u\n");
    logAndExit(logfp, outstr);
  }

  std::vector<std::string> pValToks;
  // TODO: test if this delimiter match is needed, or if \t or \n would suffice
  boost::split(pValToks, pValue, boost::is_any_of("\t\n\015"));
  if (pValToks.size() > 4) {
    std::cout << "More than 4 elements detected in string " << pValue;
    logAndExit(logfp, outstr);
  }
  // Set the output grid name
  boost::trim(outPath); // trim leading and trailing whitespaces
  opVars[pInd].gridName = outPath + pValToks.front(); // concat

  const std::filesystem::path filePath = opVars[pInd].gridName;
  if (filePath.extension() == ".tif") // Heed the dot.
  {
      std::cout << filePath.stem() << " is a valid type.";
      params[pInd].data.spType = pTif;
  }
  else
  {
      std::cout << filePath.filename() << " is an invalid filename. File extension must be '.tif'";
      logAndExit(logfp, outstr);
      // Output: "myFile.cfg is an invalid type"
  }

  // Optional second, and third tokens are recurring output start year, 
  // and recurral interval in years.  Fourth token is keyword 'monthly', or keyword 'month'.  
  // These are required as a set - ie tokens 2, 3 and 4 must be provided or none.  
  try {
    cp = pValToks.at(1);
    yearlyOutput = true;
  }
  catch (const std::out_of_range& oor) {
    std::cout << "No recurring year output detected.";
  }

  if (yearlyOutput == true) {
    // Look for start year
    if (sscanf(cp, "%d", &opVars[pInd].recurStart) != 1) {
      sprintf(outstr, "Expected an integer start year in recuring output specification on line %d\n", lineNo);
      logAndExit(logfp, outstr);
    }
    // Look for recurral interval
    cp = strtok(NULL, " \n\t\015");
    if (cp != NULL) {
      if (sscanf(cp, "%d", &opVars[pInd].recurYear) != 1) {
        sprintf(outstr, "Expected an integer interval in recuring output specification on line %d\n", lineNo); 
        logAndExit(logfp, outstr);
      }
    }
    else {
      sprintf(outstr, "Found start year but no interval in recuring output specification on line %d\n", lineNo); 
      logAndExit(logfp, outstr);
    }
    // Check non-zero interval
    if (opVars[pInd].recurYear == 0) {
      sprintf(outstr, "Found interval of zero years in recuring output specification on line %d\n", lineNo); 
      logAndExit(logfp, outstr);
    }
    // Look for 'monthly' or 'month' keywords. 
    cp = strtok(NULL, " \n\t\015");
    if (cp == NULL) {
      sprintf( outstr, "Missing either 'monthly' or 'month' keyword on line %d\n", lineNo ); 
      logAndExit( logfp, outstr ); 
    }
    else {
      if (strcmp("monthly", cp) == 0)
        opVars[pInd].recurMonthly = true;
      else if ( strcmp( "month", cp ) == 0 ) {
        opVars[pInd].recurMonthly = false; 
        // Look for the single month of output. 
        if ( ( cp = strtok(NULL, " \n\t\015") ) == NULL ) {
          sprintf( outstr, "No month number after 'month' keyword on line %d.\n", lineNo ); 
          logAndExit( logfp, outstr ); 
        }
        else {
          opVars[pInd].recurMonth = atoi( cp ); 
          if ( opVars[pInd].recurMonth == 0 ) {
            sprintf( outstr, "Bad month number after 'month' keyword on line %d.\n", lineNo );
            logAndExit( logfp, outstr ); 
          }
        }
      }
      else {
        sprintf(outstr, "Unrecognised keyword \"%s\" on line %u\n", cp, lineNo); 
        logAndExit(logfp, outstr);
      }
      monthlyOutput = true;
    }
  }

  // Mark the variable for later writing
  opVars[pInd].write = true;
  fprintf(logfp, "   variable: %-20s   grid: %-20s\n", opVars[pInd].id, pValue);
  if (opVars[pInd].recurStart) {
    fprintf(logfp, "      starting in %4d, writing every %2d years", 
      opVars[pInd].recurStart, opVars[pInd].recurYear);
    if (opVars[pInd].recurMonthly)
      fprintf(logfp, ", with monthly values");
    else if ( opVars[pInd].recurMonth != 0 )
      fprintf( logfp, ", on the %2d month", opVars[pInd].recurMonth ); 
    fprintf(logfp, ".\n");
  }
  return true;
}

//----------------------------------------------------------------------------------

bool readOtherParam(const std::string& pName, std::string pValue)
{
  // Set various miscellaneous parameters. 
  std::string cp;

  // Output directory. Allow no directory to be specified, in which case 
  // we force the current directory.  HOW CAN WE CHECK EXISTENCE!!!
  if (namesMatch("Output directory", pName)) {
    cp = strtok(pValue, " \t\n,\015");
    if (cp == NULL) 
      outPath[0] = '.';
    else
      strcpy(outPath, cp);
    int len = strlen(outPath);
    // Make sure of the trailing /
    if (strrchr(outPath, '/') != outPath+len-1)
      outPath[len] = '/';
    fprintf(logfp, "   output path: %s\n", outPath);
    return true;
  }

  // Sample points file. If the parameter is present it must have a value. 
  else if (namesMatch("sample points file", pName)) {
    cp = strtok(pValue, " \t\n,\015");
    if ( cp == NULL ) {
      sprintf(outstr, "Missing value for parameter \"%s\"\n", pName);
      logAndPrint(logfp, outstr);
      return false; 
    }
    strcpy(sampleIpFile, cp);
    // Optional keyword 'monthly'. 
    cp = strtok(NULL, " \t\n,\015");
    if (cp != NULL) {
      // strcpyTrim(str, cp);
      if (namesMatch("monthly", cp))
        samplePointsMonthly = true;
    } else
      samplePointsYearly = true;
    return true;
  }

  // Model mode (Standard 3PG or 3PGS)
  else if (namesMatch("Model mode", pName)) {
    if (strncasecmp("3PGS", pValue, 4) == 0)
      modelMode3PGS = true;
    else if (strncasecmp("3PG", pValue, 3) == 0) 
      modelMode3PGS = false;
    else {
      sprintf(outstr, "Invalid value for parameter 'Model mode': %s\n", pValue);
      logAndExit(logfp, outstr);
    }
    return true;
  }

  // Point mode output file. 
  else if (namesMatch("point mode output file", pName)) {
    cp = strtok(pValue, " \t\n,\015");
    if ( cp == NULL ) {
      sprintf(outstr, "Missing value for parameter \"%s\"\n", pName);
      logAndPrint(logfp, outstr);
      return false; 
    }
    if ((pointModeFp = fopen(cp, "w")) == NULL) {
      fprintf(stderr, "Could not open point mode output file %s\n", cp);
      exit (1);
    }
    return true;
  }
  else
    return false;
}

//----------------------------------------------------------------------------------

bool readInputManageParam(std::string *pName, FILE *paramFp, int &lineNo)
{
  // Read management table input parameters.  A table must begin on the line 
  // following the keyword identifying it.  The table has one entry per line, 
  // each entry consists of a year and a value, seperated by whitespace.  A 
  // blank line terminates the table.  Each value can be either a constant or 
  // a grid name.  
  std::string line;
  std::string tok, cp;
  PPPG_MT_PARAM *tab;
  int i, *nRead; 
  std::string tabName[30];

  // Are we reading a managment table?
  if ( namesMatch( "Management: fertility", pName ) ) {
    tab = FertMT;
    strcpy(tabName, "Fertility MT");
    nRead = &nFertility; 
  }
  else if ( namesMatch( "Management: irrigation", pName ) ) {
    tab = IrrigMT; 
    strcpy(tabName, "Irrigation MT");
    nRead = &nIrrigation; 
  }
  else if ( namesMatch( "Management: MinASW", pName ) ) {
    tab = MinAswMT;
    strcpy(tabName, "Min ASW MT");
    nRead = &nMinAvailSW;
  }
  else
    return false; 

  // Read the table
  i = 0; 
  while (fgets(line, MAXLINE, paramFp) != NULL) {
    lineNo++;
    
    // Tokenise the line, recognise space, comma, tab and newline.  
    tok = strtok(line, " ,\t\n\015");
    
    // Look for blank line to terminate series. 
    if (tok == NULL) {
      tab[i].year = -1; 
      break;
    }

    // Read the year
    if (sscanf(tok.c_str(), "%d", &tab[i].year) != 1) {
      sprintf(outstr, "Could not read year in management table at line %d\n", lineNo);
      fprintf(logfp, outstr);
      logAndExit(logfp, outstr);
    }

    tok = strtok(NULL, " ,\t\n\015"); 
    if (tok == NULL)
      return false; 

    // Is the table value constant (a float).
    if (sscanf(tok, "%lf", &tab[i].data.sval) == 1) {
      tab[i].data.spType = pScalar;
      tab[i].got = 1;
      fprintf(logfp, "   %-31s year %4d constant: %f\n", tabName, tab[i].year, tab[i].data.sval);
    }
    else {
      // Is the parameter a grid name (a string). 
      // tab[i].data.gridName = new char[GRID_NAME_LENGTH + 1];
      tab[i].data.gridName = new std::string[GRID_NAME_LENGTH + 1];
      strcpyTrim(tab[i].data.gridName, tok);
      // Type of grid. 
      cp = strrchr(tab[i].data.gridName, '.');
      if (strncasecmp(cp, ".tif", 4) == 0) 
        tab[i].data.spType = pTif;
      else {
        sprintf(outstr, "Could not determine grid type for input grid named: %s\n"
                "File extension must be '.bil' or '.flt'.\n", tab[i].data.gridName); 
        logAndExit(logfp, outstr);
      }
      fprintf(logfp, "   %-31s year %4d grid: %s\n", tabName, tab[i].year,
              tab[i].data.gridName);
    }
    i++; 
  }
  *nRead = i; 
  return true;  
}

//----------------------------------------------------------------------------------

bool readParam( PPPG_VVAL &vval, std::string *pValue )
{
  std::string *cp, cv[1000];
  double dv; 

  // Is the parameter a constant value (a float).
  if (sscanf(pValue, "%lf", &dv) == 1) {
    vval.spType = pScalar;
    vval.sval = dv; 
    return true;
  }
  else {
    // Is the parameter a grid name (a string). 
    strcpyTrim(cv, pValue);
    
    // Type of grid. 
    cp = strrchr(cv, '.');
    if (strncasecmp(cp, ".tif", 4) == 0) 
      vval.spType = pTif;
    else {
      sprintf(outstr, "Could not determine grid type for input grid named: %s\n", 
        "File extension must be '.bil' or '.flt'.\n", cv);
      logAndPrint(logfp, outstr);
      return false; 
    }
    vval.gridName = new std::string;
    if (vval.gridName == NULL) {
      sprintf(outstr, "Could not get memory for grid name in readParam.\n");
      logAndExit(logfp, outstr);
    }
    strcpy(vval.gridName, cv);
    return true; 
  }
}

//----------------------------------------------------------------------------------

bool readInputSeriesParam(std::string pName, std::string *pValue, FILE *paramFp, int &lineNo)
{
  // Read 'series' input parameters, ie climate and NDVI. 
  // Two styles of input are permitted.  Firstly, the parameter name can be 
  // followed on the same line by 12 values, in this case the 12 values will 
  // be reused for each run year.  Secondly, the parameter name can be the only 
  // thing on the line (other than a comment), on each following line must be a 
  // year followed by 12 values, until the sequence is terminated by a blank 
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
  tok = strtok(pValue, " \t\n,\015");
  if (tok != NULL) 
    // First style described above. 
    series->oneYear = true;
  else 
    // Second style above.  
    series->oneYear = false; 
  
  // Read values for one year style. 
  if ( series->oneYear ) {
    series->vlen = 1; 
    series->data = new PPPG_VVAL[series->vlen * 12]; 
    int i; 
    for ( i = 0; i < 12; i++ ) {
      if ( readParam( series->data[i], tok ) ) {
        if ( series->data[i].spType == pScalar )
          fprintf(logfp, "   %-34s month %2d constant: %12.6f\n", pName, i+1, series->data[i].sval );
        else if (series->data[i].spType == pTif)
          fprintf(logfp, "   %-34s month %2d grid: %s\n", pName, i+1, series->data[i].gridName );
      }
      else {
        sprintf( outstr, "Could not read parameter %s.\n", pName ); 
        logAndExit( logfp, outstr ); 
      }
      tok = strtok( NULL, " \t\n,\015" ); 
    }
    series->got = true; 

    if (i < 12) {
      sprintf(outstr, "Incomplete series on line %d.\n", lineNo); 
      logAndExit(logfp, outstr); 
    }
  }
  // Time series style
  else {
    // Find out how many years in the series. 
    unsigned long fpos = ftell (paramFp); 
    int ss_lineNo = lineNo; 
    prev_yr = -1; 
    while ( fgets( line, MAXLINE, paramFp ) != NULL ) {
      lineNo++; 
      tok = strtok(line, " ,\t\n\015"); 
      if ( tok == NULL )
        break;
      if ( sscanf( tok, "%d", &series_yr ) != 1 ) {
        sprintf(outstr, "Could not read year in series data at line %d\n", lineNo);
        fprintf(logfp, outstr);
        logAndExit(logfp, outstr);
      }
      if ( prev_yr < 0 ) {
        prev_yr = series_yr - 1;
        series->start = series_yr;
      }
      if (series_yr - 1 != prev_yr) {
        sprintf(outstr, "series year on line %d is not consecutive\n", lineNo);
        logAndExit(logfp, outstr);
      }
      prev_yr = series_yr; 
    }
    series->vlen = series_yr - series->start + 1; 
    
    // Allocate the space and read the series, have already checked the years. 
    fseek( paramFp, fpos, SEEK_SET ); 
    lineNo = ss_lineNo; 
    series->data = new PPPG_VVAL[series->vlen * 12]; 
    for ( int ss = 0; ss < series->vlen; ss++) {
      fgets( line, MAXLINE, paramFp );
      lineNo++; 
      tok = strtok(line, " ,\t\n\015"); 
      sscanf( tok, "%d", &series_yr );
      // Read the monthly values. 
      int si; 
      for (int mn = 0; mn < 12; mn++) {
        si = ss * 12 + mn;
        tok = strtok(NULL, " ,\t\n\015");
        if (tok == NULL) {
          sprintf(outstr, "Missing value for series on line %d\n");
          logAndExit(logfp, outstr);
        }
        if ( readParam(series->data[si], tok ) ) {
          if ( series->data[si].spType == pScalar )
            fprintf(logfp, "   %-34s  %4d/%02d constant: %12.6f\n", pName, series->start + ss, mn+1, series->data[si].sval );
          else if (series->data[si].spType == pTif)
            fprintf(logfp, "   %-34s  %4d/%02d grid: %s\n", pName, series->start + ss, mn+1, series->data[si].gridName );
        }
        else {
          sprintf( outstr, "Could not read series on line %d.\n", lineNo); 
          logAndExit( logfp, outstr ); 
        }
      }
    }
    series->got = true; 
  }
  return true;
}

//----------------------------------------------------------------------------------

void readParamFile(std::string *paramFile)
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
  const std::string line, pName;
  std::string pValue, outstr;
  std::string cp;
  int paramCount=0, lineLength=0, lineNo;
  int readingOutput=0;
  int len;

  // TODO: simplify this whole reading process using boost and std::string
  // start with: https://stackoverflow.com/questions/7868936/read-file-line-by-line-using-ifstream-in-c
  // Open file
  fprintf(logfp, "Reading input parameters from file '%s'...\n", paramFile);
  if ((paramFp = fopen(paramFile, "rb")) == NULL) {
    sprintf(outstr, "Could not open parameter file: %s\n", paramFile);
    fprintf(stderr, outstr);
    fprintf(logfp, outstr);
    exit(1);
  }
    
  // For each line. 
  lineNo=0;
  while (fgets(line, MAXLINE, paramFp) != NULL) {
    lineNo++;

    // Remove comments from end of line by inserting a null character. 
    cp = strstr(line, "//");
    if (cp != NULL)
      *cp = 0;

    // Remove whitespace from the end of the line. 
    len = strlen(line);
    if (len > 0) {
      for (cp = line + len-1; isspace(*cp) && cp >= line; cp--)
        ;
      *(cp+1) = 0;
    }

    // Consume leading whitespace. 
    cp = line + strspn(line, " \t");

    // Tokenize the line. Can't allow spaces to break as we have tokens 
    // containing spaces. 
    cp = strtok(cp, "\"\n\015");

    // Skip blank lines. 
    if (cp == NULL)
      continue;

    // Skip comments.
    if (cp[0] == '/' && cp[1] == '/')
      continue;

    // First token should be parameter name. 
    // (Seem to be overdoing the elimination of whitespace)
    strcpyTrim(pName, cp);

    // Second token should be value. Break tokens on CR for DOS files. 
    // Series parameters may have no value on this line. 
    cp = strtok(NULL, "\n\015");
    strcpyTrim(pValue, cp);

    // Get parameter index and read value. 
    if (readInputParam(pName, pValue))
      ;
    else if (readOutputParam(pName, pValue, lineNo))
      ;
    else if (readOtherParam(pName, pValue))
      ;
    else if (readInputSeriesParam(pName, pValue, paramFp, lineNo))
      ;
    else if (readInputManageParam(pName, paramFp, lineNo))
      ;
    else {
      sprintf(outstr, "Cannot read parameter in file %s, line %d: \"%s\"\n", 
        paramFile, lineNo, pName);
      logAndExit(logfp, outstr);
    }
  }
  fclose (paramFp);
}

//----------------------------------------------------------------------------------

bool haveAllParams()
{ 
  // Check that we have read a value for all parameters. 
  int i, pInd;
  int indSeed, indWRi, indWFi, indWSi;  //indexes for values to check if required variable available
  bool missing = false;

  // Parameters needed for 3PG. 
   std::string iParam3PG[] = {
    "pFS2", "pFS20", "StemConst", "StemPower", "pRx", "pRn", 
    "growthTmin", "growthTopt", "growthTmax",        // Temperature modifier (fT) 
    "kF",                                            // Frost modifier
    "gammaFx", "gammaF0", "tgammaF", "Rttover",      // Litterfall & root turnover
    "MaxCond", "CoeffCond", "BLcond",    // conductances
    "m0", "fN0", "fNn",                              // fertility effects
    "SWconst0", "SWpower0",                          // Soil water modifier (fSW)
    "wSx1000",                                       // stem numbers
    "MaxAge", "nAge", "rAge",                        // Age modifier (fAge)
    "SLA0", "SLA1", "tSLA", "k", "fullCanAge",       // Canopy structure and processes
    "alpha", "fracBB0", "fracBB1", "tBB", // Canopy structure and processes
    "y",                                  // various
    "Lat", "FR", "soilIndex", "MaxASW", "MinASW",    // 3PG site parameters. 
    "StartAge", "EndAge",              //Initial conditions
    //"WFi", "WRi", "WSi",             //Now checked along with SeedlingMass
    "StemNoi", "ASWi", "yearPlanted",  // Initial conditions. 
    "Qa", "Qb",
    "gDM_mol", "molPAR_MJ",
    "LAIgcx", "MaxIntcptn",
   "StartMonth", 
    "LAImaxIntcptn", 
    "thinPower", "mF", "mR", "mS",    //Thinning coefficients
    NULL
  };

  // Parameters needed for 3PGS.  
  std::string iParam3PGS[] = {
    "growthTmin", "growthTopt", "growthTmax",           // Temperature modifier (fT)
    "kF",                                               // Frost modifier
    "MaxCond", "CoeffCond", "BLcond",       // conductances
    "m0", "fN0",                                        // fertility effects
    "SWconst0", "SWpower0",                             // Soil water modifier (fSW)
    "SLA1", "alpha",                    // Canopy structure 
    "y",                                     // various
    "Lat", "FR", "soilIndex", "MaxASW", "MinASW",       // 3PG site parameters.
    "StartAge", "EndAge",                               // Initial conditions
    "NDVI_FPAR_intercept", "NDVI_FPAR_constant",        // FPAR from NDVI equation.
    "Qa", "Qb",
    "gDM_mol", "molPAR_MJ",
    "LAIgcx", "MaxIntcptn",
   "StartMonth", 
    "LAImaxIntcptn", 
    NULL
  };

  // Temperature series
  if ((!Tmax_vals.got) && (!userTavgSeries()))
    logAndExit(logfp, "No Tmax or Tavg data");
  if ((!Tmin_vals.got) && (!userTavgSeries()))
    logAndExit(logfp, "No Tmin or Tavg data");
  if (userTavgSeries())
    if (!userVpdSeries())
      logAndExit(logfp, "A VPD series is required if an average temperature series is used\n"); 

  if (!Rain_vals.got)
    logAndExit(logfp, "No Rain data");
  if (!SolarRad_vals.got && !NetRad_vals.got)
    logAndExit(logfp, "No Solar or Net Radiation data");
  if (!FrostDays_vals.got)
    logAndExit(logfp, "No Frost data");

  //Check for SeedlingMass and WSi, WFi, WRi

  indSeed = pNameToInd("SeedlingMass");
  indWFi = pNameToInd("WFi");
  indWRi = pNameToInd("WRi");
  indWSi = pNameToInd("WSi");

  if (!params[indSeed].got && (!params[indWFi].got || !params[indWRi].got || !params[indWSi].got))
  {
    sprintf(outstr, "Missing parameter for 3PGS mode, %s\n", params[indSeed].id);
    logAndPrint(logfp, outstr);
    missing = true;
  }

  //fill in Seedling Mass or others just in case - seems to cause errors if not there...
  if ( !haveSeedlingMass() )
  {
    sscanf("0.0", "%lf", params[indSeed].adr); 
    params[indSeed].data.spType = pScalar;
    params[indSeed].got = 0;
  } else {
    sscanf("0.0", "%lf", params[indWFi].adr); 
    params[indWFi].data.spType = pScalar;
    params[indWFi].got = 0;
    sscanf("0.0", "%lf", params[indWRi].adr); 
    params[indWRi].data.spType = pScalar;
    params[indWRi].got = 0;
    sscanf("0.0", "%lf", params[indWSi].adr); 
    params[indWSi].data.spType = pScalar;
    params[indWSi].got = 0;
  }


  


  // Check required parameters for 3PGS
  if (modelMode3PGS) {
    for (i = 0; iParam3PGS[i] != 0; i++) {
      pInd = pNameToInd(iParam3PGS[i]);
      if (pInd != 0) {
        if (!params[pInd].got) {
          sprintf(outstr, "Missing parameter for 3PGS mode, %s\n", params[pInd].id);
          logAndPrint(logfp, outstr);
          missing = true;
        }
      }
    }
    if (!NdviAvh_vals.got) 
      logAndExit(logfp, "No NDVIAHV data");
  }

  // Check required parameters for standard 3PG
  else {
    for (i = 0; iParam3PG[i] != 0; i++) {
      pInd = pNameToInd(iParam3PG[i]);
      if (pInd != 0) {
        if (!params[pInd].got) {
          sprintf(outstr, "Missing parameter for 3PG mode, %s\n", params[pInd].id);
          logAndPrint(logfp, outstr);
          missing = true;
        }
      }
    }
  }

  // Check various optional parameters
  //if ( !modelMode3PGS && NdviAvh_vals->got) {
  if ( !modelMode3PGS && NdviAvh_vals.data) {
    sprintf(outstr, "NDVI_AVH not used in 3PGS mode.\n"); 
    logAndPrint(logfp, outstr);
  }
  return !missing;
}

//----------------------------------------------------------------------------------

bool loadParamVals(std::tuple<int, int> k)
{
  // Load all model parameter values into their global variables.
  // Spatial parameters are taken from the current grid cell, and
  // non-spatial parameters are left unchanged.  k is the cell index
  // in the data array element of a FloatGrid or ByteGrid object.
  int pn;
  GDALRasterImage *fg;
  float result;
  std::string ErrorString[100];

  for (pn=1; params[pn].id != "-1"; pn++)  {
    if (params[pn].got == 0) 
    {
      //Ignore
    }
    else if (params[pn].data.spType == pTif) {
      fg = (GDALRasterImage *)params[pn].data.g;
      if (fg == NULL)
      {
        sprintf(ErrorString, "Error reading grid: %s - File not open.\n", params[pn].data.gridName);
        logAndExit(logfp, ErrorString); 
      }
      result = fg->GetVal(std::get<0>(k), std::get<1>(k));
      if (result == fg->noData) {
        return false;
      }
      *(params[pn].adr) = result;
    }
  }

    // Look for zero and do utterly bodgy things.
    //Note that any changes here should be echoed in findrunperiod.  
    //However, this is specifically for Aracruz...

    
    if ( yearPlanted < 1 ) 
          return false;  //Return nodata for years less than 1.
    if ( StartAge < 1 )
      StartAge = 1; 
    if ( EndAge < 1 )
      EndAge = 2; 
    if ( StartMonth < 1 ) 
      StartMonth = 1; 
   if ( StemNoi < 1 )
     return false;       //Return nodata for stocking values less than 1.

    return true;
}

//----------------------------------------------------------------------------------

bool openGrid( PPPG_VVAL &vval )
{
  if (vval.spType == pTif)
    fprintf(logfp, "   opening grid from %s...", vval.gridName); 
    vval.g = new GDALRasterImage(vval.gridName)
  else 
      return false;
    try {
  } catch (Exception &e) {
    sprintf(outstr, "\nException: %s\n", e.Message());
    logAndPrint(logfp, outstr);  // print as well since Exception messages don't  
    exit(1);                     // on MS VC++. 
  }
  return true; 
}

//----------------------------------------------------------------------------------
void ResetGrids(void)
{
  //Only interested in resetting the findrunperiod grids which are not series.
  for (int pn = 1; params[pn].data.spType != pNull; pn++) { // start at 1 to avoid error record. 
    if ( params[pn].data.spType == pTif ) {
      params[pn].data.g->ResetGrid();
    }
  }
}
//----------------------------------------------------------------------------------
void CloseGrids(void)
{
  for (int pn = 1; params[pn].data.spType != pNull; pn++) { // start at 1 to avoid error record. 
    if ( params[pn].data.spType == pTif ) {
      params[pn].data.g->CloseGrid();
      delete params[pn].data.g;
    }
  }
}
//----------------------------------------------------------------------------------
void PrintGrids(void)
{
  for (int pn = 1; params[pn].data.spType != pNull; pn++) { // start at 1 to avoid error record. 
    if ( params[pn].data.spType == pTif ) {
      params[pn].data.g->PrintGridState();
    }
  }
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

  fprintf(logfp, "Opening input grids...\n");
  fprintf(stdout, "Opening input grids...\r");

  for (pn = 1; params[pn].data.spType != pNull; pn++) { // start at 1 to avoid error record. 
    if (params[pn].got != 1)
    {
      //do nothing
    }
    else if ( openGrid( params[pn].data ) ) {
      spatial = true; 
      if ( first ) {
        refGrid = (GDALRasterImage *)params[pn].data.g;
        first = false; 
      }
      // TODO: Had to remove these conditions... try to find some GDAL compat way to do this
      //(fabs(refGrid->xmin - params[pn].data.g->xmin) > 0.0001)
      //    || (fabs(refGrid->ymin - params[pn].data.g->ymin) > 0.0001)
      //    || (fabs(refGrid->xmax - params[pn].data.g->xmax) > 0.0001)
      //    || (fabs(refGrid->ymax - params[pn].data.g->ymax) > 0.0001)
      else if (( refGrid->nRows != params[pn].data.g->nRows ) 
        || ( refGrid->nCols != params[pn].data.g->nCols ) ) {
        sprintf(outstr, "Grid dimensions must match, grid %s differs from first grid.\n", 
          params[pn].data.gridName ); 
        logAndExit(logfp, outstr); 
      }
    }
  }

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
        else if ( ( fabs( refGrid->xmin - ser->data[i].g->xmin ) > 0.0001 ) 
          || ( fabs( refGrid->ymin - ser->data[i].g->ymin ) > 0.0001 )
          || ( fabs( refGrid->xmax - ser->data[i].g->xmax ) > 0.0001 )
          || ( fabs( refGrid->ymax - ser->data[i].g->ymax ) > 0.0001 ) 
          || ( refGrid->nrows != ser->data[i].g->nrows ) 
          || ( refGrid->ncols != ser->data[i].g->ncols ) ) {
          sprintf(outstr, "Grid dimensions must match, grid %s differs from first grid.\n", 
            ser->data[i].gridName ); 
          logAndExit(logfp, outstr); 
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
        else if ( ( fabs( refGrid->xmin - tab[i].data.g->xmin ) > 0.0001 ) 
          || ( fabs( refGrid->ymin - tab[i].data.g->ymin ) > 0.0001 )
          || ( fabs( refGrid->xmax - tab[i].data.g->xmax ) > 0.0001 )
          || ( fabs( refGrid->ymax - tab[i].data.g->ymax ) > 0.0001 ) 
          || ( refGrid->nrows != tab[i].data.g->nrows ) 
          || ( refGrid->ncols != tab[i].data.g->ncols ) ) {
          sprintf(outstr, "Grid dimensions must match, grid %s differs from first grid.\n", 
            tab[i].data.gridName ); 
          logAndExit(logfp, outstr); 
        }
      }
    }
  }

  if (!spatial) {
    fprintf(logfp, "none\n");
    fprintf(stdout, "none\r");
    refGrid = NULL; 
  }


  return refGrid;
}

//----------------------------------------------------------------------------------

void openOutputGrids(GDALRasterImage *refGrid)
{
  // Loop through opVars array and open output grid objects for all 
  // those variables marked for output.  
  int opn; 

  // Open ordinary output grids. 
  fprintf(logfp, "Opening output grids...\n");
  fprintf(stdout, "Opening output grids...\r");
  for (opn = 0; opVars[opn].id != "-1"; opn++) {
    if (opVars[opn].write) {
      fprintf(logfp, "   float grid %s\n", opVars[opn].gridName);
      opVars[opn].g = new GDALRasterImage;
      if (opVars[opn].g->Exists(opVars[opn].gridName)) {
        sprintf(outstr, "Error, output grid named:\n"
          "   '%s'\nalready exists.\n", opVars[opn].gridName);
        fprintf(stderr, outstr);
        fprintf(logfp, outstr);
        exit(1);
      }
      // Copy the grid parameters from the reference grid, and allocate its 
      // data array.  
      *(opVars[opn].g) = *refGrid;
      opVars[opn].g->Allocate();
    }
  }
  return;
}

//----------------------------------------------------------------------------------

//bool copyHeader(GDALRasterImage *refGrid, char *fname )
//{
  // Really crappy way to generate a header file, rely on the assignment operator 
  // in the FloatGrid class.  
//  GDALRasterImage *fg;

//  fg = new GDALRasterImage;
//  *fg = *refGrid; 
//  try {
//    fg->Write(fname);
//  } catch (Exception &e) {
//    delete fg;
//    sprintf(outstr, "\nException: %s\n", e.Message()); 
//    logAndPrint(logfp, outstr);  // print as well as Exception msgs don't on Win32. 
//    return false; 
//  }
//  delete fg;
//  return true; 
//}

//----------------------------------------------------------------------------------

void openRegularOutputGrids( GDALRasterImage *refGrid, MYDate spMinMY, MYDate spMaxMY )
{
  // Open regular output grids. Unlike other grids in the program, these 
  // are managed as simple pointers to FILE.  This is because the current 
  // processing order (each cell, entire run period) means we need to 
  // append to these grids.  The pointers to FILE are in an array, being 
  // the .RO member of a opVars array entry.  In order to make the mapping from 
  // calYear, calMonth to regular output array simple, we allocate the array in 
  // whole years, and define element 0 as always being January.  THIS MEANS 
  // THAT MANY CELLS IN THE ARRAY CONTAIN NULL, in particular at the beginning 
  // and end, unless a 3PG run begins in January and finishes in December. 
 
  //int runYear, runMonth; 
  int calYear, calMonth; 
  int opn, roArrayLength; 
  FILE *fp; 
  std::string roStemName[GRID_NAME_LENGTH];
  std::string fname[GRID_NAME_LENGTH];
  std::string *cp;
  int len; 
  int maxCy, maxCm, minCy, minCm; 
  int mx;
  //int runMx; 
  int runMonths; 

  maxCy = spMaxMY.year; 
  maxCm = spMaxMY.mon; 
  minCy = spMinMY.year; 
  minCm = spMinMY.mon; 
  
  // How many months are encompassed in the run. Don't really need runMonths, 
  // but leaving it for comparison with mx below. 
  // First how many complete years. 
  runMonths = ( maxCy - minCy - 1 ) * 12; 
  // Then how many months in the (possibly) partial years at the beginning and end. 
  runMonths += 12 - minCm + 1; 
  runMonths += maxCm; 

  // In order to make the mapping from calYear, calMonth to regular output array 
  // simple, allocate the array in whole years, and define element 0 as always 
  // being January.  
  roArrayLength = ( maxCy - minCy + 1 ) * 12 + 12; // ANL last 12 is short write fix.  
  
  // for each output variable. 
  for (opn = 0; opVars[opn].id != "-1"; opn++) {
    // Is it marked for recurring output. 
    if ( !opVars[opn].recurYear )
      continue; 
    
    // Allocate and null out array of *FILE. 
    opVars[opn].RO = new FILE *[roArrayLength]; 
    for (int i = 0; i < roArrayLength; i++)
      opVars[opn].RO[i] = NULL;
    
    // Find stem name for regular output files, just lacking the 
    // year/month/extension. 
    cp = strrchr(opVars[opn].gridName, '.');
    len = cp - opVars[opn].gridName;
    strncpy( roStemName, opVars[opn].gridName, len ); 
    roStemName[len] = '\0'; 
    
    // Monthly output. mx indexes over entire ro array. 
    for ( mx = 0; mx < roArrayLength; mx++ ) {
      calYear = minCy + ( mx / 12); 
      calMonth = ( mx % 12 ) + 1; 
      
      // Skip months before the model starts. 
      if ( calYear == minCy && calMonth < minCm )
        continue; 
      
      // Skip years after the model ends. 
      if ( calYear > maxCy )
        continue; 

      // Skip months after the model ends. 
      if ( calYear == maxCy && calMonth > maxCm )
        continue; 
      
      // Reached first output year? 
      if ( calYear < opVars[opn].recurStart )
        continue; 

      // Is this an output year? 
      if ((( calYear - opVars[opn].recurStart) % opVars[opn].recurYear ) != 0 )
        continue; 
      
      // Is this an annual output?  If so only want 1 month as given in recurMonth. 
      if ( !opVars[opn].recurMonthly ) {
        if ( calMonth == opVars[opn].recurMonth ) {
          // Construct annual regular output grid name and open file. 
          sprintf( fname, "%s%4d%02d.flt", roStemName, calYear, calMonth ); 
          if ( ( fp = fopen( fname, "wb" ) ) == NULL ) {
            sprintf(outstr, "Couldn't open regular output file %s\n", fname);
            logAndExit(logfp, outstr);
          }
          if ( ! copyHeader( refGrid, fname ) ) {
            sprintf(outstr, "Couldn't open regular output header file %s\n", fname);
            logAndExit(logfp, outstr);
          }
          opVars[opn].RO[mx] = fp; 
        }
      }
      // Is a monthly output ro grid and this is the output year. 
      else {
        // Construct monthly regular output grid names and open file. 
        sprintf( fname, "%s%4d%02d.flt", roStemName, calYear, calMonth ); 
        if ( ( fp = fopen( fname, "wb" ) ) == NULL ) {
          sprintf(outstr, "Couldn't open regular output file %s\n", fname);
          logAndExit(logfp, outstr);
        }
        if ( ! copyHeader( refGrid, fname ) ) {
          sprintf(outstr, "Couldn't open regular output header file %s\n", fname);
          logAndExit(logfp, outstr);
        }
        opVars[opn].RO[mx] = fp; 
      }
    } // closes for ( mx = 0; mx < roArrayLength; mx++ )
  } // closes for each output variable. 
}
  
//----------------------------------------------------------------------------------

int writeOutputGrids(void)
{
  int opn, gridsWritten=0;

  fprintf(stdout, "Writing output grids...\r"); 
  fprintf(logfp, "Writing output grids...\n");
  for (opn=0; opVars[opn].id != "-1"; opn++) {
    if (opVars[opn].write) {
      try {
        opVars[opn].g->Write(opVars[opn].gridName);
      } catch (Exception &e) {
        sprintf(outstr, "Exception: %s\n", e.Message());
        return 0;
      }
      gridsWritten++;
    }
  }
  return gridsWritten;
}

//----------------------------------------------------------------------------------

void writeSampleFiles(std::tuple<int,int> cellIndex, int month, long calYear) 
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
    month, samplePoints[sInd].id);
  // For each variable
  for (opn = 1; opVars[opn].adr != NULL; opn++)
    fprintf(samplePoints[sInd].fp, "%f, ", *(opVars[opn].adr));
  fprintf(samplePoints[sInd].fp, "\n");
}

//----------------------------------------------------------------------------------

void writeMonthlyOutputGrids( int calYear, int calMonth, bool hitNODATA, 
                              MYDate minMY, MYDate maxMY )
{
  int opn, mx, maxInd; 
  float fval; 
  long lval; 
  FILE *fp; 
  GDALRasterImage *fg;

  int crap; 
  if ( calYear == 2008 )
    crap = 1; 

  // Number of elements in the regular output array, for sanity checking. 
  // Usually not all of these would be used.  
  maxInd = ( maxMY.year - minMY.year + 1 ) * 12 + 12;  // ANL see comment in openRegularOutputGrids

  // for each output variable. 
  for (opn = 0; opVars[opn].id != "-1"; opn++) {
    // Is it marked for recurring output. 
    if ( !opVars[opn].recurYear )
      continue; 
    
    // Skip variables not marked for recurring output. 
    if ( !opVars[opn].recurStart ) 
      continue;
    
//    // Skip variables not marked for monthly recurrence. 
//    if ( !opVars[opn].recurMonthly )
//      continue; 

    mx = ( calYear - minMY.year ) * 12 + ( calMonth - 1 );
    
    if (mx < 0)
      continue;
    
    // Sanity check. 
    if ( mx > maxInd ) {
      sprintf(outstr, "Program error, mx=%d too high in writeMonthlyOutputGrids at month/year %2d/%4d.\n", 
        mx, calYear, calMonth);
      logAndExit(logfp, outstr);
    }

    // Get file pointer for this regular output grid. 
    fp = opVars[opn].RO[mx];

    // Write it. 
    if ( !fp )
      continue;
    if ( opVars[opn].spType == pTif ) {
      fval = (float)*(opVars[opn].adr);
      fg = (GDALRasterImage *)opVars[opn].g;
      if ( hitNODATA )
        fval = fg->noData;
      fwrite( &fval, 4, 1, fp );
      /*sprintf(outstr, "Year: %d, Month: %d, Value: %f\n", calYear, calMonth, fval); 
      logOnly(logfp, outstr);*/
    }
  }
}

//----------------------------------------------------------------------------------

void writeYearlyOutputGrids( int calYear, int calMonth, bool hitNODATA, 
                            MYDate minMY, MYDate maxMY )
{
  int opn, mx, maxInd; 
  float fval; 
  long lval; 
  FILE *fp; 
  GDALRasterImage *fg;

  // Number of elements in the regular output array, for sanity checking. 
  maxInd = ( maxMY.year - minMY.year + 1 ) * 12; 

  for (opn = 0; opVars[opn].id != "-1"; opn++) {
    // Skip variables not marked for recurring output. 
    if ( !opVars[opn].recurStart ) 
      continue;

    // Skip variables marked for monthly recurrence. 
    if ( opVars[opn].recurMonthly )
      continue; 

    // Get regular output array index. 
    mx = ( calYear - minMY.year ) * 12 + ( calMonth - 1 ); 

    // Sanity check. 
    if ( mx > maxInd ) {
      sprintf(outstr, "Program error, mx=%d too high in writeYearlyOutputGrids at month/year %2d/%4d.\n", 
        mx, calYear, calMonth);
      logAndExit(logfp, outstr);
    }

    // Get file pointer for this regular output grid. 
    fp = opVars[opn].RO[mx];

    // Write it. 
    if ( !fp )
      continue;
    if ( opVars[opn].spType == pTif ) {
      fval = (float)*(opVars[opn].adr);

      
      fg = (GDALRasterImage *)opVars[opn].g;
      
      if ( hitNODATA )
        fval = fg->noData;
      fwrite( &fval, 4, 1, fp );
    }
  }
}
    
//----------------------------------------------------------------------------------

void saveVariableVals(long k, bool hitNODATA)
{
  // Copy variable values from the model to the current grid cell. 
  int opn;
  GDALRasterImage *fg;

  for (opn = 1; opVars[opn].id != "-1"; opn++) 
    if (opVars[opn].spType == pTif) {
      fg = (GDALRasterImage *)opVars[opn].g;
      fg->z[k] = (float)*(opVars[opn].adr);
      if (hitNODATA)
        fg->z[k] = fg->noData;
    }
}

//----------------------------------------------------------------------------------

FILE *openLogFile(std::string siteParamFile)
{
  // Open the site parameter file, find the output directory 
  // specification, open the logfile in that directory, close the 
  // site parameter file.  

  FILE *logFileFp, *paramFp; 
  std::string *cp, *line, *logFileName;

  if ((paramFp = fopen(siteParamFile, "rb")) == NULL) {
    fprintf(stderr, "Could not open site parameter file %s\n", siteParamFile); 
    exit(1);
  }

  line = new std::string[1000];
  logFileName = new std::string[1000];

  while (fgets(line, MAXLINE, paramFp) != NULL) {
    // Remove comments from end of line by inserting a null character. 
    cp = strstr(line, "//");
    if (cp != NULL)
      *cp = 0;

    // Consume leading whitespace. 
    cp = line + strspn(line, " \t");

    // Tokenize the line. First token ends with a closing double quote.  
    cp = strtok(cp, "\"\n\015");
    if (cp == NULL)
      continue;

    // Find the output directory specifier.  Look for ^M (ascii 13), carriage return, in 
    // DOS text files. 
    if (namesMatch("output directory", cp)) {
      cp = strtok(NULL, " \t\n\015");
      if (cp == NULL) {
        line[0] = '.';
        line[1] = '\0';
        cp = line;
      }
      int len = strlen(cp);
      if (*(cp + len - 1) != '/') {
        *(cp + len) = '/';
        *(cp + len + 1) = 0;
      }
      strcpy(logFileName, cp); 
      strcat(logFileName, "logfile.txt");

      // Open the file
      if ((logFileFp = fopen(logFileName, "w")) == NULL) {
        fprintf(stderr, "Could not open log file %s\n", logFileName);
        exit(1);
      }
      else {
        fclose(paramFp);
        delete line;
        delete logFileName;
        return logFileFp; 
      }
    }
  }
  fprintf(stderr, "Could not find \"output directory\" parameter in %s\n", siteParamFile);
  exit(1);
}

//----------------------------------------------------------------------------------

void writeStandSummaryData(FILE *outFP, int year)
{
    fprintf(outFP, "%4.0f, %4.0f, %4.0f, %5.2f, %5.2f, %5.2f, %5.1f, %4.2f, "
    "%4.2f, %3.1f, %4.2f, %4.0f, %4.0f\n", 
    StandAge + yearPlanted, StandAge, StemNo, WF, WR, WS, StandVol, LAI, MAI, avDBH, 
    cLitter, cumTransp, ASW);
}

//----------------------------------------------------------------------------------

void writeStandSummary(int year)
{
  // For point mode, write summary values to the output file pointModeFp. 

  // 3PGS
  if (modelMode3PGS) {
      // headings
      if (year == 1)
        fprintf(pointModeFp, "Year, LAI, , NPP, delWAG, cumWabv\n");

      // write initial conditions
      fprintf(pointModeFp, "%4.0f, %3.1f, %5.2f, %5.2f, %5.2f\n", 
        yearPlanted + StandAge - 1, LAI, NPP, delWAG, cumWabv);
  }

  // 3PG
  else {
    // Headings and initial conditions. 
    if (year == 0) {
        fprintf(pointModeFp, "Year, Stand age, Stem number, Wf, Wr, Ws, Stand volume, Stand LAI, "
          "Stand MAI, Average DBH, Litterfall, Total transpiration, Soil water\n");

        fprintf(pointModeFp, "%4.0f, %4.0f, %4.0f, %5.2f, %5.2f, %5.2f, %5.1f, %4.2f, %4.2f, %3.1f\n", 
          yearPlanted + StandAge, StandAge, StemNoi, WFi, WRi, WSi, 1.7 * WSi, 
          LAIi, MAIi, avDBHi);
    }

    // write annual data
    writeStandSummaryData (pointModeFp, year);

    // write peak LAI and MAI
    if (StandAge == EndAge)
      fprintf(pointModeFp, " , , , , , , ,%5.2f, %5.2f\n", LAIx, MAIx); 
  }
}

//----------------------------------------------------------------------------------

void findRunPeriod( GDALRasterImage *refGrid, MYDate &minMY, MYDate &maxMY )
{
  // Find the calendar months and years for which the model will actually run, over 
  // the whole spatial area.  Cells with NODATA on any input cause that cell to 
  // have no effect.  Cells with zero on yearPlanted get, in this version, 
  // mapped to 1993.  THIS IS A BAD IDEA AND WILL COME BACK TO BITE US.  Likewise 
  // zero in StartAge gets remapped to 1, and zero on monthPlanted gets remapped to 
  // to 1 also.  This messing with values is duplicated in loadParamVals, which 
  // is one more reason this is a bad idea.  
  int cy, cm, minCy, maxCy, minCm, maxCm; 
  float cy2, minCy2, maxCy2; 
  bool spatial, hitNoData;
  bool test;
  int nrows, ncols; 
  const int MINSEED = 999999; 
  const int MAXSEED = -1; 
  int adjAge; //Plant age and end age, taking into account endage being an actual year

  minCy = MINSEED;
  maxCy = MAXSEED;
  minCy2 = MINSEED; 
  maxCy2 = MAXSEED; 
  spatial = ( refGrid != NULL ); 
  
  // Point mode case.  If we are running in point mode yearPlanted, StartAge and 
  // EndAge are already defined. 
  if ( !haveSpatialRunYears() ) {
    minMY.year = (int)yearPlanted + (int)StartAge; 
    minMY.mon = (int)StartMonth;
    if (EndAge > yearPlanted)
      maxMY.year = (int)EndAge; //Take into account actual years rather than number of years
    else
      maxMY.year = (int)yearPlanted + (int)EndAge; 
    maxMY.mon = (int)StartMonth;
  }
  // Spatial mode case.  For each cell get yearPlanted, monthPlanted, StartAge, 
  // EndAge, endMonth, and find the earliest start and latest finish.  
  else {
    nrows = refGrid->nrows; 
    ncols = refGrid->ncols; 
    
    // For every cell
    //for (int j = 0; j < nrows-2; j++) {
    // long cellIndex = j * ncols + 1;
    
    //  for (int i = 1; i < ncols-1; i++, cellIndex++) {
    for (int cellIndex = 0; cellIndex < ncols*nrows; cellIndex++)
    {
    // Load parameter values - only really want yearPlanted, EndAge, and 
        // StartAge. 
        
        if (cellIndex == 159)
          test = true;

        hitNoData = !loadParamVals(cellIndex);

        
        // Look for NODATA, to skip it. 
        if ( hitNoData )
        {
     //     logAndPrint(logfp, "Nodata\n");
          continue; 
        }
        // Look for zero and do utterly bodgy things. 
        if ( yearPlanted < 1)
          continue;   //Treat years less than 1900 as nodata 
        if ( StartAge < 1 )
          StartAge = 1; 
        if ( EndAge < 1 )
          EndAge = 2; // Mimum end growth for a stand - just a reasonable value. AS.
        if ( StartMonth < 1 ) {
          StartMonth = 1; 
        }

        // Check reasonableness of StartMonth, must be in the range 1 - 12 inclusive. 
        if ( StartMonth < 0.999 || StartMonth > 12.0001 )
          logAndExit(logfp, "Found bad StartMonth value at XXXX\n" ); 
        
        // Minimum year
        cy2 = yearPlanted + (float)StartAge + ( StartMonth / 100 ); 
        if ( cy2 < minCy2 )
          minCy2 = cy2; 

        cy = (int)(yearPlanted + StartAge);
        cm = (int)StartMonth; 
        
        if ( cy < minCy ) {
          minCy = cy; 
          minCm = 13;  // ensure that the current cm is found as a minimum. 
        }

        // Minimum month, only update when cy is at a minimum.  
        if ( cy == minCy ) {
          if ( cm < minCm )
            minCm = cm; 
        }

        // Maximum year, month.  Currently we only allow EndAge to be a whole number 
        // of years, so the last month is always one month earlier than StartMonth. 
        // In the future might allow fractional EndAge and work out end month from that. 
        
        if (EndAge > yearPlanted)
          adjAge = (int)EndAge;
        else
          adjAge = yearPlanted + (int)EndAge; 
        cy2 = adjAge + ( ( StartMonth - 1 ) / 100 ); 
        
        if ( cy2 > maxCy2 )
          maxCy2 = cy2; 
        
        cy = (int)(adjAge); 
        cm = (int)StartMonth - 1; 
        
        if ( cm == 0 ) 
          cm = 12; 
        if ( cy > maxCy ) {
          maxCy = cy; 
          maxCm = 0;  // ensure that the current cm is found as a maximum.
        }

        // Maximum month. 
        if ( cy == maxCy )
          if ( cm > maxCm )
            maxCm = cm; 
      }
    //}

    // Check if we found useful values
    if ( minCy == MINSEED )
      logAndExit( logfp, "Failed to find valid starting year.\n" ); 
    if ( maxCy == MAXSEED )
      logAndExit( logfp, "Failed to find valid ending year.\n" ); 
    if ( minCm == 13 )
      logAndExit( logfp, "Failed to find valid starting month.\n" ); 
    if ( maxCm == 13 )
      logAndExit( logfp, "Failed to find valid ending month.\n" ); 

    minMY.year = minCy; 
    minMY.mon = minCm; 

    maxMY.year = adjAge;        //minCy + EndAge - StartAge; 
    maxMY.mon = maxCm; 
  }

  
  return;
}

//----------------------------------------------------------------------------------
//To ensure that if we don't have a variable, it doesn't get an accidental GOT value
//Initialise the params list so that everything in missing.
void InitInputParams(void)
{
  int pn;

  for (pn=0; params[pn].id != "-1"; pn++) {
    params[pn].got = 0;
  }

  //The following parameters need to be set to scalar or grids do not open.

  pn = pNameToInd("FRstart");
  params[pn].data.spType = pScalar;
  pn = pNameToInd("FRend");
  params[pn].data.spType = pScalar;
  pn = pNameToInd("FRdec");
  params[pn].data.spType = pScalar;

}

//----------------------------------------------------------------------------------

bool havePointOpFile()
{
  if (pointModeFp == NULL) {
    sprintf(outstr, "Using point mode but parameter \"point mode output file\" not set\n");
    logAndPrint(logfp, outstr);
    return false;
  }
  else
    return true;
}
//----------------------------------------------------------------------------------
bool haveSeedlingMass()
{
  int pInd;
  pInd = pNameToInd("SeedlingMass");
  if (params[pInd].got)
    return true;
  else
    return false;
}

//----------------------------------------------------------------------------------
bool haveSpatialRunYears()
{
  int pInd;

  pInd = pNameToInd("StartAge");
  if (!(params[pInd].data.spType == pScalar)) return true;

  pInd = pNameToInd("EndAge");
  if (!(params[pInd].data.spType == pScalar)) return true;

  pInd = pNameToInd("StartMonth");
  if (!(params[pInd].data.spType == pScalar)) return true;

  pInd = pNameToInd("yearPlanted");
  if (!(params[pInd].data.spType == pScalar)) return true;

  return false;

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

//----------------------------------------------------------------------------------

bool haveMinASWTG(void)
{
  int pInd;
  bool bResult;
  pInd = pNameToInd("MinASWTG");
  
  if (params[pInd].got)
    bResult = true;
  else
    bResult = false;
  
  return bResult;
}
//----------------------------------------------------------------------------------

bool haveRhoMin(void)
{
  int pInd;
  bool bResult;
  pInd = pNameToInd("rhoMin");
  
  if (params[pInd].got)
    bResult = true;
  else
    bResult = false;
  
  return bResult;  
}
//----------------------------------------------------------------------------------

bool haveRhoMax(void)
{
  int pInd;
  bool bResult;
  pInd = pNameToInd("rhoMax");
  
  if (params[pInd].got)
    bResult = true;
  else
    bResult = false;
  
  return bResult;  
}
//----------------------------------------------------------------------------------

bool haveTRho(void)
{
  int pInd;
  bool bResult;
  pInd = pNameToInd("tRho");
  
  if (params[pInd].got)
    bResult = true;
  else
    bResult = false;
  
  return bResult;  
}
//----------------------------------------------------------------------------------
//To use an age dependant modifier for fertility, all three paramters need to be 
//defined:
// FRstart - When the age dependant modifier starts effecting the FR rating
// FRend - When the age dependent modifier stops effecting the FR rating
// FRdec - The decrement per month

bool haveAgeDepFert(void)
{
  int pInd;
  
  pInd = pNameToInd("FRstart");
  if (params[pInd].got == 0)
    return false;
  
  pInd = pNameToInd("FRend");
  if (params[pInd].got == 0)
    return false;
  
  pInd = pNameToInd("FRdec");
  if (params[pInd].got == 0)
    return false;

  return true;
}