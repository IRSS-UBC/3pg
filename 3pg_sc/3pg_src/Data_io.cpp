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
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include "GDALRasterImage.hpp"
#include "util.hpp"
#include "Data_io.hpp"
#include "Params.hpp"

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
// char * outstr;

#define GRID_NAME_LENGTH 300
#define PPPG_MAX_SERIES_YEARS 150
#define PPPG_MAX_SERIES_LENGTH PPPG_MAX_SERIES_YEARS*12

// Possible types of parameter - null or not yet set, scalar (constant), 
// ByteGrid for a BIL, and FloatGrid for a floating point file.  The last 
// two are handled by the JG Grid class.  




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
//extern double StartAge, EndAge;                 // age of trees at start/end of run
//extern double StartMonth;                       // month of year to start run
//extern double yearPlanted;                      // year trees planted
//// int DaysInMonth[13];                         // array for days in months
//extern double DaysInMonth[13];                  // array for days in months
extern bool showDetailedResults;              // TRUE ==> show monthly results
extern bool showStandSummary;                 // TRUE ==> show stand summary
extern bool modelMode3PGS;
//
//// Site characteristics, site specific parameters
//extern std::string siteName;                      // name of site
//extern double Lat;                              // site latitude
//extern double MaxASW, MinASWp;                  // maximum & minimum available soil water
//extern double FRp, FR;                              // current site fertility rating
//extern double FRstart, FRend, FRdec;            // Start, end and decrement % for fertility decrease with time
//extern double soilIndex;                        // soil class index
//extern double SWconst, SWpower;                 // soil parameters for soil class
//
//// Time variant management factors
//extern int nFertility;                          // size of site fertility array
////extern MANAGE_TABLE Fertility[1000];            // time-variant site fertility
//extern int nMinAvailSW;                         // size of MinAvailSW array
////extern MANAGE_TABLE MinAvailSW[1000];           // time-variant MinAvailSW (mm)
//extern int nIrrigation;                         // size of irrigation array
////extern MANAGE_TABLE Irrigation[1000];           // time-variant irrigation (ML/y)
//extern double Irrig;                            // current annual irrigation (ML/y)
//
//// Mean monthly weather data
////int mYears;                                   // years of met data available
//// ANL changed this from int to double
//extern double mYears;                           // years of met data available
//extern double mDayLength[13];                   // day length
////int mFrostDays[13];                           // frost days/month
//// ANL changed this from int to double
//extern double mFrostDays[13];                   // frost days/month
//extern double mSolarRad[13];                    // solar radiation (MJ/m2/day)
//extern double mTx[13];                          // maximum temperature
//extern double mTn[13];                          // minimum temperature
//extern double mTav[13];                         // mean daily temperature
//extern double mVPD[13];                         // mean daily VPD
//extern double mRain[13];                        // total monthly rain + irrigation
//extern double mNDVI[13];                        // ANL monthly NDVI for 3PGS mode
//extern double mNetRad[13];                      // ANL can use net instead of short wave
//
//// Stand data
//// extern char SpeciesName[100];                // name of species
//// int StandAge;                                // stand age
//// ANL changed StandAge from int to double
//extern double StandAge;                         // stand age
//extern double ASW, ASWi;                        // available soil water
//extern double MinASWTG;                         // soil water modifier corrector
//extern double StemNoi, StemNo;                  // stem numbers
//extern double SeedlingMass;
//extern double WFi, WF;                          // foliage biomass
//extern double WRi, WR;                          // root biomass
//extern double WSi, WS;                          // stem biomass
//extern double LAIi, LAI;                        // canopy leaf area index
//extern double MAIi, MAI;                        // mean annual volume increment
//extern double avDBHi, avDBH;                    // average stem DBH
//extern double TotalW;                           // total biomass
//extern double BasArea;                          // basal area
//extern double StandVol;                         // stem volume
//extern double LAIx, ageLAIx;                    // peak LAI and age at peak LAI
//extern double MAIx, ageMAIx;                    // peak MAI and age at peak MAI
//extern double cumTransp;                        // annual stand transporation
//extern double cumIrrig;                         // annual irrig. to maintain MinASW
//
//// Stand factors that are specifically age dependent
//extern double SLA;
//extern double Littfall;
//extern double fracBB;
//extern double CanCover;
//
//// Parameter values
//// int MaxAge;
//// ANL changed MaxAge from int to double
//extern double MaxAge;
//extern double gammaFx, gammaF0, tgammaF;
//extern double Rttover;
//extern double SLA0, SLA1, tSLA;
//extern double fullCanAge;
//extern double k;
//extern double pFS2, pFS20;
//extern double StemConst, StemPower;
//extern double SWconst0, SWpower0;
//extern double Interception;
//extern double BLcond;
//extern double MaxCond, CoeffCond;
//extern double y;
//extern double growthTmax, growthTmin, growthTopt;
//extern double thinPower;           //Added 29-07-02
//extern double mF, mR, mS;          //Added 29-07-02
//extern double wSx1000;
//extern double m0, fN0, fNn;
//extern double alpha, alphaC;  //alphaC added 11/07/02
//extern double pRx, pRn;
//extern double nAge, rAge;
//extern double kF;
//extern double fracBB0, fracBB1, tBB;
//extern double fracBB; //fracBB added 11/07/02
//extern double Density;
//extern double pfsConst, pfsPower;                     // derived from pFS2, pFS20
//extern double rhoMin, rhoMax, tRho;             // Standage varying density 3-06-02 
//extern double PhysMod;
//extern double WUE;                              //Added 16/07/02
//extern double CVI;                              //Added 16/07/02
//extern double TotalLitter;                      //Added 16/07/02
//
////Conversion factors
//extern double Qa, Qb; 
//extern double gDM_mol; 
//extern double molPAR_MJ; 
//
////Additional factors (conductance)
//extern double LAIgcx;
//extern double MaxIntcptn;
//extern double LAImaxIntcptn;
//
//// Intermediate monthly results
//extern double m, alphaC;
//extern double RAD, PAR;
//extern double lightIntcptn;
//extern double fAge, fT, fFrost;
//extern double fVPD, fSW, fNutr;
//extern double CanCond;
//extern double Transp, EvapTransp;
//extern double AvStemMass;
//extern double APAR, APARu;
//extern double GPPmolc, GPPdm, NPP;
//extern double pR, pS, pF, pFS;
//extern double delWF, delWR, delWS;
//extern double delFloss, delRloss;
//extern double monthlyIrrig;
//
//// Annual results
//extern double cLAI, cGPP, cNPP, cCVI, cRainInt, cEvapTransp, cTransp, cWUE;
//extern double cumGPP, cumWabv;
//extern double abvgrndEpsilon, totalEpsilon;
//extern double StemGrthRate;
//extern double cLitter;
//extern double CumdelWF, CumdelWR, CumdelWS;
//extern double CumAPARU, cumARAD;
//extern double CumStemLoss;
//extern double CutStemMass1, CutStemMass2, CutStemMass3;
//
////----------------------------------------------------------------------------------
//
//// 3PGS variables
//extern double NDVI_FPAR_intercept, NDVI_FPAR_constant; 
//extern double delWAG;

// ANL - other globals.
bool yearlyOutput, monthlyOutput; 
bool samplePointsYearly = false, samplePointsMonthly = false;
FILE *pointModeFp;   // Output file for point mode only.  
std::string outPath = "./";

//----------------------------------------------------------------------------------

// Initialisation of parameter array. This sets up the mapping between the variable 
// and its name, which is used in parsing the parameter files.

//const std::string paramError = "paramError";


//PPPG_PARAM params[] =
//{
//  {"paramError", NULL},
//  {"pFS2",         pFS2},
//  {"pFS20",        pFS20},
//  {"StemConst",    StemConst},
//  {"StemPower",    StemPower},
//  {"pRx",          pRx},
//  {"pRn",          pRn},
//
//  // Temperature modifier (fT) | cardinal temperatures
//  // ANL - these have been renamed from just Tmax etc, to avoid confusion with the 
//  // climate variables. 
//  {"growthTmin",   growthTmin},
//  {"growthTopt",   growthTopt},
//  {"growthTmax",   growthTmax},
//
//  // Frost modifier
//  {"kF",           kF},
//
//  // Litterfall & root turnover
//  {"gammaFx",      gammaFx},
//  {"gammaF0",      gammaF0},
//  {"tgammaF",      tgammaF},
//  {"Rttover",      Rttover},
//
//  // conductances
//  {"MaxCond",      MaxCond},
//  {"CoeffCond",    CoeffCond},
//  {"BLcond",       BLcond},
//
//  // fertility effects
//  {"m0",           m0},
//  {"fN0",          fN0},
//  {"fNn",          fNn},
//
//  //Thinning effects
//  {"thinPower",    thinPower},
//  {"mF",           mF},
//  {"mR",           mR},
//  {"mS",           mS},
//
//  // Soil water modifier (fSW) | soil characteristics
//  {"SWconst0",     SWconst0},
//  {"SWpower0",     SWpower0},
//
//  // stem numbers
//  {"wSx1000",      wSx1000},
//
//  // Age modifier (fAge)
//  {"MaxAge",       MaxAge},
//  {"nAge",         nAge},
//  {"rAge",         rAge},
//
//  // Canopy structure and processes | specific leaf area
//  {"SLA0",         SLA0},
//  {"SLA1",         SLA1},
//  {"tSLA",         tSLA},
//  {"k",            k},
//  {"fullCanAge",   fullCanAge},
//  {"alpha",        alpha},
//  {"fracBB0",      fracBB0},
//  {"fracBB1",      fracBB1},
//  {"tBB",          tBB},
//
//  // various
//  {"y",            y},
//  {"rhoMin",       rhoMin},
//  {"rhoMax",       rhoMax},
//  {"tRho",         tRho},             // Standage varying density 3-06-02 
//
//  //Conversions
//  {"Qa",           Qa},
//  {"Qb",           Qb},
//  {"gDM_mol",      gDM_mol},
//  {"molPAR_MJ",    molPAR_MJ},
//
//  //Additional conversion factors 
//  {"LAIgcx",           LAIgcx},
//  {"MaxIntcptn",       MaxIntcptn},
//  {"LAImaxIntcptn",    LAImaxIntcptn},
//
//  // 3PG site parameters. 
//  {"Lat",          Lat},
//  {"FRp",          FRp},
//  {"FRstart",      FRstart},  //These three variables relate to fertility decrease with age
//  {"FRend",        FRend},
//  {"FRdec",        FRdec},
//  {"soilIndex",    soilIndex},
//  {"MaxASW",       MaxASW},
//  {"MinASWp",      MinASWp},
//
//  // Initial conditions. 
//  {"StartAge",     StartAge},
//  {"EndAge",       EndAge},
//  {"StartMonth",   StartMonth},
//  {"yearPlanted",  yearPlanted},  /* CHECK! do we still use this?*/
//  {"SeedlingMass", SeedlingMass},
//  {"WFi",          WFi},
//  {"WRi",          WRi},
//  {"WSi",          WSi},
//  {"StemNoi",      StemNoi},
//  {"ASWi",         ASWi},
//  {"MinASWTG",     MinASWTG},
//  //  {"yearPlanted",  &yearPlanted},  /* This has been moved as strange errors were occuring with grids here*/
//
//    // ANL - extras for 3PGS mode
//    {"NDVI_FPAR_intercept", NDVI_FPAR_intercept},
//    {"NDVI_FPAR_constant",  NDVI_FPAR_constant},
//
//    {"", NULL}  // NULL entries used to mark array ends. 
//};
//
////----------------------------------------------------------------------------------
//
//// Initialisation of output variable array. This lists all possible output variables 
//// and sets up the mapping of the output variable to its name, which is used in 
//// parsing the parameter file.  
//PPPG_OP_VAR opVars[] = {
//  {"opVarError", NULL},
//  {"StemNo",     StemNo},
//  {"WF",         WF},
//  {"WR",         WR},
//  {"WS",         WS},
//  {"TotalW",     TotalW},
//  {"LAI",        LAI},
//  {"cLAI",       cLAI},
//  {"MAI",        MAI},
//  {"avDBH",      avDBH},
//  {"BasArea",    BasArea},
//  {"StandVol",   StandVol},
//  {"GPP",        GPPdm},
//  {"cGPP",       cGPP},
//  {"NPP",        NPP},
//  {"cNPP",       cNPP},
//  {"delWAG",     delWAG},
//  {"cumWabv",    cumWabv},
//  {"Transp",     Transp},
//  {"cTransp",    cTransp},
//  {"ASW",        ASW},
//  {"fSW",        fSW},
//  {"fVPD",       fVPD},
//  {"fT",         fT},
//  {"fNutr",      fNutr},
//  {"fFrost",     fAge},
//  {"APAR",       APAR},
//  {"APARu",      APARu},
//  {"EvapTransp", EvapTransp},
//  {"cEvapTransp",cEvapTransp}, //Added 08/11/02
//  {"LAIx",       LAIx},
//  {"ageLAIx",    ageLAIx},
//  {"MAIx",       MAIx},    //Added 29/07/2002
//  {"ageMAIx",    ageMAIx}, //Added 29/07/2002
//  {"FR",         FR},     //Added 11/07/2002
//  {"PhysMod",    PhysMod}, //Added 11/07/2002
//  {"alphaC",     alphaC},  //Added 11/07/2002
//  {"fAge",       fAge},    //Added 11/07/2002
//  {"fracBB",     fracBB},
//  {"WUE",        WUE},     //Added 16/07/02
//  {"cWUE",       cWUE},    //Added 08/11/02
//  {"CVI",        CVI},     //Added 16/07/02
//  {"cCVI",      cCVI},    //Added 08/11/02
//  {"TotalLitter", TotalLitter}, //Added 16/07/02
//  {"cLitter",    cLitter},
//  {"",         NULL}
//};
 
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

//// 3PG management table parameters. At most one value per year. 
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

int pNameToInd(const std::string& id, const std::vector<PPPG_PARAM>& params)
{
  // TODO: this can be deprecated if the param arrays is just replaced with a map.

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
  for (pn=0; params[pn].id != ""; pn++) {
    w2 = params[pn].id.length();
    if (id.compare(params[pn].id) == 0)
      if (w1 == w2)
        return pn;
  }
 /* std::cout << "Warning, lookup of non-existent parameter name: " << id << std::endl;*/
  //fprintf(stderr, "Warning, lookup of non-existent parameter name: %s\n", 
  //        id);
  return 0;
}

//----------------------------------------------------------------------------------

int opNameToInd(const std::string& id, const std::vector<PPPG_OP_VAR>& opVars)
{
  // Exactly the same as pNameToInd, but looking at the opVars array. 
  std::size_t pn;
  std::size_t w1, w2;

  w1 = id.length();
  for (pn = 0; opVars[pn].id != ""; pn++) {
    w2 = opVars[pn].id.length();
    if (id.compare(opVars[pn].id) == 0)
      if (w1 == w2)
        return pn;
  }
  std::cout << "Warning, lookup of non-existent output variable name: " << id << std::endl;
  //fprintf(stderr, "Warning, lookup of non-existent output variable name: %s\n", 
  //        id);
  return 0;
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
    std::cout << "Program error: called lookupManageTable with invalid table\n";
    exit(EXIT_FAILURE);
    // sprintf(outstr, "Program error: called lookupManageTable with invalid table\n");
    // logAndExit(logfp, outstr);
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

bool getSeriesVal(double& val, PPPG_SERIES_PARAM& series, int calMonth, int calYear, int k)
{
    // Long run or 1 year data.
    int i;
    if (series.oneYear) {
        // calYear is irrelevant, calMonth will range from 1 to 12, series elements will 
        // be indexed from 0 to 11. 
        i = calMonth - 1;
    }
    else {
        // Long run data. 
        i = (calYear - series.start) * 12 + calMonth - 1;
        if (i > series.vlen * 12 - 1) {
            // Should not happen as we will sanity check series before running. 
            // sprintf(outstr.c_str(), "Attempted lookup of series element %d in series %d, only %d entries in series.\nCheck Start/End Ages.",
            //   i, ser, (series->vlen*12-1)); 
            std::cout << "Attempted lookup of series element " << i << " in series " << series.id << ", only " << (series.vlen * 12 - 1) << " entries in series.\nCheck Start/End Ages." << std::endl;
            exit(EXIT_FAILURE);
            // logAndExit(logfp, outstr);
        }
        if (i < 0)
        {
            // sprintf(outstr, "Attempted lookup of year before start year\n");
            std::cout << "Attempted lookup of year before start year" << std::endl;
            exit(EXIT_FAILURE);
            // logAndExit(logfp, outstr);
        }
    }
    val = series.data[i].scanline[k];
    if (series.data[i].g->IsNoData(val)) {
        return false;
    }
    else {
        return true;
    }  
    return false; 

}
//----------------------------------------------------------------------------------

void readSampleFile(GDALRasterImage *refGrid, const std::vector<PPPG_OP_VAR>& opVars)
{
  // Read a text file of sample points, one per line, in the format idstring, 
  // xcoord, ycoord; find the index number of the cell the points fall in. 
  char *line, *fname;
  char *id, *xstr, *ystr, *cp;
  int ind=0;
  int opn;
  double lat, lon;
  int cellIndex;

  if (!samplePointsMonthly && !samplePointsYearly)
    return;

    if ((sampleIpFp = fopen(sampleIpFile.c_str(), "r")) == NULL) {
      std::cout << "Could not open sample point file " << sampleIpFile << std::endl;
      exit(EXIT_FAILURE);
      // sprintf(outstr.c_str(), "Could not open sample point file %s\n", sampleIpFile);
      //fprintf(logfp, outstr.c_str());
      //fprintf(stderr, outstr.c_str());
      // exit(1);
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
    if ((samplePoints[ind].fp = fopen(fname, "w")) == NULL) {
      std::cout << "Error opening output sample file " << fname << std::endl;
      exit(EXIT_FAILURE);
      // sprintf(outstr, "Error opening output sample file %s\n", fname);
      // logAndExit(logfp, outstr);
    }

    // Write header line for each sample file. 
    // For each output sample file
    fprintf(samplePoints[ind].fp, "year, month, id, ");
    for (opn = 1; opVars[opn].val != NULL; opn++)
    //for (opn = 1; opVars[opn].id != "-1"; opn++)
      fprintf(samplePoints[ind].fp, "%s, ", opVars[opn].id.c_str());
    fprintf(samplePoints[ind].fp, "\n");

    ind++;
  }
  // Make sure end of sample is marked. 
  samplePoints[ind].id[0] = 0;
}

//----------------------------------------------------------------------------------

bool readInputParam(const std::string& pName, std::vector<std::string> pValue, std::vector<PPPG_PARAM>& params) 
{
  // For a parameter name pName and a parameter value, pValue, both as strings, read 
  // the value into an appropriate variable. The parameter name can be either the 
  // same as the variable name, or it can be a long descriptive name, which matches 
  // the description of the parameter given in the VB version. In some cases other 
  // synonymns are allowed also. If pName does not match with any of the defined 
  // parameter names, return false. 
  int pInd=0;
  const std::string cp;

  if (pValue.size() > 1) {
    // Not an InputParam because input params are of format '"pName", pVal'
    // The pVal vector must only have one value
    return false;
  }

  // Find the index within parameters array of n.  
  // Allometric relationships & partitioning
  // TODO: this can be simplified by having global list of names attached to keys (e.g "pFS2":"pFS2")
  if (namesMatch("pFS2", pName) ||
           namesMatch("Foliage:stem partitioning ratio @ D=2 cm", pName)) pInd = pNameToInd("pFS2", params);
  else if (namesMatch("pFS20", pName) ||
           namesMatch("Foliage:stem partitioning ratio @ D=20 cm", pName)) pInd = pNameToInd("pFS20", params);
  else if (namesMatch("StemConst", pName) ||
           namesMatch("Constant in the stem mass v. diam. relationship", pName)) 
    pInd = pNameToInd("StemConst", params);
  else if (namesMatch("StemPower", pName) ||
           namesMatch("Power in the stem mass v. diam. relationship", pName)) pInd = pNameToInd("StemPower", params);
  else if (namesMatch("pRx", pName) ||
           namesMatch("Maximum fraction of NPP to roots", pName)) pInd = pNameToInd("pRx", params);
  else if (namesMatch("pRn", pName) ||
           namesMatch("Minimum fraction of NPP to roots", pName)) pInd = pNameToInd("pRn", params);

  // Temperature modifier (fT) | cardinal temperatures
  else if (namesMatch("growthTmin", pName) ||
           namesMatch("Minimum temperature for growth", pName)) pInd = pNameToInd("growthTmin", params);
  else if (namesMatch("growthTopt", pName) ||
           namesMatch("Optimum temperature for growth", pName)) pInd = pNameToInd("growthTopt", params);
  else if (namesMatch("growthTmax", pName) ||
           namesMatch("Maximum temperature for growth", pName)) pInd = pNameToInd("growthTmax", params);
    
  // Frost modifier
  else if (namesMatch("kF", pName) ||
           namesMatch("Days production lost per frost day", pName)) pInd = pNameToInd("kF", params);

  // Litterfall & root turnover
  else if (namesMatch("gammaFx", pName) ||
           namesMatch("Maximum litterfall rate", pName)) pInd = pNameToInd("gammaFx", params);
  else if (namesMatch("gammaF0", pName) ||
           namesMatch("Litterfall rate at t = 0", pName)) pInd = pNameToInd("gammaF0", params);
  else if (namesMatch("tgammaF", pName) ||
           namesMatch("Age at which litterfall rate has median value", pName)) 
    pInd = pNameToInd("tgammaF", params);
  else if (namesMatch("Rttover", pName) ||
           namesMatch("Average monthly root turnover rate", pName)) pInd = pNameToInd("Rttover", params);

  // conductances
  else if (namesMatch("MaxCond", pName) ||
           namesMatch("Maximum canopy conductance", pName)) pInd = pNameToInd("MaxCond", params);
  else if (namesMatch("CoeffCond", pName) ||
           namesMatch("Defines stomatal response to VPD", pName)) pInd = pNameToInd("CoeffCond", params);
  else if (namesMatch("BLcond", pName) ||
           namesMatch("Canopy boundary layer conductance", pName)) pInd = pNameToInd("BLcond", params);

  // fertility effects
  else if (namesMatch("m0", pName) ||
           namesMatch("Value of 'm' when FR = 0", pName)) pInd = pNameToInd("m0", params);
  else if (namesMatch("fN0", pName) ||
           namesMatch("Value of 'fNutr' when FR = 0", pName)) pInd = pNameToInd("fN0", params);
  else if (namesMatch("fNn", pName) ||
           namesMatch("Power of (1-FR) in 'fNutr'", pName)) pInd = pNameToInd("fNn", params);  //added 22-07-02
  
  // Soil water modifier (fSW) | soil characteristics
  else if (namesMatch("SWconst0", pName) ||
           namesMatch("Moisture ratio deficit for fq = 0.5", pName)) pInd = pNameToInd("SWconst0", params);
  else if (namesMatch("SWpower0", pName) ||
           namesMatch("Power of moisture ratio deficit", pName)) pInd = pNameToInd("SWpower0", params);

  // stem numbers
  else if (namesMatch("wSx1000", pName) ||
           namesMatch("Max. stem mass per tree @ 1000 trees/hectare", pName)) 
    pInd = pNameToInd("wSx1000", params);

  //Thinning Parameters 29-07-02
  else if (namesMatch("thinPower", pName) ||
           namesMatch("Power in self-thinning rule", pName)) pInd = pNameToInd("thinPower", params);
  else if (namesMatch("mF", pName) ||
           namesMatch("Fraction mean single-tree foliage biomass lost per dead tree", pName)) pInd = pNameToInd("mF", params);
  else if (namesMatch("mR", pName) ||
           namesMatch("Fraction mean single-tree root biomass lost per dead tree", pName)) pInd = pNameToInd("mR", params);
  else if (namesMatch("mS", pName) ||
           namesMatch("Fraction mean single-tree stem biomass lost per dead tree", pName)) pInd = pNameToInd("mS", params);
           
  // Age modifier (fAge)
  else if (namesMatch("MaxAge", pName) ||
           namesMatch("Maximum stand age used in age modifier", pName)) pInd = pNameToInd("MaxAge", params);
  else if (namesMatch("nAge", pName) ||
           namesMatch("Power of relative age in function for fAge", pName)) pInd = pNameToInd("nAge", params);
  else if (namesMatch("rAge", pName) ||
           namesMatch("Relative age to give fAge = 0.5", pName)) pInd = pNameToInd("rAge", params);

  // Canopy structure and processes | specific leaf area
  else if (namesMatch("SLA0", pName) ||
           namesMatch("Specific leaf area at age 0", pName)) pInd = pNameToInd("SLA0", params);
  else if (namesMatch("SLA1", pName) ||
           namesMatch("Specific leaf area for mature leaves", pName)) pInd = pNameToInd("SLA1", params);
  else if (namesMatch("tSLA", pName) ||
           namesMatch("Age at which specific leaf area = (SLA0+SLA1)/2", pName)) pInd = pNameToInd("tSLA", params);
  else if (namesMatch("k", pName) ||
           namesMatch("Extinction coefficient for absorption of PAR by canopy", pName)) pInd = pNameToInd("k", params);
  else if (namesMatch("fullCanAge", pName) ||
           namesMatch("Age at canopy cover", pName)) pInd = pNameToInd("fullCanAge", params);
  else if (namesMatch("alpha", pName) ||
           namesMatch("Canopy quantum efficiency", pName)) pInd = pNameToInd("alpha", params);

  // Branch and bark fraction (fracBB)
  else if (namesMatch("fracBB0", pName) ||
           namesMatch("Branch and bark fraction at age 0", pName)) pInd = pNameToInd("fracBB0", params);
  else if (namesMatch("fracBB1", pName) ||
           namesMatch("Branch and bark fraction for mature stands", pName)) pInd = pNameToInd("fracBB1", params);
  else if (namesMatch("tBB", pName) ||
           namesMatch("Age at which fracBB = (fracBB0+fracBB1)/2", pName)) pInd = pNameToInd("tBB", params);

  // various
  else if (namesMatch("y", pName) ||
           namesMatch("Ratio NPP/GPP", pName)) pInd = pNameToInd("y", params);
  else if (namesMatch("Density", pName) ||
           namesMatch("Basic density", pName)) pInd = pNameToInd("Density", params);

  //conversion factors - Addition 29th November Anders Siggins
  else if (namesMatch("Qa", pName) ||
           namesMatch("Intercept of net v. solar radiation relationship", pName)) pInd = pNameToInd("Qa", params);
  else if (namesMatch("Qb", pName) ||
           namesMatch("Slope of net v. solar radiation relationship", pName)) pInd = pNameToInd("Qb", params);
  else if (namesMatch("gDM_mol", pName) ||
           namesMatch("Molecular weight of dry matter", pName)) pInd = pNameToInd("gDM_mol", params);
  else if (namesMatch("molPAR_MJ", pName) ||
           namesMatch("Conversion of solar radiation to PAR", pName)) pInd = pNameToInd("molPAR_MJ", params);

  //Additional conversion factors
  else if (namesMatch("LAIgcx", pName) ||
           namesMatch("LAI for maximum canopy conductance", pName)) pInd = pNameToInd("LAIgcx", params);
  else if (namesMatch("MaxIntcptn", pName) ||
           namesMatch("Maximum proportion of rainfall evaporated from canopy", pName)) pInd = pNameToInd("MaxIntcptn", params);
  else if (namesMatch("LAImaxIntcptn", pName) ||
           namesMatch("LAI for maximum rainfall interception", pName)) pInd = pNameToInd("LAImaxIntcptn", params);

  // 3PG site parameters. 
  else if (namesMatch("Lat", pName) ||
           namesMatch("Latitude", pName)) pInd = pNameToInd("Lat", params);
  else if (namesMatch("FR", pName) ||
           namesMatch("Fertility rating", pName)) pInd = pNameToInd("FRp", params);
  else if (namesMatch("soilIndex", pName) ||
           namesMatch("Soil Index", pName) ||
           namesMatch("Soil class", pName)) pInd = pNameToInd("soilIndex", params);
  else if (namesMatch("MaxASW", pName) ||
           namesMatch("Maximum ASW", pName)) pInd = pNameToInd("MaxASW", params);
  else if (namesMatch("MinASW", pName) ||
           namesMatch("Minimum ASW", pName)) pInd = pNameToInd("MinASWp", params);

  // Initial conditions. 
  else if (namesMatch("StartAge", pName) || namesMatch("Initial age", pName) || 
           namesMatch("Start age", pName)) pInd = pNameToInd("StartAge", params);
  else if (namesMatch("EndAge", pName) ||
           namesMatch("End age", pName)) pInd = pNameToInd("EndAge", params);
  else if (namesMatch("StartMonth", pName) || namesMatch("Start Month", pName) || 
           namesMatch("Start month", pName)) pInd = pNameToInd("StartMonth", params);
  else if (namesMatch("SeedlingMass", pName) || namesMatch("Seedling Mass", pName) || 
           namesMatch("Seedling mass", pName)) pInd = pNameToInd("SeedlingMass", params);
  else if (namesMatch("FRstart", pName)) pInd = pNameToInd("FRstart", params);        //Fertility modifiers that depend on age
  else if (namesMatch("FRend", pName)) pInd = pNameToInd("FRend", params);
  else if (namesMatch("FRdec", pName)) pInd = pNameToInd("FRdec", params);
  else if (namesMatch("WFi", pName) ||
           namesMatch("W foliage", pName)) pInd = pNameToInd("WFi", params);
  else if (namesMatch("WRi", pName) ||
           namesMatch("W root", pName)) pInd = pNameToInd("WRi", params);
  else if (namesMatch("WSi", pName) ||
           namesMatch("W stem", pName)) pInd = pNameToInd("WSi", params);
  else if (namesMatch("StemNoi", pName) ||
           namesMatch("Stem no", pName)) pInd = pNameToInd("StemNoi", params);
  else if (namesMatch("ASWi", pName) ||
           namesMatch("Initial soil water", pName)) pInd = pNameToInd("ASWi", params);
  else if (namesMatch("MinASWTG", pName)) pInd = pNameToInd("MinASWTG", params);
  else if (namesMatch("rhoMin", pName) ||
           namesMatch("Minimum basic density - for young trees", pName)) pInd = pNameToInd("rhoMin", params);  //Standage varying density 15/07/2002
  else if (namesMatch("rhoMax", pName) ||   //Standage varying density 15/07/2002
           namesMatch("Maximum basic density - for older trees", pName)) pInd = pNameToInd("rhoMax", params);
  else if (namesMatch("tRho", pName)   ||
           namesMatch("Age at which rho = (rhoMin+rhoMax)/2", pName)) pInd = pNameToInd("tRho", params);      //Standage varying density 15/07/2002 
  else if (namesMatch("yearPlanted", pName) ||
           namesMatch("Year Planted", pName)) pInd = pNameToInd("yearPlanted", params);

  // 3PGS mode
  else if (namesMatch("NDVI_FPAR_intercept", pName)) pInd = pNameToInd("NDVI_FPAR_intercept", params);
  else if (namesMatch("NDVI_FPAR_constant", pName)) pInd = pNameToInd("NDVI_FPAR_constant", params);
    
  // If no index was found its not a basic input parameter. 
  if (pInd == 0)
    return false;

  // Soil Index/Class is a special case, if its specified with a character code, rewrite 
  // it as an integer.  Must be carefull not to match a grid name.  
  if ( pInd == pNameToInd( "soilIndex", params) ) {
      if (namesMatch("S", pValue.front()))
          pValue[0] = "1";
      else if (namesMatch("SL", pValue.front()))
          pValue[0] = "2";
      else if (namesMatch("CL", pValue.front()))
          pValue[0] = "3";
      else if (namesMatch("C", pValue.front()))
          pValue[0] = "4";
  }
  // if (sscanf(pValue.front().c_str(), "%lf", params[pInd].adr) == 1) {
  //   // fprintf(logfp, "   %-40s constant:  % 9.3f\n", params[pInd].id, 
  //   //   *(params[pInd].adr)); 
    
  // }
  // Is the first value a number?  If so, its a constant.  If not, its a grid name.
  try {

    double f = std::stod(pValue.front());
    params[pInd].val = f;
    params[pInd].data.spType = pScalar;
    params[pInd].got = 1;
    std::cout << "   " << params[pInd].id << "         constant: " << params[pInd].val << std::endl;
    //std::cout << "Found scalar input. Set spType: " << params[pInd].data.spType << std::endl;
    return true;
  }
  catch (std::invalid_argument const& e) {
    // Is the parameter a grid name (a string). 
    // REFERENCES: https://stackoverflow.com/questions/43114174/convert-a-string-to-std-filesystem-path
    // and https://stackoverflow.com/questions/51949/
    params[pInd].data.gridName = pValue.front();
    // catch filePath exceptions
    try {
      const std::filesystem::path filePath = params[pInd].data.gridName;
      if (filePath.extension() == ".tif") // Heed the dot.
      {
          params[pInd].data.spType = pTif;
          std::cout << "   " << params[pInd].id << "         grid: " << params[pInd].data.gridName << std::endl;
          params[pInd].got = 1;
          return true;
      }
      else
      {
          std::cout << filePath.filename() << " is an invalid type. File extension must be '.tif'" << std::endl;
          exit(EXIT_FAILURE);
          // logAndExit(logfp, outstr);
          // Output: "myFile.cfg is an invalid type"
      }
    }
    catch (std::filesystem::filesystem_error const& e) {
      std::cout << " " << pValue.front() << " could not be interpreted as a scalar or grid name" << std::endl;
      std::cout << "Error: " << e.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

//----------------------------------------------------------------------------------

bool readOutputParam(const std::string& pName, const std::vector<std::string>& pValue, int lineNo, std::vector<PPPG_OP_VAR>& opVars)
{
  // For a parameter name pName and a parameter value, pValue, both as strings, read 
  // the value into an appropriate variable. The parameter name can be either the 
  // same as the variable name, or it can be a long descriptive name.  
   int pInd, pInd1, pInd2;
  std::string cp;

  pInd = pInd1 = pInd2 = 0;

  // Find the index for the p1aram in the opVars array. 
  // 3PG only output variables. 
  if (namesMatch("StemNo", pName) ||
    namesMatch("Stocking density", pName)) pInd1 = opNameToInd("StemNo", opVars);
  else if (namesMatch("WF", pName) ||
    namesMatch("Weight of foliage", pName)) pInd1 = opNameToInd("WF", opVars);
  else if (namesMatch("WR", pName) ||
    namesMatch("Weight of roots", pName)) pInd1 = opNameToInd("WR", opVars);
  else if (namesMatch("WS", pName) ||
    namesMatch("Weight of stems", pName)) pInd1 = opNameToInd("WS", opVars);
  else if (namesMatch("TotalW", pName) ||
      namesMatch("Total weight", pName)) pInd1 = opNameToInd("TotalW", opVars);
  else if (namesMatch("MAI", pName) ||
    namesMatch("Mean Annual Increment", pName)) pInd1 = opNameToInd("MAI", opVars);
  else if (namesMatch("avDBH", pName) ||
    namesMatch("Average DBH", pName)) pInd1 = opNameToInd("avDBH", opVars);
  else if (namesMatch("BasArea", pName) ||
    namesMatch("Basal Area", pName)) pInd1 = opNameToInd("BasArea", opVars);
  else if (namesMatch("StandVol", pName) ||
    namesMatch("Stand volume", pName)) pInd1 = opNameToInd("StandVol", opVars);
  else if (namesMatch("Transp", pName) ||
           namesMatch("Transpiration", pName)) pInd1 = opNameToInd("Transp", opVars);
  else if (namesMatch("cTransp", pName))  pInd1 = opNameToInd("cTransp", opVars);
  else if (namesMatch("ASW", pName) ||
           namesMatch("Available Soil Water", pName)) pInd1 = opNameToInd("ASW", opVars);
  else if (namesMatch("EvapTransp", pName)) pInd1 = opNameToInd("EvapTransp", opVars);
  else if (namesMatch("cEvapTransp", pName)) pInd1 = opNameToInd("cEvapTransp", opVars);
  else if (namesMatch("LAIx", pName)) pInd1 = opNameToInd("LAIx", opVars);
  else if (namesMatch("ageLAIx", pName)) pInd1 = opNameToInd("ageLAIx", opVars);
  else if (namesMatch("GPP", pName) ||
    namesMatch("Gross Primary Production (tDM/ha)", pName)) pInd1 = opNameToInd("GPP", opVars);
  else if (namesMatch("cGPP", pName))  pInd1 = opNameToInd("cGPP", opVars);
  else if (namesMatch("TotalLitter", pName)) pInd1 = opNameToInd("TotalLitter", opVars);
  else if (namesMatch("cLitter", pName)) pInd1 = opNameToInd("cLitter", opVars);
  
  // 3PGS only output variables
  if (namesMatch("delWAG", pName) ||
    namesMatch("change in aboveground biomass (tDM/ha)", pName)) pInd2 = opNameToInd("delWAG", opVars);
  else if (namesMatch("cumWabv", pName) ||
    namesMatch("accumulated aboveground biomass (tDM/ha)", pName)) pInd2 = opNameToInd("cumWabv", opVars);

  // Check we only matched from the appropriate set. 
  if (modelMode3PGS) {
    if (pInd1 != 0) {
      //sprintf(outstr, "The output variable %s is not supported in 3PGS mode\n", pName);
      std::cout << "The output variable " << pName << " is not supported in 3PGS mode" << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
  }
  else {
    if (pInd2 != 0) {
      //sprintf(outstr, "The output variable %s is not supported in 3PG mode\n", pName);
      std::cout << "The output variable " << pName << " is not supported in 3PG mode" << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
  }
  if (modelMode3PGS)
    pInd = pInd2; 
  else
    pInd = pInd1; 

  // Output variables common to both modes. 
  if (namesMatch("NPP", pName) ||
    namesMatch("Net Primary Production (tDM/ha)", pName)) pInd = opNameToInd("NPP", opVars); //Modified 26/07/02
  else if (namesMatch("cNPP", pName))  pInd = opNameToInd("cNPP", opVars);
  else if (namesMatch("LAI", pName) ||
    namesMatch("Leaf Area Index", pName)) pInd = opNameToInd("LAI", opVars);
  else if (namesMatch("cLAI", pName)) pInd = opNameToInd("cLAI", opVars); //Added 7 November 2002
  else if (namesMatch("FRout", pName))    pInd = opNameToInd("FR", opVars); //Added 11/07/2002
  else if (namesMatch("PhysMod", pName))  pInd = opNameToInd("PhysMod", opVars);//Added 11/07/2002
  else if (namesMatch("alphaC", pName))  pInd = opNameToInd("alphaC", opVars);//Added 11/07/2002
  else if (namesMatch("fAge", pName))  pInd = opNameToInd("fAge", opVars); //Added 11/07/2002
  else if (namesMatch("fracBB", pName))  pInd = opNameToInd("fracBB", opVars); //Added 11/07/2002
  else if (namesMatch("WUE", pName))  pInd = opNameToInd("WUE", opVars); //Added 16/07/2002
  else if (namesMatch("cWUE", pName))  pInd = opNameToInd("cWUE", opVars);
  else if (namesMatch("CVI", pName))  pInd = opNameToInd("CVI", opVars); //Added 16/07/2002
  else if (namesMatch("cCVI", pName))  pInd = opNameToInd("cCVI", opVars);
  else if (namesMatch("fSW", pName))      pInd = opNameToInd("fSW", opVars);
  else if (namesMatch("fVPD", pName))     pInd = opNameToInd("fVPD", opVars);
  else if (namesMatch("fT", pName))       pInd = opNameToInd("fT", opVars);
  else if (namesMatch("fNutr", pName))    pInd = opNameToInd("fNutr", opVars);
  else if (namesMatch("fFrost", pName))   pInd = opNameToInd("fFrost", opVars);
  else if (namesMatch("APAR", pName))     pInd = opNameToInd("APAR", opVars);
  else if (namesMatch("APARu", pName))    pInd = opNameToInd("APARu", opVars);
  else if (namesMatch("MAIx", pName))    pInd = opNameToInd("MAIx", opVars);
  else if (namesMatch("ageMAIx", pName))    pInd = opNameToInd("ageMAIx", opVars);
  
  // Did we match a name?  
  if (pInd == 0)
    return false;

  if (pValue.empty()) {
    std::cout << "No grid name for param " << pName << " on line: " << lineNo << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, outstr);

  }
  if (pValue.size() > 5) {
      std::cout << "More than 5 value elements detected for param " << pName << " on line: " << lineNo << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
  }
  // First token in the pValue is the output grid filename, outPath and filename are concatenated for the full path
  opVars[pInd].gridName = outPath + pValue.front(); 

  const std::filesystem::path filePath = opVars[pInd].gridName;
  if (filePath.extension() == ".tif") // Heed the dot.
  {
      opVars[pInd].spType = pTif;
  }
  else
  {
      std::cout << filePath.filename() << " is an invalid filename. Found " << filePath.extension() << " but must be '.tif'" << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
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
    // TODO: send this to logfile?
    std::cout << "No recurring year output detected." << std::endl;
  }
  if (yearlyOutput == true) {
    // Look for start year
    try {
      opVars[pInd].recurStart = std::stoi(cp);
    }
    catch (std::invalid_argument const& e) {
      std::cout << "Expected an integer start year in recuring output specification on line " << lineNo << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);

    }
    // Look for interval
    try {
      cp = pValue.at(2);
      try {
        const int interval = std::stoi(cp);
        if (interval == 0)
        {
          std::cout << "Found interval of zero years in recuring output specification on line " << lineNo << ". Expected non-zero" << std::endl;
          exit(EXIT_FAILURE);
          // logAndExit(logfp, outstr);

        }
        opVars[pInd].recurYear = interval;
      }
      catch (std::invalid_argument const& e) {
        std::cout << "Expected an integer interval in recuring output specification on line " << lineNo << std::endl;
        exit(EXIT_FAILURE);
        // logAndExit(logfp, outstr);
      }
    }
    catch (const std::out_of_range& oor) {
      std::cout << "Found start year but no interval in recuring output specification on line " << lineNo << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
    // Look for 'monthly' or 'month' keywords. 
    try {
      cp = pValue.at(3);
      if (cp == "monthly")
        opVars[pInd].recurMonthly = true;
      else if (cp == "month") {
        opVars[pInd].recurMonthly = false;
        // If 'month', look for the month interger
        try {
          cp = pValue.at(4);
          try {
            opVars[pInd].recurMonth = std::stoi(cp);
            if (opVars[pInd].recurMonth == 0)
            {
              std::cout << "Found month of zero in recuring output specification on line " << lineNo << ". Expected non-zero" << std::endl;
              exit(EXIT_FAILURE);
              // logAndExit(logfp, outstr);
            }
          }
          catch (std::invalid_argument const& e) {
            std::cout << "Expected an integer month in recuring output specification on line " << lineNo << std::endl;
            exit(EXIT_FAILURE);
            // logAndExit(logfp, outstr);
          }
        }
        catch (const std::out_of_range& oor) {
          std::cout << "Found 'month' keyword but no month in recuring output specification on line " << lineNo << std::endl;
          exit(EXIT_FAILURE);
          // logAndExit(logfp, outstr);
        }
      }
      else {
        std::cout << "Unrecognised keyword \"" << cp << "\" on line " << lineNo << std::endl;
        exit(EXIT_FAILURE);
        // logAndExit(logfp, outstr);
      }
    }
    catch (const std::out_of_range& oor) {
      std::cout << "Found start year and interval but no keyword in recuring output specification on line " << lineNo << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
  }
  // Mark the variable for later writing
  opVars[pInd].write = true;
  std::cout << "   variable: " << opVars[pInd].id << "   grid: " << opVars[pInd].gridName << std::endl;
  // fprintf(logfp, "   variable: %-20s   grid: %-20s\n", opVars[pInd].id, opVars[pInd].gridName.c_str());
  if (opVars[pInd].recurStart) {
    std::cout << "      starting in " << opVars[pInd].recurStart << ", writing every " << opVars[pInd].recurYear << " years";
    // fprintf(logfp, "      starting in %4d, writing every %2d years", 
    //   opVars[pInd].recurStart, opVars[pInd].recurYear);
    if (opVars[pInd].recurMonthly)
      std::cout << ", with monthly values";
      // fprintf(logfp, ", with monthly values");
    else if ( opVars[pInd].recurMonth != 0 )
      std::cout << ", on the " << opVars[pInd].recurMonth << " month";
      // fprintf( logfp, ", on the %2d month", opVars[pInd].recurMonth ); 
    // fprintf(logfp, ".\n");
  }
  std::cout << std::endl;
  return true;
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
    if (pValue.empty()) {
      std::cout << "No output directory specified." << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
    else if (pValue.size() > 1) {
      std::cout << "More than one value element detected in output directory specification." << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
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
          exit(EXIT_FAILURE);
          // logAndExit(logfp, outstr);
        }
      }
      std::cout << "   output path: " << outPath << std::endl; // "   output path: %s\n
      // fprintf(logfp, "   output path: %s\n", outPath);
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
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
    else if (pValue.size() > 2) {
      std::cout << "More than two value elements detected in sample points file specification." << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
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
        exit(EXIT_FAILURE);
        // logAndExit(logfp, outstr);
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
      // logAndPrint(logfp, outstr);
      return false;
    }
    if (pValue.size() > 1) {
      std::cout << "More than one value element detected in model mode specification." << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
    if ("3PGS" == pValue.front())
      modelMode3PGS = true;
    else if ("3PG" == pValue.front())
      modelMode3PGS = false;
    else {
      std::cout << "Invalid value for parameter 'Model mode': " << pValue.front() << std::endl;
      exit(EXIT_FAILURE);
    //  sprintf(outstr, "Invalid value for parameter 'Model mode': %s\n", pValue);
      // logAndExit(logfp, outstr);
    }
    return true;
  }

  // Look for Point mode output file. 
  else if (namesMatch("point mode output file", pName)) {
    if (pValue.empty()) {
      std::cout << "No point mode output file specified." << std::endl;
      // logAndPrint(logfp, outstr);
      return false;
    }
    if (pValue.size() > 1) {
      std::cout << "More than one value element detected in point mode output file specification." << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
    cp = pValue.front();
    if ((pointModeFp = fopen(cp.c_str(), "w")) == NULL) {
      std::cout << "Could not open point mode output file " << cp << std::endl;
      exit(EXIT_FAILURE);
      // fprintf(stderr, "Could not open point mode output file %s\n", pValue.front().c_str());
      // exit (1);
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
        std::cout << filePath.filename() << " is an invalid filetype " << filePath.extension() << std::endl;
        exit(EXIT_FAILURE);
        // logAndExit(logfp, outstr);
    }
    return false;
  }
}

//----------------------------------------------------------------------------------

bool readInputManageParam(const std::string pName, std::ifstream& inFile, int &lineNo, std::vector<PPPG_MT_PARAM> managVars)
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
    //nRead = &nFertility; 
  }
  else if ( namesMatch( "Management: irrigation", pName ) ) {
    tab = IrrigMT; 
    tabName = "Irrigation MT";
    //nRead = &nIrrigation; 
  }
  else if ( namesMatch( "Management: MinASW", pName ) ) {
    tab = MinAswMT;
    tabName = "Min ASW MT";
    //nRead = &nMinAvailSW;
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
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
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
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
    // Read the second token, which is either a constant or a grid name.
    if( !readParam( tab[i].data, tTokens.back() )) {
      std::cout << "Could not read management table value at line " << lineNo << std::endl;
      std::cout << "   " << tTokens.back() << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, outstr);
    }
    else {
      tab[i].got = 1;
    }
    if (tab[i].data.spType == pScalar) {
      std::cout << "   " << tabName << " year: " << tab[i].year << "   value: " << tab[i].data.sval << std::endl;
    }
    else {
      std::cout << "   " << tabName << " year: " << tab[i].year << "   grid: " << tab[i].data.gridName << std::endl;
    } 
    i++; 
  }
  *nRead = i; 
  return true;  
}

//----------------------------------------------------------------------------------

bool readInputSeriesParam(std::string pName, std::vector<std::string> pValue, std::ifstream& paramFp, int &lineNo, std::vector<PPPG_SERIES_PARAM>& seriesV)
{
  // Read 'series' input parameters, ie climate and NDVI. 
  // Two styles of input are permitted. 
  // Firstly, the parameter name can be followed on the same line by 12 values, in this case the 12 values will 
  // be reused for each run year.
  // Secondly, the parameter name can be the only thing on the line (other than a comment),
  // on each following line must be a year followed by 12 values, until the sequence is terminated by a blank 
  // line.  The year values must be in ascending order.  
  int ser; 
  int series_yr, prev_yr; 
  std::string tok, line;

  for (auto &series : seriesV) {

      if (pName == "Tmax") {
          if (series.id != "Tmax_vals") {
              continue;
          }
      }
      else if (pName == "Tmin") {
          if (series.id != "Tmin_vals") {
              continue;
          }
      }
      else if (pName == "Tavg") {
          if (series.id != "Tavg_vals") {
              continue;
          }
      }
      else if (pName == "Rain") {
          if (series.id != "Rain_vals") {
              continue;
          }

      }
      else if (pName == "Solar Radtn") {
          if (series.id != "SolarRad_vals") {
              continue;
          }
      }
      else if (pName == "Frost days") {
          if (series.id != "FrostDays_vals") {
              continue;
          }
      }
      else if (pName == "NDVI_AVH") {
          if (series.id != "NdviAvh_vals") {
              continue;
          }
      }
      else if (pName == "Net radtn") {
          if (series.id != "NetRad_vals") {
              continue;
          }
      }
      else if (pName == "VPD") {
          if (series.id != "Vpd_vals") {
              continue;
          }
	  }
      else {
		  return false;
	  }

      //// Which series. 
      //if (      namesMatch ("Tmax",        pName ) ) ser = SS_TMAX; 
      //else if ( namesMatch ("Tmin",        pName ) ) ser = SS_TMIN;
      //else if ( namesMatch ("Tavg",        pName ) ) ser = SS_TAVG;
      //else if ( namesMatch ("Rain",        pName ) ) ser = SS_RAIN; 
      //else if ( namesMatch ("Solar Radtn", pName ) ) ser = SS_SOLARRAD; 
      //else if ( namesMatch ("Frost days",  pName ) ) ser = SS_FROSTDAYS; 
      //else if ( namesMatch ("NDVI_AVH",    pName ) ) ser = SS_NDVI_AVH; 
      //else if ( namesMatch ("Net radtn",   pName ) ) ser = SS_NETRAD; 
      //else if ( namesMatch ("VPD",         pName ) ) ser = SS_VPD; 
      //else return false; 

      //switch (ser) {
      //case SS_TMAX:      series = &Tmax_vals; break; 
      //case SS_TMIN:      series = &Tmin_vals; break; 
      //case SS_TAVG:      series = &Tavg_vals; break;
      //case SS_RAIN:      series = &Rain_vals; break; 
      //case SS_SOLARRAD:  series = &SolarRad_vals; break; 
      //case SS_FROSTDAYS: series = &FrostDays_vals; break; 
      //case SS_NDVI_AVH:  series = &NdviAvh_vals; break; 
      //case SS_NETRAD:    series = &NetRad_vals; break; 
      //case SS_VPD:       series = &Vpd_vals; break; 
      //default: break;
      //}

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
          series.oneYear = false;
      else
          // Second style above.  
          series.oneYear = true;

      // Read values for one year style. 
      if (series.oneYear) {
          series.vlen = 1;
          series.data = new PPPG_VVAL[series.vlen * 12];
          int i;

          for (i = 0; i < 12; i++) {
              try {
                  std::string mValue = pValue.at(i);
                  if (readParam(series.data[i], pValue[i])) {
                      if (series.data[i].spType == pScalar) {
                          std::cout << "   " << pName << " month " << i + 1 << " constant: " << series.data[i].sval << std::endl;
                          // fprintf(logfp, "   %-34s month %2d constant: %12.6f\n", pName, i+1, series->data[i].sval );
                      }
                      else if (series.data[i].spType == pTif) {
                          std::cout << "   " << pName << " month " << i + 1 << " grid: " << series.data[i].gridName << std::endl;
                          // fprintf(logfp, "   %-34s month %2d grid: %s\n", pName, i+1, series->data[i].gridName );
                      }
                  }
                  else {
                      std::cout << "Could not read parameter " << pName << " at month " << i + 1 << std::endl;
                      exit(EXIT_FAILURE);
                      // sprintf( outstr, "Could not read parameter %s.\n", pName ); 
                      // logAndExit( logfp, outstr ); 
                  }
              }
              catch (const std::out_of_range& oor) {
                  std::cout << "No value for " << pName << " at month " << i + 1 << std::endl;
                  exit(EXIT_FAILURE);
                  // sprintf(outstr, "Incomplete series on line %d.\n", lineNo); 
                  // logAndExit(logfp, outstr); 
              }
          }
          series.got = true;

          if (i < 12) {
              std::cout << "Incomplete series on line " << lineNo << std::endl;
              exit(EXIT_FAILURE);
          }
      }
      // Time series style
      else {
          // Find out how many years in the series. 
          int place = paramFp.tellg();
          int ss_lineNo = lineNo;
          prev_yr = -1;
          while (std::getline(paramFp, line)) {
              lineNo++;
              std::vector<std::string> sTokens;
              sTokens = boost::split(sTokens, line, boost::is_any_of(", \n\t"));
              if (sTokens.empty())
                  break;
              try {
                  series_yr = std::stoi(sTokens.front());
              }
              catch (const std::out_of_range& oor) {
                  std::cout << "Could not read year in series data at line " << lineNo << std::endl;
                  exit(EXIT_FAILURE);
                  // logAndExit(logfp, outstr); 
              }
              if (prev_yr < 0) {
                  prev_yr = series_yr - 1;
                  series.start = series_yr;
              }
              if (series_yr - 1 != prev_yr) {
                  std::cout << "Series year on line " << lineNo << " is not consecutive." << std::endl;
                  exit(EXIT_FAILURE);
                  // sprintf(outstr, "series year on line %d is not consecutive\n", lineNo);
                  // logAndExit(logfp, outstr);
              }
              prev_yr = series_yr;
          }
          series.vlen = series_yr - series.start + 1;

          // Allocate the space and read the series, have already checked the years. 
          paramFp.seekg(place);
          lineNo = ss_lineNo;
          series.data = new PPPG_VVAL[series.vlen * 12];
          for (int ss = 0; ss < series.vlen; ss++) {
              std::getline(paramFp, line);
              lineNo++;
              std::vector<std::string> sTokens;
              sTokens = boost::split(sTokens, line, boost::is_any_of(", \n\t"));
              if (sTokens.size() != 13) {
                  std::cout << "Could not read series on line " << lineNo << std::endl;
                  exit(EXIT_FAILURE);
                  // sprintf(outstr, "Cannot parse Year and 12 months on %d.\n", lineNo); 
                  // logAndExit(logfp, outstr); 
              }
              try {
                  series_yr = std::stoi(sTokens.front());
              }
              catch (std::invalid_argument const& inv) {
                  std::cout << "Could not read series year on line " << lineNo << std::endl;
                  exit(EXIT_FAILURE);
              }
              // Read the monthly values. 
              int si;
              for (int mn = 0; mn < 12; mn++) {
                  si = ss * 12 + mn;
                  tok = pValue.at(mn);
                  if (readParam(series.data[si], tok)) {
                      if (series.data[si].spType == pScalar)
                          std::cout << "   " << pName << " year " << series.start + ss << " month " << mn + 1 << " constant: " << series.data[si].sval << std::endl;
                      // fprintf(logfp, "   %-34s  %4d/%02d constant: %12.6f\n", pName, series->start + ss, mn+1, series->data[si].sval );
                      else if (series.data[si].spType == pTif)
                          std::cout << "   " << pName << " year " << series.start + ss << " month " << mn + 1 << " grid: " << series.data[si].gridName << std::endl;
                      // fprintf(logfp, "   %-34s  %4d/%02d grid: %s\n", pName, series->start + ss, mn+1, series->data[si].gridName );
                  }
                  else {
                      std::cout << "Could not read parameter " << pName << " at line " << lineNo << std::endl;
                      exit(EXIT_FAILURE);
                      // sprintf( outstr, "Could not read series on line %d.\n", lineNo); 
                      // logAndExit( logfp, outstr ); 
                  }
              }
          }
          series.got = true;
      }
  }

  return true;
}

//----------------------------------------------------------------------------------

void readParamFile(const std::string& paramFile, std::vector<PPPG_PARAM>& params, std::vector<PPPG_OP_VAR>& opVars, std::vector<PPPG_SERIES_PARAM>& series, std::vector<PPPG_MT_PARAM>& managment)
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
  int paramCount=0, lineLength=0, lineNo=0;
  int readingOutput=0;
  int len;

  // Parse file based on these rules:
  //    1. Skip blank lines and comments indicted by '//'
  //    2. First token is the parameter name
  //    3. Second and subsequent tokens are the parameter values
  std::ifstream inFile(paramFile);
  std::cout << "Reading input parameters from file '" << paramFile << "'..." << std::endl;
  while (std::getline(inFile, line))
  {
    lineNo++;
    if (line.empty())
        continue;
    if (line[0] == '/' && line[1] == '/')
        continue;
    std::vector<std::string> tokens;
    boost::split(tokens, line, boost::is_any_of(","), boost::token_compress_on);

    // trim leading whitespace from each token using boost::trim
    for (int i = 0; i < tokens.size(); i++)
        boost::trim(tokens[i]);
    std::string pName = tokens.front();
    boost::trim_if(pName, boost::is_any_of("\""));

    // Second and subsequent tokens are the parameter values, put them all into a vector
    std::vector<std::string> pValues;
    boost::split(pValues, tokens.at(1), boost::is_any_of(" \t"), boost::token_compress_on);

    // Check input against each parameter type.
    if (readInputParam(pName, pValues, params)) { continue; }
    else if (readOutputParam(pName, pValues, lineNo, opVars)) { continue; }
    else if (readOtherParam(pName, pValues)) { continue; }
    else if (readInputSeriesParam(pName, pValues, inFile, lineNo, series)) { continue; }
    //else if (readInputManageParam(pName, inFile, lineNo, managment)) { continue; } // SVZ: axed for now
    else {
        std::cout << "Cannot read parameter in file " << paramFile << ", line: " << lineNo << ": " << pName << std::endl;
        exit(EXIT_FAILURE);
    }
  }
}

//----------------------------------------------------------------------------------

bool haveAllParams(std::vector<PPPG_PARAM>& params, std::vector<PPPG_SERIES_PARAM>& series)
{ 
  // Check that we have read a value for all parameters. 
  int i, pInd;
  int indSeed, indWRi, indWFi, indWSi;  //indexes for values to check if required variable available
  bool missing = false;

  std::cout << "Checking all required parameters have been set.." << std::endl;

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
    "Lat", "FRp", "soilIndex", "MaxASW", "MinASWp",    // 3PG site parameters. 
    "StartAge", "EndAge",              //Initial conditions
    //"WFi", "WRi", "WSi",             //Now checked along with SeedlingMass
    "StemNoi", "ASWi", "yearPlanted",  // Initial conditions. 
    "Qa", "Qb",
    "gDM_mol", "molPAR_MJ",
    "LAIgcx", "MaxIntcptn",
   "StartMonth", 
    "LAImaxIntcptn", 
    "thinPower", "mF", "mR", "mS",    //Thinning coefficients
    ""
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
    "Lat", "FRp", "soilIndex", "MaxASW", "MinASWp",       // 3PG site parameters.
    "StartAge", "EndAge",                               // Initial conditions
    "NDVI_FPAR_intercept", "NDVI_FPAR_constant",        // FPAR from NDVI equation.
    "Qa", "Qb",
    "gDM_mol", "molPAR_MJ",
    "LAIgcx", "MaxIntcptn",
   "StartMonth", 
    "LAImaxIntcptn", 
    ""
  };

 /* {"Tmax_vals"},
  { "Tmin_vals" },
  { "Tavg_vals" },
  { "Rain_vals" },
  { "SolarRad_vals" },
  { "FrostDays_vals" },
  { "NdviAvh_vals" },
  { "NetRad_vals" },
  { "Vpd_vals" },*/

      // Temperature series
  bool tmax, tavg, tmin, frost, rain, vpd, solrad, netrad;
  tmax = tavg = tmin = frost = rain = vpd = false;
  for (const PPPG_SERIES_PARAM ser : series) {
      if( ser.id == "Tmax_vals") {
          if (ser.got) {
              tmax = true;
          }
	  }
	  else if (ser.id == "Tmin_vals") {
          if (ser.got) {
              tmin = true;
          }
      }
      else if (ser.id == "Tavg_vals") {
          if (ser.got) {
              tavg = true;
          }
	  }
      else if (ser.id == "Rain_vals") {
          if (ser.got) {
              rain = true;
          }
      }
      else if (ser.id == "FrostDays_vals") {
          if (ser.got) {
              frost = true;
          }
	  }
      else if (ser.id == "Vpd_vals") {
          if (ser.got) {
              vpd = true;
          }
	  }
      else if (ser.id == "NetRad_vals") {
          if (ser.got) {
              netrad = true;
          }
	  }
      else if (ser.id == "SolarRad_vals") {
          if (ser.got) {
              solrad = true;
          }
	  }
  }

  if (!tmax && !tavg) {
    std::cout << "No Tmax or Tavg data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Tmax or Tavg data");
  }
  if (!tmin && !tavg) {
    std::cout << "No Tmin or Tavg data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Tmin or Tavg data");
  }
  if (tavg) {
    std::cout << "Using Tavg series" << std::endl;
    if (!vpd) {
      std::cout << "No VPD series but VPD series is required if Tavg series is used" << std::endl;
      exit(EXIT_FAILURE);

    }
      // logAndExit(logfp, "A VPD series is required if an average temperature series is used\n"); 
  }
  if (!rain) {
    std::cout << "No Rain data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Rain data");
  }
  if (!solrad && !netrad) {
    std::cout << "No Solar or Net Radiation data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Solar or Net Radiation data");
  }
  if (!frost) {
    std::cout << "No Frost data" << std::endl;
    exit(EXIT_FAILURE);
    // logAndExit(logfp, "No Frost data");
  }

  //Check for SeedlingMass and WSi, WFi, WRi

  indSeed = pNameToInd("SeedlingMass", params);
  indWFi = pNameToInd("WFi", params);
  indWRi = pNameToInd("WRi", params);
  indWSi = pNameToInd("WSi", params);

  if (!params[indSeed].got && (!params[indWFi].got || !params[indWRi].got || !params[indWSi].got))
  {
    std::cout << "Missing parameter for 3PGS mode, " << params[indSeed].id << std::endl;
    // sprintf(outstr, "Missing parameter for 3PGS mode, %s\n", params[indSeed].id);
    // logAndPrint(logfp, outstr);
    missing = true;
  }

  //fill in Seedling Mass or others just in case - seems to cause errors if not there...
  if ( !haveSeedlingMass(params) )
  {
    // set the values of SeedlingMass ti
    params[indSeed].val = 0.0;
    params[indSeed].data.spType = pScalar;
    params[indSeed].got = 0;
  } else {
    params[indWFi].val = 0.0; 
    params[indWFi].data.spType = pScalar;
    params[indWFi].got = 0;
    params[indWRi].val = 0.0; 
    params[indWRi].data.spType = pScalar;
    params[indWRi].got = 0;
    params[indWSi].val = 0.0; 
    params[indWSi].data.spType = pScalar;
    params[indWSi].got = 0;
  }


  


  // Check required parameters for 3PGS
  if (modelMode3PGS) {
    for (const std::string &text : iParam3PGS) {
      pInd = pNameToInd(text, params);
      if (pInd != 0) {
        if (!params[pInd].got) {
          std::cout << "Missing parameter for 3PGS mode, " << params[pInd].id << std::endl;
          // sprintf(outstr, "Missing parameter for 3PGS mode, %s\n", params[pInd].id);
          // logAndPrint(logfp, outstr);
          missing = true;
        }
      }
    }
    if (!NdviAvh_vals.got) {
      std::cout << "No NDVI_AVH data" << std::endl;
      exit(EXIT_FAILURE);
      // logAndExit(logfp, "No NDVIAHV data");
    }
  }

  // Check required parameters for standard 3PG
  else {
    for (const std::string &text : iParam3PG) {
      pInd = pNameToInd(text, params);
      if (pInd != 0) {
        if (!params[pInd].got) {
          std::cout << "Missing parameter for 3PG mode, " << params[pInd].id << std::endl;
          // sprintf(outstr, "Missing parameter for 3PG mode, %s\n", params[pInd].id);
          // logAndPrint(logfp, outstr);
          missing = true;
        }
      }
    }
  }

  // Check various optional parameters
  //if ( !modelMode3PGS && NdviAvh_vals->got) {
  if ( !modelMode3PGS && NdviAvh_vals.got) {
    std::cout << "NdviAvh_vals not used in 3PGS mode." << std::endl;
    // sprintf(outstr, "NDVI_AVH not used in 3PGS mode.\n"); 
    // logAndPrint(logfp, outstr);
  }
  return !missing;
  
}

//----------------------------------------------------------------------------------

bool loadParamVals(int k, std::vector<PPPG_PARAM>& params)
{
  // Load all model parameter values into their global variables.
  // Spatial parameters are taken from the current grid cell, and
  // non-spatial parameters are left unchanged.  k is the cell index
  // in the data array element of a FloatGrid or ByteGrid object.
  int pn;
  GDALRasterImage *fg;
  float result;
  //char ErrorString[100];

  for (pn=1; params[pn].id != ""; pn++)  {
    if (params[pn].got == true) 
    {
        if (params[pn].data.spType == pTif) {
            fg = params[pn].data.g;
            if (fg == NULL)
            {
                std::cout << "Scanline not read: " << params[pn].data.gridName << std::endl;
                exit(EXIT_FAILURE);
                // sprintf(ErrorString, "Error reading grid: %s - File not open.\n", params[pn].data.gridName);
                // logAndExit(logfp, ErrorString); 
            }
            result = fg->GetVal(k);
            if (fg->IsNoData(result)) {
                return false;
            }
            else {
                params[pn].val = result;
            }
        }
    }
  }

    // Look for zero and do utterly bodgy things.
    //Note that any changes here should be echoed in findrunperiod.  
    //However, this is specifically for Aracruz...
    double yearPlanted, StartAge, EndAge, StartMonth, StemNoi;
    yearPlanted = params[pNameToInd("yearPlanted", params)].val;
    StartAge = params[pNameToInd("StartAge", params)].val;
    EndAge = params[pNameToInd("EndAge", params)].val;
    StartMonth = params[pNameToInd("StartMonth", params)].val;
    StemNoi = params[pNameToInd("StemNoi", params)].val;

    
    if (yearPlanted < 1 || isnan(yearPlanted) )
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
  if (vval.spType == pTif) {
    //std::cout << "   opening grid from " << vval.gridName << "..." << std::endl;
    // fprintf(logfp, "   opening grid from %s...", vval.gridName); 
    try {
      vval.g = new GDALRasterImage(vval.gridName, 0);
    }
    catch (const std::exception& e) {
      std::cout << "Could not open grid " << vval.gridName << std::endl;
      // fprintf(logfp, "Could not open grid %s\n", vval.gridName);
      exit(EXIT_FAILURE);
    }
  }
  else {
    return false;
  }
  return true; 
}

//----------------------------------------------------------------------------------
// void ResetGrids(void)
// {
//   //Only interested in resetting the findrunperiod grids which are not series.
//   for (int pn = 1; params[pn].data.spType != pNull; pn++) { // start at 1 to avoid error record. 
//     if ( params[pn].data.spType == pTif ) {
//       params[pn].data.g->ResetGrid();
//     }
//   }
// } 
//----------------------------------------------------------------------------------
void CloseGrids(std::vector<PPPG_PARAM>& params, std::vector<PPPG_OP_VAR>& opVars)
{
    //for (int pn = 1; params[pn].id != ""; pn++) { // start at 1 to avoid error record. 
    //    if ( params[pn].data.spType == pTif && params[pn].got == 1) {
    //      params[pn].data.g->Close();
    //    }
    //}
    //for (int op = 1; opVars[op].id != ""; op++) {  
    //    if (opVars[op].spType == pTif) {
    //      for (PPPG_VVAL* v : opVars[op].RO) {
    //          if (v->g != NULL) {
			 //   v->g->Close();
    //          }
		  //}
	   // }
    //}
}

GDALRasterImage* openInputGrids(std::vector<PPPG_PARAM>& params, std::vector<PPPG_SERIES_PARAM>& series, std::vector<PPPG_MT_PARAM>&  mgmnt)
{
  // Open any grids in the params array, the climate and NDVI series arrays, and 
  // the management tables. 
  // Copy the grid parameters of the first grid opened to refGrid. 
  int pn, j; 
  bool spatial = false, first = true;
  PPPG_SERIES_PARAM *ser; 
  GDALRasterImage *refGrid;

  //std::cout << "Opening input grids..." << std::endl;
  // fprintf(logfp, "Opening input grids...\n");
  // fprintf(stdout, "Opening input grids...\r");

  for (pn = 1; params[pn].id != ""; pn++) { // start at 1 to avoid error record. 
    if (params[pn].data.spType == pNull || params[pn].data.spType == pScalar){
        //do nothing
    }
    else if (params[pn].got != 1)
    {
      //do nothing
    }

    else if ( openGrid( params[pn].data ) ) {
      spatial = true; 
      if ( first ) {
        refGrid = (GDALRasterImage *)params[pn].data.g;
        first = false; 
      }
      else if ((fabs(refGrid->xMin - params[pn].data.g->xMin) > 0.0001)
         || (fabs(refGrid->yMin - params[pn].data.g->yMin) > 0.0001)
         || (fabs(refGrid->xMax - params[pn].data.g->xMax) > 0.0001)
         || (fabs(refGrid->yMax - params[pn].data.g->yMax) > 0.0001)
        || ( refGrid->nRows != params[pn].data.g->nRows ) 
        || ( refGrid->nCols != params[pn].data.g->nCols ) ) {
          std::cout << "Grid dimensions must match, grid " << params[pn].data.gridName << " differs from first grid." << std::endl;
          exit(EXIT_FAILURE);
        // sprintf(outstr, "Grid dimensions must match, grid %s differs from first grid.\n", 
        //   params[pn].data.gridName ); 
        // logAndExit(logfp, outstr); 
      }
    }
  }

  // Open all series grids. 
  // For each series. 
  for (PPPG_SERIES_PARAM ser: series) {
    for (int i = 0; i < ser.vlen * 12; i++) {
      if ( openGrid( ser.data[i] ) ) {
        spatial = true; 
        if ( first ) {
          refGrid = (GDALRasterImage *)(ser.data[i].g);
          first = false;
        }
        else if ( ( fabs( refGrid->xMin - ser.data[i].g->xMin ) > 0.0001 ) 
          || ( fabs( refGrid->yMin - ser.data[i].g->yMin ) > 0.0001 )
          || ( fabs( refGrid->xMax - ser.data[i].g->xMax ) > 0.0001 )
          || ( fabs( refGrid->yMax - ser.data[i].g->yMax ) > 0.0001 ) 
          || ( refGrid->nRows != ser.data[i].g->nRows ) 
          || ( refGrid->nCols != ser.data[i].g->nCols ) ) {
            std::cout << "Grid dimensions must match, grid " << ser.data[i].gridName << " differs from first grid." << std::endl;
            exit(EXIT_FAILURE);
          // sprintf(outstr, "Grid dimensions must match, grid %s differs from first grid.\n", 
          //   ser->data[i].gridName ); 
          // logAndExit(logfp, outstr); 
        }
      }
    }
  }

  // Open all management table grids. 
  //PPPG_MT_PARAM *tab;
  //PPPG_MT_PARAM *tablist[] = { FertMT, IrrigMT, MinAswMT, NULL }; 
  // for (PPPG_SERIES_PARAM ser: series) {
  // SVZ: axed for not, re-introduce later
  //for (PPPG_MT_PARAM tab: mgmnt) {
  //  for (int i = 0; tab[i].year > 0; i++) {
  //    if ( openGrid( tab[i].data ) ) {
  //      spatial = true; 
  //      if ( first ) {
  //        refGrid = (GDALRasterImage *)tab[i].data.g;
  //        first = false; 
  //      }
  //      else if ( ( fabs( refGrid->xMin - tab[i].data.g->xMin ) > 0.0001 ) 
  //        || ( fabs( refGrid->yMin - tab[i].data.g->yMin ) > 0.0001 )
  //        || ( fabs( refGrid->xMax - tab[i].data.g->xMax ) > 0.0001 )
  //        || ( fabs( refGrid->yMax - tab[i].data.g->yMax ) > 0.0001 ) 
  //        || ( refGrid->nRows != tab[i].data.g->nRows ) 
  //        || ( refGrid->nCols != tab[i].data.g->nCols ) ) {
  //          std::cout << "Grid dimensions must match, grid " << tab[i].data.gridName << " differs from first grid." << std::endl;
  //          exit(EXIT_FAILURE);
  //        // sprintf(outstr, "Grid dimensions must match, grid %s differs from first grid.\n", 
  //        //   tab[i].data.gridName ); 
  //        // logAndExit(logfp, outstr); 
  //      }
  //    }
  //  }
  //}

  if (!spatial) {
    std::cout << "None" << std::endl;
    // fprintf(logfp, "none\n");
    refGrid = NULL; 
  }


  return refGrid;
}

//----------------------------------------------------------------------------------

int createOutputGrids(GDALRasterImage *refGrid, std::vector<PPPG_OP_VAR>& opVars)
{
  // Loop through opVars array and open output grid objects for all 
  // those variables marked for output.  
  int opn; 

  // Open ordinary output grids. 
  std::cout << "Creating output grids..." << std::endl;
  for (opn = 0; opVars[opn].id != ""; opn++) {
    if (opVars[opn].write) {
      std::cout << "   tif grid " << opVars[opn].gridName << std::endl;
      // Open output grid at gridName with same class attributes as refGrid.
      GDALRasterImage *g = new GDALRasterImage(opVars[opn].gridName, refGrid);

      if (g->Exists(opVars[opn].gridName)) {
        std::cout << "Error, output grid named: " << opVars[opn].gridName << " already exists." << std::endl;
        exit(EXIT_FAILURE);
      }
      g->Close();
      // *(opVars[opn].g) = *refGrid;
      // opVars[opn].g->Allocate();
    }
  }
  return EXIT_SUCCESS;
}

int threadOpenOutputTIFs(std::vector<PPPG_OP_VAR>& opVars)
{
    int opn, j;

    std::cout << "Opening output grids..." << std::endl;
    for (opn = 0; opVars[opn].id != ""; opn++) {
        if (opVars[opn].write) {

            try {
                opVars[opn].g = new GDALRasterImage(opVars[opn].gridName, 1);
            }
            catch (const std::exception& e) {
                std::cout << "Could not open grid " << opVars[opn].gridName << std::endl;
                // fprintf(logfp, "Could not open grid %s\n", vval.gridName);
                exit(EXIT_FAILURE);
            }
        }
    }

    return EXIT_SUCCESS;
}

int threadOpenInputTIFs(std::vector<PPPG_PARAM>& params, std::vector<PPPG_SERIES_PARAM>& series)
{
    int pn, j;

    // Open all TIFs found in Inputs
    for (pn = 1; params[pn].id != ""; pn++) { // start at 1 to avoid error record. 
        if (params[pn].data.spType == pNull || params[pn].data.spType == pScalar) {
            //do nothing
        }
        // Else if opened, got, and a TIF type -- open it for the threads
        else if (params[pn].data.spType == pTif && params[pn].got == 1)
        {
            try {
                params[pn].data.g = new GDALRasterImage(params[pn].data.gridName, 0);
            }
            catch (const std::exception& e) {
                std::cout << "Failed to input TIF " << params[pn].data.gridName << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    // Open all Series TIFs
    for (PPPG_SERIES_PARAM ser : series) {
        for (int i = 0; i < ser.vlen * 12; i++) {
            try {
                ser.data[i].g = new GDALRasterImage(ser.data[i].gridName, 0);
            }
            catch (const std::exception& e) {
                std::cout << "Failed to series TIF " << params[pn].data.gridName << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    return EXIT_SUCCESS;
}

int writeOutputGrids(bool hitNODATA, long cellIndex, const std::vector<PPPG_OP_VAR>& opVars) {
  /*
    Write the final output grids. 
  */
  int opn;
  GDALRasterImage* fp;
  float fval;

  for (opn = 0; opVars[opn].id != ""; opn++) {
      if (opVars[opn].write) {
          fval = (float)(opVars[opn].val);
          fp = opVars[opn].g;
          if (fp != NULL) {
              if (hitNODATA) {
                  fval = fp->noData;
              }
              //std::cout << "Setting value " << fval << std::endl;
              fp->SetVal(cellIndex, fval);

          }
          else {
              std::cout << "ERROR: Output grid " << opVars[opn].id << " is NULL... could not write." << std::endl;
              exit(EXIT_FAILURE);
          }
      }
  }
  return EXIT_SUCCESS;
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
int threadOpenRegularOutputTIFs(MYDate spMinMY, MYDate spMaxMY, std::vector<PPPG_OP_VAR>& opVars) {
    int calYear, calMonth;
    int opn, roArrayLength;
    //FILE *fp; 
    std::string roStemName;
    std::string fname;
    std::string* cp;
    int len;
    int maxCy, maxCm, minCy, minCm;
    int mx;
    //int runMx; 
    int runMonths;

    maxCy = spMaxMY.year;
    maxCm = spMaxMY.mon;
    minCy = spMinMY.year;
    minCm = spMinMY.mon;

    runMonths = (maxCy - minCy - 1) * 12;
    // Then how many months in the (possibly) partial years at the beginning and end. 
    runMonths += 12 - minCm + 1;
    runMonths += maxCm;

    roArrayLength = (maxCy - minCy + 1) * 12 + 12;

    for (opn = 0; opVars[opn].id != ""; opn++) {
        // Is it marked for recurring output. 
        if (opVars[opn].recurYear == -1)
            continue;

        // Find stem name for regular output files, just lacking the 
        // year/month/extension. 
        // cp = strrchr(opVars[opn].gridName, '.');
        // len = cp - opVars[opn].gridName;
        // strncpy( roStemName, opVars[opn].gridName, len ); 
        // roStemName[len] = '\0'; 
        roStemName = opVars[opn].gridName.substr(0, opVars[opn].gridName.find_last_of("."));
        // Monthly output. mx indexes over entire ro array. 
        for (mx = 0; mx < roArrayLength; mx++) {
            calYear = minCy + (mx / 12);
            calMonth = (mx % 12) + 1;

            // Skip months before the model starts. 
            if (calYear == minCy && calMonth < minCm)
                continue;

            // Skip years after the model ends. 
            if (calYear > maxCy)
                continue;

            // Skip months after the model ends. 
            if (calYear == maxCy && calMonth > maxCm)
                continue;

            // Reached first output year? 
            if (calYear < opVars[opn].recurStart)
                continue;

            // Is this an output year? 
            if (((calYear - opVars[opn].recurStart) % opVars[opn].recurYear) != 0)
                continue;

            // Is this an annual output?  If so only want 1 month as given in recurMonth. 
            // Construct the filename by concatenating the stem name, year, month and extension.
            // Open a GDALRasterImage object for the grid
            fname = roStemName + std::to_string(calYear) + std::to_string(calMonth) + ".tif";
            GDALRasterImage* fg;
            if (!opVars[opn].recurMonthly) {
                if (calMonth == opVars[opn].recurMonth) {
                    //std::cout << "   opening regular out grid from " << fname << "..." << std::endl;
                    // fprintf(logfp, "   opening grid from %s...", fname.c_str());
                    GDALRasterImage og = GDALRasterImage(fname, 1);
                    opVars[opn].RO[mx].g = &og;
                }
            }
            else {
                std::cout << "   opening regular output grid from " << fname << "..." << std::endl;
                GDALRasterImage ogm = GDALRasterImage(fname, 1);
                opVars[opn].RO[mx].g = &ogm;
            }
        }

    }
    return EXIT_SUCCESS;
}

int openRegularOutputGrids( GDALRasterImage *refGrid, MYDate spMinMY, MYDate spMaxMY, std::vector<PPPG_OP_VAR>& opVars)
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
  //FILE *fp; 
  std::string roStemName;
  std::string fname;
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
  for (opn = 0; opVars[opn].id != ""; opn++) {
    // Is it marked for recurring output. 
    if ( opVars[opn].recurYear == -1)
      continue; 
    
    // Allocate and null out array of *FILE. 
     //std::vector<GDALRasterImage*> regOut;
    opVars[opn].RO = new PPPG_VVAL[roArrayLength];
    for (int i = 0; i < roArrayLength; i++) {
        opVars[opn].RO[i].g = NULL;
    }
    
    // Find stem name for regular output files, just lacking the 
    // year/month/extension. 
    // cp = strrchr(opVars[opn].gridName, '.');
    // len = cp - opVars[opn].gridName;
    // strncpy( roStemName, opVars[opn].gridName, len ); 
    // roStemName[len] = '\0'; 
    roStemName = opVars[opn].gridName.substr(0, opVars[opn].gridName.find_last_of("."));
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
      // Construct the filename by concatenating the stem name, year, month and extension.
      // Open a GDALRasterImage object for the grid
      fname = roStemName + std::to_string(calYear) + std::to_string(calMonth) + ".tif";
      // if there is a leading slash, remove it.
  //    if (fname[0] == '/')
		//fname = fname.substr(1);
  //    if (fname[0] == '\\')
  //      fname = fname.substr(1);
      GDALRasterImage* fg;
      if ( !opVars[opn].recurMonthly ) {
        if ( calMonth == opVars[opn].recurMonth ) {
          //std::cout << "   opening regular out grid from " << fname << "..." << std::endl;
          // fprintf(logfp, "   opening grid from %s...", fname.c_str());
            try {

                fg = new GDALRasterImage(fname, refGrid);
                opVars[opn].RO[mx].g = fg;
                fg->Close();

            }
            catch (const std::exception& e) {
                std::cout << "Could not open grid " << fname << std::endl;
                // fprintf(logfp, "Could not open grid %s\n", fname.c_str());
                exit(EXIT_FAILURE);
            }
          // pointer to GDALRasterImage object for this year/month combo
          // Construct annual regular output grid name and open file. 
          // Open a GDALRasterImage object for the grid
          // sprintf( fname, "%s%4d%02d.flt", roStemName, calYear, calMonth ); 
          // if (())
          // if ( ( fp = fopen( fname, "wb" ) ) == NULL ) {
          //   sprintf(outstr, "Couldn't open regular output file %s\n", fname);
          //   logAndExit(logfp, outstr);
          // }
          // // if ( ! copyHeader( refGrid, fname ) ) {
          // //   sprintf(outstr, "Couldn't open regular output header file %s\n", fname);
          // //   logAndExit(logfp, outstr);
          // // }
          // opVars[opn].RO[mx] = fp; // pointer to start of month in array.
        }
      }
      // Is a monthly output ro grid and this is the output year. 
      else {
        std::cout << "   opening regular output grid from " << fname << "..." << std::endl;
        try {

            fg = new GDALRasterImage(fname, refGrid);
            opVars[opn].RO[mx].g = fg;
            fg->Close();

        }
        catch (const std::exception& e) {
            std::cout << "Could not open grid " << fname << std::endl;
            // fprintf(logfp, "Could not open grid %s\n", fname.c_str());
            exit(EXIT_FAILURE);
        }
        // pointer to GDALRasterImage object for this year/month combo
        // Construct monthly regular output grid names and open file. 
        // sprintf( fname, "%s%4d%02d.flt", roStemName, calYear, calMonth ); 
        // if ( ( fp = fopen( fname, "wb" ) ) == NULL ) {
        //   sprintf(outstr, "Couldn't open regular output file %s\n", fname);
        //   logAndExit(logfp, outstr);
        // }
        // // if ( ! copyHeader( refGrid, fname ) ) {
        // //   sprintf(outstr, "Couldn't open regular output header file %s\n", fname);
        // //   logAndExit(logfp, outstr);
        // // }
        // opVars[opn].RO[mx] = fp; 
      }
    } // closes for ( mx = 0; mx < roArrayLength; mx++ )

  } // closes for each output variable.

  return EXIT_SUCCESS;
}
  
//----------------------------------------------------------------------------------

// int writeOutputGrids(void)
// {
//   int opn, gridsWritten=0;

//   fprintf(stdout, "Writing output grids...\r"); 
//   fprintf(logfp, "Writing output grids...\n");
//   for (opn=0; opVars[opn].id != "-1"; opn++) {
//     if (opVars[opn].write) {
//       try {
//         opVars[opn].g->Write(opVars[opn].gridName);
//       } catch (Exception &e) {
//         sprintf(outstr, "Exception: %s\n", e.Message());
//         return 0;
//       }
//       gridsWritten++;
//     }
//   }
//   return gridsWritten;
// }

//----------------------------------------------------------------------------------

//void writeSampleFiles(int cellIndex, int month, long calYear) 
//{
//  int opn, sInd, i;
//  static bool firstTime = true;
//
//  // Do we want to sample this point? 
//  sInd = -1;
//  for (i = 0; samplePoints[i].id[0] != 0; i++) 
//    if (samplePoints[i].cellIndex == cellIndex) {
//      sInd = i;
//      break;
//    }
//  if (sInd == -1)
//    return;
//
//  fprintf(samplePoints[sInd].fp, "%d, %d, %s, ", calYear, 
//    month, samplePoints[sInd].id.c_str());
//  // For each variable
//  for (opn = 1; opVars[opn].val != NULL; opn++)
//    fprintf(samplePoints[sInd].fp, "%f, ", (opVars[opn].val));
//  fprintf(samplePoints[sInd].fp, "\n");
//}

//----------------------------------------------------------------------------------

void writeMonthlyOutputGrids(const std::vector<PPPG_OP_VAR> opVars, int calYear, int calMonth, bool hitNODATA,
                              MYDate minMY, MYDate maxMY, long cellIndex )
{
  int opn, mx, maxInd; 
  float fval; 
  long lval; 
  GDALRasterImage* fp;
  GDALRasterImage *fg;

  int crap; 
  if ( calYear == 2008  || calYear == 2020 )
    crap = 1; 

  // Number of elements in the regular output array, for sanity checking. 
  // Usually not all of these would be used.  
  maxInd = ( maxMY.year - minMY.year + 1 ) * 12 + 12;  // ANL see comment in openRegularOutputGrids

  // for each output variable. 
  for (opn = 0; opVars[opn].id != ""; opn++) {
    // Is it marked for recurring output. 
    if ( opVars[opn].recurYear == -1 )
      continue; 
    
    // Skip variables not marked for recurring output. 
    if ( !opVars[opn].recurStart ) 
      continue;
    
    // Skip variables not marked for monthly recurrence. 
  /*  if ( !opVars[opn].recurMonthly )
      continue; */

    mx = ( calYear - minMY.year ) * 12 + ( calMonth - 1 );
    
    if (mx < 0)
      continue;
    
    // Sanity check. 
    if ( mx > maxInd ) {
      std::cout << "Program error, mx=" << mx << " too high in writeMonthlyOutputGrids at month/year " << calMonth << "/" << calYear << std::endl;
      exit(EXIT_FAILURE);
      // sprintf(outstr, "Program error, mx=%d too high in writeMonthlyOutputGrids at month/year %2d/%4d.\n", 
      //   mx, calYear, calMonth);
      // logAndExit(logfp, outstr);
    }
    fp = opVars[opn].RO[mx].g;
    if (fp != NULL) {
        if (opVars[opn].spType == pTif) {

            fval = (float)(opVars[opn].val);
            fg = opVars[opn].RO[mx].g;

            if (hitNODATA) {
                fval = fg->noData;
            }
            //std::cout <<d "Writing value " << std::endl;
            fp->SetVal(cellIndex, fval);
        }
    }
      /*sprintf(outstr, "Year: %d, Month: %d, Value: %f\n", calYear, calMonth, fval); 
      logOnly(logfp, outstr);*/
    }
 }

//----------------------------------------------------------------------------------

void writeYearlyOutputGrids( const std::vector<PPPG_OP_VAR> opVars, int calYear, int calMonth, bool hitNODATA,
                            MYDate minMY, MYDate maxMY, long cellIndex )
{
  int opn, mx, maxInd; 
  float fval; 
  long lval; 
  GDALRasterImage* fp;
  GDALRasterImage *fg;

  // Number of elements in the regular output array, for sanity checking. 
  maxInd = ( maxMY.year - minMY.year + 1 ) * 12; 

  for (opn = 0; opVars[opn].id != ""; opn++) {
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
      std::cout << "Program error, mx=" << mx << " too high in writeYearlyOutputGrids at month/year " 
        << calMonth << "/" << calYear << std::endl;
      exit(EXIT_FAILURE);
      // sprintf(outstr, "Program error, mx=%d too high in writeYearlyOutputGrids at month/year %2d/%4d.\n", 
      //   mx, calYear, calMonth);
      // logAndExit(logfp, outstr);
    }

    // Get file pointer for this regular output grid. 
    fp = opVars[opn].RO[mx].g;
    if (fp != NULL) {
        if (opVars[opn].spType == pTif) {
            fval = (float)opVars[opn].val;
            if (hitNODATA) {
                fval = fg->noData;
            }
            fg = (GDALRasterImage*)opVars[opn].g;
            fp->SetVal(cellIndex, fval);
        }
    }
  }
}
    
//----------------------------------------------------------------------------------
// NOTE: I _think_ that this is covered by the writeYearlyOutputGrids function above.
//      If so, this function can be deleted.

// void saveVariableVals(long k, bool hitNODATA)
// {
//   // Copy variable values from the model to the current grid cell. 
//   int opn;
//   GDALRasterImage *fg;

//   for (opn = 1; opVars[opn].id != "-1"; opn++) 
//     if (opVars[opn].spType == pTif) {
//       fg = (GDALRasterImage *)opVars[opn].g;
//       fg->z[k] = (float)*(opVars[opn].adr);
//       if (hitNODATA)
//         fg->z[k] = fg->noData;
//     }
// }

//----------------------------------------------------------------------------------

//void writeStandSummaryData(FILE *outFP, int year)
//{
//    fprintf(outFP, "%4.0f, %4.0f, %4.0f, %5.2f, %5.2f, %5.2f, %5.1f, %4.2f, "
//    "%4.2f, %3.1f, %4.2f, %4.0f, %4.0f\n", 
//    StandAge + yearPlanted, StandAge, StemNo, WF, WR, WS, StandVol, LAI, MAI, avDBH, 
//    cLitter, cumTransp, ASW);
//}

//----------------------------------------------------------------------------------

//void writeStandSummary(int year)
//{
//  // For point mode, write summary values to the output file pointModeFp. 
//
//  // 3PGS
//  if (modelMode3PGS) {
//      // headings
//      if (year == 1)
//        fprintf(pointModeFp, "Year, LAI, , NPP, delWAG, cumWabv\n");
//
//      // write initial conditions
//      fprintf(pointModeFp, "%4.0f, %3.1f, %5.2f, %5.2f, %5.2f\n", 
//        yearPlanted + StandAge - 1, LAI, NPP, delWAG, cumWabv);
//  }
//
//  // 3PG
//  else {
//    // Headings and initial conditions. 
//    if (year == 0) {
//        fprintf(pointModeFp, "Year, Stand age, Stem number, Wf, Wr, Ws, Stand volume, Stand LAI, "
//          "Stand MAI, Average DBH, Litterfall, Total transpiration, Soil water\n");
//
//        fprintf(pointModeFp, "%4.0f, %4.0f, %4.0f, %5.2f, %5.2f, %5.2f, %5.1f, %4.2f, %4.2f, %3.1f\n", 
//          yearPlanted + StandAge, StandAge, StemNoi, WFi, WRi, WSi, 1.7 * WSi, 
//          LAIi, MAIi, avDBHi);
//    }
//
//    // write annual data
//    writeStandSummaryData (pointModeFp, year);
//
//    // write peak LAI and MAI
//    if (StandAge == EndAge)
//      fprintf(pointModeFp, " , , , , , , ,%5.2f, %5.2f\n", LAIx, MAIx); 
//  }
//}

//----------------------------------------------------------------------------------

// Could possibly replace this logic with use of min/max? Would avoid having to read every pixel.
int findRunPeriod( GDALRasterImage *refGrid, MYDate &minMY, MYDate &maxMY, std::vector<PPPG_PARAM>& params)
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
    double ypI = pNameToInd("yearPlanted", params);
    double saI = pNameToInd("StartAge", params);
    double eaI = pNameToInd("EndAge", params);
    double smI = pNameToInd("StartMonth", params);
    double yearPlanted = params[ypI].val;
    double StartAge = params[saI].val;
    double EndAge = params[eaI].val;
    double StartMonth = params[smI].val;

    if ( !haveSpatialRunYears(params) ) {
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
      nrows = refGrid->nRows; 
      ncols = refGrid->nCols; 
      
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

          hitNoData = !loadParamVals(cellIndex, params);
          yearPlanted = params[ypI].val;
          StartAge = params[saI].val;
          EndAge = params[eaI].val;
          StartMonth = params[smI].val;

          
          // Look for NODATA, to skip it. 
          if ( hitNoData )
          {
      //     logAndPrint(logfp, "Nodata\n");
            continue; 
          }
          // Look for zero and do utterly bodgy things. 
          if ( params[ypI].val < 1)
            continue;   //Treat years less than 1900 as nodata 
          if ( StartAge < 1 )
            StartAge = 1; 
          if ( EndAge < 1 )
            EndAge = 2; // Mimum end growth for a stand - just a reasonable value. AS.
          if ( StartMonth < 1 ) {
            StartMonth = 1; 
          }

          // Check reasonableness of StartMonth, must be in the range 1 - 12 inclusive. 
          if ( StartMonth < 0.999 || StartMonth > 12.0001 ) {
            std::cout << "Found bad StartMonth value at XXXX" << std::endl;
            exit(EXIT_FAILURE);
            // logAndExit(logfp, "Found bad StartMonth value at XXXX\n" ); 
          }
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
      if ( minCy == MINSEED ) {
        std::cout << "Failed to find valid starting year. Exiting..." << std::endl;
        // logAndExit( logfp, "Failed to find valid starting year.\n" );
        exit(EXIT_FAILURE);
      }
      if ( maxCy == MAXSEED ) {
        std::cout << "Failed to find valid ending year. Exiting..." << std::endl;
        // logAndExit( logfp, "Failed to find valid ending year.\n" );
        exit(EXIT_FAILURE);
      }
      if ( minCm == 13 ) {
        std::cout << "Failed to find valid starting month. Exiting..." << std::endl;
        // logAndExit( logfp, "Failed to find valid starting month.\n" );
        exit(EXIT_FAILURE); 
      }
      if ( maxCm == 13 ) {
        std::cout << "Failed to find valid ending month. Exiting..." << std::endl;
        // logAndExit( logfp, "Failed to find valid ending month.\n" ); 
        exit(EXIT_FAILURE);
      }

      minMY.year = minCy; 
      minMY.mon = minCm; 

      maxMY.year = adjAge;        //minCy + EndAge - StartAge; 
      maxMY.mon = maxCm; 
    }
    return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------------
//To ensure that if we don't have a variable, it doesn't get an accidental GOT value
//Initialise the params list so that everything in missinIg.
std::vector<PPPG_PARAM> InitInputParams(void)
{
    // Initialize the params vector with id values.
    std::vector<PPPG_PARAM> inputParams
    {
      {"paramError", NULL},
      {"pFS2"},
      {"pFS20"},
      {"StemConst"},
      {"StemPower"},
      {"pRx"},
      {"pRn"},

      // Temperature modifier (fT) | cardinal temperatures
      // ANL - these have been renamed from just Tmax etc, to avoid confusion with the 
      // climate variables. 
      {"growthTmin"},
      {"growthTopt"},
      {"growthTmax"},

      // Frost modifier
      {"kF"},

      // Litterfall & root turnover
      {"gammaFx"},
      {"gammaF0"},
      {"tgammaF"},
      {"Rttover"},

      // conductances
      {"MaxCond"},
      {"CoeffCond"},
      {"BLcond"},

      // fertility effects
      {"m0"},
      {"fN0"},
      {"fNn"},

      //Thinning effects
      {"thinPower"},
      {"mF"},
      {"mR"},
      {"mS"},

      // Soil water modifier (fSW) | soil characteristics
      {"SWconst0"},
      {"SWpower0"},

      // stem numbers
      {"wSx1000"},

      // Age modifier (fAge)
      {"MaxAge"},
      {"nAge"},
      {"rAge"},

      // Canopy structure and processes | specific leaf area
      {"SLA0"},
      {"SLA1"},
      {"tSLA"},
      {"k"},
      {"fullCanAge"},
      {"alpha"},
      {"fracBB0"},
      {"fracBB1"},
      {"tBB"},

      // various
      {"y"},
      {"rhoMin"},
      {"rhoMax"},
      {"tRho"},             // Standage varying density 3-06-02 

      //Conversions
      {"Qa"},
      {"Qb"},
      {"gDM_mol"},
      {"molPAR_MJ"},

      //Additional conversion factors 
      {"LAIgcx"},
      {"MaxIntcptn"},
      {"LAImaxIntcptn"},

      // 3PG site parameters. 
      {"Lat"},
      {"FRp"},
      {"FRstart"},  //These three variables relate to fertility decrease with age
      {"FRend"},
      {"FRdec"},
      {"soilIndex"},
      {"MaxASW"},
      {"MinASWp"},

      // Initial conditions. 
      {"StartAge"},
      {"EndAge"},
      {"StartMonth"},
      {"yearPlanted"},  /* CHECK! do we still use this?*/
      {"SeedlingMass"},
      {"WFi"},
      {"WRi"},
      {"WSi"},
      {"StemNoi"},
      {"ASWi"},
      {"MinASWTG"},
      //  {"yearPlanted",  &yearPlanted},  /* This has been moved as strange errors were occuring with grids here*/

        // ANL - extras for 3PGS mode
        {"NDVI_FPAR_intercept"},
        {"NDVI_FPAR_constant"},
        {"", NULL}  // NULL entries used to mark array ends. 
    };

  //The following parameters need to be set to scalar or grids do not open.
  int pn;
  pn = pNameToInd("FRstart", inputParams);
  inputParams[pn].data.spType = pScalar;
  pn = pNameToInd("FRend", inputParams);
  inputParams[pn].data.spType = pScalar;
  pn = pNameToInd("FRdec", inputParams);
  inputParams[pn].data.spType = pScalar;

  return inputParams;
}

std::vector<PPPG_OP_VAR> initOutputVars()
{
    std::vector<PPPG_OP_VAR> outputVars {
          {"opVarError", NULL},
          {"StemNo"},
          {"WF"},
          {"WR"},
          {"WS"},
          {"TotalW"},
          {"LAI"},
          {"cLAI"},
          {"MAI"},
          {"avDBH"},
          {"BasArea"},
          {"StandVol"},
          {"GPP"},
          {"cGPP"},
          {"NPP"},
          {"cNPP"},
          {"delWAG"},
          {"cumWabv"},
          {"Transp"},
          {"cTransp"},
          {"ASW"},
          {"fSW"},
          {"fVPD"},
          {"fT"},
          {"fNutr"},
          {"fFrost"},
          {"APAR"},
          {"APARu"},
          {"EvapTransp"},
          {"cEvapTransp"}, //Added 08/11/02
          {"LAIx"},
          {"ageLAIx"},
          {"MAIx"},    //Added 29/07/2002
          {"ageMAIx"}, //Added 29/07/2002
          {"FR"},     //Added 11/07/2002
          {"PhysMod"}, //Added 11/07/2002
          {"alphaC"},  //Added 11/07/2002
          {"fAge"},    //Added 11/07/2002
          {"fracBB"},
          {"WUE"},     //Added 16/07/02
          {"cWUE"},    //Added 08/11/02
          {"CVI"},     //Added 16/07/02
          {"cCVI"},    //Added 08/11/02
          {"TotalLitter"}, //Added 16/07/02
          {"cLitter"},
          {"",         NULL}
    };
    return outputVars;


}

std::vector<PPPG_SERIES_PARAM> initSeriesParams()
{
    std::vector<PPPG_SERIES_PARAM> seriesParams{
		{"Tmax_vals"},
		{"Tmin_vals"},
		{"Rain_vals"},
        {"SolarRad_vals"},
        {"FrostDays_vals"},
        {"NdviAvh_vals"},
        {"NetRad_vals"},
        {"Vpd_vals"},
		{"Tavg_vals"},
	};
	return seriesParams;
}

std::vector<PPPG_MT_PARAM>  initMTParams()
{
    std::vector<PPPG_MT_PARAM> mtParams {
        {"FertMT"},
        {"IrrigMT"},
        {"MinAswMT"}
	};
	return mtParams;
}

//----------------------------------------------------------------------------------

bool havePointOpFile()
{
  if (pointModeFp == NULL) {
    // sprintf(outstr, "Using point mode but parameter \"point mode output file\" not set\n");
    std::cout << "Using point mode but parameter point mode output file not set" << std::endl;
    // logAndPrint(logfp, outstr);
    return false;
  }
  else
    return true;
}
//----------------------------------------------------------------------------------
bool haveSeedlingMass(const std::vector<PPPG_PARAM>& params)
{
  int pInd;
  pInd = pNameToInd("SeedlingMass", params);
  if (params[pInd].got)
    return true;
  else
    return false;
}

//----------------------------------------------------------------------------------
bool haveSpatialRunYears(const std::vector<PPPG_PARAM>& params)
{
  int pInd;

  pInd = pNameToInd("StartAge", params);
  if (!(params[pInd].data.spType == pScalar)) return true;

  pInd = pNameToInd("EndAge", params);
  if (!(params[pInd].data.spType == pScalar)) return true;

  pInd = pNameToInd("StartMonth", params);
  if (!(params[pInd].data.spType == pScalar)) return true;

  pInd = pNameToInd("yearPlanted", params);
  if (!(params[pInd].data.spType == pScalar)) return true;

  return false;

}

//----------------------------------------------------------------------------------
bool seriesNotNull(PPPG_SERIES_PARAM& series) {
    if (series.data != NULL)
		return true;
	else
		return false;
}
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

bool haveMinASWTG(const std::vector<PPPG_PARAM>& params)
{
  int pInd;
  bool bResult;
  pInd = pNameToInd("MinASWTG", params);
  
  if (params[pInd].got)
    bResult = true;
  else
    bResult = false;
  
  return bResult;
}
//----------------------------------------------------------------------------------

bool haveRhoMin(const std::vector<PPPG_PARAM>& params)
{
  int pInd;
  bool bResult;
  pInd = pNameToInd("rhoMin", params);
  
  if (params[pInd].got)
    bResult = true;
  else
    bResult = false;
  
  return bResult;  
}
//----------------------------------------------------------------------------------

bool haveRhoMax(const std::vector<PPPG_PARAM>& params)
{
  int pInd;
  bool bResult;
  pInd = pNameToInd("rhoMax", params);
  
  if (params[pInd].got)
    bResult = true;
  else
    bResult = false;
  
  return bResult;  
}
//----------------------------------------------------------------------------------

bool haveTRho(const std::vector<PPPG_PARAM>& params)
{
  int pInd;
  bool bResult;
  pInd = pNameToInd("tRho", params);
  
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

bool haveAgeDepFert(const std::vector<PPPG_PARAM>& params)
{
  int pInd;
  
  pInd = pNameToInd("FRstart", params);
  if (params[pInd].got == 0)
    return false;
  
  pInd = pNameToInd("FRend", params);
  if (params[pInd].got == 0)
    return false;
  
  pInd = pNameToInd("FRdec", params);
  if (params[pInd].got == 0)
    return false;

  return true;
}


void readInputScanlines(std::vector<PPPG_PARAM>& inputs, std::vector<std::vector<float>>& inputScanlines, int row, int ncols) {

    int inputsCounter = 0;
    for (PPPG_PARAM& i : inputs) {
        if (i.got == 1) {
            if (i.data.spType == pTif) {
                i.data.g->GetVal(0, row, ncols, 1, inputScanlines.at(inputsCounter).data());
                i.data.scanline = inputScanlines.at(inputsCounter);
                inputsCounter++;
            }
        }
    }
}

void readSeriesScanlines(std::vector<PPPG_SERIES_PARAM>& series, std::vector<std::vector<float>>& seriesScanlines, int row, int ncols) {
    int seriesCounter = 0;
    for (PPPG_SERIES_PARAM s : series) {
        if (s.got == 1) {
            if (s.oneYear) {
                int months = 12;
                for (int m = 0; m < months; ++m) {
                    s.data[m].g->GetVal(0, row, ncols, 1, seriesScanlines.at(seriesCounter).data());
                    s.data[m].scanline = seriesScanlines.at(seriesCounter);
                    seriesCounter++;
                }
            }
            else {
                int length = (s.vlen * 12);
                for (int l = 0; l < length; ++l) {
                    s.data[l].g->GetVal(0, row, ncols, 1, seriesScanlines.at(seriesCounter).data());
                    s.data[l].scanline = seriesScanlines.at(seriesCounter);
                    seriesCounter++;
                }
            }
        }
    }

}