// static char rcsid[] = "$Id: The_3PG_Model.cpp,v 1.6 2001/08/02 06:51:42 lou026 Exp $";

/*
All source code remains the property and copyright of CSIRO. 

CSIRO accepts no responsibility for the use of 3PG(S) or of the model 3-PG in
the form supplied or as subsequently modified by third parties. CSIRO disclaims
liability for all losses, damages and costs incurred by any person as a result
of relying on this software. 
Use of this software assumes agreement to this condition of use
*/


// Model routine for 3PG.  This is basically a translation of the VB 
// version to C, particularly in the routine runTreeModel.  Most 
// of the comments are from the VB version.  
// 
// Because of this VB basis, there are a few non-standard aspects to 
// this program.  Some arrays are indexed from 1, which means that 
// the first entry (element 0) is unused.  There are a *lot* of global 
// variables, in particular the model parameters, declared in this file. 

#include <math.h>
#include <iostream>
#include <stdio.h>
#include "Data_io.hpp"

//____________________________________
//
// The code for 3PG - March 24th, 2000
//____________________________________


//Changes from August 1999 version:
//
//   1) Accumulating annual stand transpiration
//   2) Introduced minimum avail soil water, with difference made up by
//      irrigation. Can output monthly and annual irrigation.
//   3) Start mid year in southern hemisphere
//   4) Recast alpha(FR) as alpha*fNutr
//   5) Some change in how functions are parameterised to make parameters
//      more meaningful
//   6) Allometric relationships based on DBH in cm
//   7) Partioning parameterised by pFS for DBH = 2 and 20 cm
//   8) Model made strictly state-determined by ensuring that LAI,
//      partitioning, etc determined from current state not a lagged state.
//   9) SLA made stand-age dependent
//  10) Non-closed canopy allowed for (not good!)
//  11) Manner in which modifiers taken into account corrected

//NOTE: The following conversion factors are used:
//
//    1 MJ  = 2.3 mol PAR
//    1 mol = 24 gDM


#define ModelVsn "3-PG March2000.24"

#define Pi 3.141592654
#define ln2 0.693147181

#define eps 0.0001

// Controls and counters
//int StartAge, EndAge;                  // age of trees at start/end of run
//int StartMonth;                        // month of year to start run
//int yearPlanted;                       // year trees planted
// ANL changed these three from int to double
double StartAge = 0.0, EndAge = 0.0;
double StartMonth = 0.0;
double yearPlanted = 0.0;
int DaysInMonth[13] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
bool showDetailedResults = false;
bool showStandSummary = false;
bool modelMode3PGS = false;

char siteName[100] = { 0 }; // Initialize with null characters
double Lat = 1000.0;
double MaxASW = 0.0, MinASW = 0.0, MinASWp = 0.0;
double FR = 0.0, FRp = 0.0;
double FRstart = 0.0, FRend = 0.0, FRdec = 0.0;
double soilIndex = 0.0;
double SWconst = 0.0, SWpower = 0.0;

int nFertility = 0;
int nMinAvailSW = 0;
int nIrrigation = 0;
double Irrig = 0.0;

double mYears = 1.0; // Default value 1.0 as per your code
double mDayLength[13] = { 0.0 };
double mFrostDays[13] = { 0.0 };
double mSolarRad[13] = { 0.0 };
double mTx[13] = { 0.0 };
double mTn[13] = { 0.0 };
double mTav[13] = { 0.0 };
double mVPD[13] = { 0.0 };
double mRain[13] = { 0.0 };
double mNetRad[13] = { 0.0 };

char SpeciesName[100] = { 0 }; // Initialize with null characters
double SeedlingMass = 0.0;
double StandAge = 0.0;
double ASW = 0.0, ASWi = 0.0;
double MinASWTG = 0.0;
double StemNoi = 0.0, StemNo = 0.0;
double WFi = 0.0, WF = 0.0;
double WRi = 0.0, WR = 0.0;
double WSi = 0.0, WS = 0.0;
double LAIi = 0.0, LAI = 0.0;
double MAIi = 0.0, MAI = 0.0;
double avDBHi = 0.0, avDBH = 0.0;
double TotalW = 0.0;
double BasArea = 0.0;
double StandVol = 0.0;
double TotalLitter = 0.0;
double LAIx = 0.0, ageLAIx = 0.0;
double MAIx = 0.0, ageMAIx = 0.0;
double cumTransp = 0.0;
double cumIrrig = 0.0;

double SLA = 0.0;
double gammaF = 0.0;
double fracBB = 0.0;
double CanCover = 0.0;

double MaxAge = 0.0;
double gammaFx = 0.0, gammaF0 = 0.0, tgammaF = 0.0;
double Rttover = 0.0;
double SLA0 = 0.0, SLA1 = 0.0, tSLA = 0.0;
double fullCanAge = 0.0;
double k = 0.0;
double pFS2 = 0.0, pFS20 = 0.0;
double StemConst = 0.0, StemPower = 0.0;
double SWconst0 = 0.0, SWpower0 = 0.0;
double Interception = 0.0;
double BLcond = 0.0;
double MaxCond = 0.0, CoeffCond = 0.0;
double y = 0.0;
double growthTmax = 0.0, growthTmin = 0.0, growthTopt = 0.0;
double wSx1000 = 0.0;
double thinPower = 0.0;
double mF = 0.0, mR = 0.0, mS = 0.0;
double m0 = 0.0, fN0 = 0.0, fNn = 0.0;
double alpha = 0.0;
double pRx = 0.0, pRn = 0.0;
double nAge = 0.0, rAge = 0.0;
double kF = 0.0;
double fracBB0 = 0.0, fracBB1 = 0.0, tBB = 0.0;
double Density = 0.0;
double rhoMin = 0.0, rhoMax = 0.0, tRho = 0.0;
double pfsConst = 0.0, pfsPower = 0.0;
double PhysMod = 0.0;

double Qa = 0.0, Qb = 0.0;
double gDM_mol = 0.0;
double molPAR_MJ = 0.0;

double LAIgcx = 0.0;
double MaxIntcptn = 0.0;
double LAImaxIntcptn = 0.0;

double m = 0.0, alphaC = 0.0, epsilon = 0.0;
double RAD = 0.0, PAR = 0.0, RADint = 0.0;
double lightIntcptn = 0.0;
double fAge = 0.0, fT = 0.0, fFrost = 0.0;
double fVPD = 0.0, fSW = 0.0, fNutr = 0.0;
double CanCond = 0.0;
double Transp = 0.0, EvapTransp = 0.0, RainIntcptn = 0.0, WUE = 0.0;
double AvStemMass = 0.0;
double APAR = 0.0, APARu = 0.0;
double GPPmolc = 0.0, GPPdm = 0.0, NPP = 0.0;
double pR = 0.0, pS = 0.0, pF = 0.0, pFS = 0.0;
double delWF = 0.0, delWR = 0.0, delWS = 0.0, delStems = 0.0;
double delLitter = 0.0, delRloss = 0.0;
double monthlyIrrig = 0.0;
double CVI = 0.0;

double cumGPP = 0.0, cumWabv = 0.0;
double abvgrndEpsilon = 0.0, totalEpsilon = 0.0;
double StemGrthRate = 0.0;
double cumEvapTransp = 0.0;
double CumdelWF = 0.0, CumdelWR = 0.0, CumdelWS = 0.0;
double CumAPARU = 0.0, cumARAD = 0.0;
double CumStemLoss = 0.0;
double CutStemMass1 = 0.0, CutStemMass2 = 0.0, CutStemMass3 = 0.0;

double cLAI = 0.0;
double cRADint = 0.0;
double aRADint = 0.0;
double cGPP = 0.0;
double aGPP = 0.0;
double cNPP = 0.0;
double aNPP = 0.0;
double cStemDM = 0.0;
double aStemDM = 0.0;
double cCVI = 0.0;
double cLitter = 0.0;
double cWUE = 0.0;
double aWUE = 0.0;
double aSupIrrig = 0.0;
double cSupIrrig = 0.0;
double cRainInt = 0.0;
double aTransp = 0.0;
double cTransp = 0.0;
double aEvapTransp = 0.0;
double cEvapTransp = 0.0;
double cEpsilonGross = 0.0;
double aEpsilonGross = 0.0;
double cEpsilonStem = 0.0;
double aEpsilonStem = 0.0;

double mNDVI[13];
double delWAG = 0.0;
double NDVI_FPAR_intercept = 0.0, NDVI_FPAR_constant = 0.0;

extern bool samplePointsMonthly;
extern bool samplePointsYearly;


// ANL - globals defined in 3pg.cpp
//extern FILE *logfp;
extern char outstr[];

//-----------------------------------------------------------------------------

double getDayLength(double Lat, int dayOfYear)
{ 
  // gets fraction of day when sun is "up"
  double sLat, cLat, sinDec, cosH0;
  
  sLat = sin(Pi * Lat / 180);
  cLat = cos(Pi * Lat / 180);
  
  sinDec = 0.4 * sin(0.0172 * (dayOfYear - 80)); 
  cosH0 = -sinDec * sLat / (cLat * sqrt(1 - pow(sinDec,2))); 
  if (cosH0 > 1)
    return 0;
  else if (cosH0 < -1) 
    return 1;
  else {
    return acos(cosH0) / Pi;
  }
}

//-----------------------------------------------------------------------------
//Anders Siggins 29/22/01
double Minimum(double X, double Y)
{
  if (X > Y) 
    return Y;
  else
    return X;
}
//And one that is probably not used...
double Maximum(double X, double Y)
{
  if (X < Y) 
    return Y;
  else
    return X;
}
//-----------------------------------------------------------------------------

double getVPD(double Tx, double Tn)
{
  //gets daily "mean" VPD - based on daily max and min temperatures only
  double VPDx, VPDn;
  VPDx = 6.1078 * exp(17.269 * Tx / (237.3 + Tx));
  VPDn = 6.1078 * exp(17.269 * Tn / (237.3 + Tn));
  return (VPDx - VPDn) / 2.0;  
}

//-----------------------------------------------------------------------------

double getMortality(double oldN, double oldW) 
{
//This function determines the number of stems to remove to ensure the
//self-thinning rule is satisfied. It applies the Newton-Rhapson method
//to solve for N to an accuracy of 1 stem or less. To change this,
//change the value of "accuracy".
//This was the old mortality function:
//  getMortality = oldN - 1000 * (wSx1000 * oldN / oldW / 1000) ^ (1 / thinPower)
//which has been superceded by the following ...
 
  int i;
  double fN, dfN;
  double dN, n, x1, x2;
  double result;

  n = oldN / 1000;
  x1 = 1000 * mS * oldW / oldN;
  i = 0;
  while (true)
  {  
    i = i + 1;
    x2 = wSx1000 * pow(n, (1 - thinPower));
    fN = x2 - x1 * n - (1 - mS) * oldW;
    dfN = (1 - thinPower) * x2 / n - x1;
    dN = -fN / dfN;
    n = n + dN;
    if ((fabs(dN) <= eps) || (i >= 5))
      break;
  }
  
  result = oldN - 1000 * n;
  return result;
}
//-----------------------------------------------------------------------------

double CanopyTranspiration(double Q, double VPD, double h, 
                           double gBL, double gC, double NR, bool hNRS)
//Penman-Monteith equation for computing canopy transpiration
//in kg/m2/day, which is conmverted to mm/day.
//The following are constants in the PM formula (Landsberg & Gower, 1997)
{
  double const e20 = 2.2;          // rate of change of saturated VP with T at 20C
  double const rhoAir = 1.2;       // density of air, kg/m3
  double const lambda = 2460000;   // latent heat of vapourisation of H2O (J/kg)
  double const VPDconv = 0.000622; // convert VPD to saturation deficit = 18/29/1000
  double netRad;
  double defTerm;
  double div; 
  double Etransp;
  double CT;                       //Canopy Transpiration return value
  
  if (hNRS)
    netRad = NR;
  else
    netRad = Qa + Qb * (Q * pow(10,  6)) / h;                // Q in MJ/m2/day --> W/m2
  
  defTerm = rhoAir * lambda * (VPDconv * VPD) * gBL;
  div = (1 + e20 + gBL / gC);
  Etransp = (e20 * netRad + defTerm) / div;           // in J/m2/s
  CT = Etransp / lambda * h;         // converted to kg/m2/day
  
  return CT;
}


//-----------------------------------------------------------------------------

void assignDefaultParameters(void)
{ 
  // We don't actually use this currently.  
  // Public
  MaxAge = 50;          // Determines rate of "physiological decline" of forest
  SLA0 = 4;             // specific leaf area at age 0 (m^2/kg)
  SLA1 = 4;             // specific leaf area for mature trees (m^2/kg)
  tSLA = 2.5;           // stand age (years) for SLA = (SLA0+SLA1)/2
  fullCanAge = 0;       // Age at full canopy cover
  k = 0.5;              // Radiation extinction coefficient
  gammaFx = 0.03;       // Coefficients in monthly litterfall rate
  gammaF0 = 0.001;
  tgammaF = 24;
  Rttover = 0.015;      // Root turnover rate per month
  SWconst0 = 0.7;       // SW constants are 0.7 for sand,0.6 for sandy-loam,
                       //   0.5 for clay-loam, 0.4 for clay
  SWpower0 = 9;         // Powers in the eqn for SW modifiers are 9 for sand,
                       //   7 for sandy-loam, 5 for clay-loam and 3 for clay
  Interception = 0.15;  // Proportion of rainfall intercepted by canopy
  MaxCond = 0.02;       // Maximum canopy conductance (gc, m/s)
  BLcond = 0.2;         // Canopy boundary layer conductance, assumed constant
  CoeffCond = 0.05;     // Determines response of canopy conductance to VPD
  y = 0.47;             // Assimilate use efficiency
  growthTmax = 32;            // "Critical" biological temperatures: max, min
  growthTmin = 2;             //   and optimum. Reset if necessary/appropriate
  growthTopt = 20;
  kF = 1;               // Number of days production lost per frost days
  pFS2 = 1;             // Foliage:stem partitioning ratios for D = 2cm
  pFS20 = 0.15;         //   and D = 20cm
  StemConst = 0.095;    // Stem allometric parameters
  StemPower = 2.4;
  pRx = 0.8;            // maximum root biomass partitioning
  pRn = 0.25;           // minimum root biomass partitioning
  m0 = 0;               // Value of m when FR = 0
  fN0 = 1;              // Value of fN when FR = 0
  alpha = 0.055;        // Canopy quantum efficiency
  wSx1000 = 300;     // Max tree stem mass (kg) likely in mature stands of 1000 trees/ha
  nAge = 4;             // Parameters in age-modifier
  rAge = 0.95;
  fracBB0 = 0.15;       // branch & bark fraction at age 0 (m^2/kg)
  fracBB1 = 0.15;       // branch & bark fraction for mature trees (m^2/kg)
  tBB = 1.5;            // stand age (years) for fracBB = (fracBB0+fracBB1)/2
  Density = 0.5;        // Basic density (t/m3)
}

//-----------------------------------------------------------------------------

void Initialisation(void)
{
  // Private
  //double Tanav;
  //int i;

  //Tanav = 0;
  //for (i=1; i<=12; i++)
  //  Tanav = Tanav + (mTx[i] + mTn[i]) / 24;   
  // ANL - does this actually work, ie calculate annual average temp?  
  // Is that what its supposed to do?  As we don't actually use Tanav 
  // anywhre its a moot point, and I've commented it out. 

  // At initialisation param file has only possible value to use.  
  MinASW = MinASWp; 

// If Tanav <= 13.4 Then gammaFx = gammaFx * (1 + 0.04 * (Tanav - 13.4) ^ 4)
  // Assign the SWconst and SWpower parameters for this soil class
  if (soilIndex != 0) {
    SWconst = 0.8 - 0.1 * soilIndex;
    SWpower = 11 - 2 * soilIndex;
  }
  else {
    SWconst = SWconst0;
    SWpower = SWpower0;
  }
  // Derive some parameters
  pfsPower = log(pFS20 / pFS2) / log(10);
  pfsConst = pFS2 / pow(2, pfsPower);
  // Initial ASW must be between min and max ASW
  if (ASWi <= MinASW) 
    ASWi = MinASW;
  else if (ASWi >= MaxASW) 
    ASWi = MaxASW;
  // Initialise ages
  MAIx = 0;
  LAIx = 0;

  //Basic density
  if (haveRhoMax() == false)
    rhoMax = 0.5;         //basic density for older trees (t/m3)
  if (haveRhoMin() == false)
    rhoMin = 0.5;         //ratio of basic density of young to old trees
  if (haveTRho() == false)
    tRho = 4;             //age at which density = average of old and young values
}

//-----------------------------------------------------------------------------

//Standage function translated from March beta of Excel 3-PG
//StartAge and StandAge are considered global variables...
void GetStandAge(void)
{

  //Assign initial stand age
  if (StartAge < yearPlanted) 
    StartAge = yearPlanted + StartAge;
  StandAge = (StartAge + StartMonth / 12) - (yearPlanted + StartMonth / 12);
  //Get and check StartAge
  StartAge = int(StandAge);
  if (StartAge < 0)
      std::cout << "Invalid Age Limits: StartAge must be greater than 0" << std::endl;
    //fprintf(logfp, "Invalid Age Limits: StartAge must be greater than 0");
  else if (StartAge > EndAge) 
      std::cout << "Invalid Age Limits: StartAge is greater than EndAge" << std::endl;
    //fprintf(logfp, "Invalid Age Limits: StartAge is greater than EndAge");

}
//-----------------------------------------------------------------------------

bool AssignMonthlyMetData(int calMonth, int calYear, int cellIndex,
                          double &SolarRad, double &FrostDays, double &Rain,
                          double &NetRad, double &Tav, double &Tx, double &Tn,
                          double &VPD, double &NDVI_AVH)
{
  bool hitNODATA;

  hitNODATA = false;

  if ( !getSeriesVal( SolarRad,  SS_SOLARRAD,  calMonth, calYear, cellIndex ) ) hitNODATA = true; 
  if ( !getSeriesVal( FrostDays, SS_FROSTDAYS, calMonth, calYear, cellIndex ) ) hitNODATA = true; 
  if ( !getSeriesVal( Rain,      SS_RAIN,      calMonth, calYear, cellIndex ) ) hitNODATA = true;  
  if ( userNetRadSeries() ) {
    if ( !getSeriesVal( NetRad,    SS_NETRAD,    calMonth, calYear, cellIndex ) ) 
      hitNODATA = true; 
  }
  
  if(userTavgSeries())
  {
    if ( !getSeriesVal( Tav,      SS_TAVG,       calMonth, calYear, cellIndex ) ) 
      hitNODATA = true; 
  }
  else
  {
    if ( !getSeriesVal( Tx,        SS_TMAX,      calMonth, calYear, cellIndex ) ) hitNODATA = true; 
    if ( !getSeriesVal( Tn,        SS_TMIN,      calMonth, calYear, cellIndex ) ) hitNODATA = true;
    Tav = (Tx + Tn) / 2;
  }
  
  if (userVpdSeries()) {
    if ( !getSeriesVal( VPD,       SS_VPD,       calMonth, calYear, cellIndex ) ) 
      hitNODATA = true; 
  }
  else
    VPD = getVPD(Tx, Tn);
  if (modelMode3PGS) {
    if ( !getSeriesVal( NDVI_AVH, SS_NDVI_AVH, calMonth, calYear, cellIndex ) ) 
      hitNODATA = true; 
  }
  return hitNODATA;
}

//-----------------------------------------------------------------------------

// This is the main routine for the 3PG model

//void runTreeModel(int minCy, int maxCy, bool spatial, long cellIndex)
void runTreeModel( MYDate minMY, MYDate maxMY, bool spatial, long cellIndex )
{
//  int minCy, maxCy; 
  
  // ANL - Note that in a spatial run, almost any input parameter having a
  // NODATA value will result in NODATA outputs.  To support this,
  // once we hit NODATA, we use a sequence of goto's to ensure that we
  // still go through all of the output hoops, but don't really do any
  // calculations.  The exception to this is that grids that are part of 
  // management tables may have NODATA cells, these cells cause that 
  // row in the management table to be ignored, for just that cell. 

  // The following variables probabaly could be Public so they can be
  // printed as part of the monthly output ...

  double RelAge;
  int dayofyr;
  double MoistRatio, NetRad; // PhysMod has been moved
  double wSmax;
  double delStemNo;

  double cumLAI;
  double oldVol;    //Added 16/07/02 as part of CVI
  
  bool hitNODATA = false, yrPreStart = false, yrPstEnd = false; 
  
  //New Soilwater modifier adjuster

  bool useMinASWTG = false;  
  double ASWmod;

  // monthly met data
  double SolarRad, Tx, Tn, Tav, VPD, FrostDays, dayLength, Rain; 

  // year and month counters, etc
  int year, calYear, calMonth, runYear, cm, cy; 

  int thinEventNo, defoltnEventNo;


  bool haveAvgTempSeries = false;  //Change: - Needs to be available iff 
								   //VPD is available

  // 3PGS - variables for 3PGS
  double FPAR_AVH;
  double NDVI_AVH;
  bool haveVpdSeries = false; 
  bool haveNetRadSeries = false; 

  // month counters
  //int firstRunMonth, lastRunMonth;
  int runMonth; 


  // Compute daylengths
  // ANL - only do dayLength here, as Tav and VPD potentially need recalculation each year. 
  dayofyr = -15;
  for (int mn = 1; mn <= 12; mn++) {
    dayofyr = dayofyr + DaysInMonth[mn];
    mDayLength[mn] = 86400 * getDayLength(Lat, dayofyr);
  }
  
  if (haveMinASWTG())
    useMinASWTG = true;
  else
    useMinASWTG = false;

  // ANL - Load the parameter values.  On NODATA write NODATA to output 
  // grids. 
  hitNODATA = !loadParamVals(cellIndex); 

  // For this cell, what is the first run month and the last run month, with reference 
  // to the January in minMY.  Use this later to spot periods that we skip calculations. 
  // calMonth and calYear relate to these values like this:
  //   ( calYear - minCY.year ) * 12 + calMonth - 1; 
  //firstRunMonth = ( ( yearPlanted ) - minMY.year ) * 12 + StartMonth - 1; // (yearPlanted + 1)
  //lastRunMonth  =  ( yearPlanted + EndAge - minMY.year ) * 12 + StartMonth - 2; 

  // May have hit nodata in StartMonth, yearPlanted and EndAge, in which case 
  // firstRunMonth and LastRunMonth will be meaningless.  In any case, if we aren't at 
  // NODATA, we have to do all the pre-year stuff below.  
  if ( hitNODATA )
    goto skipPreYearCalcs;
  
  Initialisation();

  // VPD and NetRad from internal model, or user specified series? 
  haveVpdSeries = userVpdSeries(); 
  haveNetRadSeries = userNetRadSeries(); 
  haveAvgTempSeries = userTavgSeries();

  // Assign initial state of the stand
  
  //StandAge = StartAge; //(StartAge-1) + 1.0/12.0; //18-01-02
  //StandAge = 0; //(minMY.year - yearPlanted); // + (cm - StartMonth) / 12.0;

  //New StandAge function
  GetStandAge();
  StemNo = StemNoi;
  //StartMonth++; //Synchronise with vb version 20-01-02

  //Fix for aracruz work.  Implements SeedlingMass distribution
  //that Peter Sands uses for multisite data.

  if ( haveSeedlingMass() )
  {
    WFi = (0.5 * StemNoi * SeedlingMass)/pow(10,6);
    WRi = (0.25 * StemNoi * SeedlingMass)/pow(10,6);
    WSi = WRi;
  }

  WS = WSi;
  WF = WFi;
  WR = WRi;
  
  ASW = ASWi;
  TotalLitter = 0;
  thinEventNo = 1;
  defoltnEventNo = 1;

  AvStemMass = WS * 1000 / StemNo;                             //  kg/tree
  avDBH = pow((AvStemMass / StemConst), (1 / StemPower));
  BasArea = ((pow((avDBH / 200), 2)) * Pi) * StemNo;
  SLA = SLA1 + (SLA0 - SLA1) * exp(-ln2 * pow((StandAge / tSLA), 2)); //Modified StandAge
  LAI = WF * SLA * 0.1;
  cLAI = LAI;
  
  fracBB = fracBB1 + (fracBB0 - fracBB1) * exp(-ln2 * (StandAge  / tBB)); //Modified StandAge
  Density = rhoMax + (rhoMin - rhoMax) * exp(-ln2 * (StandAge / tRho));

  StandVol = WS * (1 - fracBB) / Density;
  oldVol = StandVol;

  if (StandAge > 0) 
      MAI = StandVol / StandAge ;    //UnModified StandAge
  else MAI = 0;
  
  avDBHi = avDBH;
  LAIi = LAI;
  CumStemLoss = 0;
  
skipPreYearCalcs:

  // Do annual calculations.  The year loop here is controlled by minMY and maxMY, 
  // which refer to the overall run start and end, across all cells.
  

  //Print first month results

  calYear = minMY.year;
  calMonth = (int)StartMonth;
  
  //Find out if there is supposed to be any data here in the first place...

  hitNODATA = AssignMonthlyMetData(calMonth, calYear, cellIndex, 
                                       SolarRad, FrostDays, Rain, NetRad, Tav, Tx, Tn, VPD, NDVI_AVH) || hitNODATA;
      
  if (FrostDays > 30)
    FrostDays = 30;

  FR = FRp;

  if (spatial) {
      // 3PGS. Monthly output of some grids.  Note that yrPstEnd is not in this check, to ensure
      //previous calculated values are written instead of nodata
      
      writeMonthlyOutputGrids( calYear, calMonth, hitNODATA || yrPreStart, minMY, maxMY, cellIndex ); 

      // Monthly sample point output
      if (samplePointsMonthly)
        writeSampleFiles(cellIndex, calMonth, calYear);
  }


  //Start processing loop

  for ( cy = minMY.year; cy <= maxMY.year; cy++) {
    runYear = cy; 
    calYear = cy; 
    calMonth = (int)StartMonth; 
    runMonth = ( calYear - minMY.year ) * 12 + calMonth - 1; 

    //if ( calYear == 2008 )
    //  runYear = runYear; 

    // If we've already encountered NODATA we don't care about any annual variable 
    // except runYear. 
    if (hitNODATA)
      goto skipYearStartCalcs; 
    /*
    if ( runMonth < firstRunMonth )
      yrPreStart = true; 
    else 
      yrPreStart = false; 
    if ( runMonth > lastRunMonth )
      yrPstEnd = true; 
    else 
      yrPstEnd = false; 
    */
    year = cy - (int)yearPlanted;   // seem to still need year for point mode output. 
        
    // Once we've encountered nodata just cycle through as quickly as possible.  
    if ( hitNODATA || yrPreStart || yrPstEnd )
      goto skipYearStartCalcs; 

    // Initialise cumulative variables
    cLitter = 0;
    CumdelWF = 0;
    CumdelWR = 0;
    CumdelWS = 0;
    CumAPARU = 0;
    cumARAD = 0;
    cumLAI = 0;
    cumGPP = 0;
    cumWabv = 0;            //Now known as cumWabvgrnd
    cumTransp = 0;
    cumEvapTransp = 0;
    cumIrrig = 0;

    //Initialise annual cumulative variables
    aStemDM = 0;
    aRADint = 0;
    aGPP = 0;
    aNPP = 0;
    aEvapTransp = 0;
    aTransp = 0;
    aSupIrrig = 0;

    
    // Get management-related options for current year and cell. 
    // First load param file values, then possibly override them with management table values. 
    if (!haveAgeDepFert())
      FR = FRp; 
    if (nFertility > 0) 
      FR = lookupManageTable( runYear, MT_FERTILITY, FRp, cellIndex ); 

    MinASW = MinASWp; 
    if (nMinAvailSW > 0) 
      MinASW = lookupManageTable( runYear, MT_MINASW, MinASWp, cellIndex ); 

    Irrig = 0; 
    if (nIrrigation > 0) {
      Irrig = lookupManageTable( runYear, MT_IRRIGATION, 0, cellIndex ); 
    }
    else Irrig = 0;
    

skipYearStartCalcs:


    //Fill in noData values for first year. AS 20/01/02

    if (spatial) {

      if (calYear == minMY.year)

        for(int beforeCalcMonth = 1; beforeCalcMonth < StartMonth; beforeCalcMonth++)
          writeMonthlyOutputGrids( calYear, beforeCalcMonth, true, minMY, maxMY, cellIndex );
    }

    //Initialise output step cumulative variables

    delStemNo = 0;
    cRADint = 0;
    cLAI = 0;
    cCVI = 0;
    cGPP = 0;
    cNPP = 0;
    cStemDM = 0;
    cTransp = 0;
    cEvapTransp = 0;
    cRainInt = 0;
    cSupIrrig = 0;
    cLitter = 0;

    // Do monthly calculations

    for ( cm = (int)StartMonth + 1; cm < (int)StartMonth + 13; cm++ ) {
      //Note that the added one is to sync in with the VB code, which always
      //incrememt to the next month before starting...
       if ( cm >= 13 ) {
        calYear = cy + 1;
        if (calYear > maxMY.year) 
          goto skipYearEndCalcs;
        calMonth = cm - 12; 
      }
      else {
        calYear = cy;  
        calMonth = cm; 
      }

      //Check to see the year we are currently in is before the plant year
      if(cm == StartMonth)
      {
        //Initialise output step cumulative variables
        delStemNo = 0;
        cRADint = 0;
        cLAI = 0;
        cCVI = 0;
        cGPP = 0;
        cNPP = 0;
        cStemDM = 0;
        cTransp = 0;
        cEvapTransp = 0;
        cRainInt = 0;
        cSupIrrig = 0;
        cLitter = 0;
      }

      yrPreStart = false;
      yrPstEnd = false;
      if (calYear < yearPlanted)
        yrPreStart = true;
      if ((calYear == yearPlanted) && (calMonth < StartMonth))
        yrPreStart = true;
      if (calYear > (yearPlanted + EndAge)) 
        yrPstEnd = true;
      if ((calYear == (yearPlanted + EndAge)) && (calMonth > StartMonth))
        yrPstEnd = true;

 
      if ( hitNODATA || yrPreStart || yrPstEnd )
        goto skipMonthCalcs; 
      
      hitNODATA = AssignMonthlyMetData(calMonth, calYear, cellIndex, 
                                       SolarRad, FrostDays, Rain, NetRad, Tav, Tx, Tn, VPD, NDVI_AVH) || hitNODATA;
      if ( hitNODATA )
        goto skipMonthCalcs; 
      
      dayLength = mDayLength[calMonth]; 
      
      // Determine the various environmental modifiers

      //Fertility.
      if (nFertility > 0)
      {
        //Do nothing
      } else {
        //If we are in a period where we wish FR to decay, make it so.
        if ( haveAgeDepFert() && (FRstart <= StandAge) && (FRend > StandAge))
        {
          FR = FR - FR*FRdec;
        }
      }

      // calculate temperature response function to apply to alpha
      if ((Tav <= growthTmin) || (Tav >= growthTmax)) 
        fT = 0;
      else 
        fT = ((Tav - growthTmin) / (growthTopt - growthTmin)) * 
          pow(((growthTmax - Tav) / (growthTmax - growthTopt)), 
              ((growthTmax - growthTopt) / (growthTopt - growthTmin)));
      
      // calculate VPD modifier
      fVPD = exp(-CoeffCond * VPD);
      
      // calculate soil water modifier
      if (useMinASWTG)
      { 
        double dAdjMod;
        if (MaxASW <= MinASWTG)
        {
          dAdjMod = (MaxASW - MinASWTG)/MinASWTG;
          ASWmod = pow(2.718281828459045235, dAdjMod);
        }
         else
          ASWmod = 1;          
      }
      else
        ASWmod = 1;

      MoistRatio = ASWmod * ASW / MaxASW;
      fSW = 1 / (1 + pow(((1 - MoistRatio) / SWconst), SWpower));
      fSW = 0.5;

      if (fSW == 1)
        bool test = true;

      if (fNn == 0)
        fNutr = 1;
      else
         fNutr = 1 - (1 - fN0) * pow((1 - FR), fNn);
      
      // calculate frost modifier
      fFrost = 1 - kF * (FrostDays / 30.0);
      
      // calculate age modifier
      RelAge = StandAge  / MaxAge;  //Modified StandAge
      if (modelMode3PGS)
        fAge = 1;
      else
        fAge = (1 / (1 + pow((RelAge / rAge), nAge)));
         
      // PhysMod is the physiological modifier to be applied to canopy conductance
      // and APARu. It is the lesser of the soil-water and VPD modifier, times the
      // age modifier:

      PhysMod = Minimum(fVPD, fSW) * fAge;
      
      // Determine gross and net biomass production
      
      // canopy cover and light interception.
      CanCover = 1;
      if ((fullCanAge > 0) && (StandAge  < fullCanAge))  //Modified StandAge
        CanCover = StandAge / fullCanAge; //Modified StandAge
      lightIntcptn = (1 - (exp(-k * LAI)));
      
      
      // 3PGS. 
      // Calculate FPAR_AVH and LAI from NDVI data. 
      if (modelMode3PGS) {
        // Initial value of FPAR_AVH from linear fit. 
        FPAR_AVH = (NDVI_AVH * NDVI_FPAR_constant) + NDVI_FPAR_intercept;
        // Constrain FPAR_AVH to within threshhold values. 
        if (FPAR_AVH > 0.98) 
          FPAR_AVH = 0.98;
        else if (FPAR_AVH < 0) 
          FPAR_AVH = 0;
        // LAI
        if (FPAR_AVH < 0.05)
          LAI = 0.2;
        else
          LAI = -2.0 * log(1 - FPAR_AVH);
      }

      // Calculate PAR, APAR, APARu and GPP
      //   APARu, "utilisable PAR", is intercepted PAR multiplied by the physiological
      //     modifier applied to conductance (PhysMod).
      //   alphaC is alpha, the nominal "canopy qunatum efficiency", multiplied by
      //     modifiers to take into account effects of nutrition, temperature and
      //     frost on photosynthetic rate
      RAD = SolarRad * DaysInMonth[calMonth];        // MJ/m^2
      PAR = RAD * molPAR_MJ;                      // mol/m^2
      // 3PGS
      if (modelMode3PGS)
        APAR = PAR * FPAR_AVH;
      else 
        APAR = PAR * lightIntcptn * CanCover;
      APARu = APAR * PhysMod;

      alphaC = alpha * fNutr * fT * fFrost * PhysMod;   //22-07-02 for Excel March beta consis.
      epsilon = gDM_mol * molPAR_MJ * alphaC;
      RADint = RAD * lightIntcptn * CanCover;
      GPPmolc = APARu * alphaC;                   // mol/m^2
      GPPdm = epsilon * RADint / 100;               // tDM/ha
      NPP = GPPdm * y;                            // assumes respiratory rate is constant
     
     // Determine biomass increments and losses
     
      // calculate partitioning coefficients
      m = m0 + (1 - m0) * FR;
      pFS = pfsConst * pow(avDBH, pfsPower);
      if (fabs(APAR) < 0.000001) APAR = 0.000001;
      pR = pRx * pRn / (pRn + (pRx - pRn) * (APARu / APAR) * m);
      pS = (1 - pR) / (1 + pFS);
      pF = 1 - pR - pS;

      // calculate biomass increments
      delWF = NPP * pF;
      delWR = NPP * pR;
      delWS = NPP * pS;

      // calculate litterfall & root turnover -
      gammaF = gammaFx * gammaF0 / 
        (gammaF0 + (gammaFx - gammaF0) * 
         exp(-12 * log(1 + gammaFx / gammaF0) * StandAge / tgammaF)); 
      
      delLitter = gammaF * WF;
      delRloss = Rttover * WR;
      
      // Calculate end-of-month biomass
      
      if (!modelMode3PGS) {
        WF = WF + delWF - delLitter;
        WR = WR + delWR - delRloss;
        WS = WS + delWS;
        TotalW = WF + WR + WS;
        TotalLitter = TotalLitter + delLitter;
      }
      
      
      // Now do the water balance ...
      
      // calculate canopy conductance from stomatal conductance

      CanCond = MaxCond * PhysMod * Minimum(1.0, LAI / LAIgcx);
      //if (fabs(0 - CanCond) < eps)
      if (CanCond == 0)
        CanCond = 0.0001;

      //transpiration from Penman-Monteith (mm/day converted to mm/month)
      Transp = CanopyTranspiration(SolarRad, VPD, dayLength, BLcond, 
                                   CanCond, NetRad, haveNetRadSeries);
      Transp = DaysInMonth[calMonth] * Transp;
      
      // do soil water balance

      if (LAImaxIntcptn <= 0)
        Interception = MaxIntcptn;
      else 
        Interception = MaxIntcptn * Minimum(1, LAI / LAImaxIntcptn);
      EvapTransp = Transp + Interception * Rain;
      ASW = ASW + Rain + (100 * Irrig / 12) - EvapTransp;        //Irrig is Ml/ha/year
      monthlyIrrig = 0;
      if (ASW < MinASW) {
        if (MinASW > 0) {               // make up deficit with irrigation
          monthlyIrrig = MinASW - ASW;
          cumIrrig = cumIrrig + monthlyIrrig;
        }
        ASW = MinASW;
      }
      else if (ASW > MaxASW) { 
        ASW = MaxASW;
      }
      
      WUE = 100 * NPP / EvapTransp;

      //StandAge = (cy - yearPlanted) + (cm - StartMonth + 1) / 12.0; //OG position
      StandAge = StandAge + 1.0/12.0;
     
      if (StandAge < 0)
      {
          std::cout << "Negative StandAge" << std::endl;
        //fprintf(logfp, "Negative StandAge");
      }
      

      if (!modelMode3PGS) {
        
        // Update tree and stand data
        
        //Calculate mortality

        wSmax = wSx1000 * pow((1000/StemNo), thinPower); 
        AvStemMass = WS * 1000 / StemNo;
        delStems = 0;
        if (wSmax < AvStemMass)
        {
          delStems = getMortality(StemNo, WS);
          WF = WF - mF * delStems * (WF / StemNo);
          WR = WR - mR * delStems * (WR / StemNo);
          WS = WS - mS * delStems * (WS / StemNo);
          StemNo = StemNo - delStems;
          wSmax = wSx1000 * pow((1000 / StemNo), thinPower);
          AvStemMass = WS * 1000 / StemNo;
          delStemNo = delStemNo + delStems;
        }

        //update age-dependent factors
        SLA = SLA1 + (SLA0 - SLA1) * exp(-ln2 * pow((StandAge  / tSLA), 2));  //Modified StandAge
        fracBB = fracBB1 + (fracBB0 - fracBB1) * exp(-ln2 * (StandAge  / tBB));  //Modified StandAge
        Density = rhoMax + (rhoMin - rhoMax) * exp(-ln2 * (StandAge / tRho));

        //update stsand characteristics
        LAI = WF * SLA * 0.1;
        avDBH = pow((AvStemMass / StemConst), (1 / StemPower));
        BasArea = (pow((avDBH / 200), 2) * Pi) * StemNo;
        StandVol = WS * (1 - fracBB) / Density;
        
        CVI = StandVol - oldVol;       //Added 16/07/02 
        oldVol = StandVol;

        if (StandAge  > 0)             //Modified StandAge
          MAI =  StandVol / StandAge ;  //UnModified StandAge
        else 
          MAI = 0;
        
        // Update accumulated totals
        
        cRADint = cRADint + RADint;
        aRADint = aRADint + RADint;
        cGPP = cGPP + GPPdm;
        aGPP = aGPP + GPPdm;
        cNPP = cNPP + NPP;
        aNPP = aNPP + NPP;
        cCVI = cCVI + CVI;
        cLitter = cLitter + delLitter;
        cStemDM = cStemDM + delWS;
        aStemDM = aStemDM + delWS;
        cRainInt = cRainInt + RainIntcptn;
        cTransp = cTransp + Transp;
        aTransp = aTransp + Transp;
        cEvapTransp = cEvapTransp + EvapTransp;
        aEvapTransp = aEvapTransp + EvapTransp;
        aSupIrrig = aSupIrrig + monthlyIrrig;
        cSupIrrig = cSupIrrig + monthlyIrrig;
        cWUE = 100 * cNPP / cEvapTransp;
        cLAI = cLAI + LAI / 12.0;

        // Accumulate biomass increments and LAI
        //cumTransp = cumTransp + Transp;
        //cumEvapTransp = cumEvapTransp + EvapTransp ;  //unknown what or why
        //CumdelWF = CumdelWF + delWF;
        //CumdelWR = CumdelWR + delWR;
        //CumdelWS = CumdelWS + delWS;
        //cLitter = cLitter + delLitter;
        cumWabv = cumWabv + delWF + delWS - delLitter;  // ANL - PROBLEM?  
        //cumGPP = cumGPP + GPPdm;
        //cumLAI = cumLAI + LAI;
        
        // Accumulate intercepted radiation (MJ/m2) and production (t/ha)
        cumARAD = cumARAD + RAD * lightIntcptn * CanCover;
        CumAPARU = CumAPARU + APARu;
      }
      
      // 3PGS
      if (modelMode3PGS) {
        delWAG = NPP * (1 - pR);
        cumWabv += delWAG;
      } 

skipMonthCalcs:


      if (spatial) {
        // 3PGS. Monthly output of some grids.  Note that yrPstEnd is not in this check, to ensure
        //previous calculated values are written instead of nodata
        writeMonthlyOutputGrids( calYear, calMonth, hitNODATA || yrPreStart, minMY, maxMY, cellIndex ); 

        // Monthly sample point output
        if (samplePointsMonthly)
          writeSampleFiles(cellIndex, calMonth, calYear);
      }
      // if (showDetailedResults) writeMonthlySummary(lastMonthFile, monthCounter, year);

    } 
    
    if ( hitNODATA || yrPreStart || yrPstEnd )
      goto skipYearEndCalcs;

    // Calculate above ground and total Epsilon
    if (!modelMode3PGS) {
      if (aRADint == 0) {
//        sprintf(outstr, 
//          "Warning: No growth occurred in year with Standage = %4.0f "
//          "at cell index %d\n", StandAge, cellIndex);
//        fprintf(logfp, outstr);
//        fprintf(stderr, outstr);
      }
      else {
        aEpsilonStem = 100 * aStemDM / aRADint;    //100 converts to gDM/MJ
        aEpsilonGross = 100 * aGPP / aRADint;
      }
    }
    
    // Update some stand characteristics
    LAI = cumLAI / 12.0;
    fracBB = fracBB1 + (fracBB0 - fracBB1) * exp(-ln2 * (StandAge  / tBB));  //Modified StandAge
    StandVol = WS * (1 - fracBB) / Density;
    if (StandAge  > 0)              //Modified StandAge
      MAI = StandVol / StandAge ;   //Modified StandAge
    else 
      MAI = 0;

    // Determine peak LAI & MAI and age at peaks
    if (LAI > LAIx) {
      LAIx = LAI;
      ageLAIx = StandAge ;  //Modified StandAge
    }
    if (MAI > MAIx) {
      MAIx = MAI;
      ageMAIx = StandAge ;  //Modified StandAge
    }
    
skipYearEndCalcs:

    // ANL if (showDetailedResults) writeAnnualResults(year);
    // ANL if (showStandSummary) writeStandSummary(year);if(calMonth == 1)
      
    if (!spatial) {
      writeStandSummary(year);
    }
    else {
      // ANL - Annual sample point output. 
      if (samplePointsYearly)
        writeSampleFiles(cellIndex, 12, calYear);
    }

    // Restore LAI
    LAI = WF * SLA * 0.1; 
    
  } 

  // ANL - Save wanted variables into grids. 
  // if (spatial)
  //   saveVariableVals(cellIndex, hitNODATA);
}

