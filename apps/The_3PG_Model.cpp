//static char rcsid[] = "$Id: The_3PG_Model.cpp,v 1.6 2001/08/02 06:51:42 lou026 Exp $";

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
#include <iostream>
#include <math.h>
#include <stdio.h>
#include "Data_io.hpp"
#include "The_3PG_Model.hpp"
//____________________________________
//
// The code for 3PG - March 24th, 2000
//____________________________________

//NOTE: The following conversion factors are used:
//
//    1 MJ  = 2.3 mol PAR
//    1 mol = 24 gDM

#define Pi 3.141592654
#define ln2 0.693147181
#define eps 0.0001

// Controls and counters
double StartAge, EndYear;                 // age of trees at start/end of run
double StartMonth;                       // month of year to start run
double yearPlanted;                      // year trees planted
int DaysInMonth[13] = {                  // array for days in months 
  0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};
bool modelMode3PGS = false;

bool yrPreStart = false;
bool yrPstEnd = false;

// Site characteristics, site specific parameters
double Lat = 1000;                       // site latitude
double MaxASW, MinASW, MinASWp;          // maximum & minimum available soil water, current, param file. 
double FRp;                          // site fertility rating, current, param file. 
double FRstart, FRend, FRdec;            // Start, end and decrement % for fertility decrease with time
double soilIndex;                        // soil class index
double SWconst, SWpower;                 // soil parameters for soil class

// Time variant management factors
int nFertility;                          // size of site fertility array
int nMinAvailSW;                         // size of MinAvailSW array
int nIrrigation;                         // size of irrigation array
double Irrig;                            // current annual irrigation (ML/y)

// day length
double mDayLength[13];                   

double SeedlingMass;                     // Alternative way of deriving initial distribution 
                                         // of mass using seedling mass constant
double StandAge;                         // stand age
double ASWi;                        // available soil water
double MinASWTG;
double StemNoi;                  // stem numbers
double WFi;                          // foliage biomass
double WRi;                          // root biomass
double WSi;                          // stem biomass
double LAIi;                        // canopy leaf area index
double MAIi;                        // mean annual volume increment                                                                        

// Stand factors that are specifically age dependent
double SLA;
double gammaF;
double CanCover;

// Parameter values
// int MaxAge;
// ANL changed MaxAge from int to double
double MaxAge;
double gammaFx, gammaF0, tgammaF;
double Rttover;
double SLA0, SLA1, tSLA;
double fullCanAge;
double k;
double pFS2, pFS20;
double StemConst, StemPower;
double SWconst0, SWpower0;
double Interception;
double BLcond;
double MaxCond, CoeffCond;
double y;
double growthTmax, growthTmin, growthTopt;
double wSx1000;
double thinPower;
double mF, mR, mS;
double m0, fN0, fNn;                      //added 22-07-02
double alpha;
double pRx, pRn;
double nAge, rAge;
double kF;
double fracBB0, fracBB1, tBB;
double Density;
double rhoMin, rhoMax, tRho;             // Standage varying density 3-06-02 
double pfsConst, pfsPower;               // derived from pFS2, pFS20

//Conversion factors
double Qa, Qb;
double gDM_mol;
double molPAR_MJ;

//Additional factors (conductance)
double LAIgcx;
double MaxIntcptn;
double LAImaxIntcptn;

// Intermediate monthly results
double m, epsilon;
double RAD, PAR, RADint;
double lightIntcptn;
double CanCond;
double RainIntcptn; //Added 16/07/02

double AvStemMass;
double GPPdm;
double pR, pS, pF, pFS;
double delWF, delWR, delWS, delStems;
double delLitter, delRloss;

// ANL - other globals
double NDVI_FPAR_intercept, NDVI_FPAR_constant;

extern bool samplePointsMonthly;
extern bool samplePointsYearly;

// ANL - globals defined in 3pg.cpp
//extern FILE* logfp;
extern char outstr[];

//-----------------------------------------------------------------------------

double getDayLength(double Lat, int dayOfYear)
{
    // gets fraction of day when sun is "up"
    double sLat, cLat, sinDec, cosH0;

    sLat = sin(Pi * Lat / 180);
    cLat = cos(Pi * Lat / 180);

    sinDec = 0.4 * sin(0.0172 * (dayOfYear - 80));
    cosH0 = -sinDec * sLat / (cLat * sqrt(1 - pow(sinDec, 2)));
    if (cosH0 > 1)
        return 0;
    else if (cosH0 < -1)
        return 1;
    else {
        return acos(cosH0) / Pi;
    }
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
        netRad = Qa + Qb * (Q * pow(10, 6)) / h;                // Q in MJ/m2/day --> W/m2

    defTerm = rhoAir * lambda * (VPDconv * VPD) * gBL;
    div = (1 + e20 + gBL / gC);
    Etransp = (e20 * netRad + defTerm) / div;           // in J/m2/s
    CT = Etransp / lambda * h;         // converted to kg/m2/day

    return CT;
}

//-----------------------------------------------------------------------------

void Initialisation()
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
        std::cout << "Invalid StartAge: StartAge must be greater than 0" << std::endl;
    //fprintf(logfp, "Invalid Age Limits: StartAge must be greater than 0");
    else if (StartAge > EndYear)
        std::cout << "Invalid Age Limits: StartAge is greater than EndYear" << std::endl;
        //fprintf(logfp, "Invalid Age Limits: StartAge is greater than EndYear");

}
//-----------------------------------------------------------------------------

bool AssignMonthlyMetData(int calMonth, int calYear, long cellIndex,
    double& SolarRad, double& FrostDays, double& Rain,
    double& NetRad, double& Tav, double& Tx, double& Tn,
    double& VPD, double& NDVI_AVH)
{
    bool hitNODATA;

    hitNODATA = false;

    if (!getSeriesVal(SolarRad, SS_SOLARRAD, calMonth, calYear, cellIndex)) hitNODATA = true;
    if (!getSeriesVal(FrostDays, SS_FROSTDAYS, calMonth, calYear, cellIndex)) hitNODATA = true;
    if (!getSeriesVal(Rain, SS_RAIN, calMonth, calYear, cellIndex)) hitNODATA = true;
    if (userNetRadSeries()) {
        if (!getSeriesVal(NetRad, SS_NETRAD, calMonth, calYear, cellIndex))
            hitNODATA = true;
    }

    if (userTavgSeries())
    {
        if (!getSeriesVal(Tav, SS_TAVG, calMonth, calYear, cellIndex))
            hitNODATA = true;
    }
    else
    {
        if (!getSeriesVal(Tx, SS_TMAX, calMonth, calYear, cellIndex)) hitNODATA = true;
        if (!getSeriesVal(Tn, SS_TMIN, calMonth, calYear, cellIndex)) hitNODATA = true;
        Tav = (Tx + Tn) / 2;
    }

    if (userVpdSeries()) {
        if (!getSeriesVal(VPD, SS_VPD, calMonth, calYear, cellIndex))
            hitNODATA = true;
    }
    else
        VPD = getVPD(Tx, Tn);
    if (modelMode3PGS) {
        if (!getSeriesVal(NDVI_AVH, SS_NDVI_AVH, calMonth, calYear, cellIndex))
            hitNODATA = true;
    }
    return hitNODATA;
}

//-----------------------------------------------------------------------------

struct Vars {
    double WF;
    double LAIx = 0;
    double StemNo;
    double WS;
    double BasArea;
    double WR;
    double MAIx = 0;
    double ASW;
    double fracBB;
    double TotalLitter;
    double WUE;
    double cNPP;
    double avDBH;
    double LAI;
    double alphaC;
    double cLAI;
    double cTransp;
    double StandVol;
    double cEvapTransp;
    double MAI;
    double FR;
    double cLitter;
    double cumWabv;
    double cCVI;
    double cGPP;
    double fT;
    double fVPD;
    double fSW;
    double fNutr;
    double fFrost;
    double fAge;
    double PhysMod;
    double APAR;
    double APARu;
    double NPP;
    double TotalW;
    double Transp;
    double EvapTransp;
    double ageMAIx;
    double ageLAIx;
    double CVI;
    double cWUE;
    double delWAG;
};

void copyVars(Vars vars, std::unordered_map<std::string, PPPG_OP_VAR>& opVars) {
    opVars["WF"].v = vars.WF;
    opVars["LAIx"].v = vars.LAIx;
    opVars["StemNo"].v = vars.StemNo;
    opVars["WS"].v = vars.WS;
    opVars["BasArea"].v = vars.BasArea;
    opVars["WR"].v = vars.WR;
    opVars["MAIx"].v = vars.MAIx;
    opVars["ASW"].v = vars.ASW;
    opVars["fracBB"].v = vars.fracBB;
    opVars["TotalLitter"].v = vars.TotalLitter;
    opVars["WUE"].v = vars.WUE;
    opVars["cNPP"].v = vars.cNPP;
    opVars["avDBH"].v = vars.avDBH;
    opVars["LAI"].v = vars.LAI;
    opVars["alphaC"].v = vars.alphaC;
    opVars["cLAI"].v = vars.cLAI;
    opVars["cTransp"].v = vars.cTransp;
    opVars["StandVol"].v = vars.StandVol;
    opVars["cEvapTransp"].v = vars.cEvapTransp;
    opVars["MAI"].v = vars.MAI;
    opVars["FR"].v = vars.FR;
    opVars["cLitter"].v = vars.cLitter;
    opVars["cumWabv"].v = vars.cumWabv;
    opVars["cCVI"].v = vars.cCVI;
    opVars["cGPP"].v = vars.cGPP;
    opVars["fT"].v = vars.fT;
    opVars["fVPD"].v = vars.fVPD;
    opVars["fSW"].v = vars.fSW;
    opVars["fNutr"].v = vars.fNutr;
    opVars["fFrost"].v = vars.fFrost;
    opVars["fAge"].v = vars.fAge;
    opVars["PhysMod"].v = vars.PhysMod;
    opVars["APAR"].v = vars.APAR;
    opVars["APARu"].v = vars.APARu;
    opVars["NPP"].v = vars.NPP;
    opVars["TotalW"].v = vars.TotalW;
    opVars["Transp"].v = vars.Transp;
    opVars["EvapTransp"].v = vars.EvapTransp;
    opVars["ageMAIx"].v = vars.ageMAIx;
    opVars["ageLAIx"].v = vars.ageLAIx;
    opVars["CVI"].v = vars.CVI;
    opVars["cWUE"].v = vars.cWUE;
    opVars["delWAG"].v = vars.delWAG;
}

// This is the main routine for the 3PG model
void runTreeModel(std::unordered_map<std::string, PPPG_OP_VAR> opVars, MYDate spMinMY, MYDate spMaxMY, long cellIndex)
{
    Vars vars;
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

    //std::cout << "\n\nCELL INDEX " << cellIndex << "!!!!!!!!!!!!!!!!!!!!!\n\n" << std::endl;
    // Compute daylengths
    // ANL - only do dayLength here, as Tav and VPD potentially need recalculation each year. 
    dayofyr = -15;
    for (int mn = 1; mn <= 12; mn++) {
        dayofyr = dayofyr + DaysInMonth[mn];
        mDayLength[mn] = 86400 * getDayLength(Lat, dayofyr);
    }

    useMinASWTG = haveMinASWTG();

    if (!loadParamVals(cellIndex)) {
        return;
    }

    Initialisation();

    // VPD and NetRad from internal model, or user specified series? 
    haveVpdSeries = userVpdSeries();
    haveNetRadSeries = userNetRadSeries();
    haveAvgTempSeries = userTavgSeries();

    // Assign initial state of the stand
    GetStandAge();
    vars.StemNo = StemNoi;
    //StartMonth++; //Synchronise with vb version 20-01-02

    //Fix for aracruz work.  Implements SeedlingMass distribution
    //that Peter Sands uses for multisite data.
    if (haveSeedlingMass())
    {
        WFi = (0.5 * StemNoi * SeedlingMass) / pow(10, 6);
        WRi = (0.25 * StemNoi * SeedlingMass) / pow(10, 6);
        WSi = WRi;
    }

    vars.WS = WSi;
    vars.WF = WFi;
    vars.WR = WRi;

    vars.ASW = ASWi;
    vars.TotalLitter = 0;
    thinEventNo = 1;
    defoltnEventNo = 1;

    AvStemMass = vars.WS * 1000 / vars.StemNo;                             //  kg/tree
    vars.avDBH = pow((AvStemMass / StemConst), (1 / StemPower));
    vars.BasArea = ((pow((vars.avDBH / 200), 2)) * Pi) * vars.StemNo;
    SLA = SLA1 + (SLA0 - SLA1) * exp(-ln2 * pow((StandAge / tSLA), 2)); //Modified StandAge
    vars.LAI = vars.WF * SLA * 0.1;
    vars.cLAI = vars.LAI;

    vars.fracBB = fracBB1 + (fracBB0 - fracBB1) * exp(-ln2 * (StandAge / tBB)); //Modified StandAge
    Density = rhoMax + (rhoMin - rhoMax) * exp(-ln2 * (StandAge / tRho));

    vars.StandVol = vars.WS * (1 - vars.fracBB) / Density;
    oldVol = vars.StandVol;

    if (StandAge > 0)
        vars.MAI = vars.StandVol / StandAge;    //UnModified StandAge
    else vars.MAI = 0;

    LAIi = vars.LAI;

    // Do annual calculations.  The year loop here is controlled by minMY and maxMY, 
    // which refer to the overall run start and end, across all cells.

    //Print first month results
    calYear = spMinMY.year;
    calMonth = (int)StartMonth;

    //Find out if there is supposed to be any data here in the first place...
    if (AssignMonthlyMetData(calMonth, calYear, cellIndex, SolarRad, FrostDays, Rain, NetRad, Tav, Tx, Tn, VPD, NDVI_AVH)) {
        return;
    }

    if (FrostDays > 30)
        FrostDays = 30;

    vars.FR = FRp;

    //write initial state of output variables
    copyVars(vars, opVars);
    writeMonthlyOutputGrids(opVars, calYear, calMonth, spMinMY, spMaxMY, cellIndex);

    // Monthly sample point output
    if (samplePointsMonthly)
        writeSampleFiles(opVars, cellIndex, calMonth, calYear);

    //Start processing loop
    for (cy = spMinMY.year; cy <= spMaxMY.year; cy++) {
        runYear = cy;
        calYear = cy;
        calMonth = (int)StartMonth;

        year = cy - (int)yearPlanted;   // seem to still need year for point mode output. 

        // Initialise cumulative variables
        vars.cLitter = 0;
        cumLAI = 0;
        vars.cumWabv = 0;            //Now known as cumWabvgrnd

        // Get management-related options for current year and cell. 
        // First load param file values, then possibly override them with management table values. 
        if (!haveAgeDepFert())
            vars.FR = FRp;

        if (nFertility > 0)
            vars.FR = lookupManageTable(runYear, MT_FERTILITY, FRp, cellIndex);

        MinASW = MinASWp;
        if (nMinAvailSW > 0)
            MinASW = lookupManageTable(runYear, MT_MINASW, MinASWp, cellIndex);

        Irrig = 0;
        if (nIrrigation > 0) {
            Irrig = lookupManageTable(runYear, MT_IRRIGATION, 0, cellIndex);
        }
        else Irrig = 0;

        //Initialise output step cumulative variables
        delStemNo = 0;
        vars.cLAI = 0;
        vars.cCVI = 0;
        vars.cNPP = 0;
        vars.cGPP = 0;
        vars.cTransp = 0;
        vars.cEvapTransp = 0;
        vars.cLitter = 0;

        // Do monthly calculations
        for (cm = (int)StartMonth + 1; cm < (int)StartMonth + 13; cm++) {
            //Note that the added one is to sync in with the VB code, which always
            //incrememt to the next month before starting...
            if (cm >= 13) {
                calYear = cy + 1;
                if (calYear > spMaxMY.year)
                    goto skipYearEndCalcs;
                calMonth = cm - 12;
            }
            else {
                calYear = cy;
                calMonth = cm;
            }

            //Check to see the year we are currently in is before the plant year
            if (cm == StartMonth)
            {
                //Initialise output step cumulative variables
                delStemNo = 0;
                vars.cLAI = 0;
                vars.cCVI = 0;
                vars.cGPP = 0;
                vars.cNPP = 0;
                vars.cTransp = 0;
                vars.cEvapTransp = 0;
                vars.cLitter = 0;
            }

            yrPreStart = false;
            yrPstEnd = false;
            if (calYear < yearPlanted)
                yrPreStart = true;
            if ((calYear == yearPlanted) && (calMonth < StartMonth))
                yrPreStart = true;
            if (calYear > (EndYear))
                yrPstEnd = true;
            if ((calYear == (EndYear)) && (calMonth > StartMonth))
                yrPstEnd = true;


            if (yrPreStart || yrPstEnd)
                goto skipMonthCalcs;

            if (AssignMonthlyMetData(calMonth, calYear, cellIndex, SolarRad, FrostDays, Rain, NetRad, Tav, Tx, Tn, VPD, NDVI_AVH)) {
                return;
            }

            dayLength = mDayLength[calMonth];

            // Determine the various environmental modifiers

            //Fertility.
            if (nFertility > 0)
            {
                //Do nothing
            }
            else {
                //If we are in a period where we wish FR to decay, make it so.
                if (haveAgeDepFert() && (FRstart <= StandAge) && (FRend > StandAge))
                {
                    vars.FR = vars.FR - vars.FR * FRdec;
                }
            }

            // calculate temperature response function to apply to alpha
            if ((Tav <= growthTmin) || (Tav >= growthTmax))
                vars.fT = 0;
            else
                vars.fT = ((Tav - growthTmin) / (growthTopt - growthTmin)) *
                pow(((growthTmax - Tav) / (growthTmax - growthTopt)),
                    ((growthTmax - growthTopt) / (growthTopt - growthTmin)));

            // calculate VPD modifier
            vars.fVPD = exp(-CoeffCond * VPD);

            // calculate soil water modifier
            if (useMinASWTG)
            {
                double dAdjMod;
                if (MaxASW <= MinASWTG)
                {
                    dAdjMod = (MaxASW - MinASWTG) / MinASWTG;
                    ASWmod = pow(2.718281828459045235, dAdjMod);
                }
                else
                    ASWmod = 1;
            }
            else
                ASWmod = 1;

            MoistRatio = ASWmod * vars.ASW / MaxASW;
            vars.fSW = 1 / (1 + pow(((1 - MoistRatio) / SWconst), SWpower));

            if (vars.fSW == 1)
                bool test = true;

            if (fNn == 0)
                vars.fNutr = 1;
            else
                vars.fNutr = 1 - (1 - fN0) * pow((1 - vars.FR), fNn);

            // calculate frost modifier
            vars.fFrost = 1 - kF * (FrostDays / 30.0);

            // calculate age modifier
            RelAge = StandAge / MaxAge;  //Modified StandAge
            if (modelMode3PGS)
                vars.fAge = 1;
            else
                vars.fAge = (1 / (1 + pow((RelAge / rAge), nAge)));

            // PhysMod is the physiological modifier to be applied to canopy conductance
            // and APARu. It is the lesser of the soil-water and VPD modifier, times the
            // age modifier:

            vars.PhysMod = std::min(vars.fVPD, vars.fSW) * vars.fAge;

            // Determine gross and net biomass production

            // canopy cover and light interception.
            CanCover = 1;
            if ((fullCanAge > 0) && (StandAge < fullCanAge))  //Modified StandAge
                CanCover = (StandAge) / fullCanAge; //Modified StandAge
            lightIntcptn = (1 - (exp(-k * vars.LAI)));


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
                    vars.LAI = 0.2;
                else
                    vars.LAI = -2.0 * log(1 - FPAR_AVH);
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
                vars.APAR = PAR * FPAR_AVH;
            else
                vars.APAR = PAR * lightIntcptn * CanCover;
            vars.APARu = vars.APAR * vars.PhysMod;


            vars.alphaC = alpha * vars.fNutr * vars.fT * vars.fFrost * vars.PhysMod;   //22-07-02 for Excel March beta consis.
            epsilon = gDM_mol * molPAR_MJ * vars.alphaC;
            RADint = RAD * lightIntcptn * CanCover;
            GPPdm = epsilon * RADint / 100;               // tDM/ha
            vars.NPP = GPPdm * y;                            // assumes respiratory rate is constant

            // Determine biomass increments and losses

             // calculate partitioning coefficients
            m = m0 + (1 - m0) * vars.FR;
            pFS = pfsConst * pow(vars.avDBH, pfsPower);
            if (fabs(vars.APAR) < 0.000001) vars.APAR = 0.000001;
            pR = pRx * pRn / (pRn + (pRx - pRn) * (vars.APARu / vars.APAR) * m);
            pS = (1 - pR) / (1 + pFS);
            pF = 1 - pR - pS;

            // calculate biomass increments
            delWF = vars.NPP * pF;
            delWR = vars.NPP * pR;
            delWS = vars.NPP * pS;

            // calculate litterfall & root turnover -
            // print out each gamma variable value before computing gammaF
            //std::cout << "StandAge = " << StandAge << std::endl;
            gammaF = gammaFx * gammaF0 /
                (gammaF0 + (gammaFx - gammaF0) *
                    exp(-12 * log(1 + gammaFx / gammaF0) * StandAge / tgammaF));
            //std::cout << "gammaF = " << gammaF << std::endl;
            delLitter = gammaF * vars.WF;
            delRloss = Rttover * vars.WR;

            // Calculate end-of-month biomass

            if (!modelMode3PGS) {
                vars.WF = vars.WF + delWF - delLitter;
                vars.WR = vars.WR + delWR - delRloss;
                vars.WS = vars.WS + delWS;
                vars.TotalW = vars.WF + vars.WR + vars.WS;
                vars.TotalLitter = vars.TotalLitter + delLitter;
            }

            // Now do the water balance ...
            // calculate canopy conductance from stomatal conductance

            CanCond = MaxCond * vars.PhysMod * std::min(1.0, vars.LAI / LAIgcx);
            //if (fabs(0 - CanCond) < eps)
            if (CanCond == 0)
                CanCond = 0.0001;

            //transpiration from Penman-Monteith (mm/day converted to mm/month)
            vars.Transp = CanopyTranspiration(SolarRad, VPD, dayLength, BLcond,
                CanCond, NetRad, haveNetRadSeries);
            vars.Transp = DaysInMonth[calMonth] * vars.Transp;

            // do soil water balance

            if (LAImaxIntcptn <= 0)
                Interception = MaxIntcptn;
            else
                Interception = MaxIntcptn * std::min((double)1, vars.LAI / LAImaxIntcptn);
            vars.EvapTransp = vars.Transp + Interception * Rain;
            vars.ASW = vars.ASW + Rain + (100 * Irrig / 12) - vars.EvapTransp;        //Irrig is Ml/ha/year
            if (vars.ASW < MinASW) {
                vars.ASW = MinASW;
            }
            else if (vars.ASW > MaxASW) {
                vars.ASW = MaxASW;
            }

            vars.WUE = 100 * vars.NPP / vars.EvapTransp;

            //StandAge = (cy - yearPlanted) + (cm - StartMonth + 1) / 12.0; //OG position
            StandAge = StandAge + 1.0 / 12.0;

            if (StandAge < 0)
            {
                std::cout << "Negative StandAge" << std::endl;
                //fprintf(logfp, "Negative StandAge");
            }

            if (!modelMode3PGS) {

                // Update tree and stand data

                //Calculate mortality

                wSmax = wSx1000 * pow((1000 / vars.StemNo), thinPower);
                AvStemMass = vars.WS * 1000 / vars.StemNo;
                delStems = 0;
                if (wSmax < AvStemMass)
                {
                    delStems = getMortality(vars.StemNo, vars.WS);
                    vars.WF = vars.WF - mF * delStems * (vars.WF / vars.StemNo);
                    vars.WR = vars.WR - mR * delStems * (vars.WR / vars.StemNo);
                    vars.WS = vars.WS - mS * delStems * (vars.WS / vars.StemNo);
                    vars.StemNo = vars.StemNo - delStems;
                    wSmax = wSx1000 * pow((1000 / vars.StemNo), thinPower);
                    AvStemMass = vars.WS * 1000 /  vars.StemNo;
                    delStemNo = delStemNo + delStems;
                }

                //update age-dependent factors
                SLA = SLA1 + (SLA0 - SLA1) * exp(-ln2 * pow((StandAge / tSLA), 2));  //Modified StandAge
                vars.fracBB = fracBB1 + (fracBB0 - fracBB1) * exp(-ln2 * (StandAge / tBB));  //Modified StandAge
                Density = rhoMax + (rhoMin - rhoMax) * exp(-ln2 * (StandAge / tRho));

                //update stsand characteristics
                vars.LAI = vars.WF * SLA * 0.1;
                vars.avDBH = pow((AvStemMass / StemConst), (1 / StemPower));
                vars.BasArea = (pow((vars.avDBH / 200), 2) * Pi) *  vars.StemNo;
                vars.StandVol =  vars.WS * (1 - vars.fracBB) / Density;

                vars.CVI =  vars.StandVol - oldVol;       //Added 16/07/02 
                oldVol = vars.StandVol;

                if (StandAge > 0)             //Modified StandAge
                    vars.MAI = vars.StandVol / StandAge;  //UnModified StandAge
                else
                    vars.MAI = 0;

                // Update accumulated totals
                vars.cGPP = vars.cGPP + GPPdm;
                vars.cNPP = vars.cNPP + vars.NPP;
                vars.cCVI = vars.cCVI + vars.CVI;
                vars.cLitter = vars.cLitter + delLitter;
                vars.cTransp = vars.cTransp + vars.Transp;
                vars.cEvapTransp = vars.cEvapTransp + vars.EvapTransp;
                vars.cWUE = 100 * vars.cNPP / vars.cEvapTransp;
                vars.cLAI = vars.cLAI +  vars.LAI / 12.0;

                // Accumulate biomass increments and LAI
                //cumTransp = cumTransp + Transp;
                //cumEvapTransp = cumEvapTransp + EvapTransp ;  //unknown what or why
                //CumdelWF = CumdelWF + delWF;
                //CumdelWR = CumdelWR + delWR;
                //CumdelWS = CumdelWS + delWS;
                //cLitter = cLitter + delLitter;
                vars.cumWabv = vars.cumWabv + delWF + delWS - delLitter;  // ANL - PROBLEM?  
                //cumGPP = cumGPP + GPPdm;
                //cumLAI = cumLAI + LAI;
            }

            // 3PGS
            if (modelMode3PGS) {
                vars.delWAG = vars.NPP * (1 - pR);
                vars.cumWabv += vars.delWAG;
            }

        skipMonthCalcs:

            //Joe is not sure what the following comment means, but figures it might be important so he's leaving it.
            // 3PGS. Monthly output of some grids.  Note that yrPstEnd is not in this check, to ensure
            //previous calculated values are written instead of nodata

            //Joe: the startMonth, cy, cm, yrPreStart, yrPstEnd code has to be some of the most unecessarily complicated
            //code I've ever seen... Add on to that that even if we are PAST THE END, we continue iterating through a for
            //loop, just not changing any values (???). AND the month counter (cm) starts iterating at 2 and finishes at 14...
            //I'm working right now on changing how data output works to make future parallelization possible while refactoring 
            //as little as possible (to keep changes as more manageable chunks), but in the future, we MUST change how this works, 
            //it *should* be pretty simple to just use minMY and maxMY.
            if (!yrPreStart && !yrPstEnd && !(calYear == spMaxMY.year && calMonth == spMaxMY.mon)) {
                

                //the !(calYear == maxMY.year && calMonth == maxMY.mon) is so that at the last iteration we don't write to a monthly
                //output, but rather skip it and the values are eventually written via writeOutputGrids().
                copyVars(vars, opVars);
                writeMonthlyOutputGrids(opVars, calYear, calMonth, spMinMY, spMaxMY, cellIndex);
            }

            // Monthly sample point output
            if (samplePointsMonthly)
                writeSampleFiles(opVars, cellIndex, calMonth, calYear);
            // if (showDetailedResults) writeMonthlySummary(lastMonthFile, monthCounter, year);

        }

        if (yrPreStart || yrPstEnd)
            goto skipYearEndCalcs;

        // Update some stand characteristics
        vars.LAI = cumLAI / 12.0;
        vars.fracBB = fracBB1 + (fracBB0 - fracBB1) * exp(-ln2 * (StandAge / tBB));  //Modified StandAge
        vars.StandVol =  vars.WS * (1 - vars.fracBB) / Density;
        if (StandAge > 0)              //Modified StandAge
            vars.MAI = vars.StandVol / StandAge;   //Modified StandAge
        else
            vars.MAI = 0;

        // Determine peak LAI & MAI and age at peaks
        if (vars.LAI > vars.LAIx) {
            vars.LAIx = vars.LAI;
            vars.ageLAIx = StandAge;  //Modified StandAge
        }
        if (vars.MAI > vars.MAIx) {
            vars.MAIx = vars.MAI;
            vars.ageMAIx = StandAge;  //Modified StandAge
        }

    skipYearEndCalcs:

        // ANL if (showDetailedResults) writeAnnualResults(year);
        // ANL if (showStandSummary) writeStandSummary(year);if(calMonth == 1)

      
        // ANL - Annual sample point output. 
        if (samplePointsYearly) {
            writeSampleFiles(opVars, cellIndex, 12, calYear);
        }

        // Restore LAI
        vars.LAI = vars.WF * SLA * 0.1;

    }
    copyVars(vars, opVars);
    writeOutputGrids(opVars, cellIndex);
}
