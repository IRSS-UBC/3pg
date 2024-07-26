/*
All source code remains the property and copyright of CSIRO.

CSIRO accepts no responsibility for the use of 3PG(S) or of the model 3-PG in
the form supplied or as subsequently modified by third parties. CSIRO disclaims
liability for all losses, damages and costs incurred by any person as a result
of relying on this software.
Use of this software assumes agreement to this condition of use
*/

#include <iostream>
#include <math.h>
#include <stdio.h>
#include "Data_io.hpp"
#include "The_3PG_Model.hpp"
//____________________________________
//
// The code for 3PG - July, 2024
//____________________________________

//NOTE: The following conversion factors are used:
//
//    1 MJ  = 2.3 mol PAR
//    1 mol = 24 gDM

#define Pi 3.141592654
#define ln2 0.693147181
#define eps 0.0001

// Controls and counters
int DaysInMonth[13] = {                  // array for days in months 
  0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};
bool showDetailedResults;                // TRUE ==> show monthly results
bool showStandSummary;                   // TRUE ==> show stand summary
bool modelMode3PGS = false;                         

// Time variant management factors
int nFertility;                          // size of site fertility array
int nMinAvailSW;                         // size of MinAvailSW array
int nIrrigation;                         // size of irrigation array
double Irrig;                            // current annual irrigation (ML/y)

extern bool samplePointsMonthly;
extern bool samplePointsYearly;

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

double getMortality(double oldN, double oldW, InputParams& params)
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
    x1 = 1000 * params.mS * oldW / oldN;
    i = 0;
    while (true)
    {
        i = i + 1;
        x2 = params.wSx1000 * pow(n, (1 - params.thinPower));
        fN = x2 - x1 * n - (1 - params.mS) * oldW;
        dfN = (1 - params.thinPower) * x2 / n - x1;
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
    double gBL, double gC, double NR, bool hNRS, InputParams& params)
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
        netRad = params.Qa + params.Qb * (Q * pow(10, 6)) / h;                // Q in MJ/m2/day --> W/m2

    defTerm = rhoAir * lambda * (VPDconv * VPD) * gBL;
    div = gC * (1 + e20) + gBL;
    Etransp = gC * (e20 * netRad + defTerm) / div;           // in J/m2/s
    CT = Etransp / lambda * h;         // converted to kg/m2/day

    return CT;
}

//-----------------------------------------------------------------------------

//Standage function translated from March beta of Excel 3-PG
//StartAge and StandAge are considered global variables...
void GetStandAge(double& StandAge, InputParams & params)
{
    //Assign initial stand age
    if (params.StartAge < params.yearPlanted)
        params.StartAge = params.yearPlanted + params.StartAge;
    StandAge = (params.StartAge + params.StartMonth / 12) - (params.yearPlanted + params.StartMonth / 12);
    //Get and check StartAge
    params.StartAge = int(StandAge);
    if (params.StartAge < 0)
        std::cout << "Invalid StartAge: StartAge must be greater than 0" << std::endl;
    //fprintf(logfp, "Invalid Age Limits: StartAge must be greater than 0");
    else if (params.StartAge > params.EndYear)
        std::cout << "Invalid Age Limits: StartAge is greater than EndYear" << std::endl;
        //fprintf(logfp, "Invalid Age Limits: StartAge is greater than EndYear");

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
void runTreeModel(std::unordered_map<std::string, PPPG_OP_VAR> opVars, MYDate spMinMY, MYDate spMaxMY, long cellIndex, DataInput *dataInput)
{
    Vars vars;

    // ANL - Load the parameter values.  On NODATA return, because all pixels are initialized to nodata
    InputParams params;
    if (!dataInput->getInputParams(cellIndex, params)) {
        return;
    }

    //various variables that used to be global but I have no idea how necessary they are to the modol
    double aEpsilonGross;
    double aEpsilonStem;
    double aEvapTransp;
    double aGPP;
    double aNPP;
    double aRADint;
    double aStemDM;
    double aSupIrrig;
    double aTransp;
    double avDBHi;
    double cRADint;
    double cRainInt;
    double cStemDM;
    double cSupIrrig;
    double CumAPARU;
    double cumARAD;
    double CumdelWF;
    double CumdelWR;
    double CumdelWS;
    double cumEvapTransp;
    double cumGPP;
    double cumIrrig;
    double CumStemLoss;
    double cumTransp;
    double GPPmolc;
    double LAIi;

    double TranspScaleFactor;
    double RunOff;
    double monthlyIrrig;
    double poolFractn = 0;
    double excessSW;
    double pooledSW = 0;
    double RainIntcptn = 0;

    //before start or after end indication
    bool yrPreStart = false;
    bool yrPstEnd = false;

    // At initialisation param file has only possible value to use.  
    double MinASW = params.MinASWp;

    //soil parameters for soil class
    double SWconst;
    double SWpower;

    //day length (by month)
    double mDayLength[13];

    //stand age
    double StandAge;

    // Stand factors that are specifically age dependent
    double SLA;
    double gammaF;
    double CanCover;

    // Derive some parameters
    double pfsPower = log(params.pFS20 / params.pFS2) / log(10);
    double pfsConst = params.pFS2 / pow(2, pfsPower);

    double Interception;         // Proportion of rainfall intercepted by canopy (used to be assigned 0.15 in assignDefaultParameters)
    double Density;              // Basic density (t/m3) (used to be assigned 0.5 in assignDefaultParameters)

    // Intermediate monthly results
    double m;
    double epsilon;
    double RAD;
    double PAR;
    double RADint;
    double lightIntcptn;
    double CanCond;
    double AvStemMass;
    double GPPdm;
    double pR;
    double pS;
    double pF;
    double pFS;
    double delWF;
    double delWR;
    double delWS;
    double delStems;
    double delLitter;
    double delRloss;

    double RelAge;
    int dayofyr;
    double MoistRatio; // PhysMod has been moved
    double wSmax;
    double delStemNo;

    double cumLAI;
    double oldVol;    //Added 16/07/02 as part of CVI

    //New Soilwater modifier adjuster
    bool useMinASWTG = false;
    double ASWmod;

    // monthly met data
    double dayLength;

    // year and month counters, etc
    int year, calYear, calMonth, runYear, cm, cy;

    int thinEventNo, defoltnEventNo;

    bool haveAvgTempSeries = false;  //Change: - Needs to be available iff 
    //VPD is available

    // 3PGS - variables for 3PGS
    double FPAR_AVH;

    //std::cout << "\n\nCELL INDEX " << cellIndex << "!!!!!!!!!!!!!!!!!!!!!\n\n" << std::endl;
    // Compute daylengths
    // ANL - only do dayLength here, as Tav and VPD potentially need recalculation each year. 
    dayofyr = -15;
    for (int mn = 1; mn <= 12; mn++) {
        dayofyr = dayofyr + DaysInMonth[mn];
        mDayLength[mn] = 86400 * getDayLength(params.Lat, dayofyr);
    }

    useMinASWTG = dataInput->haveMinASWTG;

    // Assign the SWconst and SWpower parameters for this soil class
    if (params.soilIndex != 0) {
        SWconst = 0.8 - 0.1 * params.soilIndex;
        SWpower = 11 - 2 * params.soilIndex;
    }
    else {
        SWconst = params.SWconst0;
        SWpower = params.SWpower0;
    }

    // Initial ASW must be between min and max ASW
    if (params.ASWi <= MinASW) {
        params.ASWi = MinASW;
    }
    else if (params.ASWi >= params.MaxASW) {
        params.ASWi = params.MaxASW;
    }

    // Initialise ages
    vars.MAIx = 0;
    vars.LAIx = 0;

    //assign initial age of stand
    GetStandAge(StandAge, params);
    vars.StemNo = params.StemNoi;
    //StartMonth++; //Synchronise with vb version 20-01-02

    if (dataInput->haveSeedlingMass)
    {
        params.WFi = (0.5 * params.StemNoi * params.SeedlingMass) / pow(10, 6);
        params.WRi = (0.25 * params.StemNoi * params.SeedlingMass) / pow(10, 6);
        params.WSi = params.WRi;
    }

    vars.WS = params.WSi;
    vars.WF = params.WFi;
    vars.WR = params.WRi;

    vars.ASW = params.ASWi;
    vars.TotalLitter = 0;
    thinEventNo = 1;
    defoltnEventNo = 1;

    AvStemMass = vars.WS * 1000 / vars.StemNo;                             //  kg/tree
    vars.avDBH = pow((AvStemMass / params.StemConst), (1 / params.StemPower));
    vars.BasArea = ((pow((vars.avDBH / 200), 2)) * Pi) * vars.StemNo;
    SLA = params.SLA1 + (params.SLA0 - params.SLA1) * exp(-ln2 * pow((StandAge / params.tSLA), 2)); //Modified StandAge
    vars.LAI = vars.WF * SLA * 0.1;
    vars.cLAI = vars.LAI;

    vars.fracBB = params.fracBB1 + (params.fracBB0 - params.fracBB1) * exp(-ln2 * (StandAge / params.tBB)); //Modified StandAge
    Density = params.rhoMax + (params.rhoMin - params.rhoMax) * exp(-ln2 * (StandAge / params.tRho));
    gammaF = params.gammaFx * params.gammaF0 /
        (params.gammaF0 + (params.gammaFx - params.gammaF0) *
            exp(-12 * log(1 + params.gammaFx / params.gammaF0) * StandAge / params.tgammaF));

    vars.StandVol = vars.WS * (1 - vars.fracBB) / Density;
    oldVol = vars.StandVol;

    if (StandAge > 0)
        vars.MAI = vars.StandVol / StandAge;    //UnModified StandAge
    else vars.MAI = 0;

    avDBHi = vars.avDBH;
    LAIi = vars.LAI;
    CumStemLoss = 0;
    

    // Do annual calculations.  The year loop here is controlled by minMY and maxMY, 
    // which refer to the overall run start and end, across all cells.

    //Print first month results
    calYear = spMinMY.year;
    calMonth = (int)params.StartMonth;

    //Find out if there is supposed to be any data here in the first place...

    SeriesParams sParams;
    if (!dataInput->getSeriesParams(cellIndex, calYear, calMonth, sParams)) {
        return;
    }

    if (sParams.FrostDays > 30) {
        sParams.FrostDays = 30;
    }

    // init pool fraction
    poolFractn = std::max(0.0, std::min(1.0, poolFractn));

    vars.FR = params.FRp;

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
        calMonth = (int)params.StartMonth;

        year = cy - (int)params.yearPlanted;   // seem to still need year for point mode output. 

        // Initialise cumulative variables
        vars.cLitter = 0;
        CumdelWF = 0;
        CumdelWR = 0;
        CumdelWS = 0;
        CumAPARU = 0;
        cumARAD = 0;
        cumLAI = 0;
        cumGPP = 0;
        vars.cumWabv = 0;            //Now known as cumWabvgrnd
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
        if (!dataInput->haveAgeDepFert)
            vars.FR = params.FRp;
        if (nFertility > 0)
            vars.FR = lookupManageTable(runYear, MT_FERTILITY, params.FRp, cellIndex);

        MinASW = params.MinASWp;
        if (nMinAvailSW > 0)
            MinASW = lookupManageTable(runYear, MT_MINASW, params.MinASWp, cellIndex);

        Irrig = 0;
        if (nIrrigation > 0) {
            Irrig = lookupManageTable(runYear, MT_IRRIGATION, 0, cellIndex);
        }
        else Irrig = 0;

        //Initialise output step cumulative variables
        delStemNo = 0;
        cRADint = 0;
        cRainInt = 0;
        cStemDM = 0;
        cSupIrrig = 0;
        vars.cLAI = 0;
        vars.cCVI = 0;
        vars.cNPP = 0;
        vars.cGPP = 0;
        vars.cTransp = 0;
        vars.cEvapTransp = 0;
        vars.cLitter = 0;

        // Do monthly calculations
        for (cm = (int)params.StartMonth + 1; cm < (int)params.StartMonth + 13; cm++) {
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
            if (cm == params.StartMonth) {
                //Initialise output step cumulative variables
                delStemNo = 0;
                cRADint = 0;
                cStemDM = 0;
                cRainInt = 0;
                cSupIrrig = 0;
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
            if (calYear < params.yearPlanted)
                yrPreStart = true;
            if ((calYear == params.yearPlanted) && (calMonth < params.StartMonth))
                yrPreStart = true;
            if (calYear > (params.EndYear))
                yrPstEnd = true;
            if ((calYear == (params.EndYear)) && (calMonth > params.StartMonth))
                yrPstEnd = true;


            if (yrPreStart || yrPstEnd)
                goto skipMonthCalcs;

            if (!dataInput->getSeriesParams(cellIndex, calYear, calMonth, sParams)) {
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
                if (dataInput->haveAgeDepFert && (params.FRstart <= StandAge) && (params.FRend > StandAge)) {
                    vars.FR = vars.FR - vars.FR * params.FRdec;
                }
            }

            // calculate temperature response function to apply to alpha
            if ((sParams.Tavg <= params.growthTmin) || (sParams.Tavg >= params.growthTmax))
                vars.fT = 0;
            else
                vars.fT = ((sParams.Tavg - params.growthTmin) / (params.growthTopt - params.growthTmin)) *
                pow(((params.growthTmax - sParams.Tavg) / (params.growthTmax - params.growthTopt)),
                    ((params.growthTmax - params.growthTopt) / (params.growthTopt - params.growthTmin)));

            // calculate VPD modifier
            vars.fVPD = exp(-params.CoeffCond * sParams.VPD);

            // calculate soil water modifier
            if (useMinASWTG)
            {
                double dAdjMod;
                if (params.MaxASW <= params.MinASWTG)
                {
                    dAdjMod = (params.MaxASW - params.MinASWTG) / params.MinASWTG;
                    ASWmod = pow(2.718281828459045235, dAdjMod);
                }
                else
                    ASWmod = 1;
            }
            else
                ASWmod = 1;

            MoistRatio = ASWmod * vars.ASW / params.MaxASW;
            vars.fSW = 1 / (1 + pow(((1 - MoistRatio) / SWconst), SWpower));

            if (vars.fSW == 1)
                bool test = true;
            // calculate soil nutrition
            if (params.fNn == 0)
                vars.fNutr = 1;
            else
                vars.fNutr = 1 - (1 - params.fN0) * pow((1 - vars.FR), params.fNn);

            // calculate frost modifier
            vars.fFrost = 1 - params.kF * (sParams.FrostDays / 30.0);

            // calculate age modifier
            RelAge = StandAge / params.MaxAge;  //Modified StandAge
            if (modelMode3PGS)
                vars.fAge = 1;
            else
                vars.fAge = (1 / (1 + pow((RelAge / params.rAge), params.nAge)));

            // PhysMod is the physiological modifier to be applied to canopy conductance
            // and APARu. It is the lesser of the soil-water and VPD modifier, times the
            // age modifier:

            vars.PhysMod = std::min(vars.fVPD, vars.fSW) * vars.fAge;

            // Determine gross and net biomass production

            // canopy cover and light interception.
            CanCover = 1;
            if ((params.fullCanAge > 0) && (StandAge < params.fullCanAge))  //Modified StandAge
                CanCover = (StandAge) / params.fullCanAge; //Modified StandAge
            lightIntcptn = (1 - (exp(-params.k * vars.LAI / CanCover)));

            // 3PGS. 
            // Calculate FPAR_AVH and LAI from NDVI data. 
            if (modelMode3PGS) {
                // Initial value of FPAR_AVH from linear fit. 
                FPAR_AVH = (sParams.NDVI_AVH * params.NDVI_FPAR_constant) + params.NDVI_FPAR_intercept;

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
            RAD = sParams.SolarRad * DaysInMonth[calMonth];        // MJ/m^2
            PAR = RAD * params.molPAR_MJ;                      // mol/m^2
            // 3PGS
            if (modelMode3PGS)
                vars.APAR = PAR * FPAR_AVH;
            else
                vars.APAR = PAR * lightIntcptn * CanCover;
            vars.APARu = vars.APAR * vars.PhysMod;

            vars.alphaC = params.alpha * vars.fNutr * vars.fT * vars.fFrost * vars.PhysMod;   //22-07-02 for Excel March beta consis.
            epsilon = params.gDM_mol * params.molPAR_MJ * vars.alphaC;
            RADint = RAD * lightIntcptn * CanCover;
            GPPmolc = vars.APARu * vars.alphaC;                   // mol/m^2
            GPPdm = epsilon * RADint / 100;               // tDM/ha
            vars.NPP = GPPdm * params.y;                            // assumes respiratory rate is constant

            // calculate canopy conductance
            double gC;
            if (vars.LAI <= params.LAIgcx) {
                gC = params.MaxCond * vars.LAI / params.LAIgcx;
            }
            else {
                gC = params.MaxCond;
            }
            CanCond = gC * vars.PhysMod;

            // calculate transpiration from Penman-Monteith (mm/day converted to mm/month)
            vars.Transp = CanopyTranspiration(sParams.SolarRad, sParams.VPD, dayLength, params.BLcond,
                CanCond, sParams.NetRad, dataInput->haveNetRadParam(), params);
            vars.Transp = DaysInMonth[calMonth] * vars.Transp;

            // rainfall interception
            if (params.LAImaxIntcptn <= 0)
                Interception = params.MaxIntcptn;
            else
                Interception = params.MaxIntcptn * std::min(1.0, vars.LAI / params.LAImaxIntcptn);
            RainIntcptn = sParams.Rain * Interception;
            
            // water balance
            monthlyIrrig = 0;
            RunOff = 0;
            vars.ASW = vars.ASW + sParams.Rain + (100 * Irrig / 12) + pooledSW;        //Irrig is Ml/ha/year
            vars.EvapTransp = std::min(vars.ASW, vars.Transp + RainIntcptn);
            excessSW = std::max(vars.ASW - vars.EvapTransp - params.MaxASW, 0.0);
            vars.ASW = vars.ASW - vars.EvapTransp - excessSW;
            pooledSW = poolFractn * excessSW;
            RunOff = (1 - poolFractn) * excessSW;
            if (vars.ASW < params.MinASWp) {
                monthlyIrrig = params.MinASWp - vars.ASW;
                cumIrrig = cumIrrig + monthlyIrrig;
            }

            // correct for actual ET
            TranspScaleFactor = vars.EvapTransp / (vars.Transp + RainIntcptn);      
            GPPdm = TranspScaleFactor * GPPdm;
            vars.NPP = TranspScaleFactor * vars.NPP;
            if (vars.EvapTransp != 0) {
                vars.WUE = 100 * vars.NPP / vars.EvapTransp;
            }
                   
            // calculate partitioning coefficients
            m = params.m0 + (1 - params.m0) * vars.FR;
            pFS = pfsConst * pow(vars.avDBH, pfsPower);
            if (fabs(vars.APAR) < 0.000001) vars.APAR = 0.000001;
            pR = params.pRx * params.pRn / (params.pRn + (params.pRx - params.pRn) * (vars.APARu / vars.APAR) * m);
            pS = (1 - pR) / (1 + pFS);
            pF = 1 - pR - pS;

            // calculate biomass increments
            delWF = vars.NPP * pF;
            delWR = vars.NPP * pR;
            delWS = vars.NPP * pS;

            // calculate litterfall & root turnover -
            delLitter = gammaF * vars.WF;
            delRloss = params.Rttover * vars.WR;

            // Calculate end-of-month biomass
            if (!modelMode3PGS) {
                vars.WF = vars.WF + delWF - delLitter;
                vars.WR = vars.WR + delWR - delRloss;
                vars.WS = vars.WS + delWS;
                vars.TotalW = vars.WF + vars.WR + vars.WS;
                vars.TotalLitter = vars.TotalLitter + delLitter;
            }

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

                wSmax = params.wSx1000 * pow((1000 / vars.StemNo), params.thinPower);
                AvStemMass = vars.WS * 1000 / vars.StemNo;
                delStems = 0;
                if (wSmax < AvStemMass)
                {
                    delStems = getMortality(vars.StemNo, vars.WS, params);
                    vars.WF = vars.WF - params.mF * delStems * (vars.WF / vars.StemNo);
                    vars.WR = vars.WR - params.mR * delStems * (vars.WR / vars.StemNo);
                    vars.WS = vars.WS - params.mS * delStems * (vars.WS / vars.StemNo);
                    vars.StemNo = vars.StemNo - delStems;
                    wSmax = params.wSx1000 * pow((1000 / vars.StemNo), params.thinPower);
                    AvStemMass = vars.WS * 1000 /  vars.StemNo;
                    delStemNo = delStemNo + delStems;
                }

                //update age-dependent factors
                SLA = params.SLA1 + (params.SLA0 - params.SLA1) * exp(-ln2 * pow((StandAge / params.tSLA), 2));  //Modified StandAge
                vars.fracBB = params.fracBB1 + (params.fracBB0 - params.fracBB1) * exp(-ln2 * (StandAge / params.tBB));  //Modified StandAge
                Density = params.rhoMax + (params.rhoMin - params.rhoMax) * exp(-ln2 * (StandAge / params.tRho));
                gammaF = params.gammaFx * params.gammaF0 /
                    (params.gammaF0 + (params.gammaFx - params.gammaF0) *
                        exp(-12 * log(1 + params.gammaFx / params.gammaF0) * StandAge / params.tgammaF));

                //update stsand characteristics
                vars.LAI = vars.WF * SLA * 0.1;
                vars.avDBH = pow((AvStemMass / params.StemConst), (1 / params.StemPower));
                vars.BasArea = (pow((vars.avDBH / 200), 2) * Pi) * vars.StemNo;
                vars.StandVol = vars.WS * (1 - vars.fracBB) / Density;

                vars.CVI =  vars.StandVol - oldVol;       //Added 16/07/02 
                oldVol = vars.StandVol;

                if (StandAge > 0)             //Modified StandAge
                    vars.MAI = vars.StandVol / StandAge;  //UnModified StandAge
                else
                    vars.MAI = 0;

                // Update accumulated totals

                cRADint = cRADint + RADint;
                aRADint = aRADint + RADint;
                vars.cGPP = vars.cGPP + GPPdm;
                aGPP = aGPP + GPPdm;
                vars.cNPP = vars.cNPP + vars.NPP;
                aNPP = aNPP + vars.NPP;
                vars.cCVI = vars.cCVI + vars.CVI;
                vars.cLitter = vars.cLitter + delLitter;
                cStemDM = cStemDM + delWS;
                aStemDM = aStemDM + delWS;
                cRainInt = cRainInt + RainIntcptn;
                vars.cTransp = vars.cTransp + vars.Transp;
                aTransp = aTransp + vars.Transp;
                vars.cEvapTransp = vars.cEvapTransp + vars.EvapTransp;
                aEvapTransp = aEvapTransp + vars.EvapTransp;
                aSupIrrig = aSupIrrig + monthlyIrrig;
                cSupIrrig = cSupIrrig + monthlyIrrig;
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

                // Accumulate intercepted radiation (MJ/m2) and production (t/ha)
                cumARAD = cumARAD + RAD * lightIntcptn * CanCover;
                CumAPARU = CumAPARU + vars.APARu;
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
            // if (showDetailedResults) writeMonthlySummary(lastMonthFile, monthCounter, year)
        }

        if (yrPreStart || yrPstEnd)
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
        vars.LAI = cumLAI / 12.0;
        vars.fracBB = params.fracBB1 + (params.fracBB0 - params.fracBB1) * exp(-ln2 * (StandAge / params.tBB));  //Modified StandAge
        vars.StandVol = vars.WS * (1 - vars.fracBB) / Density;
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