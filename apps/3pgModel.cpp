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
#include "3pgModel.hpp"
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

//-----------------------------------------------------------------------------

double getDayLength(double Lat, int dayOfYear)
{
    // gets fraction of day when sun is "up"
    double sLat, cLat, sinDec, cosH0;

    sLat = sin(Pi * Lat / 180);
    cLat = cos(Pi * Lat / 180);

    sinDec = 0.4 * sin(0.0172 * (dayOfYear - 80));
    cosH0 = -sinDec * sLat / (cLat * sqrt(1 - pow(sinDec, 2)));
    
    if (cosH0 > 1) {
        return 0;
    }
    else if (cosH0 < -1) {
        return 1;
    }

    return acos(cosH0) / Pi;
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

    n = oldN / 1000;
    x1 = 1000 * params.mS * oldW / oldN;
    i = 0;

    while (true) {
        i = i + 1;
        x2 = params.wSx1000 * pow(n, (1 - params.thinPower));
        fN = x2 - x1 * n - (1 - params.mS) * oldW;
        dfN = (1 - params.thinPower) * x2 / n - x1;
        dN = -fN / dfN;
        n = n + dN;
        if ((fabs(dN) <= eps) || (i >= 5)) {
            break;
        }
    }

    return oldN - 1000 * n;
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

    if (hNRS) {
        netRad = NR;
    }
    else {
        netRad = params.Qa + params.Qb * (Q * pow(10, 6)) / h;                // Q in MJ/m2/day --> W/m2
    }

    defTerm = rhoAir * lambda * (VPDconv * VPD) * gBL;
    div = gC * (1 + e20) + gBL;
    Etransp = gC * (e20 * netRad + defTerm) / div;           // in J/m2/s
    CT = Etransp / lambda * h;         // converted to kg/m2/day

    return CT;
}

double getGammaFoliage(InputParams& params, double& StandAge) {
    return params.gammaFx * params.gammaF0 / 
        (params.gammaF0 + (params.gammaFx - params.gammaF0) * exp(-12 * log(1 + params.gammaFx / params.gammaF0) * StandAge / params.tgammaF));
}

//-----------------------------------------------------------------------------

struct Vars {
    double WF = 0;
    double LAIx = 0;
    double StemNo = 0;
    double WS = 0;
    double BasArea = 0;
    double WR = 0;
    double MAIx = 0;
    double ASW = 0;
    double fracBB = 0;
    double TotalLitter = 0;
    double WUE = 0;
    double cNPP = 0;
    double avDBH = 0;
    double LAI = 0;
    double alphaC = 0;
    double cLAI = 0;
    double cTransp = 0;
    double StandVol = 0;
    double cEvapTransp = 0;
    double MAI = 0;
    double FR = 0;
    double cLitter = 0;
    double cumWabv = 0;
    double cCVI = 0;
    double GPP = 0;
    double cGPP = 0;
    double fT = 0;
    double fVPD = 0;
    double fSW = 0;
    double fNutr = 0;
    double fFrost = 0;
    double fAge = 0;
    double PhysMod = 0;
    double APAR = 0;
    double APARu = 0;
    double NPP = 0;
    double TotalW = 0;
    double Transp = 0;
    double EvapTransp = 0;
    double ageMAIx = 0;
    double ageLAIx = 0;
    double CVI = 0;
    double cWUE = 0;
    double delWAG = 0;
};

void copyVars(Vars vars, std::unordered_map<std::string, double>& opVarVals) {
    opVarVals["wf"] = vars.WF;
    opVarVals["laix"] = vars.LAIx;
    opVarVals["stemno"] = vars.StemNo;
    opVarVals["ws"] = vars.WS;
    opVarVals["basarea"] = vars.BasArea;
    opVarVals["wr"] = vars.WR;
    opVarVals["maix"] = vars.MAIx;
    opVarVals["asw"] = vars.ASW;
    opVarVals["fracbb"] = vars.fracBB;
    opVarVals["totallitter"] = vars.TotalLitter;
    opVarVals["wue"] = vars.WUE;
    opVarVals["cnpp"] = vars.cNPP;
    opVarVals["avdbh"] = vars.avDBH;
    opVarVals["lai"] = vars.LAI;
    opVarVals["alphac"] = vars.alphaC;
    opVarVals["clai"] = vars.cLAI;
    opVarVals["ctransp"] = vars.cTransp;
    opVarVals["standvol"] = vars.StandVol;
    opVarVals["cevaptransp"] = vars.cEvapTransp;
    opVarVals["mai"] = vars.MAI;
    opVarVals["fr"] = vars.FR;
    opVarVals["clitter"] = vars.cLitter;
    opVarVals["cumwabv"] = vars.cumWabv;
    opVarVals["ccvi"] = vars.cCVI;
    opVarVals["gpp"] = vars.GPP;
    opVarVals["cgpp"] = vars.cGPP;
    opVarVals["ft"] = vars.fT;
    opVarVals["fvpd"] = vars.fVPD;
    opVarVals["fsw"] = vars.fSW;
    opVarVals["fnutr"] = vars.fNutr;
    opVarVals["ffrost"] = vars.fFrost;
    opVarVals["fage"] = vars.fAge;
    opVarVals["physmod"] = vars.PhysMod;
    opVarVals["apar"] = vars.APAR;
    opVarVals["aparu"] = vars.APARu;
    opVarVals["npp"] = vars.NPP;
    opVarVals["totalw"] = vars.TotalW;
    opVarVals["transp"] = vars.Transp;
    opVarVals["evaptransp"] = vars.EvapTransp;
    opVarVals["agemaix"] = vars.ageMAIx;
    opVarVals["agelaix"] = vars.ageLAIx;
    opVarVals["cvi"] = vars.CVI;
    opVarVals["cwue"] = vars.cWUE;
    opVarVals["delwag"] = vars.delWAG;
}

// This is the main routine for the 3PG model
void runTreeModel(long cellIndex, DataInput& dataInput, DataOutput& dataOutput)
{
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    step 1: initialize variables
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    std::unordered_map<std::string, double> opVarVals;
    Vars vars;

    // load input parameters, and exit if it's a NODATA pixel
    InputParams params;
    if (!dataInput.getInputParams(cellIndex, params)) {
        return;
    }

    //get the overall run period (not just of this pixel)
    RunPeriod runPeriod = dataInput.getRunPeriod();

    double Irrig;                            // current annual irrigation (ML/y)

    //various necessary variables (hopefully better comment to come)
    double TranspScaleFactor;
    double poolFractn = 0;
    double excessSW;
    double pooledSW = 0;
    double RainIntcptn = 0;

    // At initialisation param file has only possible value to use.  
    double MinASW = params.MinASWp;

    //soil parameters for soil class
    double SWconst;
    double SWpower;

    //day length (by month)
    double mDayLength[13];

    // Stand factors that are specifically age dependent
    double SLA;
    double gammaF;
    double CanCover;

    // Derive some parameters
    double pfsPower = log(params.pFS20 / params.pFS2) / log(10);
    double pfsConst = params.pFS2 / pow(2, pfsPower);

    double fCalphax = params.fCalpha700 / (2 - params.fCalpha700);
    double fCg0 = params.fCg700 / (2 * params.fCg700 - 1);

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

    int dayofyr;
    double MoistRatio; // PhysMod has been moved
    double wSmax;

    double oldVol;    //Added 16/07/02 as part of CVI

    //New Soilwater modifier adjuster
    double ASWmod;

    // monthly met data
    double dayLength;

    // co2 modifiers
    double Co2Slope, curYearCo2, fCalpha, fCg;

    // 3PGS - variables for 3PGS
    double FPAR_AVH;

    // Compute daylengths
    dayofyr = -15;
    for (int mn = 1; mn <= 12; mn++) {
        dayofyr = dayofyr + DaysInMonth[mn];
        mDayLength[mn] = 86400 * getDayLength(params.Lat, dayofyr);
    }

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
    double StandAge = params.StartAge;
    vars.StemNo = params.StemNoi;

    if (dataInput.haveSeedlingMass) {
        params.WFi = (0.5 * params.StemNoi * params.SeedlingMass) / pow(10, 6);
        params.WRi = (0.25 * params.StemNoi * params.SeedlingMass) / pow(10, 6);
        params.WSi = params.WRi;
    }

    vars.WS = params.WSi;
    vars.WF = params.WFi;
    vars.WR = params.WRi;

    vars.ASW = params.ASWi;
    vars.TotalLitter = 0;

    AvStemMass = vars.WS * 1000 / vars.StemNo;                             //  kg/tree
    vars.avDBH = pow((AvStemMass / params.StemConst), (1 / params.StemPower));
    vars.BasArea = ((pow((vars.avDBH / 200), 2)) * Pi) * vars.StemNo;
    SLA = params.SLA1 + (params.SLA0 - params.SLA1) * exp(-ln2 * pow((StandAge / params.tSLA), 2)); //Modified StandAge
    vars.LAI = vars.WF * SLA * 0.1;
    vars.cLAI = vars.LAI; // set cLAI to LAI for initial monthly output

    vars.fracBB = params.fracBB1 + (params.fracBB0 - params.fracBB1) * exp(-ln2 * (StandAge / params.tBB)); //Modified StandAge
    Density = params.rhoMax + (params.rhoMin - params.rhoMax) * exp(-ln2 * (StandAge / params.tRho));
    gammaF = getGammaFoliage(params, StandAge);

    vars.StandVol = vars.WS * (1 - vars.fracBB) / Density;
    oldVol = vars.StandVol;

    vars.MAI = 0;
    if (StandAge > 0) {
        vars.MAI = vars.StandVol / StandAge;    //UnModified StandAge
    }
        
    // Compute slope of CO2 growth over the run period
    Co2Slope = (params.CO2End - params.CO2Start) / (runPeriod.EndYear - runPeriod.StartYear);
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    step 2: load management params for startage
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    //set fertility to management value if it exists, otherwise set to param
    //also set fertility decay flag if we don't have the management param but do have params FRStart, FREnd, and FRDec
    bool fertilityDecay = false;
    if (!dataInput.getManagementParam(ManagementIndex::FERTILITY, cellIndex, (int)params.yearPlanted + (int)params.StartAge, MinASW)) {
        fertilityDecay = dataInput.haveAgeDepFert;
        vars.FR = params.FRp;
    }

    //set MinASW to management value if it exists, otherwise set to param
    if (!dataInput.getManagementParam(ManagementIndex::MINASW, cellIndex, (int)params.yearPlanted + (int)params.StartAge, MinASW)) {
        MinASW = params.MinASWp;
    }

    //set irrigation to management value if it exists, otherwise set to 0
    if (!dataInput.getManagementParam(ManagementIndex::IRRIGATION, cellIndex, (int)params.yearPlanted + (int)params.StartAge, Irrig)) {
        Irrig = 0;
    }

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    step 3: write start month initialization values to output
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    //Find out if there is supposed to be any data here in the first place...

    SeriesParams sParams;
    if (!dataInput.getSeriesParams(cellIndex, (int)params.yearPlanted + (int)params.StartAge, (int)params.StartMonth, sParams)) {
        return;
    }

    if (sParams.FrostDays > 30) {
        sParams.FrostDays = 30;
    }

    // init pool fraction
    poolFractn = std::max(0.0, std::min(1.0, poolFractn));

    vars.FR = params.FRp;

    //write initial state of output variables
    copyVars(vars, opVarVals);
    dataOutput.writeMonthlyOutputGrids(opVarVals, (int)params.yearPlanted + (int)params.StartAge, (int)params.StartMonth, cellIndex);

    //set cumuliative yearly LAI back to 0
    vars.cLAI = 0;

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    step 4: determine start and end month/years
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    //first year might be different than yearPlanted + startAge! (This happens if startMonth == 12)
    //calculate the starting iteration month/year using the start month.
    //the first is the month following the start month, since the start month state is the initial state.

    //start monthly/yearly iteration on StartMonth + 1 unless StartMonth is 12, then start on month 1 of the following year
    int firstYear = ((int)params.StartMonth == 12) ? (int)params.yearPlanted + (int)params.StartAge + 1 : (int)params.yearPlanted + (int)params.StartAge;
    int firstMonth = ((int)params.StartMonth == 12) ? 1 : (int)params.StartMonth + 1;

    //end monthly/yearly iteration on the last month and year according to spMaxMY
    int lastYear = (int)params.EndYear;
    int lastMonth = (int)params.StartMonth;

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    step 5: start yearly processing loop
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    for (int year = firstYear; year <= lastYear; year++) {

        //find start and end month for the current year iteration
        int startMonth = (year == firstYear) ? firstMonth : 1;
        int endMonth = (year == lastYear) ? lastMonth : 12;

        //get management params for current year
        //set fertility to management value if it exists, otherwise set to param (and set fertilityDecay flag)
        fertilityDecay = false;
        if (!dataInput.getManagementParam(ManagementIndex::FERTILITY, cellIndex, year, MinASW)) {
            fertilityDecay = dataInput.haveAgeDepFert;
            vars.FR = params.FRp;
        }

        //set MinASW to management value if it exists, otherwise set to param
        if (!dataInput.getManagementParam(ManagementIndex::MINASW, cellIndex, year, MinASW)) {
            MinASW = params.MinASWp;
        }

        //set irrigation to management value if it exists, otherwise set to 0
        if (!dataInput.getManagementParam(ManagementIndex::IRRIGATION, cellIndex, year, Irrig)) {
            Irrig = 0;
        }

        //calculate the co2 for the current year
        curYearCo2 = params.CO2Start + Co2Slope * (year - runPeriod.StartYear);

        /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        step 6: start monthly processing loop
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
        for (int month = startMonth; month <= endMonth; month++) {
            //get series parameters for current year and month
            if (!dataInput.getSeriesParams(cellIndex, year, month, sParams)) {
                return;
            }

            //get current day length based on month
            dayLength = mDayLength[month];

            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            step 7: determine various environmental modifiers
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

            //calculate fertility using fertility decay if it applies to the current month on the current pixel
            //ie. If we DON'T have fertility management param, and we DO have 'FRStart' 'FREnd' and 'FRDec' input params
            if (fertilityDecay && (params.FRstart <= StandAge) && (params.FRend > StandAge)) {
                vars.FR = vars.FR - vars.FR * params.FRdec;
            }

            // calculate temperature response function to apply to alpha
            vars.fT = 0;
            if (!(sParams.Tavg <= params.growthTmin) && !(sParams.Tavg >= params.growthTmax)) {
                vars.fT = ((sParams.Tavg - params.growthTmin) / (params.growthTopt - params.growthTmin)) *
                    pow(((params.growthTmax - sParams.Tavg) / (params.growthTmax - params.growthTopt)),
                        ((params.growthTmax - params.growthTopt) / (params.growthTopt - params.growthTmin)));
            }

            // calculate VPD modifier
            vars.fVPD = exp(-params.CoeffCond * sParams.VPD);

            // calculate soil water modifier
            ASWmod = 1;
            if (dataInput.haveMinASWTG && params.MaxASW <= params.MinASWTG) {
                double dAdjMod = (params.MaxASW - params.MinASWTG) / params.MinASWTG;
                ASWmod = pow(2.718281828459045235, dAdjMod);
            }
            MoistRatio = ASWmod * vars.ASW / params.MaxASW;
            vars.fSW = 1 / (1 + pow(((1 - MoistRatio) / SWconst), SWpower));

            // calculate soil nutrition
            vars.fNutr = 1;
            if (params.fNn != 0) {
                vars.fNutr = 1 - (1 - params.fN0) * pow((1 - vars.FR), params.fNn);
            }

            // calculate frost modifier
            vars.fFrost = 1 - params.kF * (sParams.FrostDays / 30.0);

            // calculate co2 modifiers
            fCalpha = fCalphax * curYearCo2 / (350 * (fCalphax - 1) + curYearCo2);
            fCg = fCg0 / (1 + (fCg0 - 1) * curYearCo2 / 350);

            // calculate age modifier
            vars.fAge = 1;
            if (!dataInput.modelMode3PGS) {
                double RelAge = StandAge / params.MaxAge;  //Modified StandAge
                vars.fAge = (1 / (1 + pow((RelAge / params.rAge), params.nAge)));
            }

            // calculate physiological modifier applied to conductance and alpha
            vars.PhysMod = std::min(vars.fVPD, vars.fSW) * vars.fAge;

            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            step 8: determine gross and net biomass production
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

            // canopy cover and light interception.
            CanCover = 1;
            if ((params.fullCanAge > 0) && (StandAge < params.fullCanAge)) {
                CanCover = (StandAge) / params.fullCanAge; //Modified StandAge
            }
            lightIntcptn = (1 - (exp(-params.k * vars.LAI / CanCover)));

            // Calculate FPAR_AVH and LAI from NDVI data if in 3PGS mode 
            if (dataInput.modelMode3PGS) {
                // Initial value of FPAR_AVH from linear fit. 
                FPAR_AVH = (sParams.NDVI_AVH * params.NDVI_FPAR_constant) + params.NDVI_FPAR_intercept;

                // Constrain FPAR_AVH to within threshhold values. 
                if (FPAR_AVH > 0.98) {
                    FPAR_AVH = 0.98;
                }
                else if (FPAR_AVH < 0) {
                    FPAR_AVH = 0;
                }

                // LAI
                if (FPAR_AVH < 0.05) {
                    vars.LAI = 0.2;
                }
                else {
                    vars.LAI = -2.0 * log(1 - FPAR_AVH);
                }
            }

            // Calculate PAR, APAR, and APARu. Use FPAR_AVG if in 3PGS mode
            RAD = sParams.SolarRad * DaysInMonth[month];
            PAR = RAD * params.molPAR_MJ;
            if (dataInput.modelMode3PGS) {
                vars.APAR = PAR * FPAR_AVH;
            }
            else {
                vars.APAR = PAR * lightIntcptn * CanCover;
            }
            vars.APARu = vars.APAR * vars.PhysMod;

            //calculate NPP
            vars.alphaC = params.alpha * vars.fNutr * vars.fT * vars.fFrost * vars.PhysMod * fCalpha;   //22-07-02 for Excel March beta consis.
            epsilon = params.gDM_mol * params.molPAR_MJ * vars.alphaC;
            RADint = RAD * lightIntcptn * CanCover;
            vars.GPP = epsilon * RADint / 100;
            vars.NPP = vars.GPP * params.y;

            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            step 9: determine water balance
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

            // calculate canopy conductance
            double gC;
            if (vars.LAI <= params.LAIgcx) {
                gC = params.MaxCond * vars.LAI / params.LAIgcx;
            }
            else {
                gC = params.MaxCond;
            }
            CanCond = gC * vars.PhysMod * fCg;

            // calculate transpiration from Penman-Monteith (mm/day converted to mm/month)
            vars.Transp = CanopyTranspiration(
                sParams.SolarRad, 
                sParams.VPD, 
                dayLength, 
                params.BLcond,
                CanCond, 
                sParams.NetRad, 
                dataInput.haveNetRadParam(), 
                params
            );
            vars.Transp = DaysInMonth[month] * vars.Transp;

            // rainfall interception
            Interception = params.MaxIntcptn;
            if (params.LAImaxIntcptn > 0) {
                Interception = params.MaxIntcptn * std::min(1.0, vars.LAI / params.LAImaxIntcptn);
            }
            RainIntcptn = sParams.Rain * Interception;
            
            // water balance
            vars.ASW = vars.ASW + sParams.Rain + (100 * Irrig / 12) + pooledSW;        //Irrig is Ml/ha/year
            vars.EvapTransp = std::min(vars.ASW, vars.Transp + RainIntcptn);
            excessSW = std::max(vars.ASW - vars.EvapTransp - params.MaxASW, 0.0);
            vars.ASW = vars.ASW - vars.EvapTransp - excessSW;
            pooledSW = poolFractn * excessSW;

            // correct for actual ET
            TranspScaleFactor = vars.EvapTransp / (vars.Transp + RainIntcptn);      
            vars.GPP = TranspScaleFactor * vars.GPP;
            vars.NPP = TranspScaleFactor * vars.NPP;
            if (vars.EvapTransp != 0) {
                vars.WUE = 100 * vars.NPP / vars.EvapTransp;
            }
                   
            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            step 10: determine biomass increments and losses
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

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
            if (!dataInput.modelMode3PGS) {
                vars.WF = vars.WF + delWF - delLitter;
                vars.WR = vars.WR + delWR - delRloss;
                vars.WS = vars.WS + delWS;
                vars.TotalW = vars.WF + vars.WR + vars.WS;
                vars.TotalLitter = vars.TotalLitter + delLitter;
            }

            //increment stand age with yearly unit
            StandAge = StandAge + 1.0 / 12.0;

            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            step 11: update tree and stand data
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
            if (!dataInput.modelMode3PGS) {
                //Calculate mortality
                wSmax = params.wSx1000 * pow((1000 / vars.StemNo), params.thinPower);
                AvStemMass = vars.WS * 1000 / vars.StemNo;
                delStems = 0;
                if (wSmax < AvStemMass) {
                    delStems = getMortality(vars.StemNo, vars.WS, params);
                    vars.WF = vars.WF - params.mF * delStems * (vars.WF / vars.StemNo);
                    vars.WR = vars.WR - params.mR * delStems * (vars.WR / vars.StemNo);
                    vars.WS = vars.WS - params.mS * delStems * (vars.WS / vars.StemNo);
                    vars.StemNo = vars.StemNo - delStems;
                    wSmax = params.wSx1000 * pow((1000 / vars.StemNo), params.thinPower);
                    AvStemMass = vars.WS * 1000 /  vars.StemNo;
                }

                //update age-dependent factors
                SLA = params.SLA1 + (params.SLA0 - params.SLA1) * exp(-ln2 * pow((StandAge / params.tSLA), 2));  //Modified StandAge
                vars.fracBB = params.fracBB1 + (params.fracBB0 - params.fracBB1) * exp(-ln2 * (StandAge / params.tBB));  //Modified StandAge
                Density = params.rhoMax + (params.rhoMin - params.rhoMax) * exp(-ln2 * (StandAge / params.tRho));
                gammaF = getGammaFoliage(params, StandAge);

                //update stand characteristics
                vars.LAI = vars.WF * SLA * 0.1;
                vars.avDBH = pow((AvStemMass / params.StemConst), (1 / params.StemPower));
                vars.BasArea = (pow((vars.avDBH / 200), 2) * Pi) * vars.StemNo;
                vars.StandVol = vars.WS * (1 - vars.fracBB) / Density;

                vars.CVI =  vars.StandVol - oldVol;
                oldVol = vars.StandVol;

                vars.MAI = 0;
                if (StandAge > 0) {
                    vars.MAI = vars.StandVol / StandAge;  //UnModified StandAge
                }

                // Update accumulated totals
                vars.cGPP = vars.cGPP + vars.GPP;
                vars.cNPP = vars.cNPP + vars.NPP;
                vars.cCVI = vars.cCVI + vars.CVI;
                vars.cLitter = vars.cLitter + delLitter;
                vars.cTransp = vars.cTransp + vars.Transp;
                vars.cEvapTransp = vars.cEvapTransp + vars.EvapTransp;
                vars.cWUE = 100 * vars.cNPP / vars.cEvapTransp;
                vars.cLAI = vars.cLAI +  vars.LAI / 12.0;

                // Accumulate biomass increments and LAI
                vars.cumWabv = vars.cumWabv + delWF + delWS - delLitter;  // ANL - PROBLEM?  

                vars.MAI = 0;
                if (StandAge > 0) {
                    vars.MAI = vars.StandVol / StandAge;
                }

                // Update accumulated totals
                vars.cGPP += vars.GPP;
                vars.cNPP += vars.NPP;
                vars.cCVI += vars.CVI;
                vars.cLitter += delLitter;
                vars.cTransp += vars.Transp;
                vars.cEvapTransp += vars.EvapTransp;
                vars.cWUE = 100 * vars.cNPP / vars.cEvapTransp;
                vars.cLAI +=  vars.LAI / 12;
                vars.cumWabv += delWF + delWS - delLitter;
            }
            else {
                // 3PGS
                vars.delWAG = vars.NPP * (1 - pR);
                vars.cumWabv += vars.delWAG;
            }

            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            step 12: write monthly output if not in final month
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
            if (!(year == params.EndYear && month == params.StartMonth)) {
                copyVars(vars, opVarVals);
                dataOutput.writeMonthlyOutputGrids(opVarVals, year, month, cellIndex);
            }

            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            step 13: reset cumulative variables and run year-end calculations.
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
            if (month == params.StartMonth) {
                //year end calculations, run after the model has run for 12 months -- NOT on December every year

                // Update some stand characteristics
                vars.fracBB = params.fracBB1 + (params.fracBB0 - params.fracBB1) * exp(-ln2 * (StandAge / params.tBB));  //Modified StandAge
                vars.StandVol = vars.WS * (1 - vars.fracBB) / Density;
                vars.MAI = 0;
                if (StandAge > 0) {
                    vars.MAI = vars.StandVol / StandAge;
                }


                // Determine peak LAI & MAI and age at peaks
                // 
                //cLAI is the average LAI across the year. This is how the excel 
                //VBA code does it (which is our reference) so that's how I'm doing it
                // - Joe
                if (vars.cLAI > vars.LAIx) {
                    vars.LAIx = vars.cLAI;
                    vars.ageLAIx = StandAge;
                }
                if (vars.MAI > vars.MAIx) {
                    vars.MAIx = vars.MAI;
                    vars.ageMAIx = StandAge;
                }

                // Restore LAI
                vars.LAI = vars.WF * SLA * 0.1;

                // reset cumulative variables
                vars.cumWabv = 0;
                vars.cLAI = 0;
                vars.cCVI = 0;
                vars.cNPP = 0;
                vars.cGPP = 0;
                vars.cTransp = 0;
                vars.cEvapTransp = 0;
                vars.cLitter = 0;
            }
        }
    }
    copyVars(vars, opVarVals);
    dataOutput.writeOutputGrids(opVarVals, cellIndex);
}