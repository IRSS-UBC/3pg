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

double getGammaFoliage(InputParams& params, double& StandAge) {
    double gf;
    gf = params.gammaFx* params.gammaF0 /
        (params.gammaF0 + (params.gammaFx - params.gammaF0) *
            exp(-12 * log(1 + params.gammaFx / params.gammaF0) * StandAge / params.tgammaF));
    return gf;
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
void runTreeModel(MYDate spMinMY, MYDate spMaxMY, long cellIndex, DataInput& dataInput, DataOutput& dataOutput)
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

    //various necessary variables (hopefully better comment to come)
    double aEpsilonGross;
    double aEpsilonStem;
    double avDBHi;
    double CumStemLoss;
    double GPPmolc;
    double LAIi;
    double TranspScaleFactor;
    double RunOff;
    double monthlyIrrig;
    double poolFractn = 0;
    double excessSW;
    double pooledSW = 0;
    double RainIntcptn = 0;

    // Initialise cumulative variables
    double CumdelWF = 0;
    double CumdelWR = 0;
    double CumdelWS = 0;
    double CumAPARU = 0;
    double cumARAD = 0;
    double cumLAI = 0;
    double cumGPP = 0;
    double cumTransp = 0;
    double cumEvapTransp = 0;
    double cumIrrig = 0;

    //Initialise annual cumulative variables
    double aStemDM = 0;
    double aRADint = 0;
    double aGPP = 0;
    double aNPP = 0;
    double aEvapTransp = 0;
    double aTransp = 0;
    double aSupIrrig = 0;

    //Initialise output step cumulative variables
    double delStemNo = 0;
    double cRADint = 0;
    double cRainInt = 0;
    double cStemDM = 0;
    double cSupIrrig = 0;

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

    double oldVol;    //Added 16/07/02 as part of CVI

    //New Soilwater modifier adjuster
    bool useMinASWTG = false;
    double ASWmod;

    // monthly met data
    double dayLength;

    int thinEventNo, defoltnEventNo;

    bool haveAvgTempSeries = false;  //Change: - Needs to be available iff 
    //VPD is available

    // 3PGS - variables for 3PGS
    double FPAR_AVH;

    // Compute daylengths
    dayofyr = -15;
    for (int mn = 1; mn <= 12; mn++) {
        dayofyr = dayofyr + DaysInMonth[mn];
        mDayLength[mn] = 86400 * getDayLength(params.Lat, dayofyr);
    }

    useMinASWTG = dataInput.haveMinASWTG;

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

    if (dataInput.haveSeedlingMass)
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
    gammaF = getGammaFoliage(params, StandAge);

    vars.StandVol = vars.WS * (1 - vars.fracBB) / Density;
    oldVol = vars.StandVol;

    if (StandAge > 0)
        vars.MAI = vars.StandVol / StandAge;    //UnModified StandAge
    else vars.MAI = 0;

    avDBHi = vars.avDBH;
    LAIi = vars.LAI;
    CumStemLoss = 0;
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    step 2: load management params if they exist
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    //set fertility to management value if it exists, otherwise set to param
    //also set fertility decay flag if we don't have the management param but do have params FRStart, FREnd, and FRDec
    bool fertilityDecay = false;
    if (!dataInput.getManagementParam(ManagementIndex::FERTILITY, cellIndex, spMinMY.year, MinASW)) {
        fertilityDecay = dataInput.haveAgeDepFert;
        vars.FR = params.FRp;
    }

    //set MinASW to management value if it exists, otherwise set to param
    if (!dataInput.getManagementParam(ManagementIndex::MINASW, cellIndex, spMinMY.year, MinASW)) {
        MinASW = params.MinASWp;
    }

    //set irrigation to management value if it exists, otherwise set to 0
    if (!dataInput.getManagementParam(ManagementIndex::IRRIGATION, cellIndex, spMinMY.year, Irrig)) {
        Irrig = 0;
    }

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    step 3: write start month initialization values to output
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    //Find out if there is supposed to be any data here in the first place...

    SeriesParams sParams;
    if (!dataInput.getSeriesParams(cellIndex, spMinMY.year, (int)params.StartMonth, sParams)) {
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
    dataOutput.writeMonthlyOutputGrids(opVarVals, spMinMY.year, (int)params.StartMonth, spMinMY, spMaxMY, cellIndex);

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    step 4: determine start and end month/years
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    //calculate the starting iteration month/year using the start month.
    //the first is the month following the start month, since the start month state is the initial state.

    //start monthly/yearly iteration on StartMonth + 1 unless StartMonth is 12, then start on month 1 of the following year
    int firstYear = ((int)params.StartMonth == 12) ? spMinMY.year + 1 : spMinMY.year;
    int firstMonth = ((int)params.StartMonth == 12) ? 1 : params.StartMonth + 1;

    //end monthly/yearly iteration on the last month and year according to spMaxMY
    int lastYear = spMaxMY.year;
    int lastMonth = spMaxMY.mon;

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    step 5: start yearly processing loop
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    for (int year = firstYear; year <= lastYear; year++) {

        //find start and end month for the current year iteration
        int startMonth = (year == firstYear) ? firstMonth : 1;
        int endMonth = (year == lastYear) ? lastMonth : 12;

        /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        step 6: start monthly processing loop
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
        for (int month = startMonth; month <= endMonth; month++) {

            if (!dataInput.getSeriesParams(cellIndex, year, month, sParams)) {
                return;
            }

            dayLength = mDayLength[month];

            // Determine the various environmental modifiers

            //calculate fertility using fertility decay if it applies to the current month on the current pixel
            if (fertilityDecay && (params.FRstart <= StandAge) && (params.FRend > StandAge)) {
                vars.FR = vars.FR - vars.FR * params.FRdec;
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
            RAD = sParams.SolarRad * DaysInMonth[month];        // MJ/m^2
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
            vars.GPP = epsilon * RADint / 100;               // tDM/ha
            vars.NPP = vars.GPP * params.y;                            // assumes respiratory rate is constant

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
                CanCond, sParams.NetRad, dataInput.haveNetRadParam(), params);
            vars.Transp = DaysInMonth[month] * vars.Transp;

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
            vars.GPP = TranspScaleFactor * vars.GPP;
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

            StandAge = StandAge + 1.0 / 12.0;

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
                gammaF = getGammaFoliage(params, StandAge);

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
                vars.cGPP = vars.cGPP + vars.GPP;
                aGPP = aGPP + vars.GPP;
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

            if (!(year == spMaxMY.year && month == spMaxMY.mon)) {
                //the !(year == maxMY.year && calMonth == maxMY.mon) is so that at the last iteration we don't write to a monthly
                //output, but rather skip it and the values are eventually written via writeOutputGrids().
                copyVars(vars, opVarVals);
                dataOutput.writeMonthlyOutputGrids(opVarVals, year, month, spMinMY, spMaxMY, cellIndex);
            }

            //reset cumulative variables and run year-end calculations
            if (month == params.StartMonth) {
                // reset cumulative variables
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

                // reset annual cumulative variables
                aStemDM = 0;
                aRADint = 0;
                aGPP = 0;
                aNPP = 0;
                aEvapTransp = 0;
                aTransp = 0;
                aSupIrrig = 0;

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

                //get new yearly management params
                //set fertility to management value if it exists, otherwise set to param (and set fertilityDecay flag)
                fertilityDecay = false;
                if (!dataInput.getManagementParam(ManagementIndex::FERTILITY, cellIndex, spMinMY.year, MinASW)) {
                    fertilityDecay = dataInput.haveAgeDepFert;
                    vars.FR = params.FRp;
                }

                //set MinASW to management value if it exists, otherwise set to param
                if (!dataInput.getManagementParam(ManagementIndex::MINASW, cellIndex, spMinMY.year, MinASW)) {
                    MinASW = params.MinASWp;
                }

                //set irrigation to management value if it exists, otherwise set to 0
                if (!dataInput.getManagementParam(ManagementIndex::IRRIGATION, cellIndex, spMinMY.year, Irrig)) {
                    Irrig = 0;
                }

                //year end calculations, run after the model has run for 12 months -- NOT on December every year

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

                // Restore LAI
                vars.LAI = vars.WF * SLA * 0.1;
            }
        }

        // Restore LAI here too... (this is how the model worked before so I'm keeping it -- Joe)
        vars.LAI = vars.WF * SLA * 0.1;

    }
    copyVars(vars, opVarVals);
    dataOutput.writeOutputGrids(opVarVals, cellIndex);
}