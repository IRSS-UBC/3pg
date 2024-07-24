#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include "GDALRasterImage.hpp"
#include "ParamStructs.hpp"
#include "MYDate.h"

struct InputParams {
	double pFS2;
	double pFS20;
	double StemConst;
	double StemPower;
	double pRx;
	double pRn;
	double growthTmin;
	double growthTopt;
	double growthTmax;
	double kF;
	double gammaFx;
	double gammaF0;
	double tgammaF;
	double Rttover;
	double MaxCond;
	double CoeffCond;
	double BLcond;
	double m0;
	double fN0;
	double fNn;
	double thinPower;
	double mF;
	double mR;
	double mS;
	double SWconst0;
	double SWpower0;
	double wSx1000;
	double MaxAge;
	double nAge;
	double rAge;
	double SLA0;
	double SLA1;
	double tSLA;
	double k;
	double fullCanAge;
	double alpha;
	double fracBB0;
	double fracBB1;
	double tBB;
	double y;
	double rhoMin;
	double rhoMax;
	double tRho;
	double Qa;
	double Qb;
	double gDM_mol;
	double molPAR_MJ;
	double LAIgcx;
	double MaxIntcptn;
	double LAImaxIntcptn;
	double Lat;
	double FRp;
	double FRstart;
	double FRend;
	double FRdec;
	double soilIndex;
	double MaxASW;
	double MinASWp;
	double StartAge;
	double EndYear;
	double StartMonth;
	double yearPlanted;
	double SeedlingMass;
	double WFi;
	double WRi;
	double WSi;
	double StemNoi;
	double ASWi;
	double MinASWTG;
	double NDVI_FPAR_intercept;
	double NDVI_FPAR_constant;
};

struct SeriesParams {
	double Tavg;
	double Rain;
	double SolarRad;
	double FrostDays;
	double NDVI_AVH;
	double NetRad;
	double VPD;
};

typedef enum {
	NONE = -1,
	TMAX = 0,
	TMIN = 1,
	TAVG = 2,
	RAIN = 3,
	SOLAR_RAD = 4,
	FROST_DAYS = 5,
	NDVI_AVH = 6,
	NET_RAD = 7,
	VPD = 8,
} SeriesIndex;

class DataInput {
private:
	//maps and sets for dealing with input parameters
	std::unordered_map<std::string, std::string> inputParamNames = {
		{"Foliage:stem partitioning ratio @ D=2 cm", "pFS2"},
		{"Foliage:stem partitioning ratio @ D=20 cm", "pFS20"},
		{"Constant in the stem mass v. diam. relationship", "StemConst"},
		{"Power in the stem mass v. diam. relationship", "StemPower"},
		{"Maximum fraction of NPP to roots", "pRx"},
		{"Minimum fraction of NPP to roots", "pRn"},
		{"Minimum temperature for growth", "growthTmin"},
		{"Optimum temperature for growth", "growthTopt"},
		{"Maximum temperature for growth", "growthTmax"},
		{"Days production lost per frost day", "kF"},
		{"Maximum litterfall rate", "gammaFx"},
		{"Litterfall rate at t = 0", "gammaF0"},
		{"Age at which litterfall rate has median value", "tgammaF"},
		{"Average monthly root turnover rate", "Rttover"},
		{"Maximum canopy conductance", "MaxCond"},
		{"Defines stomatal response to VPD", "CoeffCond"},
		{"Canopy boundary layer conductance", "BLcond"},
		{"Value of 'm' when FR = 0", "m0"},
		{"Value of 'fNutr' when FR = 0", "fN0"},
		{"Power of (1-FR) in 'fNutr'", "fNn"},
		{"Moisture ratio deficit for fq = 0.5", "SWconst0"},
		{"Power of moisture ratio deficit", "SWpower0"},
		{"Max. stem mass per tree @ 1000 trees/hectare", "wSx1000"},
		{"Power in self-thinning rule", "thinPower"},
		{"Fraction mean single-tree foliage biomass lost per dead tree", "mF"},
		{"Fraction mean single-tree root biomass lost per dead tree", "mR"},
		{"Fraction mean single-tree stem biomass lost per dead tree", "mS"},
		{"Maximum stand age used in age modifier", "MaxAge"},
		{"Power of relative age in function for fAge", "nAge"},
		{"Relative age to give fAge = 0.5", "rAge"},
		{"Specific leaf area at age 0", "SLA0"},
		{"Specific leaf area for mature leaves", "SLA1"},
		{"Age at which specific leaf area = (SLA0+SLA1)/2", "tSLA"},
		{"Extinction coefficient for absorption of PAR by canopy", "k"},
		{"Age at canopy cover", "fullCanAge"},
		{"Canopy quantum efficiency", "alpha"},
		{"Branch and bark fraction at age 0", "fracBB0"},
		{"Branch and bark fraction for mature stands", "fracBB1"},
		{"Age at which fracBB = (fracBB0+fracBB1)/2", "tBB"},
		{"Ratio NPP/GPP", "y"},
		{"Basic density", "Density"},
		{"Intercept of net v. solar radiation relationship", "Qa"},
		{"Slope of net v. solar radiation relationship", "Qb"},
		{"Molecular weight of dry matter", "gDM_mol"},
		{"Conversion of solar radiation to PAR", "molPAR_MJ"},
		{"LAI for maximum canopy conductance", "LAIgcx"},
		{"Maximum proportion of rainfall evaporated from canopy", "MaxIntcptn"},
		{"LAI for maximum rainfall interception", "LAImaxIntcptn"},
		{"Latitude", "Lat"},
		{"FR", "FRp"},
		{"Fertility rating", "FRp"},
		{"Soil Index", "soilIndex"},
		{"Soil class", "soilIndex"},
		{"Maximum ASW", "MaxASW"},
		{"MinASW", "MinASWp"},
		{"Minimum ASW", "MinASWp"},
		{"Initial age", "StartAge"},
		{"Start age", "StartAge"},
		{"End year", "EndYear"},
		{"Start Month", "StartMonth"},
		{"Start month", "StartMonth"},
		{"Seedling Mass", "SeedlingMass"},
		{"Seedling mass", "SeedlingMass"},
		{"W foliage", "WFi"},
		{"W root", "WRi"},
		{"W stem", "WSi"},
		{"Stem no", "StemNoi"},
		{"Initial soil water", "ASWi"},
		{"Minimum basic density - for young trees", "rhoMin"},
		{"Maximum basic density - for older trees", "rhoMax"},
		{"Age at which rho = (rhoMin+rhoMax)/2", "tRho"},
		{"Year Planted", "yearPlanted"},
	};
	std::unordered_set<std::string> allInputParams = {
		"pFS2",
		"pFS20",
		"StemConst",
		"StemPower",
		"pRx",
		"pRn",
		"growthTmin",
		"growthTopt",
		"growthTmax",
		"kF",
		"gammaFx",
		"gammaF0",
		"tgammaF",
		"Rttover",
		"MaxCond",
		"CoeffCond",
		"BLcond",
		"m0",
		"fN0",
		"fNn",
		"thinPower",
		"mF",
		"mR",
		"mS",
		"SWconst0",
		"SWpower0",
		"wSx1000",
		"MaxAge",
		"nAge",
		"rAge",
		"SLA0",
		"SLA1",
		"tSLA",
		"k",
		"fullCanAge",
		"alpha",
		"fracBB0",
		"fracBB1",
		"tBB",
		"y",
		"rhoMin",
		"rhoMax",
		"tRho",
		"Qa",
		"Qb",
		"gDM_mol",
		"molPAR_MJ",
		"LAIgcx",
		"MaxIntcptn",
		"LAImaxIntcptn",
		"Lat",
		"FRp",
		"FRstart",
		"FRend",
		"FRdec",
		"soilIndex",
		"MaxASW",
		"MinASWp",
		"StartAge",
		"EndYear",
		"StartMonth",
		"yearPlanted",
		"SeedlingMass",
		"WFi",
		"WRi",
		"WSi",
		"StemNoi",
		"ASWi",
		"MinASWTG",
		"NDVI_FPAR_intercept",
		"NDVI_FPAR_constant",
	};
	std::unordered_set<std::string> requiredInputParams3PG = {
		"pFS2", "pFS20", "StemConst", "StemPower", "pRx", "pRn",
		"growthTmin", "growthTopt", "growthTmax",
		"kF",
		"gammaFx", "gammaF0", "tgammaF", "Rttover",
		"MaxCond", "CoeffCond", "BLcond",
		"m0", "fN0", "fNn",
		"SWconst0", "SWpower0",
		"wSx1000",
		"MaxAge", "nAge", "rAge",
		"SLA0", "SLA1", "tSLA", "k", "fullCanAge",
		"alpha", "fracBB0", "fracBB1", "tBB",
		"y",
		"Lat", "FRp", "soilIndex", "MaxASW", "MinASWp",
		"StartAge", "EndYear",
		"StemNoi", "ASWi", "yearPlanted", 
		"Qa", "Qb",
		"gDM_mol", "molPAR_MJ",
		"LAIgcx", "MaxIntcptn",
		"StartMonth",
		"LAImaxIntcptn",
		"thinPower", "mF", "mR", "mS",
	};
	std::unordered_set<std::string> requiredInputParams3PGS = {
		"growthTmin", "growthTopt", "growthTmax",
		"kF",
		"MaxCond", "CoeffCond", "BLcond",
		"m0", "fN0",
		"SWconst0", "SWpower0",
		"SLA1", "alpha",
		"y",
		"Lat", "FRp", "soilIndex", "MaxASW", "MinASWp",
		"StartAge","EndYear",
		"NDVI_FPAR_intercept", "NDVI_FPAR_constant",
		"Qa", "Qb",
		"gDM_mol", "molPAR_MJ",
		"LAIgcx", "MaxIntcptn",
		"StartMonth",
		"LAImaxIntcptn",
	}; //TODO do we need to add yearPlanted here???
	std::unordered_map<std::string, PPPG_PARAM> inputParams;

	//maps and sets for dealing with series parameters
	std::unordered_map<std::string, SeriesIndex> seriesParamNameMap = {
		{"Tmax", SeriesIndex::TMAX},
		{"Tmin", SeriesIndex::TMIN},
		{"Tavg", SeriesIndex::TAVG},
		{"Rain", SeriesIndex::RAIN},
		{"Solar radtn", SeriesIndex::SOLAR_RAD},
		{"Frost", SeriesIndex::FROST_DAYS},
		{"NDVI_AVH", SeriesIndex::NDVI_AVH},
		{"Net radtn", SeriesIndex::NET_RAD},
		{"VPD", SeriesIndex::VPD},
	};
	PPPG_SERIES_PARAM seriesParams[9];
	std::unordered_set<std::string> acquiredSeriesParams;
	//std::unordered_map<std::string, std::unordered_map<int, std::vector<PPPG_PARAM>>> seriesParams;

	bool haveTavg;
	bool haveVPD;
	bool haveNDVI;
	bool haveNetRad;

	GDALRasterImage* refGrid;
	bool finishedInput = false;

	bool getScalar(std::string value, PPPG_PARAM& param);
	bool getGrid(std::string value, PPPG_PARAM& param);
	double getValFromInputParam(std::string paramName, long cellIndex);
	double getValFromSeriesParam(int paramIndex, int year, int month, long cellIndex);
	bool openCheckGrid(std::string path, GDALRasterImage*& grid);
public:
	DataInput();
	~DataInput();
	bool tryAddInputParam(std::string pname, std::vector<std::string> value);
	bool tryAddSeriesParam(std::string name, std::vector<std::string> value, std::ifstream& paramFp, int& lineNo);
	bool inputFinished(bool modelMode3PGS);
	bool getInputParams(long cellIndex, InputParams& params);
	bool getSeriesParams(long cellIndex, int year, int month, SeriesParams& params);
	bool haveNetRadParam();
	GDALRasterImage* getRefGrid();
	void findRunPeriod(MYDate& minMY, MYDate& maxMY);

	bool haveSeedlingMass;
	bool haveMinASWTG;
	bool haveAgeDepFert;
};