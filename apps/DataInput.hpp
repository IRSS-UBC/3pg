#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include "GDALRasterImage.hpp"
#include "ParamStructs.hpp"

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
	double CO2;
	double fCalpha700;
	double fCg700;
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

typedef enum {
	FERTILITY = 0,
	MINASW = 1,
	IRRIGATION = 2,
} ManagementIndex;

class DataInput {
private:
	//data structures for dealing with input parameters
	std::unordered_map<std::string, std::string> inputParamNames = {
		{"foliage:stem partitioning ratio @ d=2 cm", "pfs2"},
		{"foliage:stem partitioning ratio @ d=20 cm", "pfs20"},
		{"constant in the stem mass v. diam. relationship", "stemconst"},
		{"power in the stem mass v. diam. relationship", "stempower"},
		{"maximum fraction of npp to roots", "prx"},
		{"minimum fraction of npp to roots", "prn"},
		{"minimum temperature for growth", "growthtmin"},
		{"optimum temperature for growth", "growthtopt"},
		{"maximum temperature for growth", "growthtmax"},
		{"days production lost per frost day", "kf"},
		{"maximum litterfall rate", "gammafx"},
		{"litterfall rate at t = 0", "gammaf0"},
		{"age at which litterfall rate has median value", "tgammaf"},
		{"average monthly root turnover rate", "rttover"},
		{"maximum canopy conductance", "maxcond"},
		{"defines stomatal response to vpd", "coeffcond"},
		{"canopy boundary layer conductance", "blcond"},
		{"value of 'm' when fr = 0", "m0"},
		{"value of 'fnutr' when fr = 0", "fn0"},
		{"power of (1-fr) in 'fnutr'", "fnn"},
		{"moisture ratio deficit for fq = 0.5", "swconst0"},
		{"power of moisture ratio deficit", "swpower0"},
		{"max. stem mass per tree @ 1000 trees/hectare", "wsx1000"},
		{"power in self-thinning rule", "thinpower"},
		{"fraction mean single-tree foliage biomass lost per dead tree", "mf"},
		{"fraction mean single-tree root biomass lost per dead tree", "mr"},
		{"fraction mean single-tree stem biomass lost per dead tree", "ms"},
		{"maximum stand age used in age modifier", "maxage"},
		{"power of relative age in function for fage", "nage"},
		{"relative age to give fage = 0.5", "rage"},
		{"specific leaf area at age 0", "sla0"},
		{"specific leaf area for mature leaves", "sla1"},
		{"age at which specific leaf area = (sla0+sla1)/2", "tsla"},
		{"extinction coefficient for absorption of par by canopy", "k"},
		{"age at canopy cover", "fullcanage"},
		{"canopy quantum efficiency", "alpha"},
		{"branch and bark fraction at age 0", "fracbb0"},
		{"branch and bark fraction for mature stands", "fracbb1"},
		{"age at which fracbb = (fracbb0+fracbb1)/2", "tbb"},
		{"ratio npp/gpp", "y"},
		{"basic density", "density"},
		{"intercept of net v. solar radiation relationship", "qa"},
		{"slope of net v. solar radiation relationship", "qb"},
		{"molecular weight of dry matter", "gdm_mol"},
		{"conversion of solar radiation to par", "molpar_mj"},
		{"lai for maximum canopy conductance", "laigcx"},
		{"maximum proportion of rainfall evaporated from canopy", "maxintcptn"},
		{"lai for maximum rainfall interception", "laimaxintcptn"},
		{"latitude", "lat"},
		{"fr", "frp"},
		{"fertility rating", "frp"},
		{"soil index", "soilindex"},
		{"soil class", "soilindex"},
		{"maximum asw", "maxasw"},
		{"minasw", "minaswp"},
		{"minimum asw", "minaswp"},
		{"initial age", "startage"},
		{"start age", "startage"},
		{"end year", "endyear"},
		{"start month", "startmonth"},
		{"seedling mass", "seedlingmass"},
		{"w foliage", "wfi"},
		{"w root", "wri"},
		{"w stem", "wsi"},
		{"stem no", "stemnoi"},
		{"initial soil water", "aswi"},
		{"minimum basic density - for young trees", "rhomin"},
		{"maximum basic density - for older trees", "rhomax"},
		{"age at which rho = (rhomin+rhomax)/2", "trho"},
		{"year planted", "yearplanted"},
		{"atmospheric co2", "co2"},
		{"assimialtion enhancement factor at 700 ppm", "fcalpha700"},
		{"canopy conductance enhancement factor at 700 ppm", "fcg700"},
	};
	std::unordered_set<std::string> allInputParams = {
		"pfs2",
		"pfs20",
		"stemconst",
		"stempower",
		"prx",
		"prn",
		"growthtmin",
		"growthtopt",
		"growthtmax",
		"kf",
		"gammafx",
		"gammaf0",
		"tgammaf",
		"rttover",
		"maxcond",
		"coeffcond",
		"blcond",
		"m0",
		"fn0",
		"fnn",
		"thinpower",
		"mf",
		"mr",
		"ms",
		"swconst0",
		"swpower0",
		"wsx1000",
		"maxage",
		"nage",
		"rage",
		"sla0",
		"sla1",
		"tsla",
		"k",
		"fullcanage",
		"alpha",
		"fracbb0",
		"fracbb1",
		"tbb",
		"y",
		"rhomin",
		"rhomax",
		"trho",
		"qa",
		"qb",
		"gdm_mol",
		"molpar_mj",
		"laigcx",
		"maxintcptn",
		"laimaxintcptn",
		"lat",
		"frp",
		"frstart",
		"frend",
		"frdec",
		"soilindex",
		"maxasw",
		"minaswp",
		"startage",
		"endyear",
		"startmonth",
		"yearplanted",
		"seedlingmass",
		"wfi",
		"wri",
		"wsi",
		"stemnoi",
		"aswi",
		"minaswtg",
		"ndvi_fpar_intercept",
		"ndvi_fpar_constant",
		"co2",
		"fcalpha700",
		"fcg700"
	};
	std::unordered_set<std::string> requiredInputParams3PG = {
		"pfs2", "pfs20", "stemconst", "stempower", "prx", "prn",
		"growthtmin", "growthtopt", "growthtmax",
		"kf",
		"gammafx", "gammaf0", "tgammaf", "rttover",
		"maxcond", "coeffcond", "blcond",
		"m0", "fn0", "fnn",
		"swconst0", "swpower0",
		"wsx1000",
		"maxage", "nage", "rage",
		"sla0", "sla1", "tsla", "k", "fullcanage",
		"alpha", "fracbb0", "fracbb1", "tbb",
		"y",
		"lat", "frp", "soilindex", "maxasw", "minaswp",
		"startage", "endyear",
		"stemnoi", "aswi", "yearplanted", 
		"qa", "qb",
		"gdm_mol", "molpar_mj",
		"laigcx", "maxintcptn",
		"startmonth",
		"laimaxintcptn",
		"thinpower", "mf", "mr", "ms",
		"co2", "fcalpha700", "fcg700"
	};
	std::unordered_set<std::string> requiredInputParams3PGS = {
		"growthtmin", "growthtopt", "growthtmax",
		"kf",
		"maxcond", "coeffcond", "blcond",
		"m0", "fn0",
		"swconst0", "swpower0",
		"sla1", "alpha",
		"y",
		"lat", "frp", "soilindex", "maxasw", "minaswp",
		"startage","endyear", "yearplanted",
		"ndvi_fpar_intercept", "ndvi_fpar_constant",
		"qa", "qb",
		"gdm_mol", "molpar_mj",
		"laigcx", "maxintcptn",
		"startmonth",
		"laimaxintcptn",
	};
	std::unordered_map<std::string, std::unique_ptr<PPPG_PARAM>> inputParams;

	//data structures for dealing with series parameters
	std::unordered_map<std::string, SeriesIndex> seriesParamNameMap = {
		{"tmax", SeriesIndex::TMAX},
		{"tmin", SeriesIndex::TMIN},
		{"tavg", SeriesIndex::TAVG},
		{"rain", SeriesIndex::RAIN},
		{"solar radtn", SeriesIndex::SOLAR_RAD},
		{"frost", SeriesIndex::FROST_DAYS},
		{"frost days", SeriesIndex::FROST_DAYS},
		{"ndvi_avh", SeriesIndex::NDVI_AVH},
		{"net radtn", SeriesIndex::NET_RAD},
		{"vpd", SeriesIndex::VPD},
	};
	PPPG_SERIES_PARAM seriesParams[9];
	std::unordered_set<std::string> acquiredSeriesParams;
	
	//data structures for dealing with output paramenters
	std::unordered_map<std::string, std::string> outputParamNames = {
		{"stocking density","stemno"},
		{"weight of foliage","wf" },
		{"weight of roots","wr"},
		{"weight of stems","ws"},
		{"total weight","totalw"},
		{"mean annual increment","mai"},
		{"average dbh","avdbh"},
		{"basal area","basarea"},
		{"stand volume","standvol"},
		{"transpiration","transp"},
		{"available soil water","asw"},
		{"gross primary production (tdm/ha)","gpp"},
		{"change in aboveground boimass (tdm/ha)","delwag"},
		{"accumulated aboveground biomass (tdm/ha)","cumwabv"},
		{"net primary production (tdm/ha)","npp"},
		{"leaf area index","lai"},
		{"frout","fr"},
	};
	std::unordered_set<std::string> allOutputParams = {
		"stemno",
		"wf",
		"wr",
		"ws",
		"totalw",
		"mai",
		"avdbh",
		"basarea",
		"standvol",
		"transp",
		"ctransp",
		"asw",
		"evaptransp",
		"cevaptransp",
		"laix",
		"agelaix",
		"gpp",
		"cgpp",
		"totallitter",
		"clitter",
		"delwag",
		"cumwabv",
		"npp",
		"cnpp",
		"lai",
		"clai",
		"fr",
		"physmod",
		"alphac",
		"fage",
		"fracbb",
		"wue",
		"cwue",
		"cvi",
		"ccvi",
		"fsw",
		"fvpd",
		"ft",
		"fnutr",
		"ffrost",
		"apar",
		"aparu",
		"maix",
		"agemaix"
	};
	std::unordered_set<std::string> only3PG = {
		"stemno",
		"wf",
		"wr",
		"ws",
		"totalw",
		"mai",
		"avdbh",
		"basarea",
		"standvol",
		"transp",
		"ctransp",
		"asw",
		"evaptransp",
		"cevaptransp",
		"laix",
		"agelaix",
		"gpp",
		"cgpp",
		"totallitter",
		"clitter" };
	std::unordered_set<std::string> only3PGS = {
		"delwag",
		"cumwabv" 
	};
	bool allow3PG = true;
	bool allow3PGS = true;
	std::unordered_map<std::string, PPPG_OP_VAR> outputParams;
	std::function<void(std::string)> log;

	//data structures for dealing with
	PPPG_MT_PARAM managementTables[3];

	bool haveTavg = false;
	bool haveVPD = false;
	bool haveNDVI = false;
	bool haveNetRad = false;

	RefGridProperties refGrid;
	bool finishedInput = false;

	std::string outPath;

	bool getScalar(std::string value, PPPG_PARAM* param);
	bool getGrid(std::string value, PPPG_PARAM* param);
	double getValFromInputParam(std::string paramName, long cellIndex);
	double getValFromSeriesParam(int paramIndex, int year, int month, long cellIndex);
	bool openCheckGrid(std::string path, std::unique_ptr<GDALRasterImage>& grid);

public:
	DataInput(const std::string& speciesFile, 
		const std::string& siteFile, 
		std::function<void(std::string)>& log
	);
	DataInput();//for no logger
	bool tryAddInputParam(std::string name, std::vector<std::string> value);
	bool tryAddSeriesParam(std::string name, std::vector<std::string> value, std::ifstream& paramFp, int& lineNo);
	bool tryAddOutputParam(std::string name, std::vector<std::string> value, int lineno);
	bool tryAddManagementParam(std::string name, std::ifstream& inFile, int& lineNo);
	bool inputFinished(bool modelModeSpatial);
	bool getInputParams(long cellIndex, InputParams& params);
	bool getSeriesParams(long cellIndex, int year, int month, SeriesParams& params);
	bool getManagementParam(ManagementIndex index, long cellIndex, int year, double& val);
	std::unordered_map<std::string, PPPG_OP_VAR> getOpVars();
	bool haveNetRadParam();
	RefGridProperties getRefGrid();

	bool haveSeedlingMass = false;
	bool haveMinASWTG = false;
	bool haveAgeDepFert = false;

	bool modelMode3PGS = false;
};