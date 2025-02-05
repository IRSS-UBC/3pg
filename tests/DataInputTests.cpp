#include <gtest/gtest.h>
#include <../apps/DataInput.hpp>
#include <fstream>
#include <vector>
#include <tuple>

#include <../apps/GDALRasterImage.hpp>

TEST(DataInputTests, correctParamNames) {
	DataInput* dataInput = new DataInput();

	//long versions according to legacy documentations (some lowercase)
	EXPECT_TRUE(dataInput->tryAddInputParam("Foliage:stem partitioning ratio @ D=2 cm", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("foliage:stem partitioning ratio @ d=20 cm", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Constant in the stem mass v. diam. relationship", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("power in the stem mass v. diam. relationship", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Maximum fraction of NPP to roots", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("minimum fraction of npp to roots", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Minimum temperature for growth", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("optimum temperature for growth", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Maximum temperature for growth", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("days production lost per frost day", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Moisture ratio deficit for fq = 0.5", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("power of moisture ratio deficit", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Value of 'm' when FR = 0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Value of 'fNutr' when FR = 0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Power of (1-FR) in 'fNutr'", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Maximum stand age used in age modifier", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Power of relative age in function for fAge", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Relative age to give fAge = 0.5", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Maximum litterfall rate", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Litterfall rate at t = 0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Age at which litterfall rate has median value", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Average monthly root turnover rate", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Maximum canopy conductance", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("LAI for maximum canopy conductance", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Defines stomatal response to VPD", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Canopy boundary layer conductance", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Max. stem mass per tree @ 1000 trees/hectare", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Power in self-thinning rule", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Fraction mean single-tree foliage biomass lost per dead tree", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Fraction mean single-tree root biomass lost per dead tree", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Fraction mean single-tree stem biomass lost per dead tree", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Specific leaf area at age 0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Specific leaf area for mature leaves", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Age at which specific leaf area = (SLA0+SLA1)/2", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Extinction coefficient for absorption of PAR by canopy", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Age at canopy cover", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Maximum proportion of rainfall evaporated from canopy", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("LAI for maximum rainfall interception", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Canopy quantum efficiency", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Ratio NPP/GPP", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Branch and bark fraction at age 0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Branch and bark fraction for mature stands", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Age at which fracBB = (fracBB0+fracBB1)/2", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Minimum basic density - for young trees", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Maximum basic density - for older trees", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Age at which rho = (rhoMin+rhoMax)/2", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Intercept of net v. solar radiation relationship", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Slope of net v. solar radiation relationship", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Molecular weight of dry matter", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Conversion of solar radiation to PAR", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("atmospheric co2 - start of run period", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("atmospheric co2 - end of run period", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("assimialtion enhancement factor at 700 ppm", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("canopy conductance enhancement factor at 700 ppm", { "1" }));

	//reset DataInput
	delete dataInput;
	dataInput = new DataInput();

	//short versions (some lowercase)
	EXPECT_TRUE(dataInput->tryAddInputParam("pfs2", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("pFS20", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("stemconst", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("StemPower", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("prx", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("pRn", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("growthtmin", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("growthTopt", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("growthtmax", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("kF", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("gammafx", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("gammaF0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("tgammaf", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Rttover", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("maxcond", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("CoeffCond", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("blcond", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("m0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("fN0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("fNn", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("thinPower", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("mf", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("mR", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("mS", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("SWconst0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("SWpower0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("wSx1000", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("MaxAge", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("nAge", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("rAge", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("SLA0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("SLA1", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("tSLA", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("k", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("fullCanAge", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("alpha", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("fracBB0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("fracBB1", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("tBB", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("y", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("rhoMin", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("rhoMax", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("tRho", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Qa", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Qb", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("gDM_mol", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("molPAR_MJ", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("LAIgcx", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("MaxIntcptn", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("LAImaxIntcptn", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Lat", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("FRp", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("FRstart", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("FRend", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("FRdec", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("soilIndex", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("MaxASW", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("MinASWp", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("StartAge", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("EndYear", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("StartMonth", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("yearPlanted", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("SeedlingMass", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("WFi", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("WRi", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("WSi", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("StemNoi", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("ASWi", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("MinASWTG", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("NDVI_FPAR_intercept", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("NDVI_FPAR_constant", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("CO2Start", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("CO2End", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("fCalpha700", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("fCg700", { "1" }));

	std::ifstream paramFp("");
	int lineNo;
	EXPECT_TRUE(dataInput->tryAddSeriesParam("Tmin", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_TRUE(dataInput->tryAddSeriesParam("Tmax", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_TRUE(dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_TRUE(dataInput->tryAddSeriesParam("Rain", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_TRUE(dataInput->tryAddSeriesParam("Solar radtn", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_TRUE(dataInput->tryAddSeriesParam("Frost", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_TRUE(dataInput->tryAddSeriesParam("Net radtn", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_TRUE(dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_TRUE(dataInput->tryAddSeriesParam("NDVI_AVH", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
 
	//clean up
	delete dataInput;
}

TEST(DataInputTests, incorrectParamNames) {
	DataInput* dataInput = new DataInput();

	//a smattering of incorrect input params which should return false
	EXPECT_FALSE(dataInput->tryAddInputParam("", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("	", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("1", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("3A0YUBtBlA", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Foliagestem partitioning ratio @ D=2 cm", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Foliagestem partitioning ratio @ D=20 cm", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Age at which specific leaf area = (SLA0+SLA1)2", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Powerin self-thinning rule", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Power  in self-thinning rule", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Power in selfthinning rule", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Power in self thinning rule", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Foliagestem partitioning ratio D=20 cm", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Foliagestem partitioning ratio @ D= cm", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Q", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("b", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Tmax", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("Rain", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("VPD", { "1" }));
	EXPECT_FALSE(dataInput->tryAddInputParam("C02", { "1" }));

	//incorrect series param names
	std::ifstream paramFp("");
	int lineNo;
	EXPECT_FALSE(dataInput->tryAddSeriesParam("min", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_FALSE(dataInput->tryAddSeriesParam("", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_FALSE(dataInput->tryAddSeriesParam("	", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_FALSE(dataInput->tryAddSeriesParam("1", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_FALSE(dataInput->tryAddSeriesParam("Solarradtn", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_FALSE(dataInput->tryAddSeriesParam("3A0YUBtBlA", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_FALSE(dataInput->tryAddSeriesParam("Foliage:stem partitioning ratio @ D=2 cm", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_FALSE(dataInput->tryAddSeriesParam("Foliage:stem partitioning ratio @ D=20 cm", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
	EXPECT_FALSE(dataInput->tryAddSeriesParam("Constant in the stem mass v. diam. relationship", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo));
 
	delete dataInput;
}

TEST(DataInputTests, scalarValues) {
	DataInput* dataInput = new DataInput();

	//ensure each input goes to the correct parameter
	dataInput->tryAddInputParam("pFS2", { "1" });
	dataInput->tryAddInputParam("pFS20", { "2" });
	dataInput->tryAddInputParam("StemConst", { "3" });
	dataInput->tryAddInputParam("StemPower", { "4" });
	dataInput->tryAddInputParam("pRx", { "5" });
	dataInput->tryAddInputParam("pRn", { "6" });
	dataInput->tryAddInputParam("growthTmin", { "7" });
	dataInput->tryAddInputParam("growthTopt", { "8" });
	dataInput->tryAddInputParam("growthTmax", { "9" });
	dataInput->tryAddInputParam("kF", { "10" });
	dataInput->tryAddInputParam("gammaFx", { "11" });
	dataInput->tryAddInputParam("gammaF0", { "21" });
	dataInput->tryAddInputParam("tgammaF", { "13" });
	dataInput->tryAddInputParam("Rttover", { "14" });
	dataInput->tryAddInputParam("MaxCond", { "15" });
	dataInput->tryAddInputParam("CoeffCond", { "16" });
	dataInput->tryAddInputParam("BLcond", { "17" });
	dataInput->tryAddInputParam("m0", { "18" });
	dataInput->tryAddInputParam("fN0", { "19" });
	dataInput->tryAddInputParam("fNn", { "20" });
	dataInput->tryAddInputParam("thinPower", { "21" });
	dataInput->tryAddInputParam("mF", { "22" });
	dataInput->tryAddInputParam("mR", { "23" });
	dataInput->tryAddInputParam("mS", { "24" });
	dataInput->tryAddInputParam("SWconst0", { "25" });
	dataInput->tryAddInputParam("SWpower0", { "26" });
	dataInput->tryAddInputParam("wSx1000", { "27" });
	dataInput->tryAddInputParam("MaxAge", { "28" });
	dataInput->tryAddInputParam("nAge", { "29" });
	dataInput->tryAddInputParam("rAge", { "30" });
	dataInput->tryAddInputParam("SLA0", { "31" });
	dataInput->tryAddInputParam("SLA1", { "32" });
	dataInput->tryAddInputParam("tSLA", { "33" });
	dataInput->tryAddInputParam("k", { "34" });
	dataInput->tryAddInputParam("fullCanAge", { "35" });
	dataInput->tryAddInputParam("alpha", { "36" });
	dataInput->tryAddInputParam("fracBB0", { "37" });
	dataInput->tryAddInputParam("fracBB1", { "38" });
	dataInput->tryAddInputParam("tBB", { "39" });
	dataInput->tryAddInputParam("y", { "40" });
	dataInput->tryAddInputParam("rhoMin", { "41" });
	dataInput->tryAddInputParam("rhoMax", { "42" });
	dataInput->tryAddInputParam("tRho", { "43" });
	dataInput->tryAddInputParam("Qa", { "44" });
	dataInput->tryAddInputParam("Qb", { "45" });
	dataInput->tryAddInputParam("gDM_mol", { "46" });
	dataInput->tryAddInputParam("molPAR_MJ", { "47" });
	dataInput->tryAddInputParam("LAIgcx", { "48" });
	dataInput->tryAddInputParam("MaxIntcptn", { "49" });
	dataInput->tryAddInputParam("LAImaxIntcptn", { "50" });
	dataInput->tryAddInputParam("Lat", { "51" });
	dataInput->tryAddInputParam("FRp", { "52" });
	dataInput->tryAddInputParam("FRstart", { "53" });
	dataInput->tryAddInputParam("FRend", { "54" });
	dataInput->tryAddInputParam("FRdec", { "55" });
	dataInput->tryAddInputParam("soilIndex", { "56" });
	dataInput->tryAddInputParam("MaxASW", { "57" });
	dataInput->tryAddInputParam("MinASWp", { "58" });
	dataInput->tryAddInputParam("StartAge", { "59" });
	dataInput->tryAddInputParam("EndYear", { "125" }); //121 instead of 60 because it must be larger than startAge(59) + yearPlanted(62)
	dataInput->tryAddInputParam("StartMonth", { "61" });
	dataInput->tryAddInputParam("yearPlanted", { "62" });
	dataInput->tryAddInputParam("CO2Start", { "63" });
	dataInput->tryAddInputParam("CO2End", { "64" });
	dataInput->tryAddInputParam("fCalpha700", { "65" });
	dataInput->tryAddInputParam("fCg700", { "66" });

	//and ensure doubles are fully handled
	dataInput->tryAddInputParam("SeedlingMass", { ".0000000000001"});
	dataInput->tryAddInputParam("WFi", { ".00000000000001" });
	dataInput->tryAddInputParam("WRi", { ".000000000000001" });
	dataInput->tryAddInputParam("WSi", { ".0000000000000001" });
	dataInput->tryAddInputParam("StemNoi", { ".00000000000000001" });
	dataInput->tryAddInputParam("ASWi", { std::to_string(std::numeric_limits<double>::min() + 1) });
	dataInput->tryAddInputParam("MinASWTG", { std::to_string(std::numeric_limits<double>::max() - 1) });
	dataInput->tryAddInputParam("NDVI_FPAR_intercept", { std::to_string(std::numeric_limits<double>::min() + .01) });
	dataInput->tryAddInputParam("NDVI_FPAR_constant", { std::to_string(std::numeric_limits<double>::max() - .01) });

	//add series parameters
	std::ifstream paramFp("");
	int lineNo;
	dataInput->tryAddSeriesParam("Tavg", { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("Rain", { "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("Solar radtn", { "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("Frost", { "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("NDVI_AVH", { "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("Net radtn", { "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("VPD", { "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84" }, paramFp, lineNo);

	//indicate input finished
	dataInput->inputFinished(true);

	//get the input params
	InputParams params;
	dataInput->getInputParams(0, params);

	//check the input params
	EXPECT_EQ(params.pFS2, (double)1);
	EXPECT_EQ(params.pFS20, (double)2);
	EXPECT_EQ(params.StemConst, (double)3);
	EXPECT_EQ(params.StemPower, (double)4);
	EXPECT_EQ(params.pRx, (double)5);
	EXPECT_EQ(params.pRn, (double)6);
	EXPECT_EQ(params.growthTmin, (double)7);
	EXPECT_EQ(params.growthTopt, (double)8);
	EXPECT_EQ(params.growthTmax, (double)9);
	EXPECT_EQ(params.kF, (double)10);
	EXPECT_EQ(params.gammaFx, (double)11);
	EXPECT_EQ(params.gammaF0, (double)21);
	EXPECT_EQ(params.tgammaF, (double)13);
	EXPECT_EQ(params.Rttover, (double)14);
	EXPECT_EQ(params.MaxCond, (double)15);
	EXPECT_EQ(params.CoeffCond, (double)16);
	EXPECT_EQ(params.BLcond, (double)17);
	EXPECT_EQ(params.m0, (double)18);
	EXPECT_EQ(params.fN0, (double)19);
	EXPECT_EQ(params.fNn, (double)20);
	EXPECT_EQ(params.thinPower, (double)21);
	EXPECT_EQ(params.mF, (double)22);
	EXPECT_EQ(params.mR, (double)23);
	EXPECT_EQ(params.mS, (double)24);
	EXPECT_EQ(params.SWconst0, (double)25);
	EXPECT_EQ(params.SWpower0, (double)26);
	EXPECT_EQ(params.wSx1000, (double)27);
	EXPECT_EQ(params.MaxAge, (double)28);
	EXPECT_EQ(params.nAge, (double)29);
	EXPECT_EQ(params.rAge, (double)30);
	EXPECT_EQ(params.SLA0, (double)31);
	EXPECT_EQ(params.SLA1, (double)32);
	EXPECT_EQ(params.tSLA, (double)33);
	EXPECT_EQ(params.k, (double)34);
	EXPECT_EQ(params.fullCanAge, (double)35);
	EXPECT_EQ(params.alpha, (double)36);
	EXPECT_EQ(params.fracBB0, (double)37);
	EXPECT_EQ(params.fracBB1, (double)38);
	EXPECT_EQ(params.tBB, (double)39);
	EXPECT_EQ(params.y, (double)40);
	EXPECT_EQ(params.rhoMin, (double)41);
	EXPECT_EQ(params.rhoMax, (double)42);
	EXPECT_EQ(params.tRho, (double)43);
	EXPECT_EQ(params.Qa, (double)44);
	EXPECT_EQ(params.Qb, (double)45);
	EXPECT_EQ(params.gDM_mol, (double)46);
	EXPECT_EQ(params.molPAR_MJ, (double)47);
	EXPECT_EQ(params.LAIgcx, (double)48);
	EXPECT_EQ(params.MaxIntcptn, (double)49);
	EXPECT_EQ(params.LAImaxIntcptn, (double)50);
	EXPECT_EQ(params.Lat, (double)51);
	EXPECT_EQ(params.FRp, (double)52);
	EXPECT_EQ(params.FRstart, (double)53);
	EXPECT_EQ(params.FRend, (double)54);
	EXPECT_EQ(params.FRdec, (double)55);
	EXPECT_EQ(params.soilIndex, (double)56);
	EXPECT_EQ(params.MaxASW, (double)57);
	EXPECT_EQ(params.MinASWp, (double)58);
	EXPECT_EQ(params.StartAge, (double)59);
	EXPECT_EQ(params.EndYear, (double)125);
	EXPECT_EQ(params.StartMonth, (double)61);
	EXPECT_EQ(params.yearPlanted, (double)62);
	EXPECT_EQ(params.CO2Start, (double)63);
	EXPECT_EQ(params.CO2End, (double)64);
	EXPECT_EQ(params.fCalpha700, (double)65);
	EXPECT_EQ(params.fCg700, (double)66);
	EXPECT_EQ(params.SeedlingMass, (double).0000000000001);
	EXPECT_EQ(params.WFi, (double).00000000000001);
	EXPECT_EQ(params.WRi, (double).000000000000001);
	EXPECT_EQ(params.WSi, (double).0000000000000001);
	EXPECT_EQ(params.StemNoi, (double).00000000000000001);
	EXPECT_EQ(params.ASWi, std::numeric_limits<double>::min() + 1);
	EXPECT_EQ(params.MinASWTG, std::numeric_limits<double>::max() - 1);
	EXPECT_EQ(params.NDVI_FPAR_intercept, std::numeric_limits<double>::min() + .01);
	EXPECT_EQ(params.NDVI_FPAR_constant, std::numeric_limits<double>::max() - .01);

	SeriesParams sParams;
	//series params January
	dataInput->getSeriesParams(0, 0, 1, sParams);
	EXPECT_EQ(sParams.Tavg, (double)1);
	EXPECT_EQ(sParams.Rain, (double)13);
	EXPECT_EQ(sParams.SolarRad, (double)25);
	EXPECT_EQ(sParams.FrostDays, (double)37);
	EXPECT_EQ(sParams.NDVI_AVH, (double)49);
	EXPECT_EQ(sParams.NetRad, (double)61);
	EXPECT_EQ(sParams.VPD, (double)73);

	//series params February
	dataInput->getSeriesParams(0, 0, 2, sParams);
	EXPECT_EQ(sParams.Tavg, (double)2);
	EXPECT_EQ(sParams.Rain, (double)14);
	EXPECT_EQ(sParams.SolarRad, (double)26);
	EXPECT_EQ(sParams.FrostDays, (double)38);
	EXPECT_EQ(sParams.NDVI_AVH, (double)50);
	EXPECT_EQ(sParams.NetRad, (double)62);
	EXPECT_EQ(sParams.VPD, (double)74);

	//series params March
	dataInput->getSeriesParams(0, 0, 3, sParams);
	EXPECT_EQ(sParams.Tavg, (double)3);
	EXPECT_EQ(sParams.Rain, (double)15);
	EXPECT_EQ(sParams.SolarRad, (double)27);
	EXPECT_EQ(sParams.FrostDays, (double)39);
	EXPECT_EQ(sParams.NDVI_AVH, (double)51);
	EXPECT_EQ(sParams.NetRad, (double)63);
	EXPECT_EQ(sParams.VPD, (double)75);

	//series params April
	dataInput->getSeriesParams(0, 0, 4, sParams);
	EXPECT_EQ(sParams.Tavg, (double)4);
	EXPECT_EQ(sParams.Rain, (double)16);
	EXPECT_EQ(sParams.SolarRad, (double)28);
	EXPECT_EQ(sParams.FrostDays, (double)40);
	EXPECT_EQ(sParams.NDVI_AVH, (double)52);
	EXPECT_EQ(sParams.NetRad, (double)64);
	EXPECT_EQ(sParams.VPD, (double)76);

	//series params May
	dataInput->getSeriesParams(0, 0, 5, sParams);
	EXPECT_EQ(sParams.Tavg, (double)5);
	EXPECT_EQ(sParams.Rain, (double)17);
	EXPECT_EQ(sParams.SolarRad, (double)29);
	EXPECT_EQ(sParams.FrostDays, (double)41);
	EXPECT_EQ(sParams.NDVI_AVH, (double)53);
	EXPECT_EQ(sParams.NetRad, (double)65);
	EXPECT_EQ(sParams.VPD, (double)77);

	//series params June
	dataInput->getSeriesParams(0, 0, 6, sParams);
	EXPECT_EQ(sParams.Tavg, (double)6);
	EXPECT_EQ(sParams.Rain, (double)18);
	EXPECT_EQ(sParams.SolarRad, (double)30);
	EXPECT_EQ(sParams.FrostDays, (double)42);
	EXPECT_EQ(sParams.NDVI_AVH, (double)54);
	EXPECT_EQ(sParams.NetRad, (double)66);
	EXPECT_EQ(sParams.VPD, (double)78);

	//series params July
	dataInput->getSeriesParams(0, 0, 7, sParams);
	EXPECT_EQ(sParams.Tavg, (double)7);
	EXPECT_EQ(sParams.Rain, (double)19);
	EXPECT_EQ(sParams.SolarRad, (double)31);
	EXPECT_EQ(sParams.FrostDays, (double)43);
	EXPECT_EQ(sParams.NDVI_AVH, (double)55);
	EXPECT_EQ(sParams.NetRad, (double)67);
	EXPECT_EQ(sParams.VPD, (double)79);

	//series params August
	dataInput->getSeriesParams(0, 0, 8, sParams);
	EXPECT_EQ(sParams.Tavg, (double)8);
	EXPECT_EQ(sParams.Rain, (double)20);
	EXPECT_EQ(sParams.SolarRad, (double)32);
	EXPECT_EQ(sParams.FrostDays, (double)44);
	EXPECT_EQ(sParams.NDVI_AVH, (double)56);
	EXPECT_EQ(sParams.NetRad, (double)68);
	EXPECT_EQ(sParams.VPD, (double)80);

	//series params September
	dataInput->getSeriesParams(0, 0, 9, sParams);
	EXPECT_EQ(sParams.Tavg, (double)9);
	EXPECT_EQ(sParams.Rain, (double)21);
	EXPECT_EQ(sParams.SolarRad, (double)33);
	EXPECT_EQ(sParams.FrostDays, (double)45);
	EXPECT_EQ(sParams.NDVI_AVH, (double)57);
	EXPECT_EQ(sParams.NetRad, (double)69);
	EXPECT_EQ(sParams.VPD, (double)81);

	//series params October
	dataInput->getSeriesParams(0, 0, 10, sParams);
	EXPECT_EQ(sParams.Tavg, (double)10);
	EXPECT_EQ(sParams.Rain, (double)22);
	EXPECT_EQ(sParams.SolarRad, (double)34);
	EXPECT_EQ(sParams.FrostDays, (double)46);
	EXPECT_EQ(sParams.NDVI_AVH, (double)58);
	EXPECT_EQ(sParams.NetRad, (double)70);
	EXPECT_EQ(sParams.VPD, (double)82);

	//series params November
	dataInput->getSeriesParams(0, 0, 11, sParams);
	EXPECT_EQ(sParams.Tavg, (double)11);
	EXPECT_EQ(sParams.Rain, (double)23);
	EXPECT_EQ(sParams.SolarRad, (double)35);
	EXPECT_EQ(sParams.FrostDays, (double)47);
	EXPECT_EQ(sParams.NDVI_AVH, (double)59);
	EXPECT_EQ(sParams.NetRad, (double)71);
	EXPECT_EQ(sParams.VPD, (double)83);

	//series params December
	dataInput->getSeriesParams(0, 0, 12, sParams);
	EXPECT_EQ(sParams.Tavg, (double)12);
	EXPECT_EQ(sParams.Rain, (double)24);
	EXPECT_EQ(sParams.SolarRad, (double)36);
	EXPECT_EQ(sParams.FrostDays, (double)48);
	EXPECT_EQ(sParams.NDVI_AVH, (double)60);
	EXPECT_EQ(sParams.NetRad, (double)72);
	EXPECT_EQ(sParams.VPD, (double)84);
 
	//clean up
	delete dataInput;
}


std::string frost1 = "test_files/DataInputTests/NFFD01.tif";
std::string frost2 = "test_files/DataInputTests/NFFD02.tif";
std::string frost3 = "test_files/DataInputTests/NFFD03.tif";
std::string frost4 = "test_files/DataInputTests/NFFD04.tif";
std::string frost5 = "test_files/DataInputTests/NFFD05.tif";
std::string frost6 = "test_files/DataInputTests/NFFD06.tif";
std::string frost7 = "test_files/DataInputTests/NFFD07.tif";
std::string frost8 = "test_files/DataInputTests/NFFD08.tif";
std::string frost9 = "test_files/DataInputTests/NFFD09.tif";
std::string frost10 = "test_files/DataInputTests/NFFD10.tif";
std::string frost11 = "test_files/DataInputTests/NFFD11.tif";
std::string frost12 = "test_files/DataInputTests/NFFD12.tif";

//test grid parameters
TEST(DataInputTests, gridValues) {

	std::string testFile1 = "test_files/DataInputTests/test1.tif";
	std::string testFile2 = "test_files/DataInputTests/test2.tif";
	std::string testFile3 = "test_files/DataInputTests/test3.tif";
 
	std::unique_ptr<GDALRasterImage> testTif1 = std::make_unique<GDALRasterImage>(testFile1);
	std::unique_ptr<GDALRasterImage> testTif2 = std::make_unique<GDALRasterImage>(testFile2);
	std::unique_ptr<GDALRasterImage> testTif3 = std::make_unique<GDALRasterImage>(testFile3);
	std::unique_ptr<GDALRasterImage> frostTif1 = std::make_unique<GDALRasterImage>(frost1);
	std::unique_ptr<GDALRasterImage> frostTif2 = std::make_unique<GDALRasterImage>(frost2);
	std::unique_ptr<GDALRasterImage> frostTif3 = std::make_unique<GDALRasterImage>(frost3);
	std::unique_ptr<GDALRasterImage> frostTif4 = std::make_unique<GDALRasterImage>(frost4);
	std::unique_ptr<GDALRasterImage> frostTif5 = std::make_unique<GDALRasterImage>(frost5);
	std::unique_ptr<GDALRasterImage> frostTif6 = std::make_unique<GDALRasterImage>(frost6);
	std::unique_ptr<GDALRasterImage> frostTif7 = std::make_unique<GDALRasterImage>(frost7);
	std::unique_ptr<GDALRasterImage> frostTif8 = std::make_unique<GDALRasterImage>(frost8);
	std::unique_ptr<GDALRasterImage> frostTif9 = std::make_unique<GDALRasterImage>(frost9);
	std::unique_ptr<GDALRasterImage> frostTif10 = std::make_unique<GDALRasterImage>(frost10);
	std::unique_ptr<GDALRasterImage> frostTif11 = std::make_unique<GDALRasterImage>(frost11);
	std::unique_ptr<GDALRasterImage> frostTif12 = std::make_unique<GDALRasterImage>(frost12);

	DataInput* dataInput = new DataInput();
	std::ifstream paramFp("");
	int lineNo;
	
	//ensure we are able to add grids correctly
	ASSERT_TRUE(dataInput->tryAddInputParam("growthTmin", { testFile1 })); 
	ASSERT_TRUE(dataInput->tryAddInputParam("growthTopt", { testFile2 })); 
	ASSERT_TRUE(dataInput->tryAddInputParam("growthTmax", { testFile3 }));
	ASSERT_TRUE(dataInput->tryAddSeriesParam("Frost", {frost1,frost2,frost3,frost4,frost5,frost6,frost7,frost8,frost9,frost10,frost11,frost12,}, paramFp, lineNo));

	//these are here because we would otherwise not be able to get the input params grids
	dataInput->tryAddInputParam("kF", { "1" });
	dataInput->tryAddInputParam("MaxCond", { "1" });
	dataInput->tryAddInputParam("CoeffCond", { "1" }); 
	dataInput->tryAddInputParam("BLcond", { "1" });
	dataInput->tryAddInputParam("m0", { "1" }); 
	dataInput->tryAddInputParam("fN0", { "1" });
	dataInput->tryAddInputParam("SWconst0", { "1" }); 
	dataInput->tryAddInputParam("SWpower0", { "1" });
	dataInput->tryAddInputParam("SLA1", { "1" }); 
	dataInput->tryAddInputParam("alpha", { "1" });
	dataInput->tryAddInputParam("y", { "1" });
	dataInput->tryAddInputParam("Lat", { "1" }); 
	dataInput->tryAddInputParam("FRp", { "1" }); 
	dataInput->tryAddInputParam("soilIndex", { "1" }); 
	dataInput->tryAddInputParam("MaxASW", { "1" }); 
	dataInput->tryAddInputParam("MinASWp", { "1" });
	dataInput->tryAddInputParam("StartAge", { "1" }); 
	dataInput->tryAddInputParam("EndYear", { "3" }); //3 because EndYear must be greater than StartAge(1) + yearPlanted(1)
	dataInput->tryAddInputParam("yearPlanted", { "1" });
	dataInput->tryAddInputParam("NDVI_FPAR_intercept", { "1" }); 
	dataInput->tryAddInputParam("NDVI_FPAR_constant", { "1" });
	dataInput->tryAddInputParam("Qa", { "1" }); 
	dataInput->tryAddInputParam("Qb", { "1" });
	dataInput->tryAddInputParam("gDM_mol", { "1" }); 
	dataInput->tryAddInputParam("molPAR_MJ", { "1" });
	dataInput->tryAddInputParam("LAIgcx", { "1" }); 
	dataInput->tryAddInputParam("MaxIntcptn", { "1" });
	dataInput->tryAddInputParam("StartMonth", { "1" });
	dataInput->tryAddInputParam("LAImaxIntcptn", { "1" });
	dataInput->tryAddInputParam("SeedlingMass", { "1" });
	dataInput->tryAddInputParam("CO2Start", { "1" });
	dataInput->tryAddInputParam("CO2End", { "1" });
	dataInput->tryAddInputParam("fCalpha700", { "1" });
	dataInput->tryAddInputParam("fCg700", { "1" });
	dataInput->tryAddSeriesParam("Tavg", { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("Rain", { "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("Solar radtn", { "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("NDVI_AVH", { "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("Net radtn", { "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("VPD", { "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84" }, paramFp, lineNo);

	dataInput->inputFinished(true);

	InputParams params;
	SeriesParams sParams;

	for (int i = 0; i < testTif1->nRows; i++) {
		for (int j = 0; j < testTif1->nCols; j++) {
			long cellIndex = (long)i * (long)testTif1->nCols + (long)j;
			//test input parameter grids
			if (std::isnan(testTif1->GetVal(cellIndex))) {
				ASSERT_FALSE(dataInput->getInputParams(cellIndex, params));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 1, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 2, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 3, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 4, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 5, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 6, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 7, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 8, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 9, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 10, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 11, sParams));
				ASSERT_FALSE(dataInput->getSeriesParams(cellIndex, 0, 12, sParams));
				continue;
			}

			dataInput->getInputParams(cellIndex, params);

			ASSERT_EQ(params.growthTmin, testTif1->GetVal(cellIndex));
			ASSERT_EQ(params.growthTopt, testTif2->GetVal(cellIndex));
			ASSERT_EQ(params.growthTmax, testTif3->GetVal(cellIndex));
				
			//test series parameter January
			dataInput->getSeriesParams(cellIndex, 0, 1, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif1->GetVal(cellIndex));

			//test series parameter February
			dataInput->getSeriesParams(cellIndex, 0, 2, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif2->GetVal(cellIndex));

			//test series parameter March
			dataInput->getSeriesParams(cellIndex, 0, 3, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif3->GetVal(cellIndex));

			//test series parameter April
			dataInput->getSeriesParams(cellIndex, 0, 4, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif4->GetVal(cellIndex));

			//test series parameter May
			dataInput->getSeriesParams(cellIndex, 0, 5, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif5->GetVal(cellIndex));

			//test series parameter June
			dataInput->getSeriesParams(cellIndex, 0, 6, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif6->GetVal(cellIndex));

			//test series parameter July
			dataInput->getSeriesParams(cellIndex, 0, 7, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif7->GetVal(cellIndex));

			//test series parameter August
			dataInput->getSeriesParams(cellIndex, 0, 8, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif8->GetVal(cellIndex));

			//test series parameter September
			dataInput->getSeriesParams(cellIndex, 0, 9, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif9->GetVal(cellIndex));

			//test series parameter October
			dataInput->getSeriesParams(cellIndex, 0, 10, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif10->GetVal(cellIndex));

			//test series parameter November
			dataInput->getSeriesParams(cellIndex, 0, 11, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif11->GetVal(cellIndex));

			//test series parameter December
			dataInput->getSeriesParams(cellIndex, 0, 12, sParams);
			ASSERT_EQ(sParams.FrostDays, frostTif12->GetVal(cellIndex));
		}
	}
 
	//clean up
	delete dataInput;
}

//test bad grids
TEST(DataInputTests, badGrids) {
	std::string testFile = "test_files/DataInputTests/test1.tif";
	std::string differentLocationFile = "test_files/DataInputTests/differentLocation.tif";
	std::string differentDimensionsFile = "test_files/DataInputTests/differentDimensions.tif";

	std::unique_ptr<GDALRasterImage> testTif = std::make_unique<GDALRasterImage>(testFile);
	std::unique_ptr<GDALRasterImage> difLocationTif = std::make_unique<GDALRasterImage>(differentLocationFile);
	std::unique_ptr<GDALRasterImage> difDimensionsTif = std::make_unique<GDALRasterImage>(differentDimensionsFile);

	DataInput* dataInput = new DataInput();
	std::ifstream paramFp("");
	int lineNo;

	ASSERT_TRUE(dataInput->tryAddInputParam("kF", { testFile }));
	EXPECT_EXIT(dataInput->tryAddInputParam("MaxCond", { differentLocationFile }), testing::ExitedWithCode(1), "");
	EXPECT_EXIT(dataInput->tryAddInputParam("CoeffCond", { differentDimensionsFile }), testing::ExitedWithCode(1), "");
	EXPECT_EXIT(
		dataInput->tryAddSeriesParam("Frost", {"1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", differentLocationFile}, paramFp, lineNo), 
		testing::ExitedWithCode(1), "");
	EXPECT_EXIT(
		dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", differentDimensionsFile }, paramFp, lineNo), 
		testing::ExitedWithCode(1), "");

	//clean up
	delete dataInput;
}

//test required parameters
std::vector<std::string> requiredInputParams3PG = {
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
	"Lat", 
	"FRp", 
	"soilIndex", 
	"MaxASW", 
	"MinASWp",
	"StartAge", 
	"EndYear",
	"StemNoi", 
	"ASWi", 
	"yearPlanted",
	"Qa", 
	"Qb",
	"gDM_mol", 
	"molPAR_MJ",
	"LAIgcx", 
	"MaxIntcptn",
	"StartMonth",
	"LAImaxIntcptn",
	"thinPower", 
	"mF", 
	"mR", 
	"mS", 
	"CO2Start",
	"CO2End",
	"fCalpha700",
	"fCg700",
};
std::vector<std::string> requiredSeriesParams3PG = {
	"Rain",
	"Solar radtn", 
	"Frost"
};
std::vector<std::string> requiredInputParams3PGS = {
	"growthTmin", 
	"growthTopt", 
	"growthTmax",
	"kF",
	"MaxCond", 
	"CoeffCond", 
	"BLcond",
	"m0", 
	"fN0",
	"SWconst0", 
	"SWpower0",
	"SLA1", 
	"alpha",
	"y",
	"Lat", 
	"FRp", 
	"soilIndex",
	"MaxASW", 
	"MinASWp",
	"StartAge",
	"EndYear",
	"yearPlanted",
	"NDVI_FPAR_intercept", 
	"NDVI_FPAR_constant",
	"Qa", 
	"Qb",
	"gDM_mol", 
	"molPAR_MJ",
	"LAIgcx", 
	"MaxIntcptn",
	"StartMonth",
	"LAImaxIntcptn", 
};
std::vector<std::string> requiredSeriesParams3PGS = {
	"Rain",
	"Solar radtn",
	"Frost",
	"NDVI_AVH"
};

TEST(DataInputTests, requiredParameters3PG) {
	//check that we fail when missing each required 3PG input parameter
	for (int exclude = 0; exclude < requiredInputParams3PG.size(); exclude++) {
		DataInput* dataInput = new DataInput();
		std::ifstream paramFp("");
		int lineNo;

		for (int i = 0; i < requiredInputParams3PG.size(); i++) {
			if (i == exclude) {
				continue;
			}

			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "1" });
		}
		dataInput->tryAddInputParam("SeedlingMass", { "1" });

		for (int i = 0; i < requiredSeriesParams3PG.size(); i++) {
			dataInput->tryAddSeriesParam(requiredSeriesParams3PG[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
		}
		dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
		dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);

		ASSERT_FALSE(dataInput->inputFinished(false));

		//clean up here because dataInput goes out of scope when for loop finishes an iteration
		delete dataInput;
	}

	//check that we fail when missing each required 3PG series parameter
	for (int exclude = 0; exclude < requiredSeriesParams3PG.size(); exclude++) {
		DataInput* dataInput = new DataInput();
		std::ifstream paramFp("");
		int lineNo;
	
		for (int i = 0; i < requiredInputParams3PG.size(); i++) {
			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "1" });
		}
		dataInput->tryAddInputParam("SeedlingMass", { "1" });
	
		for (int i = 0; i < requiredSeriesParams3PG.size(); i++) {
			if (i == exclude) {
				continue;
			}
	
			dataInput->tryAddSeriesParam(requiredSeriesParams3PG[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
		}
		dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
		dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	
		ASSERT_FALSE(dataInput->inputFinished(false));

		//clean up here because dataInput goes out of scope when for loop finishes an iteration
		delete dataInput;
	}

	//check that we pass with all required params
	DataInput* dataInput = new DataInput();
	std::ifstream paramFp("");
	int lineNo;

	for (int i = 0; i < requiredInputParams3PG.size(); i++) {
		if (requiredInputParams3PG[i] == "EndYear") {
			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "3" });
		}
		else {
			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "1" });
		}
	}
	for (int i = 0; i < requiredSeriesParams3PG.size(); i++) {
		dataInput->tryAddSeriesParam(requiredSeriesParams3PG[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	}
	dataInput->tryAddInputParam("SeedlingMass", { "1" });
	dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	ASSERT_TRUE(dataInput->inputFinished(false));

	//clean up
	delete dataInput;
}

TEST(DataInputTests, requiredParameters3PGS) {
	//check that we fail when missing each required 3PGS input parametergit 
	for (int exclude = 0; exclude < requiredInputParams3PGS.size(); exclude++) {
		DataInput* dataInput = new DataInput();
		std::ifstream paramFp("");
		int lineNo;

		for (int i = 0; i < requiredInputParams3PGS.size(); i++) {
			if (i == exclude) {
				continue;
			}

			if (requiredInputParams3PGS[i] == "EndYear") {
				dataInput->tryAddInputParam(requiredInputParams3PGS[i], { "3" });
			}
			else {
				dataInput->tryAddInputParam(requiredInputParams3PGS[i], { "1" });
			}
		}
		dataInput->tryAddInputParam("SeedlingMass", { "1" });

		for (int i = 0; i < requiredSeriesParams3PGS.size(); i++) {
			dataInput->tryAddSeriesParam(requiredSeriesParams3PGS[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
		}
		dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
		dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);

		ASSERT_FALSE(dataInput->inputFinished(true));

		//clean up here because dataInput goes out of scope when for loop finishes an iteration
		delete dataInput;
	}

	//check that we fail when missing each required 3PG series parameter
	for (int exclude = 0; exclude < requiredSeriesParams3PGS.size(); exclude++) {
		DataInput* dataInput = new DataInput();
		std::ifstream paramFp("");
		int lineNo;

		for (int i = 0; i < requiredInputParams3PGS.size(); i++) {
			if (requiredInputParams3PGS[i] == "EndYear") {
				dataInput->tryAddInputParam(requiredInputParams3PGS[i], { "3" });
			}
			else {
				dataInput->tryAddInputParam(requiredInputParams3PGS[i], { "1" });
			}
		}
		dataInput->tryAddInputParam("SeedlingMass", { "1" });

		for (int i = 0; i < requiredSeriesParams3PG.size(); i++) {
			if (i == exclude) {
				continue;
			}

			dataInput->tryAddSeriesParam(requiredSeriesParams3PGS[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
		}
		dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
		dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);

		ASSERT_FALSE(dataInput->inputFinished(true));

		//clean up here because dataInput goes out of scope when for loop finishes an iteration
		delete dataInput;
	}

	//check that we pass with all required params
	DataInput* dataInput = new DataInput();
	std::ifstream paramFp("");
	int lineNo;

	for (int i = 0; i < requiredInputParams3PGS.size(); i++) {
		if (requiredInputParams3PGS[i] == "EndYear") {
			dataInput->tryAddInputParam(requiredInputParams3PGS[i], { "3" });
		}
		else {
			dataInput->tryAddInputParam(requiredInputParams3PGS[i], { "1" });
		}
	}
	for (int i = 0; i < requiredSeriesParams3PGS.size(); i++) {
		dataInput->tryAddSeriesParam(requiredSeriesParams3PGS[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	}
	dataInput->tryAddInputParam("SeedlingMass", { "1" });
	dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	ASSERT_TRUE(dataInput->inputFinished(true));

	//clean up
	delete dataInput;
}

TEST(DataInputTests, requiredParametersEdgeCases) {
	//seeding mass edge case
	DataInput* dataInput = new DataInput();
	std::ifstream paramFp("");
	int lineNo;

	for (int i = 0; i < requiredInputParams3PG.size(); i++) {
		if (requiredInputParams3PG[i] == "EndYear") {
			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "3" });
		}
		else {
			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "1" });
		}
	}
	for (int i = 0; i < requiredSeriesParams3PG.size(); i++) {
		dataInput->tryAddSeriesParam(requiredSeriesParams3PG[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	}
	dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	ASSERT_FALSE(dataInput->inputFinished(false));

	//must have all of WRi, WFi, and WSi if we don't have seedlingmass
	dataInput->tryAddInputParam("WRi", { "1" });
	dataInput->tryAddInputParam("WFi", { "1" });
	dataInput->tryAddInputParam("WSi", { "1" });
	ASSERT_TRUE(dataInput->inputFinished(false));
	
	//reset dataInput
	delete dataInput;
	dataInput = new DataInput();

	//series param VPD and Tavg which require Tmin and Tmax
	for (int i = 0; i < requiredInputParams3PG.size(); i++) {
		if (requiredInputParams3PG[i] == "EndYear") {
			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "3" });
		}
		else {
			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "1" });
		}
	}
	for (int i = 0; i < requiredSeriesParams3PG.size(); i++) {
		dataInput->tryAddSeriesParam(requiredSeriesParams3PG[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	}
	dataInput->tryAddInputParam("SeedlingMass", { "1" });
	ASSERT_FALSE(dataInput->inputFinished(false));

	//must have both Tmin and Tmax if we don't have either VPD or Tavg
	dataInput->tryAddSeriesParam("Tmin", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("Tmax", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);

	//clean up
	delete dataInput;
}

//test different input method for series params
TEST(DataInputTests, seriesParamInputFormat) {
	std::unique_ptr<GDALRasterImage> frostTif1 = std::make_unique<GDALRasterImage>(frost1);
	std::unique_ptr<GDALRasterImage> frostTif2 = std::make_unique<GDALRasterImage>(frost2);
	std::unique_ptr<GDALRasterImage> frostTif3 = std::make_unique<GDALRasterImage>(frost3);
	std::unique_ptr<GDALRasterImage> frostTif4 = std::make_unique<GDALRasterImage>(frost4);
	std::unique_ptr<GDALRasterImage> frostTif5 = std::make_unique<GDALRasterImage>(frost5);
	std::unique_ptr<GDALRasterImage> frostTif6 = std::make_unique<GDALRasterImage>(frost6);
	std::unique_ptr<GDALRasterImage> frostTif7 = std::make_unique<GDALRasterImage>(frost7);
	std::unique_ptr<GDALRasterImage> frostTif8 = std::make_unique<GDALRasterImage>(frost8);
	std::unique_ptr<GDALRasterImage> frostTif9 = std::make_unique<GDALRasterImage>(frost9);
	std::unique_ptr<GDALRasterImage> frostTif10 = std::make_unique<GDALRasterImage>(frost10);
	std::unique_ptr<GDALRasterImage> frostTif11 = std::make_unique<GDALRasterImage>(frost11);
	std::unique_ptr<GDALRasterImage> frostTif12 = std::make_unique<GDALRasterImage>(frost12);

	DataInput* dataInput = new DataInput();
	std::ifstream paramFp("test_files//DataInputTests//formatTest.txt");
	std::string line;
	//read line containing: "Frost", 
	std::getline(paramFp, line); 
	int lineNo = 0;
	//don't do anything with that line

	for (int i = 0; i < requiredInputParams3PG.size(); i++) {
		std::string curParam = requiredInputParams3PG[i];
		if (
			curParam != "yearPlanted" &&
			curParam != "StartAge" &&
			curParam != "EndYear" &&
			curParam != "StartMonth") {
			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "1" });
		}
	}
	dataInput->tryAddInputParam("SeedlingMass", { "1" });
	for (int i = 0; i < requiredSeriesParams3PG.size(); i++) {
		std::string curParam = requiredSeriesParams3PG[i];
		if (curParam != "Frost") {
			dataInput->tryAddSeriesParam(requiredSeriesParams3PG[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
		}
	}
	dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);

	// add parameter
	ASSERT_TRUE(dataInput->tryAddSeriesParam("Frost", {}, paramFp, lineNo));

	//add run period parameters for specific years of test
	dataInput->tryAddInputParam("StartAge", { "1" });
	dataInput->tryAddInputParam("StartMonth", { "1" });
	dataInput->tryAddInputParam("EndYear", { "1872" });
	dataInput->tryAddInputParam("yearPlanted", { "1869" });

	ASSERT_TRUE(dataInput->inputFinished(false));

	for (int year = 1870; year <= 1872; year++) {
		//1872 has grid series parameter
		if (year == 1872) {
			for (int i = 0; i < frostTif1->nRows; i++) {
				for (int j = 0; j < frostTif1->nCols; j++) {
					long cellIndex = (long)i * (long)frostTif1->nCols + (long)j;
					SeriesParams sParams;

					if (std::isnan(frostTif1->GetVal(cellIndex))) {
						//if nan, getSeriesParam should return false for all months
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 1, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 2, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 3, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 4, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 5, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 6, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 7, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 8, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 9, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 10, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 11, sParams));
						EXPECT_FALSE(dataInput->getSeriesParams(cellIndex, year, 12, sParams));
					}
					else {
						//otherwise, check values against test tif files for all months
						//January
						dataInput->getSeriesParams(cellIndex, year, 1, sParams);
						EXPECT_EQ(frostTif1->GetVal(cellIndex), sParams.FrostDays);

						//February
						dataInput->getSeriesParams(cellIndex, year, 2, sParams);
						EXPECT_EQ(frostTif2->GetVal(cellIndex), sParams.FrostDays);

						//March
						dataInput->getSeriesParams(cellIndex, year, 3, sParams);
						EXPECT_EQ(frostTif3->GetVal(cellIndex), sParams.FrostDays);

						//April
						dataInput->getSeriesParams(cellIndex, year, 4, sParams);
						EXPECT_EQ(frostTif4->GetVal(cellIndex), sParams.FrostDays);

						//May
						dataInput->getSeriesParams(cellIndex, year, 5, sParams);
						EXPECT_EQ(frostTif5->GetVal(cellIndex), sParams.FrostDays);

						//June
						dataInput->getSeriesParams(cellIndex, year, 6, sParams);
						EXPECT_EQ(frostTif6->GetVal(cellIndex), sParams.FrostDays);

						//July
						dataInput->getSeriesParams(cellIndex, year, 7, sParams);
						EXPECT_EQ(frostTif7->GetVal(cellIndex), sParams.FrostDays);

						//August
						dataInput->getSeriesParams(cellIndex, year, 8, sParams);
						EXPECT_EQ(frostTif8->GetVal(cellIndex), sParams.FrostDays);

						//September
						dataInput->getSeriesParams(cellIndex, year, 9, sParams);
						EXPECT_EQ(frostTif9->GetVal(cellIndex), sParams.FrostDays);

						//October
						dataInput->getSeriesParams(cellIndex, year, 10, sParams);
						EXPECT_EQ(frostTif10->GetVal(cellIndex), sParams.FrostDays);

						//November
						dataInput->getSeriesParams(cellIndex, year, 11, sParams);
						EXPECT_EQ(frostTif11->GetVal(cellIndex), sParams.FrostDays);

						//December
						dataInput->getSeriesParams(cellIndex, year, 12, sParams);
						EXPECT_EQ(frostTif12->GetVal(cellIndex), sParams.FrostDays);
					}
				}
			}
		}
		else {
			//non 1872 years have scalar parameters
			for (int month = 1; month <= 12; month++) {
				//the integer put into the test .txt file follows this pattern
				int expectedVal = (year - 1870) * 12 + month;

				//get params
				SeriesParams sParams;
				dataInput->getSeriesParams(0, year, month, sParams);

				//test frost (the param we are using to test the different format)
				EXPECT_EQ(expectedVal, sParams.FrostDays);
			}
		}	
	}

	//clean up
	delete dataInput;
}

//test management table functionality
TEST(DataInputTests, managementTable) {
	DataInput *dataInput = new DataInput();

	std::filesystem::path imagePath = std::filesystem::current_path();
	imagePath /= "test_files";
	imagePath /= "DataInputTests";
	imagePath /= "irr1997.tif";
	std::cout << imagePath.string() << std::endl;

	std::ofstream mgmtTest1("test_files//DataInputTests//manageTableTest1.txt");
	std::string line1 = "1994            8\n";
	std::string line2 = "1997 " + imagePath.string() + "\n";
	std::string line3 = "2000 5\n";
	mgmtTest1.write(line1.c_str(), line1.length());
	mgmtTest1.write(line2.c_str(), line2.length());
	mgmtTest1.write(line3.c_str(), line3.length());

	std::ofstream mgmtTest2("test_files//DataInputTests//manageTableTest2.txt");
	line1 = "1994	1\n";
	line2 = "	1995	2\n";
	line3 = "1996	3\n";
	mgmtTest2.write(line1.c_str(), line1.length());
	mgmtTest2.write(line2.c_str(), line2.length());
	mgmtTest2.write(line3.c_str(), line3.length());

	std::ofstream mgmtTest3("test_files//DataInputTests//manageTableTest3.txt");
	line1 = "    2 " + imagePath.string() + "\n";
	line2 = "3000 .0005\n";
	mgmtTest3.write(line1.c_str(), line1.length());
	mgmtTest3.write(line2.c_str(), line2.length());
	
	mgmtTest1.close();
	mgmtTest2.close();
	mgmtTest3.close();

	int lineNo1 = 0;
	int lineNo2 = 0;
	int lineNo3 = 0;
	
	std::ifstream paramFp1("test_files//DataInputTests//manageTableTest1.txt");
	std::ifstream paramFp2("test_files//DataInputTests//manageTableTest2.txt");
	std::ifstream paramFp3("test_files//DataInputTests//manageTableTest3.txt");

	ASSERT_TRUE(dataInput->tryAddManagementParam("Management: irrigation", paramFp1, lineNo1));
	ASSERT_TRUE(dataInput->tryAddManagementParam("management: minasw", paramFp2, lineNo2));
	ASSERT_TRUE(dataInput->tryAddManagementParam("management: FERTILITY", paramFp3, lineNo3));
	
	double val;
	
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~
	Irrigation
	~~~~~~~~~~~~~~~~~~~~~~~~~ */
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 0, 1, val));
	EXPECT_EQ(val, 8);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 0, 1993, val));
	EXPECT_EQ(val, 8);

	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 0, 1994, val));
	EXPECT_EQ(val, 1);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 1, 1994, val));
	EXPECT_EQ(val, 2);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 2, 1994, val)); //nan
	EXPECT_EQ(val, 8); //nan should go to previous value
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 3, 1994, val));
	EXPECT_EQ(val, 3);
	
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 0, 1996, val));
	EXPECT_EQ(val, 1);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 1, 1996, val));
	EXPECT_EQ(val, 2);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 2, 1996, val)); //nan
	EXPECT_EQ(val, 8); //nan should go to previous value
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 3, 1996, val));
	EXPECT_EQ(val, 3);
	
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 0, 1997, val));
	EXPECT_EQ(val, 5);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 0, 2001, val));
	EXPECT_EQ(val, 5);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::IRRIGATION, 0, 5000, val));
	EXPECT_EQ(val, 5);
	
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~
	MinASW
	~~~~~~~~~~~~~~~~~~~~~~~~~ */
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::MINASW, 999, 1, val));
	EXPECT_EQ(val, 1);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::MINASW, 999, 1993, val));
	EXPECT_EQ(val, 1);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::MINASW, 999, 1994, val));
	EXPECT_EQ(val, 2);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::MINASW, 999, 1995, val));
	EXPECT_EQ(val, 3);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::MINASW, 999, 1996, val));
	EXPECT_EQ(val, 3);
	
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~
	Fertility
	~~~~~~~~~~~~~~~~~~~~~~~~~ */
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::FERTILITY, 0, 1, val));
	EXPECT_EQ(val, 1);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::FERTILITY, 1, 1, val));
	EXPECT_EQ(val, 2);
	EXPECT_FALSE(dataInput->getManagementParam(ManagementIndex::FERTILITY, 2, 1, val)); //false because nan
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::FERTILITY, 3, 1, val));
	EXPECT_EQ(val, 3);
	
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::FERTILITY, 0, 2, val));
	EXPECT_EQ(val, .0005);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::FERTILITY, 0, 3000, val));
	EXPECT_EQ(val, .0005);
	EXPECT_TRUE(dataInput->getManagementParam(ManagementIndex::FERTILITY, 0, std::numeric_limits<int>::max(), val));
	EXPECT_EQ(val, .0005);
}

std::string startAge1 = "test_files//DataInputTests//runPeriod//StartAge1.tif";
std::string startAgeVariable = "test_files//DataInputTests//runPeriod//StartAgeVariable.tif";
std::string endYear1872 = "test_files//DataInputTests//runPeriod//EndYear1872.tif";
std::string endYearVariable = "test_files//DataInputTests//runPeriod//EndYearVariable.tif";
std::string yearPlanted1869 = "test_files//DataInputTests//runPeriod//yearPlanted1869.tif";
std::string yearPlantedVariable = "test_files//DataInputTests//runPeriod//yearPlantedVariable.tif";

//test run period
void addRequiredParamsExcludingRunPeriod(DataInput* dataInput) {
	std::ifstream paramFp("");
	int lineNo;

	for (int i = 0; i < requiredInputParams3PG.size(); i++) {
		std::string curParam = requiredInputParams3PG[i];
		if (
			curParam != "yearPlanted" &&
			curParam != "StartAge" &&
			curParam != "EndYear") {
			dataInput->tryAddInputParam(requiredInputParams3PG[i], { "1" });
		}
	}
	dataInput->tryAddInputParam("SeedlingMass", { "1" });
	for (int i = 0; i < requiredSeriesParams3PG.size(); i++) {
		dataInput->tryAddSeriesParam(requiredSeriesParams3PG[i], { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	}
	dataInput->tryAddSeriesParam("Tavg", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	dataInput->tryAddSeriesParam("VPD", { "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1" }, paramFp, lineNo);
	ASSERT_FALSE(dataInput->inputFinished(false));
}

TEST(DataInputTests, runPeriod1) {
	//add required params except run period params
	DataInput* dataInput = new DataInput();
	addRequiredParamsExcludingRunPeriod(dataInput);

	//add run period params
	dataInput->tryAddInputParam("StartAge", { "1" });
	dataInput->tryAddInputParam("EndYear", { "1872" });
	dataInput->tryAddInputParam("yearPlanted", { "1869" });

	//finish input
	ASSERT_TRUE(dataInput->inputFinished(false));

	//get run period
	RunPeriod runPeriod = dataInput->getRunPeriod();

	//test run period outputs
	EXPECT_EQ(runPeriod.StartYear, 1870);
	EXPECT_EQ(runPeriod.EndYear, 1872);

	//clean up
	delete dataInput;
}

TEST(DataInputTests, runPeriod2) {
	//add required params except run period params
	DataInput* dataInput = new DataInput();
	addRequiredParamsExcludingRunPeriod(dataInput);

	//add run period params
	dataInput->tryAddInputParam("StartAge", { startAge1 });
	dataInput->tryAddInputParam("EndYear", { "1872" });
	dataInput->tryAddInputParam("yearPlanted", { "1869" });

	//finish input
	ASSERT_TRUE(dataInput->inputFinished(false));

	//get run period
	RunPeriod runPeriod = dataInput->getRunPeriod();

	//test run period outputs
	EXPECT_EQ(runPeriod.StartYear, 1870);
	EXPECT_EQ(runPeriod.EndYear, 1872);

	//clean up
	delete dataInput;
}

TEST(DataInputTests, runPeriod3) {
	//add required params except run period params
	DataInput* dataInput = new DataInput();
	addRequiredParamsExcludingRunPeriod(dataInput);

	//add run period params
	dataInput->tryAddInputParam("StartAge", { startAge1 });
	dataInput->tryAddInputParam("EndYear", { endYear1872 });
	dataInput->tryAddInputParam("yearPlanted", { yearPlanted1869 });

	//finish input
	ASSERT_TRUE(dataInput->inputFinished(false));

	//get run period
	RunPeriod runPeriod = dataInput->getRunPeriod();

	//test run period outputs
	EXPECT_EQ(runPeriod.StartYear, 1870);
	EXPECT_EQ(runPeriod.EndYear, 1872);

	//clean up
	delete dataInput;
}

TEST(DataInputTests, runPeriod4) {
	//add required params except run period params
	DataInput* dataInput = new DataInput();
	addRequiredParamsExcludingRunPeriod(dataInput);

	//add run period params
	dataInput->tryAddInputParam("StartAge", { "1" });
	dataInput->tryAddInputParam("EndYear", { endYear1872 });
	dataInput->tryAddInputParam("yearPlanted", { yearPlanted1869 });

	//finish input
	ASSERT_TRUE(dataInput->inputFinished(false));

	//get run period
	RunPeriod runPeriod = dataInput->getRunPeriod();

	//test run period outputs
	EXPECT_EQ(runPeriod.StartYear, 1870);
	EXPECT_EQ(runPeriod.EndYear, 1872);

	//clean up
	delete dataInput;
}

TEST(DataInputTests, runPeriod5) {
	//add required params except run period params
	DataInput* dataInput = new DataInput();
	addRequiredParamsExcludingRunPeriod(dataInput);

	//add run period params
	dataInput->tryAddInputParam("StartAge", { startAgeVariable });
	dataInput->tryAddInputParam("EndYear", { endYearVariable });
	dataInput->tryAddInputParam("yearPlanted", { yearPlantedVariable });

	//finish input
	ASSERT_TRUE(dataInput->inputFinished(false));

	//get run period
	RunPeriod runPeriod = dataInput->getRunPeriod();

	//test run period outputs
	EXPECT_EQ(runPeriod.StartYear, 1870);
	EXPECT_EQ(runPeriod.EndYear, 1874);

	//clean up
	delete dataInput;
}

TEST(DataInputTests, runPeriod6) {
	//add required params except run period params
	DataInput* dataInput = new DataInput();
	addRequiredParamsExcludingRunPeriod(dataInput);

	//add run period params
	dataInput->tryAddInputParam("StartAge", { startAge1 });
	dataInput->tryAddInputParam("EndYear", { endYearVariable });
	dataInput->tryAddInputParam("yearPlanted", { yearPlantedVariable });

	//finish input
	ASSERT_TRUE(dataInput->inputFinished(false));

	//get run period
	RunPeriod runPeriod = dataInput->getRunPeriod();

	//test run period outputs
	EXPECT_EQ(runPeriod.StartYear, 1870);
	EXPECT_EQ(runPeriod.EndYear, 1874);

	//clean up
	delete dataInput;
}

TEST(DataInputTests, runPeriod7) {
	//add required params except run period params
	DataInput* dataInput = new DataInput();
	addRequiredParamsExcludingRunPeriod(dataInput);

	//add run period params
	dataInput->tryAddInputParam("StartAge", { startAgeVariable });
	dataInput->tryAddInputParam("EndYear", { endYear1872 });
	dataInput->tryAddInputParam("yearPlanted", { yearPlanted1869 });

	//finish input
	ASSERT_TRUE(dataInput->inputFinished(false));

	//get run period
	RunPeriod runPeriod = dataInput->getRunPeriod();

	//test run period outputs
	EXPECT_EQ(runPeriod.StartYear, 1870);
	EXPECT_EQ(runPeriod.EndYear, 1872);

	//clean up
	delete dataInput;
}