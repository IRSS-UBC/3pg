#include <gtest/gtest.h>
#include <../apps/DataInput.hpp>
#include <fstream>

#include <../apps/GDALRasterImage.hpp>

TEST(DataInputTests, correctParamNames) {
	DataInput* dataInput = new DataInput();

	//long versions according to legacy documentations
	EXPECT_TRUE(dataInput->tryAddInputParam("Foliage:stem partitioning ratio @ D=2 cm", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Foliage:stem partitioning ratio @ D=20 cm", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Constant in the stem mass v. diam. relationship", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Power in the stem mass v. diam. relationship", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Maximum fraction of NPP to roots", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Minimum fraction of NPP to roots", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Minimum temperature for growth", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Optimum temperature for growth", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Maximum temperature for growth", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Days production lost per frost day", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Moisture ratio deficit for fq = 0.5", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Power of moisture ratio deficit", { "1" }));
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

	//reset DataInput
	delete dataInput;
	dataInput = new DataInput();

	//short versions
	EXPECT_TRUE(dataInput->tryAddInputParam("pFS2", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("pFS20", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("StemConst", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("StemPower", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("pRx", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("pRn", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("growthTmin", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("growthTopt", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("growthTmax", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("kF", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("gammaFx", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("gammaF0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("tgammaF", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("Rttover", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("MaxCond", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("CoeffCond", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("BLcond", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("m0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("fN0", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("fNn", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("thinPower", { "1" }));
	EXPECT_TRUE(dataInput->tryAddInputParam("mF", { "1" }));
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
	dataInput->tryAddInputParam("EndYear", { "60" });
	dataInput->tryAddInputParam("StartMonth", { "61" });
	dataInput->tryAddInputParam("yearPlanted", { "62" });

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
	ASSERT_EQ(params.pFS2, (double)1);
	ASSERT_EQ(params.pFS20, (double)2);
	ASSERT_EQ(params.StemConst, (double)3);
	ASSERT_EQ(params.StemPower, (double)4);
	ASSERT_EQ(params.pRx, (double)5);
	ASSERT_EQ(params.pRn, (double)6);
	ASSERT_EQ(params.growthTmin, (double)7);
	ASSERT_EQ(params.growthTopt, (double)8);
	ASSERT_EQ(params.growthTmax, (double)9);
	ASSERT_EQ(params.kF, (double)10);
	ASSERT_EQ(params.gammaFx, (double)11);
	ASSERT_EQ(params.gammaF0, (double)21);
	ASSERT_EQ(params.tgammaF, (double)13);
	ASSERT_EQ(params.Rttover, (double)14);
	ASSERT_EQ(params.MaxCond, (double)15);
	ASSERT_EQ(params.CoeffCond, (double)16);
	ASSERT_EQ(params.BLcond, (double)17);
	ASSERT_EQ(params.m0, (double)18);
	ASSERT_EQ(params.fN0, (double)19);
	ASSERT_EQ(params.fNn, (double)20);
	ASSERT_EQ(params.thinPower, (double)21);
	ASSERT_EQ(params.mF, (double)22);
	ASSERT_EQ(params.mR, (double)23);
	ASSERT_EQ(params.mS, (double)24);
	ASSERT_EQ(params.SWconst0, (double)25);
	ASSERT_EQ(params.SWpower0, (double)26);
	ASSERT_EQ(params.wSx1000, (double)27);
	ASSERT_EQ(params.MaxAge, (double)28);
	ASSERT_EQ(params.nAge, (double)29);
	ASSERT_EQ(params.rAge, (double)30);
	ASSERT_EQ(params.SLA0, (double)31);
	ASSERT_EQ(params.SLA1, (double)32);
	ASSERT_EQ(params.tSLA, (double)33);
	ASSERT_EQ(params.k, (double)34);
	ASSERT_EQ(params.fullCanAge, (double)35);
	ASSERT_EQ(params.alpha, (double)36);
	ASSERT_EQ(params.fracBB0, (double)37);
	ASSERT_EQ(params.fracBB1, (double)38);
	ASSERT_EQ(params.tBB, (double)39);
	ASSERT_EQ(params.y, (double)40);
	ASSERT_EQ(params.rhoMin, (double)41);
	ASSERT_EQ(params.rhoMax, (double)42);
	ASSERT_EQ(params.tRho, (double)43);
	ASSERT_EQ(params.Qa, (double)44);
	ASSERT_EQ(params.Qb, (double)45);
	ASSERT_EQ(params.gDM_mol, (double)46);
	ASSERT_EQ(params.molPAR_MJ, (double)47);
	ASSERT_EQ(params.LAIgcx, (double)48);
	ASSERT_EQ(params.MaxIntcptn, (double)49);
	ASSERT_EQ(params.LAImaxIntcptn, (double)50);
	ASSERT_EQ(params.Lat, (double)51);
	ASSERT_EQ(params.FRp, (double)52);
	ASSERT_EQ(params.FRstart, (double)53);
	ASSERT_EQ(params.FRend, (double)54);
	ASSERT_EQ(params.FRdec, (double)55);
	ASSERT_EQ(params.soilIndex, (double)56);
	ASSERT_EQ(params.MaxASW, (double)57);
	ASSERT_EQ(params.MinASWp, (double)58);
	ASSERT_EQ(params.StartAge, (double)59);
	ASSERT_EQ(params.EndYear, (double)60);
	ASSERT_EQ(params.StartMonth, (double)61);
	ASSERT_EQ(params.yearPlanted, (double)62);
	ASSERT_EQ(params.SeedlingMass, (double).0000000000001);
	ASSERT_EQ(params.WFi, (double).00000000000001);
	ASSERT_EQ(params.WRi, (double).000000000000001);
	ASSERT_EQ(params.WSi, (double).0000000000000001);
	ASSERT_EQ(params.StemNoi, (double).00000000000000001);
	ASSERT_EQ(params.ASWi, std::numeric_limits<double>::min() + 1);
	ASSERT_EQ(params.MinASWTG, std::numeric_limits<double>::max() - 1);
	ASSERT_EQ(params.NDVI_FPAR_intercept, std::numeric_limits<double>::min() + .01);
	ASSERT_EQ(params.NDVI_FPAR_constant, std::numeric_limits<double>::max() - .01);

	SeriesParams sParams;
	//series params January
	dataInput->getSeriesParams(0, 0, 1, sParams);
	ASSERT_EQ(sParams.Tavg, (double)1);
	ASSERT_EQ(sParams.Rain, (double)13);
	ASSERT_EQ(sParams.SolarRad, (double)25);
	ASSERT_EQ(sParams.FrostDays, (double)37);
	ASSERT_EQ(sParams.NDVI_AVH, (double)49);
	ASSERT_EQ(sParams.NetRad, (double)61);
	ASSERT_EQ(sParams.VPD, (double)73);

	//series params February
	dataInput->getSeriesParams(0, 0, 2, sParams);
	ASSERT_EQ(sParams.Tavg, (double)2);
	ASSERT_EQ(sParams.Rain, (double)14);
	ASSERT_EQ(sParams.SolarRad, (double)26);
	ASSERT_EQ(sParams.FrostDays, (double)38);
	ASSERT_EQ(sParams.NDVI_AVH, (double)50);
	ASSERT_EQ(sParams.NetRad, (double)62);
	ASSERT_EQ(sParams.VPD, (double)74);

	//series params March
	dataInput->getSeriesParams(0, 0, 3, sParams);
	ASSERT_EQ(sParams.Tavg, (double)3);
	ASSERT_EQ(sParams.Rain, (double)15);
	ASSERT_EQ(sParams.SolarRad, (double)27);
	ASSERT_EQ(sParams.FrostDays, (double)39);
	ASSERT_EQ(sParams.NDVI_AVH, (double)51);
	ASSERT_EQ(sParams.NetRad, (double)63);
	ASSERT_EQ(sParams.VPD, (double)75);

	//series params April
	dataInput->getSeriesParams(0, 0, 4, sParams);
	ASSERT_EQ(sParams.Tavg, (double)4);
	ASSERT_EQ(sParams.Rain, (double)16);
	ASSERT_EQ(sParams.SolarRad, (double)28);
	ASSERT_EQ(sParams.FrostDays, (double)40);
	ASSERT_EQ(sParams.NDVI_AVH, (double)52);
	ASSERT_EQ(sParams.NetRad, (double)64);
	ASSERT_EQ(sParams.VPD, (double)76);

	//series params May
	dataInput->getSeriesParams(0, 0, 5, sParams);
	ASSERT_EQ(sParams.Tavg, (double)5);
	ASSERT_EQ(sParams.Rain, (double)17);
	ASSERT_EQ(sParams.SolarRad, (double)29);
	ASSERT_EQ(sParams.FrostDays, (double)41);
	ASSERT_EQ(sParams.NDVI_AVH, (double)53);
	ASSERT_EQ(sParams.NetRad, (double)65);
	ASSERT_EQ(sParams.VPD, (double)77);

	//series params June
	dataInput->getSeriesParams(0, 0, 6, sParams);
	ASSERT_EQ(sParams.Tavg, (double)6);
	ASSERT_EQ(sParams.Rain, (double)18);
	ASSERT_EQ(sParams.SolarRad, (double)30);
	ASSERT_EQ(sParams.FrostDays, (double)42);
	ASSERT_EQ(sParams.NDVI_AVH, (double)54);
	ASSERT_EQ(sParams.NetRad, (double)66);
	ASSERT_EQ(sParams.VPD, (double)78);

	//series params July
	dataInput->getSeriesParams(0, 0, 7, sParams);
	ASSERT_EQ(sParams.Tavg, (double)7);
	ASSERT_EQ(sParams.Rain, (double)19);
	ASSERT_EQ(sParams.SolarRad, (double)31);
	ASSERT_EQ(sParams.FrostDays, (double)43);
	ASSERT_EQ(sParams.NDVI_AVH, (double)55);
	ASSERT_EQ(sParams.NetRad, (double)67);
	ASSERT_EQ(sParams.VPD, (double)79);

	//series params August
	dataInput->getSeriesParams(0, 0, 8, sParams);
	ASSERT_EQ(sParams.Tavg, (double)8);
	ASSERT_EQ(sParams.Rain, (double)20);
	ASSERT_EQ(sParams.SolarRad, (double)32);
	ASSERT_EQ(sParams.FrostDays, (double)44);
	ASSERT_EQ(sParams.NDVI_AVH, (double)56);
	ASSERT_EQ(sParams.NetRad, (double)68);
	ASSERT_EQ(sParams.VPD, (double)80);

	//series params September
	dataInput->getSeriesParams(0, 0, 9, sParams);
	ASSERT_EQ(sParams.Tavg, (double)9);
	ASSERT_EQ(sParams.Rain, (double)21);
	ASSERT_EQ(sParams.SolarRad, (double)33);
	ASSERT_EQ(sParams.FrostDays, (double)45);
	ASSERT_EQ(sParams.NDVI_AVH, (double)57);
	ASSERT_EQ(sParams.NetRad, (double)69);
	ASSERT_EQ(sParams.VPD, (double)81);

	//series params October
	dataInput->getSeriesParams(0, 0, 10, sParams);
	ASSERT_EQ(sParams.Tavg, (double)10);
	ASSERT_EQ(sParams.Rain, (double)22);
	ASSERT_EQ(sParams.SolarRad, (double)34);
	ASSERT_EQ(sParams.FrostDays, (double)46);
	ASSERT_EQ(sParams.NDVI_AVH, (double)58);
	ASSERT_EQ(sParams.NetRad, (double)70);
	ASSERT_EQ(sParams.VPD, (double)82);

	//series params November
	dataInput->getSeriesParams(0, 0, 11, sParams);
	ASSERT_EQ(sParams.Tavg, (double)11);
	ASSERT_EQ(sParams.Rain, (double)23);
	ASSERT_EQ(sParams.SolarRad, (double)35);
	ASSERT_EQ(sParams.FrostDays, (double)47);
	ASSERT_EQ(sParams.NDVI_AVH, (double)59);
	ASSERT_EQ(sParams.NetRad, (double)71);
	ASSERT_EQ(sParams.VPD, (double)83);

	//series params December
	dataInput->getSeriesParams(0, 0, 12, sParams);
	ASSERT_EQ(sParams.Tavg, (double)12);
	ASSERT_EQ(sParams.Rain, (double)24);
	ASSERT_EQ(sParams.SolarRad, (double)36);
	ASSERT_EQ(sParams.FrostDays, (double)48);
	ASSERT_EQ(sParams.NDVI_AVH, (double)60);
	ASSERT_EQ(sParams.NetRad, (double)72);
	ASSERT_EQ(sParams.VPD, (double)84);
}

//test grid parameters
TEST(DataInputTests, gridValues) {
	std::string testFile1 = "test_files/DataInputTests/test1.tif";
	std::string testFile2 = "test_files/DataInputTests/test2.tif";
	std::string testFile3 = "test_files/DataInputTests/test3.tif";
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

	GDALRasterImage* testTif1 = new GDALRasterImage(testFile1);
	GDALRasterImage* testTif2 = new GDALRasterImage(testFile2);
	GDALRasterImage* testTif3 = new GDALRasterImage(testFile3);
	GDALRasterImage* frostTif1 = new GDALRasterImage(frost1);
	GDALRasterImage* frostTif2 = new GDALRasterImage(frost2);
	GDALRasterImage* frostTif3 = new GDALRasterImage(frost3);
	GDALRasterImage* frostTif4 = new GDALRasterImage(frost4);
	GDALRasterImage* frostTif5 = new GDALRasterImage(frost5);
	GDALRasterImage* frostTif6 = new GDALRasterImage(frost6);
	GDALRasterImage* frostTif7 = new GDALRasterImage(frost7);
	GDALRasterImage* frostTif8 = new GDALRasterImage(frost8);
	GDALRasterImage* frostTif9 = new GDALRasterImage(frost9);
	GDALRasterImage* frostTif10 = new GDALRasterImage(frost10);
	GDALRasterImage* frostTif11 = new GDALRasterImage(frost11);
	GDALRasterImage* frostTif12 = new GDALRasterImage(frost12);

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
	dataInput->tryAddInputParam("EndYear", { "1" });
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
}

//test bad grids
TEST(DataInputTests, badGrids) {
	std::string testFile = "test_files/DataInputTests/test1.tif";
	std::string differentLocationFile = "test_files/DataInputTests/differentLocation.tif";
	std::string differentDimensionsFile = "test_files/DataInputTests/differentDimensions.tif";

	GDALRasterImage* testTif = new GDALRasterImage(testFile);
	GDALRasterImage* difLocationTif = new GDALRasterImage(differentLocationFile);
	GDALRasterImage* difDimensionsTif = new GDALRasterImage(differentDimensionsFile);

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
}

//test run period

//test required parameters

//test different input method for series params

//test series param calculations (using Tmin and Tmax);