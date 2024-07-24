#include <gtest/gtest.h>
#include <../apps/DataInput.hpp>
#include <fstream>

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
