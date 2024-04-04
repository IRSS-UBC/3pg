// 3PG Input routines. 


#include <string>
#include <vector>
#include "Params.hpp"
#include "GDALRasterImage.hpp"
// Must include FloatGrid.hpp prior to this. 
#include "MYDate.h" 

// Management tables
#define MT_FERTILITY  1
#define MT_MINASW     2
#define MT_IRRIGATION 3

// Time series input data. 
#define SS_TMAX      0
#define SS_TMIN      1
#define SS_RAIN      2
#define SS_SOLARRAD  3
#define SS_FROSTDAYS 4
#define SS_NDVI_AVH  5
#define SS_NETRAD    6
#define SS_VPD       7
#define SS_TAVG      8

bool loadParamVals(int k, std::vector<PPPG_PARAM>& params);
// table must be one of MT_FERTILITY, MT_MINASW, MT_MINASW. 
double lookupManageTable( int year, int table, double def, int cellIndex ); 
void writeMonthlyOutputGrids(const std::vector<PPPG_OP_VAR> opVars, int calYear, int calMonth, bool hitNODATA, MYDate minMY, MYDate maxMY, long cellIndex );
void writeYearlyOutputGrids(const std::vector<PPPG_OP_VAR> opVars, int calYear, int calMonth, bool hitNODATA, MYDate minMY, MYDate maxMY, long cellIndex );
void writeSampleFiles(int cellIndex, int month, long calYear);
void saveVariableVals(int k, bool hitNODATA);
// FILE *openLogFile(const std::string& siteParamFile);
void readParamFile(const std::string& paramFile, std::vector<PPPG_PARAM>& params, std::vector<PPPG_OP_VAR>& opVars, std::vector<PPPG_SERIES_PARAM>& series, std::vector<PPPG_MT_PARAM>& managment);
GDALRasterImage* openInputGrids(std::vector<PPPG_PARAM>& params, std::vector<PPPG_SERIES_PARAM>& series, std::vector<PPPG_MT_PARAM>& mgmnt);
int threadOpenInputTIFs(std::vector<PPPG_PARAM>& params, std::vector<PPPG_SERIES_PARAM>& series);
bool haveAllParams(std::vector<PPPG_PARAM>& params, std::vector<PPPG_SERIES_PARAM>& series);
bool havePointOpFile();
int findRunPeriod( GDALRasterImage*refGrid, MYDate &minMY, MYDate &maxMY, std::vector<PPPG_PARAM>& params);
int createOutputGrids( GDALRasterImage*refGrid, std::vector<PPPG_OP_VAR>& opVars);
int threadOpenOutputTIFs(std::vector<PPPG_OP_VAR>& opVars);
void ResetGrids(void);
void CloseGrids(std::vector<PPPG_PARAM>& params, std::vector<PPPG_OP_VAR>& opVars);
void PrintGrids(void);
int pNameToInd(const std::string& id, const std::vector<PPPG_PARAM>& params);
int opNameToInd(const std::string& id, const std::vector<PPPG_OP_VAR>& opVars);
int openRegularOutputGrids( GDALRasterImage*refGrid, MYDate spMinMY, MYDate spMaxMY, std::vector<PPPG_OP_VAR>& opVars);
int threadOpenRegularOutputTIFs(MYDate spMinMY, MYDate spMaxMY, std::vector<PPPG_OP_VAR>& opVars);
void readSampleFile( GDALRasterImage*refGrid, const std::vector<PPPG_OP_VAR>& opVars);
int writeOutputGrids(bool hitNODATA, long cellIndex, const std::vector<PPPG_OP_VAR>& opVars);
void writeStandSummary(int year);
std::vector<PPPG_PARAM> InitInputParams(void);
std::vector<PPPG_OP_VAR> initOutputVars(void);
std::vector<PPPG_SERIES_PARAM> initSeriesParams(void);
std::vector<PPPG_MT_PARAM> initMTParams(void);
bool userVpdSeries(void);
bool userNetRadSeries(void);
bool userTavgSeries(void);
bool haveAgeDepFert(const std::vector<PPPG_PARAM>& params);
bool haveMinASWTG(const std::vector<PPPG_PARAM>& params);
bool haveSeedlingMass(const std::vector<PPPG_PARAM>& params);
bool haveSpatialRunYears(const std::vector<PPPG_PARAM>& params);
bool haveRhoMin(const std::vector<PPPG_PARAM>& params);  //Standage dependant Density 15/07/2002
bool haveRhoMax(const std::vector<PPPG_PARAM>& params);  //Standage dependant Density 15/07/2002
bool haveTRho(const std::vector<PPPG_PARAM>& params);    //Standage dependant Density 15/07/2002
bool getSeriesVal(double &val, PPPG_SERIES_PARAM& series, int calMonth, int calYear, int k);
bool seriesNotNull(PPPG_SERIES_PARAM& series);
void readInputScanlines(std::vector<PPPG_PARAM>& inputs, std::vector<std::vector<float>>& inputScanlines, int row, int ncols);
void readSeriesScanlines(std::vector<PPPG_SERIES_PARAM>& series, std::vector<std::vector<float>>& seriesScanlines, int row, int ncols);
