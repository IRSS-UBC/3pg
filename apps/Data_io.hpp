// 3PG Input routines. 


#include <string>
#include "GDALRasterImage.hpp"
#include "DataOutput.hpp"
#include "DataInput.hpp"
#include "ParamStructs.hpp"
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

GDALRasterImage* openInputGrids();
void initDataOutput(GDALRasterImage* refGrid);

double lookupManageTable( int year, int table, double def, int cellIndex ); 

void readSpeciesParamFile(const std::string& speciesFile, DataInput *dataInput);
std::unordered_map<std::string, PPPG_OP_VAR> readSiteParamFile(const std::string& paramFile, DataInput *dataInput);
void readSampleFile(std::unordered_map<std::string, PPPG_OP_VAR> &opVars, GDALRasterImage* refGrid);
PPPG_OP_VAR readOutputParam(const std::string& pName, const std::vector<std::string>& pValue, int lineNo);

void writeSampleFiles(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, int cellIndex, int month, long calYear);
void writeMonthlyOutputGrids(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, int calYear, int calMonth, bool hitNODATA, MYDate minMY, MYDate maxMY, long cellIndex );
int writeOutputGrids(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, bool hitNODATA, long cellIndex);

bool userVpdSeries(void);
bool userNetRadSeries(void);
bool userTavgSeries(void);

bool haveAllParams();

bool getSeriesVal(double &val, int ser, int calMonth, int calYear, int k);
std::string getOutPathTMP(const std::string& siteParamFile);
void deleteDataOutput();