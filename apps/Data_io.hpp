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
void InitInputParams(void);
void initDataOutput(GDALRasterImage* refGrid);

bool loadParamVals(int k);
// table must be one of MT_FERTILITY, MT_MINASW, MT_MINASW. 
double lookupManageTable( int year, int table, double def, int cellIndex ); 
/**
 * Find the month/year period over which the model will run.
 *
 * @param[out] minMY the minimum (start) year
 * @param[out] maxMY the maximum (end) month and year
 */
int findRunPeriod( MYDate &minMY, MYDate &maxMY );
// Check that min and max month/years are valid.
bool validRunPeriod(const MYDate& minMY, const MYDate& maxMY);

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
bool haveAgeDepFert(void);
bool haveMinASWTG(void);
bool haveSeedlingMass(void);
bool haveSpatialRunYears(void);
bool haveRhoMin(void);  //Standage dependant Density 15/07/2002
bool haveRhoMax(void);  //Standage dependant Density 15/07/2002
bool haveTRho(void);    //Standage dependant Density 15/07/2002

bool getSeriesVal(double &val, int ser, int calMonth, int calYear, int k);
std::string getOutPathTMP(const std::string& siteParamFile);
void deleteDataOutput();