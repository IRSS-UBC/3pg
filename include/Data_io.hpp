// 3PG Input routines. 


#include "GDALRasterImage.hpp"
#include <string>
#include "gdal.h"
#include "gdal_priv.h"
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

bool loadParamVals(int k);
// table must be one of MT_FERTILITY, MT_MINASW, MT_MINASW. 
double lookupManageTable( int year, int table, double def, int cellIndex ); 
void writeMonthlyOutputGrids( int calYear, int calMonth, bool hitNODATA, MYDate minMY, MYDate maxMY );
void writeYearlyOutputGrids( int calYear, int calMonth, bool hitNODATA, MYDate minMY, MYDate maxMY );
void writeSampleFiles(int cellIndex, int month, long calYear);
void saveVariableVals(int k, bool hitNODATA);
// FILE *openLogFile(const std::string& siteParamFile);
void readParamFile(const std::string& paramFile);
GDALRasterImage* openInputGrids();
bool haveAllParams();
bool havePointOpFile();
int findRunPeriod( GDALRasterImage*refGrid, MYDate &minMY, MYDate &maxMY );
int openOutputGrids( GDALRasterImage*refGrid);
void ResetGrids(void);
void CloseGrids(void);
void PrintGrids(void);
int openRegularOutputGrids( GDALRasterImage*refGrid, MYDate spMinMY, MYDate spMaxMY );
void readSampleFile( GDALRasterImage*refGrid );
int writeOutputGrids(void);
void writeStandSummary(int year);
void InitInputParams(void);
bool userVpdSeries(void);
bool userNetRadSeries(void);
bool userTavgSeries(void);
bool haveAgeDepFert(void);
bool haveMinASWTG(void);
bool haveSeedlingMass(void);
bool haveSpatialRunYears(void);
bool haveRhoMin(void);  //Standage dependant Density 15/07/2002
bool haveRhoMax(void);  //Standage dependant Density 15/07/2002
bool haveTRho(void);    //Standage dependant Density 15/07/2002
bool getSeriesVal(double &val, int ser, int calMonth, int calYear, int k);
