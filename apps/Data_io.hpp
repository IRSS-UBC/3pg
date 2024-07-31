#include <string>
#include "GDALRasterImage.hpp"
#include "DataOutput.hpp"
#include "DataInput.hpp"
#include "ParamStructs.hpp"
#include "MYDate.h" 

// Management tables
#define MT_FERTILITY  1
#define MT_MINASW     2
#define MT_IRRIGATION 3

GDALRasterImage* openInputGrids();
void initDataOutput(std::shared_ptr<GDALRasterImage> refGrid);

double lookupManageTable( int year, int table, double def, int cellIndex ); 

void readSpeciesParamFile(const std::string& speciesFile, DataInput *dataInput);
std::unordered_map<std::string, PPPG_OP_VAR> readSiteParamFile(const std::string& paramFile, DataInput *dataInput);
PPPG_OP_VAR readOutputParam(const std::string& pName, const std::vector<std::string>& pValue, int lineNo);

void writeMonthlyOutputGrids(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, int calYear, int calMonth, MYDate minMY, MYDate maxMY, long cellIndex );
int writeOutputGrids(const std::unordered_map<std::string, PPPG_OP_VAR>& opVars, long cellIndex);

std::string getOutPathTMP(const std::string& siteParamFile);
void deleteDataOutput();
CPLErr writeRowDataOutput(int row);