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

void setLogFunc(std::function<void(std::string)>& log);

void readSpeciesParamFile(const std::string& speciesFile, DataInput&dataInput);
void readSiteParamFile(const std::string& paramFile, DataInput& dataInput);

std::string getOutPathTMP(const std::string& siteParamFile);