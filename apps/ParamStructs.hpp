#pragma once 

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "GDALRasterImage.hpp"

typedef enum { 
	pNull = 0, 
	pScalar = 1, 
	pTif = 2 
} ParamSpatial;

typedef struct PPPG_PARAM {
	//variable name
	std::string id = "";

	//indication of scalar or grid parameter
	ParamSpatial spType = pNull;

	//scalar value
	double val;
	
	//grid reference
	std::unique_ptr<GDALRasterImage> g;
} PPPG_PARAM;

// 3PG output variables. In spatial mode output variables may be written 
// repeatedly, on a time step defined by recurStart, recurYear, and recurMonthly. 
typedef struct PPPG_OP_VAR {
	//variable name
	std::string id; 

	//file name if given file
	std::string gridName;

	//first year to write regular output
	int recurStart;

	// Interval on which to write regular output. 
	int recurYear = -1;

	// Single month number we want output in. 
	int recurMonth;    

	// Whether to write every month in an output year. 
	bool recurMonthly;              
} PPPG_OP_VAR;
 
typedef struct PPPG_SERIES_PARAM {
	//monthlyParams will always contain a multiple of 12 params, the first year is the
	//first year in the vector (monthlyparams[0] will be January of the first year)
	//first year will be set to -1 if the same monthly params will be used in every year.
	int firstYear = 0;

	//end year, used to ensure we don't try to get an element in the vector that doesn't exist
	//there will always be (lastYear + 1 - firstYear) * 12 elements in the vector
	//unless first year is set to -1 indicating there are only 12 elements.
	int lastYear = 0;

	//vector containing monthly params
	std::vector<std::unique_ptr<PPPG_PARAM>> monthlyParams;
} PPPG_SERIES_PARAM;

// 3PG 'management table' parameters. 
typedef struct PPPG_MT_PARAM {
	//start year of the yearlyParams vector. The year that a yearlyParams entry represents is the vector index + firstYear
	int firstYear = -1;

	//converts a year to an index in the yearlyParams vector
	std::vector<int> yearToIndex;

	//management params by year
	std::vector<std::unique_ptr<PPPG_PARAM>> yearlyParams;
} PPPG_MT_PARAM;