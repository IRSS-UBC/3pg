#pragma once 

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "GDALRasterImage.hpp"

typedef enum { pNull, pScalar, pTif } ParamSpatial;

struct PPPG_VVAL {
	ParamSpatial spType = pNull;                   // Scalar, grid, or null
	double sval;                           // Scalar value. 
	std::string gridName;                        // ptr to grid file name
	GDALRasterImage* g;                               // ptr to grid value
};

typedef struct PPPG_PARAM {
	//variable name
	std::string id = "";

	//indication of scalar or grid parameter
	ParamSpatial spType = pNull;

	//scalar value
	double val;
	
	//grid reference
	std::shared_ptr<GDALRasterImage> g;
} PPPG_PARAM;

// 3PG output variables. In spatial mode output variables may be written 
// repeatedly, on a time step defined by recurStart, recurYear, and recurMonthly. 
typedef struct PPPG_OP_VAR {
	std::string id;                       // String version of the variable name. 
	double v;                    // The address of the model variable. 
	ParamSpatial spType;            // If its a spatial parameter and what grid type. 
	std::string gridName;  // The gridname, in spatial mode. 
	//GDALRasterImage *g;                        // The final output grid, in spatial mode. 
	bool write; // Whether the variable is wanted. 
	int recurStart;                 // First year to write regular output. 
	int recurYear = -1;                  // Interval on which to write regular output. 
	int recurMonth;                 // Single month number we want output in. 
	bool recurMonthly;              // Whether to write every month in an output year. 
	//std::vector<GDALRasterImage*> RO;                      // The output tifs for regular output. 
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

	//vector containing monthl params
	std::vector<PPPG_PARAM> monthlyParams;
} PPPG_SERIES_PARAM;

// 3PG 'management table' parameters.  Only one value per year is allowed.  
typedef struct PPPG_MT_PARAM {
	int year;                                   // Calendar year
	bool got = 0;
	PPPG_VVAL data;
} PPPG_MT_PARAM;