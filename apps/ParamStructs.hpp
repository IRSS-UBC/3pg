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

// 3PG 'parameters'. These are all stored as double.  The 'id' string field 
// is set to the name of the variable in the initialisation below.  Within 
// the model itself we don't reference the parameters via this type.  Its 
// used to help identify parameter lines, and to help get grid values into 
// the model.  
typedef struct PPPG_PARAM {
	std::string id = "";                        // String version of the variable name. 
	double* adr;                     // The address of the model variable. 
	bool got = 0;                        // Has the parameter been set? 
	PPPG_VVAL data;                  // Variant value
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

// 3PG 'series' parameters.  This is any parameter with a time series for value, 
// in particular the climate parameters, and NDVI.  
typedef struct PPPG_SERIES_PARAM {
	int start;                             // Calendar year of first entry. 
	PPPG_VVAL* data;                       // Array of variant values.
	int vlen;                              // Number of entries in array. 
	bool oneYear;                          // Array is of a single 'average' year (eg esoclim). 
	bool got = 0;                              // Have read the series. 
} PPPG_SERIES_PARAM;

// 3PG 'management table' parameters.  Only one value per year is allowed.  
typedef struct PPPG_MT_PARAM {
	int year;                                   // Calendar year
	bool got = 0;
	PPPG_VVAL data;
} PPPG_MT_PARAM;