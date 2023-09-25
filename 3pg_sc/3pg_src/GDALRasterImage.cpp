#include <stdexcept>
#include <tuple>
#include <gdal.h>
#include <gdal_priv.h>
#include <math.h>
#include <iostream>
#include <string>
#include "GDALRasterImage.hpp"

GDALRasterImage::GDALRasterImage(std::string filename) {
	// Open the raster source located at `filename`
	GDALAllRegister();
	const GDALAccess eAccess = GA_ReadOnly;
	dataset = GDALDataset::FromHandle(GDALOpen(filename.c_str(), eAccess));
	if (!dataset) {
		throw std::invalid_argument("File cannot be opened.");
	}

	// Assume there is only one band in the raster source and use that
	if (dataset->GetRasterCount() != 1) {
		throw std::invalid_argument("TIF must have only one band.");
	}
	band = dataset->GetRasterBand(1);

	// Get the inverse geo transform to map from geo location -> pixel location
	double datasetTransform[6] = {};
	double inverseTransform[6] = {};
	dataset->GetGeoTransform(datasetTransform);
	if (GDALInvGeoTransform(datasetTransform, inverseTransform) != CE_None) {
		throw std::invalid_argument("Cannot get inverse transform.");
	}

	// get crs
	crs = dataset->GetProjectionRef();
	if (crs == NULL) {
		throw std::invalid_argument("CRS is not defined.");
	}

	noData = band->GetNoDataValue();
	nRows = band->GetYSize();
	nCols = band->GetXSize();

	// get extent of raster
	xMin = datasetTransform[0];
	yMax = datasetTransform[3];
	xMax = xMin + datasetTransform[1] * nCols;
	yMin = yMax + datasetTransform[5] * nRows;

};

GDALRasterImage::GDALRasterImage(std::string filename, GDALRasterImage* refGrid) {
	// Create a new GDALRasterImage dataset with one band with the same extent, transform, and crs as refGrid
	if (Exists(filename)) {
		throw std::invalid_argument("File already exists.");
	};
	GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* dataset = driver->Create(filename.c_str(), refGrid->nCols, refGrid->nRows, 1, GDT_Float32, NULL);
	dataset->SetGeoTransform(refGrid->datasetTransform);
	dataset->SetProjection(refGrid->crs);
	GDALRasterBand* band = dataset->GetRasterBand(1);
	band->SetNoDataValue(refGrid->noData);

	noData = band->GetNoDataValue();
	nRows = band->GetYSize();
	nCols = band->GetXSize();

	// set class variables
	dataset->GetGeoTransform(datasetTransform);
	if (GDALInvGeoTransform(datasetTransform, inverseTransform) != CE_None) {
		throw std::invalid_argument("Cannot get inverse transform.");
	}
	crs = dataset->GetProjectionRef();
	noData = band->GetNoDataValue();
	xMin = datasetTransform[0];
	yMax = datasetTransform[3];
	xMax = xMin + datasetTransform[1] * nCols;
	yMin = yMax + datasetTransform[5] * nRows;

};

GDALRasterImage::~GDALRasterImage() {
	// From API docs, why `GDALClose()` and not `~GDALDataset()`: 
	// Equivalent of the C callable GDALClose(). Except that GDALClose() first decrements the reference count, and then closes only if it has dropped to zero.
	// For Windows users, it is not recommended to use the delete operator on the dataset object because of known issues when allocating and freeing memory
	// across module boundaries. Calling GDALClose() is then a better option.
	GDALClose(GDALDataset::ToHandle(dataset));
};

void GDALRasterImage::Close() {
	// From API docs, why `GDALClose()` and not `~GDALDataset()`: 
	// Equivalent of the C callable GDALClose(). Except that GDALClose() first decrements the reference count, and then closes only if it has dropped to zero.
	// For Windows users, it is not recommended to use the delete operator on the dataset object because of known issues when allocating and freeing memory
	// across module boundaries. Calling GDALClose() is then a better option.
	GDALClose(GDALDataset::ToHandle(dataset));
};

void GDALRasterImage::Create(std::string fname) {
	// Create a new dataset with the same extent, transform, and crs as the source
	// TODO: possibly have this create a new GDALRasterImage object?

	GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* outDataset = driver->Create(fname.c_str(), nCols, nRows, 1, GDT_Float32, NULL);
	outDataset->SetGeoTransform(datasetTransform);
	outDataset->SetProjection(crs);
	GDALRasterBand* outBand = outDataset->GetRasterBand(1);
	outBand->SetNoDataValue(noData);
	GByte abyRaster[nCols*nRows];
	outBand->RasterIO(GF_Write, 0, 0, nCols, nRows, abyRaster, nCols, nRows, GDT_Float32, 0, 0);
	GDALClose(GDALDataset::ToHandle(outDataset));
};

float GDALRasterImage::GetVal(int x, int y) {
	float pixelValue;
	band->RasterIO(GF_Read, x, y, 1, 1, &pixelValue, 1, 1, GDT_Float32, 0, 0);
	// check that dataset isn't CE_None
	return pixelValue;

};

std::tuple<int, int> GDALRasterImage::XYfrom(double lat, double lon) {
	int x = static_cast<int>(floor(inverseTransform[0] + inverseTransform[1] * lon + inverseTransform[2] * lat));
	int y = static_cast<int>(floor(inverseTransform[3] + inverseTransform[4] * lon + inverseTransform[5] * lat));
	if ( x < 0 || x > nCols || y < 0 || y > nRows) {
		throw std::invalid_argument("Lat/lon is outside of raster extent.");
	}

	// int32_t pixelValue;
	// assert(GDALRasterIO(band, GF_Read, x, y, 1, 1, &pixelValue, 1, 1, GDT_Int32, 0, 0) == CE_None);
	return std::make_tuple(x, y);
};

bool GDALRasterImage::Exists(std::string fname) {
	// Check if a GDALDataset already exists with the given filename
	const GDALAccess eAccess = GA_ReadOnly;
	GDALDataset* dataset = GDALDataset::FromHandle(GDALOpen(fname.c_str(), eAccess));
	if (dataset) {
		GDALClose(GDALDataset::ToHandle(dataset));
		return true;
	}
	return false;
};

float GDALRasterImage::GetMin() {
	float min;
	float max;
	band->ComputeRasterMinMax(&min, &max);
	return min;
};

float GDALRasterImage::GetMax() {
	float min;
	float max;
	band->ComputeRasterMinMax(&min, &max);
	return max;
};

