#include <stdexcept>
#include <tuple>
#include <gdal_priv.h>
#include <gdal.h>
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

	name = filename;
	// Assume there is only one band in the raster source and use that
	if (dataset->GetRasterCount() != 1) {
		GDALClose(GDALDataset::ToHandle(dataset));
		throw std::invalid_argument("TIF must have only one band.");
	}
	band = dataset->GetRasterBand(1);

	// Get the inverse geo transform to map from geo location -> pixel location
	if (dataset->GetGeoTransform(datasetTransform) == CE_Failure) {
		GDALClose(GDALDataset::ToHandle(dataset));
		throw std::runtime_error("Cannot get transform.");
	}
	if (GDALInvGeoTransform(datasetTransform, inverseTransform) == false) {
		GDALClose(GDALDataset::ToHandle(dataset));
		throw std::runtime_error("Cannot get inverse transform.");
	}

	// get crs
	crs = dataset->GetProjectionRef();
	if (crs == NULL) {
		GDALClose(GDALDataset::ToHandle(dataset));
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
	GDALAllRegister();
	CPLPushErrorHandler(CPLQuietErrorHandler); // suppress error messages that Exists() throws for non-existent files
	if (Exists(filename)) {
		throw std::invalid_argument("File already exists.");
	};
	GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	dataset = driver->Create(filename.c_str(), refGrid->nCols, refGrid->nRows, 1, GDT_Float32, NULL);
	name = filename;
	dataset->SetGeoTransform(refGrid->datasetTransform);
	dataset->SetProjection(refGrid->crs);
	band = dataset->GetRasterBand(1);
	band->SetNoDataValue(refGrid->noData);

	noData = band->GetNoDataValue();
	nRows = band->GetYSize();
	nCols = band->GetXSize();

	// set class variables
	dataset->GetGeoTransform(datasetTransform);
	if (GDALInvGeoTransform(datasetTransform, inverseTransform) == false) {
		GDALClose(GDALDataset::ToHandle(dataset));
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

float GDALRasterImage::GetVal(int x, int y) {
	float pixelValue;
	if (band->RasterIO(GF_Read, x, y, 1, 1, &pixelValue, 1, 1, GDT_Float32, 0, 0) != CE_None) {
		throw std::invalid_argument("Cannot read pixel value.");
	}
	return pixelValue;
};

float GDALRasterImage::GetVal(int index) {
	// Get the value of the pixel at the given index
	std::tuple<int, int> xy = IndexToXY(index);
	float pixelValue;

	if (band->RasterIO(GF_Read, std::get<0>(xy), std::get<1>(xy), 1, 1, &pixelValue, 1, 1, GDT_Float32, 0, 0) != CE_None) {
		throw std::invalid_argument("Cannot read pixel value.");
	}
	return pixelValue;
};

CPLErr GDALRasterImage::SetVal(int x, int y, float val) {
	// Set the value of the pixel at the given x,y coordinates
	return band->RasterIO(GF_Write, x, y, 1, 1, &val, 1, 1, GDT_Float32, 0, 0);
};

CPLErr GDALRasterImage::SetVal(int index, float val) {
	// Set the value of the pixel at the given index
	std::tuple<int, int> xy = IndexToXY(index);
	//std::cout << "From " << index << " Setting pixel value at " << std::get<0>(xy) << ", " << std::get<1>(xy) << " to " << val << std::endl;
	return band->RasterIO(GF_Write, std::get<0>(xy), std::get<1>(xy), 1, 1, &val, 1, 1, GDT_Float32, 0, 0);
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

int GDALRasterImage::IndexFrom(double lat, double lon) {
	// Get the index of the pixel at the given lat/lon
	int x = static_cast<int>(floor(inverseTransform[0] + inverseTransform[1] * lon + inverseTransform[2] * lat));
	int y = static_cast<int>(floor(inverseTransform[3] + inverseTransform[4] * lon + inverseTransform[5] * lat));
	if ( x < 0 || x > nCols || y < 0 || y > nRows) {
		throw std::invalid_argument("Lat/lon is outside of raster extent.");
	}
	//std::cout << "From " << lat << ", " << lon << " got " << x << ", " << y << std::endl;
	return y * nCols + x;
};

std::tuple<int,int> GDALRasterImage::IndexToXY(int index) {
	// Get the x,y coordinates of the pixel at the given index
	int x = index % nCols;
	int y = index / nCols;
	return std::make_tuple(x, y);
};

bool GDALRasterImage::Exists(std::string fname) {
	// Check if a GDALDataset already exists with the given filename
	const GDALAccess eAccess = GA_ReadOnly;
	try {
		GDALDataset* dataset = GDALDataset::FromHandle(GDALOpen(fname.c_str(), eAccess));
		if (dataset) {
			GDALClose(GDALDataset::ToHandle(dataset));
			return true;
		}
	}
	catch (std::exception& e) {
		return false;
	}
	return false;
};

bool GDALRasterImage::IsNoData(float val) {
	if (std::isnan(noData)) {
		return (std::isnan(val));
	}
	else {
		return (val == noData);
	}
};

// float GDALRasterImage::GetMin() {
// 	float min;
// 	float max;
// 	band->ComputeRasterMinMax(&min, &max);
// 	return min;
// };

// float GDALRasterImage::GetMax() {
// 	float min;
// 	float max;
// 	band->ComputeRasterMinMax(&min, &max);
// 	return max;
// };

