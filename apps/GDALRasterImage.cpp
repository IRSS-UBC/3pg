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

GDALRasterImage::GDALRasterImage(std::string filename, RefGridProperties& refGrid) {
	// Create a new GDALRasterImage dataset with one band with the same extent, transform, and crs as refGrid
	GDALAllRegister();
	CPLPushErrorHandler(CPLQuietErrorHandler); // suppress error messages that Exists() throws for non-existent files
	if (Exists(filename)) {
		throw std::invalid_argument("File already exists.");
	};
	GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	dataset = driver->Create(filename.c_str(), refGrid.nCols, refGrid.nRows, 1, GDT_Float32, NULL);
	name = filename;
	dataset->SetGeoTransform(refGrid.datasetTransform);
	dataset->SetProjection(refGrid.crs);
	band = dataset->GetRasterBand(1);
	band->SetNoDataValue(refGrid.noData);

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
	//exception thrown when this isn't locked. I *assume* that's because there's potential race conditions and GDAL is smart
	//enough to warn us of them in debug mode.
	this->mutex.lock();
	if (band->RasterIO(GF_Read, x, y, 1, 1, &pixelValue, 1, 1, GDT_Float32, 0, 0) != CE_None) {
		throw std::invalid_argument("Cannot read pixel value.");
	}
	this->mutex.unlock();
	return pixelValue;
};

float GDALRasterImage::GetVal(int index) {
	// Get the value of the pixel at the given index
	std::tuple<int, int> xy = IndexToXY(index);
	return GDALRasterImage::GetVal(std::get<0>(xy), std::get<1>(xy));
};

std::tuple<int,int> GDALRasterImage::IndexToXY(int index) {
	// Get the x,y coordinates of the pixel at the given index
	int x = static_cast<int>(index % nCols);
	int y = static_cast<int>(index / nCols);
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
	catch (std::exception&) {
		return false;
	}
	return false;
};

bool GDALRasterImage::IsNoData(float val) {
	return std::isnan(val) || val == this->noData;
};

 float GDALRasterImage::GetMin() {
	float min = std::numeric_limits<float>::max();
	for (int i = 0; i < this->nRows; i++) {
		for (int j = 0; j < this->nCols; j++) {
			float check = this->GetVal(i, j);
			if (!this->IsNoData(check)) {
				min = (min > check) ? check : min;
			}
		}
	}
 	return min;
 };

 float GDALRasterImage::GetMax() {
	 float max = std::numeric_limits<float>::min();
	 for (int i = 0; i < this->nRows; i++) {
		 for (int j = 0; j < this->nCols; j++) {
			 float check = this->GetVal(i, j);
			 if (!this->IsNoData(check)) {
				 max = (max < check) ? check : max;
			 }
		 }
	 }
	 return max;
 };

 CPLErr GDALRasterImage::writeRow(int row, float* buffer) {
	 //for info on how RasterIO works see:
	 //https://gdal.org/api/gdalrasterband_cpp.html#_CPPv4N14GDALRasterBand8RasterIOE10GDALRWFlagiiiiPvii12GDALDataType8GSpacing8GSpacingP20GDALRasterIOExtraArg

	 this->mutex.lock();
	 CPLErr retval = band->RasterIO(
		 GF_Write,		//eRWFlag: Either GF_Read or GF_Write
		 0,				//column index
		 row,			//row index
		 nCols,			//number of columns we're writing (all of them)
		 1,				//number of rows we're writing
		 buffer,		//data buffer
		 nCols,			//number of columns in the data buffer
		 1,				//number of rows in the data buffer
		 GDT_Float32,	//buffer type
		 0,				//byte offset between scanlines in buffer. 0 automatically sets to default of eBufType * nBufXSize
		 0				//extra arguments
	 );

	 this->mutex.unlock();
	 return retval;
 }

 RefGridProperties GDALRasterImage::getRefGrid() {
	 RefGridProperties retval;

	 retval.name = this->name;
	 retval.nRows = this->nRows;
	 retval.nCols = this->nCols;
	 retval.xMin = this->xMin;
	 retval.xMax = this->xMax;
	 retval.yMin = this->yMin;
	 retval.yMax = this->yMax;
	 retval.noData = this->noData;
	 retval.crs = this->crs;
	 retval.name = this->name;

	 retval.datasetTransform[0] = this->datasetTransform[0];
	 retval.datasetTransform[1] = this->datasetTransform[1];
	 retval.datasetTransform[2] = this->datasetTransform[2];
	 retval.datasetTransform[3] = this->datasetTransform[3];
	 retval.datasetTransform[4] = this->datasetTransform[4];
	 retval.datasetTransform[5] = this->datasetTransform[5];

	 return retval;
 }
