#include <stdexcept>
#include <tuple>
#include <gdal.h>
#include <gdal_priv.h>
#include <math.h>
#include "GDALRasterImage.hpp"


GDALRasterImage::GDALRasterImage(const char* filename) {

	GDALAllRegister();
	const GDALAccess eAccess = GA_ReadOnly;
	dataset = GDALDataset::FromHandle(GDALOpen(filename, eAccess));
	//if (!dataset)
	//{
		//throw std::invalid_argument("File cannot be opened.");
	//}

	// Assume there is only one band in the raster source and use that
	// assert(dataset.GetRasterCount() == 1);
	band = dataset->GetRasterBand(1);

	// Get the inverse geo transform to map from geo location -> pixel location
	double datasetTransform[6] = {};
	// assert(dataset.GetGeoTransform(dataset, datasetTransform) == CE_None);
	// assert(GDALInvGeoTransform(datasetTransform, inverseTransform));

	noData = band->GetNoDataValue();
	nRows = band->GetYSize();
	nCols = band->GetXSize();
};

GDALRasterImage::~GDALRasterImage() {
	// From API docs, why `GDALClose()` and not `~GDALDataset()`: 
	// Equivalent of the C callable GDALClose(). Except that GDALClose() first decrements the reference count, and then closes only if it has dropped to zero.
	// For Windows users, it is not recommended to use the delete operator on the dataset object because of known issues when allocating and freeing memory
	// across module boundaries. Calling GDALClose() is then a better option.
	GDALClose(GDALDataset::ToHandle(dataset));
};

void GDALRasterImage::Write(char* fname) {
	// GDALDataset* newDataset;
	// char** papszOptions = NULL;
	// poDstDS = poDriver->Create(pszDstFilename, 512, 512, 1, GDT_Byte,
		// papszOptions);
	//
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

	// int32_t pixelValue;
	// assert(GDALRasterIO(band, GF_Read, x, y, 1, 1, &pixelValue, 1, 1, GDT_Int32, 0, 0) == CE_None);
	return std::make_tuple(x, y);
};

