#include <tuple>
#include <gdal_priv.h>
#pragma once

// from: https://gis.stackexchange.com/questions/393549
class GDALRasterImage {

	GDALDataset* dataset;
	GDALRasterBand* band;
	double inverseTransform[6];

public:

	double noData{ 0 };
	int nRows { 0 };
	int nCols { 0 };
	GDALRasterImage(const char* filename);
	~GDALRasterImage();
	std::tuple<int, int> XYfrom(double lat, double lon);
	float GetVal(int x, int y);
	void Write(char* fname);

};

