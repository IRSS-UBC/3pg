#include <tuple>
#include <string>
#include <gdal_priv.h>
#pragma once

// from: https://gis.stackexchange.com/questions/393549
class GDALRasterImage {

	GDALDataset* dataset;
	GDALRasterBand* band;
	double datasetTransform[6];
	double inverseTransform[6];
	const char* crs;

public:

	double noData{ 0 };
	int nRows { 0 };
	int nCols { 0 };
	int xMin { 0 };
	int xMax { 0 };
	int yMin { 0 };
	int yMax { 0 };

	GDALRasterImage(std::string filename);
	~GDALRasterImage();
	std::tuple<int, int> XYfrom(double lat, double lon);
	float GetVal(int x, int y);
	void Create(std::string fname);
	bool Exists(std::string fname);

};

