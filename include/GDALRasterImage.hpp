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
	double xMin { 0 };
	double xMax { 0 };
	double yMin { 0 };
	double yMax { 0 };

	GDALRasterImage(std::string filename);
	GDALRasterImage(std::string filename, GDALRasterImage* refGrid);
	~GDALRasterImage();
	std::tuple<int, int> XYfrom(double lat, double lon);
	float GetVal(int x, int y);
	void Create(std::string fname);
	bool Exists(std::string fname);
	float GetMin();
	float GetMax();
	void Close();

};

