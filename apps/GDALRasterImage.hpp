#include <tuple>
#include <string>
#include <mutex>
#include <gdal_priv.h>
#pragma once

// from: https://gis.stackexchange.com/questions/393549
class GDALRasterImage {

	GDALDataset* dataset;
	GDALRasterBand* band;
	double datasetTransform[6];
	double inverseTransform[6];
	const char* crs;
	std::mutex mutex;

public:

	std::string name;
	double noData{ 0 };
	int nRows { 0 };
	int nCols { 0 };
	double xMin { 0 };
	double xMax { 0 };
	double yMin { 0 };
	double yMax { 0 };

	GDALRasterImage(std::string filename);
	GDALRasterImage(std::string filename, std::shared_ptr<GDALRasterImage> refGrid);
	~GDALRasterImage();
	std::tuple<int, int> XYfrom(double lat, double lon);
	int IndexFrom(double lat, double lon);
	std::tuple<int, int> IndexToXY(int index);
	float GetVal(int x, int y);
	float GetVal(int index);
	bool IsNoData(float val);
	void Create(std::string fname);
	bool Exists(std::string fname);
	float GetMin();
	float GetMax();
	void Close();
	std::vector<std::pair<int, int>> getIndicesWhere(const double& value);
	double minFromIndices(const std::vector<std::pair<int, int>>& indices);
	double maxFromIndices(const std::vector<std::pair<int, int>>& indices);
	CPLErr writeRow(int row, float* buffer);

};

