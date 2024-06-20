#include <string>
#include <unordered_map>
#include <mutex>
#include "GDALRasterImage.hpp"

#pragma once
class DataOutput {
private:
	GDALRasterImage* refGrid;
	std::string outpath;

	//thread safe wrapper of GDALRasterImage which can:
	// - create a GDALRasterImage
	// - write a value to a specific index
	// - close the GDALRasterImage
	class ImageWrapper {
	private:
		GDALRasterImage* image;
		std::mutex mutex;
	public:
		ImageWrapper(std::string filename, GDALRasterImage* refGrid);
		~ImageWrapper();

		CPLErr setVal(int index, float val, bool hitNODATA);
		void close();
	};

	//map for storing GDALRasterImage wrappers:
	//	we need a lock for our usage of unordered_map could cause race conditions
	//	and unordered_map does not guarantee thread safety.
	std::unordered_map<std::string, ImageWrapper*> images;
	std::mutex imagesMutex;
	ImageWrapper* getImageWrapper(std::string filename);
public:
	DataOutput(GDALRasterImage* refGrid, std::string outpath);
	~DataOutput();
	CPLErr write(int year, int month, std::string name, int index, float val, bool hitNODATA);
};

