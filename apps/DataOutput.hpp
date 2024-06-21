#pragma once

#include <string>
#include <unordered_map>
#include <mutex>
#include "GDALRasterImage.hpp"

//a class which should only be defined once that synchronizes and abstracts away the writing
//of pixel values to GDALRasterImage output files.
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
	//	we need a lock since our usage of unordered_map (potential synchronous writes) 
	//	could cause race conditions and unordered_map does not guarantee thread safety.
	std::unordered_map<std::string, ImageWrapper*> images;
	std::mutex imagesMutex;

	//check the images map. Return the image wrapper associated with the filename key if it exists.
	//Otherwise, create a new GDALRasterImage at that filepath and return the image wrapper. 
	ImageWrapper* getImageWrapper(std::string filename);
public:
	DataOutput(GDALRasterImage* refGrid, std::string outpath);
	~DataOutput();

	//determine the filepath of the output given year, month, name.
	//call getImageWrapper() to get the associated wrapper.
	//write the correct val at the correct index using index, val, and hitNODATA.
	CPLErr write(int year, int month, std::string name, int index, float val, bool hitNODATA);
};

