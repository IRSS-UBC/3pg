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
	class ImageBuffer {
	private:
		GDALRasterImage* image;
		std::vector<std::vector<float>*> rows;
	public:
		ImageBuffer(std::string filename, GDALRasterImage* refGrid);
		~ImageBuffer();

		void setVal(int index, float val);
		CPLErr writeRow(int row);
		void close();
	};

	//map for storing GDALRasterImage wrappers:
	//	we need a lock since our usage of unordered_map (potential synchronous writes) 
	//	could cause race conditions and unordered_map does not guarantee thread safety.
	std::unordered_map<std::string, std::unordered_map<int, ImageBuffer*>*> imageBuffers;
	std::mutex imageBuffersMutex;

	//check the images map. Return the image wrapper associated with the filename key if it exists.
	//Otherwise, create a new GDALRasterImage at that filepath and return the image wrapper. 
	ImageBuffer* getImageBuffer(std::string filename, int year, int month);
public:
	DataOutput(GDALRasterImage* refGrid, std::string outpath);
	~DataOutput();

	//determine the filepath of the output given year, month, name.
	//call getImageBuffer() to get the associated wrapper.
	//write the correct val at the correct index using index, val, and hitNODATA.
	void setVal(int year, int month, std::string name, int index, float val);
	CPLErr writeRow(int row);
};

