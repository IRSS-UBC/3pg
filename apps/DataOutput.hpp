#pragma once

#include <string>
#include <unordered_map>
#include <mutex>
#include "GDALRasterImage.hpp"

//a class which should only be defined once that synchronizes and abstracts away the writing
//of pixel values to GDALRasterImage output files.
class DataOutput {
private:
	std::shared_ptr<GDALRasterImage> refGrid;
	std::string outpath;

	//thread safe wrapper of GDALRasterImage which can:
	// - create a GDALRasterImage
	// - write a value to a specific index
	// - close the GDALRasterImage
	class ImageBuffer {
	private:
		std::shared_ptr<GDALRasterImage> image;
		std::vector<std::unique_ptr<std::vector<float>>> rows;
	public:
		ImageBuffer(std::string filename, std::shared_ptr<GDALRasterImage> refGrid);

		void setVal(int index, float val);
		CPLErr writeRow(int row);
	};

	//map for storing GDALRasterImage wrappers:
	//	we need a lock since our usage of unordered_map (potential synchronous writes) 
	//	could cause race conditions and unordered_map does not guarantee thread safety.
	typedef std::unordered_map<int, std::unique_ptr<ImageBuffer>> varMap;
	std::unordered_map<std::string, std::unique_ptr<varMap>> imageBuffers;
	std::mutex imageBuffersMutex;

public:
	DataOutput(std::shared_ptr<GDALRasterImage> refGrid, std::string outpath);

	//determine the filepath of the output given year, month, name.
	//call getImageBuffer() to get the associated wrapper.
	//write the correct val at the correct index using index, val, and hitNODATA.
	void setVal(int year, int month, std::string name, int index, float val);
	CPLErr writeRow(int row);
};

