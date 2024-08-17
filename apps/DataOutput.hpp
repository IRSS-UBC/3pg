#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <mutex>
#include "GDALRasterImage.hpp"
#include "ParamStructs.hpp"
#include "MYDate.h"
#include "util.hpp"

//a class which should only be defined once that synchronizes and abstracts away the writing
//of pixel values to GDALRasterImage output files.
class DataOutput {
private:
	RefGridProperties refGrid;
	std::string outpath;
	std::unordered_map<std::string, PPPG_OP_VAR> vars;

	//thread safe wrapper of GDALRasterImage which can:
	// - create a GDALRasterImage
	// - write a value to a specific index
	// - close the GDALRasterImage
	class ImageBuffer {
	private:
		std::unique_ptr<GDALRasterImage> image;
		std::vector<std::unique_ptr<std::vector<float>>> rows;
	public:
		ImageBuffer(std::string filename, RefGridProperties refGrid);

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
	DataOutput(RefGridProperties& refGrid, std::string outpath, std::unordered_map<std::string, PPPG_OP_VAR> vars);

	//determine the filepath of the output given year, month, name.
	//call getImageBuffer() to get the associated wrapper.
	//write the correct val at the correct index using index, val.
	void setVal(int year, int month, std::string name, int index, float val);

	int writeOutputGrids(const std::unordered_map<std::string, double>& opVarVals, long cellIndex);
	void writeMonthlyOutputGrids(const std::unordered_map<std::string, double>& opVarVals, int calYear, int calMonth, MYDate minMY, MYDate maxMY, long cellIndex);

	CPLErr writeRow(int row);
};

