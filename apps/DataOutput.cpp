#include "DataOutput.hpp"

DataOutput::ImageBuffer::ImageBuffer(std::string filepath, std::shared_ptr<GDALRasterImage> refGrid) {
	//no synchronization required for the constructor
	this->image = std::make_shared<GDALRasterImage>(filepath, refGrid);

	//for every row index, add an empty std::vector<float> to the rows vector
	for (int i = 0; i < this->image->nRows; i++) {
		this->rows.push_back(nullptr);
	}
}

void DataOutput::ImageBuffer::setVal(int index, float val) {
	std::tuple<int, int> indices = this->image->IndexToXY(index);
	int x = std::get<0>(indices); //x is column index
	int y = std::get<1>(indices); //y is row index

	//initialize the row to the correct size, with nodata
	if (!this->rows[y]) {
		this->rows[y] = std::make_unique<std::vector<float>>(this->image->nCols, static_cast<float>(this->image->noData));
	}

	//write new value to row 
	// this is not synchronized, as only 1 thread should have access to any given row.
	(*this->rows[y])[x] = val;
}

CPLErr DataOutput::ImageBuffer::writeRow(int row) {
	CPLErr retval = CE_None;

	if (!this->rows[row]) {
		//in rare cases, a row may exist that has never had a value set by this class
		//for example, if a startdate raster is given and there is supposed to be an
		//output file for before some of the pixels even start. Those rows would still 
		// need to be written as nodata despite never being set.
		
		//generate a nodata row with std::vector
		std::vector<float> noDataRow(this->image->nCols, static_cast<float>(this->image->noData));

		//write the row (memory will be cleaned up automatically with std::vector)
		retval = this->image->writeRow(row, noDataRow.data());
	}
	else {
		//write the row
		retval = this->image->writeRow(row, this->rows[row]->data());

		//clean up memory that will no longer be used
		this->rows[row].reset();
	}

	return retval;
}

DataOutput::DataOutput(std::shared_ptr<GDALRasterImage> refGrid, std::string outpath) {
	this->refGrid = refGrid;
	this->outpath = outpath;
}

void DataOutput::setVal(int year, int month, std::string name, int index, float val) {

	if (!this->imageBuffers.contains(name)) {
		//acquire mutex since we may be writing to map
		this->imageBuffersMutex.lock();

		//search again, since another thread may have got the mutex first and already created
		//the new entry in the map
		if (!this->imageBuffers.contains(name)) {
			//create a varMap if it doesn't exist
			this->imageBuffers.emplace(name, std::make_unique<varMap>());
		}

		//release lock before continuing
		this->imageBuffersMutex.unlock();
	}

	//get varImages map
	varMap* varImages = this->imageBuffers.at(name).get();

	//search through second layer of map (integer representing year and month)
	//month will never be larger than 4 bits, so left shift year and add month.
	//searchInt will be unique (unless year is astronomically learge)
	int searchInt = (year << 4) + month;

	if (!varImages->contains(searchInt)) {
		//acquire lock since we may be writing to map
		this->imageBuffersMutex.lock();

		//search again, since another thread may have got the mutex first and already created
		//the new entry in the map
		if (!varImages->contains(searchInt)) {
			//create a new ImageBuffer if it doesn't exist
			std::string filepath;
			if (month == -1 && year == -1) {
				filepath = this->outpath + name + ".tif";
			}
			else {
				filepath = this->outpath + name + std::to_string(year) + std::to_string(month) + ".tif";
			}
			varImages->emplace(searchInt, std::make_unique<ImageBuffer>(filepath, this->refGrid));
		}

		//release lock before continuing
		this->imageBuffersMutex.unlock();
	}

	//set image buffer value
	varImages->at(searchInt)->setVal(index, val);
}

CPLErr DataOutput::writeRow(int row) {
	//iterate through every image
	for (auto varImages = this->imageBuffers.begin(); varImages != this->imageBuffers.end(); varImages++) {
		for (auto image = varImages->second->begin(); image != varImages->second->end(); image++) {
			//try to write current images row
			CPLErr error = image->second->writeRow(row);

			//return with error if there is one
			if (error != CE_None) {
				return error;
			}
		}
	}

	//return with no error if we've made it this far
	return CE_None;
}