#include "DataOutput.hpp"

DataOutput::ImageBuffer::ImageBuffer(std::string filepath, std::shared_ptr<GDALRasterImage> refGrid) {
	//no synchronization required for the constructor
	this->image = std::make_shared<GDALRasterImage>(filepath, refGrid);

	//for every row index, add an empty std::vector<float> to the rows vector
	for (int i = 0; i < this->image->nRows; i++) {
		this->rows.push_back(nullptr);
	}
}

DataOutput::ImageBuffer::~ImageBuffer() {
	this->close();
}

void DataOutput::ImageBuffer::setVal(int index, float val) {
	std::tuple<int, int> indices = this->image->IndexToXY(index);
	int x = std::get<0>(indices); //x is column index
	int y = std::get<1>(indices); //y is row index
	std::vector<float>* row = this->rows[y];

	//initialize the row to the correct size, with nodata
	if (row == nullptr) {
		row = new std::vector<float>(this->image->nCols, this->image->noData);
		this->rows[y] = row;
	}

	//write new value to row 
	// this is not synchronized, as only 1 thread should have access to any given row.
	(*row)[x] = val;
}

CPLErr DataOutput::ImageBuffer::writeRow(int row) {
	CPLErr retval = CE_None;
	std::vector<float>* imageRow = this->rows[row];

	//in rare cases, a row may exist that has never had a value set by this class
	//for example, if a startdate raster is given and there is supposed to be an
	//output file for before some of the pixels even start. Those rows would still 
	// need to be written as nodata despite never being set.
	if (imageRow != nullptr) {
		//write the row
		retval = this->image->writeRow(row, imageRow->data());

		//clean up memory that will no longer be used
		delete this->rows[row];
		this->rows[row] = nullptr;
	}
	else {
		//generate a nodata row with std::vector
		std::vector<float> noDataRow(this->image->nCols, this->image->noData);

		//write the row (memory will be cleaned up automatically with std::vector)
		retval = this->image->writeRow(row, noDataRow.data());
	}

	return retval;
}

void DataOutput::ImageBuffer::close() {	
	//remove all dynamically allocated vectors
	for (int i = 0; i < this->rows.size(); i++) {
		if (this->rows[i] != nullptr) {
			delete this->rows[i];
			this->rows[i] = nullptr;
		}
	}
}

DataOutput::DataOutput(std::shared_ptr<GDALRasterImage> refGrid, std::string outpath) {
	this->refGrid = refGrid;
	this->outpath = outpath;
}

DataOutput::~DataOutput() {
	//delete all of the image allocations we made
	for (auto imageMap = this->imageBuffers.begin(); imageMap != this->imageBuffers.end(); imageMap++) {
		for (auto image = imageMap->second->begin(); image != imageMap->second->end(); image++) {
			image->second->close();
			delete image->second;
			image->second = nullptr;
		}
	}
}

DataOutput::ImageBuffer* DataOutput::getImageBuffer(std::string var, int year, int month) {

	ImageBuffer* retval = nullptr;
	std::unordered_map<int, ImageBuffer*>* varImages;

	if (this->imageBuffers.contains(var)) {
		//set varImages if we've found an entry
		varImages = this->imageBuffers.at(var);
	}
	else {
		//acquire mutex since we may be writing to map
		this->imageBuffersMutex.lock();

		//search again, since another thread may have got the mutex first and already created
		//the new entry in the map
		if (this->imageBuffers.contains(var)) {
			//set varImages if we've found an entry
			varImages = this->imageBuffers.at(var);
		}
		else {
			//create a varImages map if it doesn't exist
			varImages = new std::unordered_map<int, ImageBuffer*>;
			this->imageBuffers.emplace(var, varImages);
		}

		//release lock before continuing
		this->imageBuffersMutex.unlock();
	}

	//search through second layer of map (integer representing year and month)
	//month will never be larger than 4 bits, so left shift year and add month.
	//searchInt will be unique (unless year is astronomically learge)
	int searchInt = (year << 4) + month;

	if (varImages->contains(searchInt)) {
		//set Images if we've found an entry
		retval = varImages->at(searchInt);
	}
	else {
		//acquire lock since we may be writing to map
		this->imageBuffersMutex.lock();

		//search again, since another thread may have got the mutex first and already created
		//the new entry in the map
		if (varImages->contains(searchInt)) {
			//set Images if we've found an entry
			retval = varImages->at(searchInt);
		}
		else {
			//create a new ImageBuffer if it doesn't exist
			std::string filepath;
			if (month == -1 && year == -1) {
				filepath = this->outpath + var + ".tif";
			}
			else {
				filepath = this->outpath + var + std::to_string(year) + std::to_string(month) + ".tif";
			}
			retval = new ImageBuffer(filepath, this->refGrid);
			varImages->emplace(searchInt, retval);
		}

		//release lock before continuing
		this->imageBuffersMutex.unlock();
	}

	return retval;
}

void DataOutput::setVal(int year, int month, std::string name, int index, float val) {

	//retrieve wrapper image from map
	ImageBuffer* image = this->getImageBuffer(name, year, month);

	//set the value of the particular image
	image->setVal(index, val);
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