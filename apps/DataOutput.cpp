#include "DataOutput.hpp"

DataOutput::ImageWrapper::ImageWrapper(std::string filename, GDALRasterImage* refGrid) {
	//no synchronization required for the constructor
	this->image = new GDALRasterImage(filename, refGrid);
}

DataOutput::ImageWrapper::~ImageWrapper() {
	if (this->image != nullptr) {
		this->close();
	}
}

CPLErr DataOutput::ImageWrapper::setVal(int index, float val, bool hitNODATA) {
	//acquire lock
	this->mutex.lock();
	CPLErr retval;

	//write value
	float inVal = hitNODATA ? this->image->noData : val;
	retval = this->image->SetVal(index, inVal);

	//release lock
	this->mutex.unlock();
	return retval;
}

void DataOutput::ImageWrapper::close() {
	//acquire lock
	this->mutex.lock();

	//close and delete image
	this->image->Close();
	delete this->image;
	this->image = nullptr;

	//release lock
	this->mutex.unlock();
}

DataOutput::DataOutput(GDALRasterImage* refGrid, std::string outpath) {
	this->refGrid = refGrid;
	this->outpath = outpath;
}

DataOutput::~DataOutput() {
	//delete all the allocations we made with new
	for (auto image = this->images.begin(); image != this->images.end(); image++) {
		image->second->close();
		delete image->second;
		image->second = nullptr;
	}
}

DataOutput::ImageWrapper* DataOutput::getImageWrapper(std::string filename) {
	this->imagesMutex.lock();
	ImageWrapper* retval = nullptr;

	//search for the filename in the map
	auto search = this->images.find(filename);

	if (search != this->images.end()) {
		//if an image exists with that filename, return the image
		retval = search->second;
	} 
	else {
		//otherwise, create an image with that filename and add it to the images map
		retval = new ImageWrapper(filename, this->refGrid);
		this->images.emplace(filename, retval);
	}

	//release the lock then return the image
	this->imagesMutex.unlock();
	return retval;
}

CPLErr DataOutput::write(int year, int month, std::string name, int index, float val, bool hitNODATA) {
	std::string filepath;
	std::string filename;

	//determine filepath from year, month, and variable name
	if (year == -1 && month == -1) {
		filename = name;
	}
	else {
		filename = name + std::to_string(year) + std::to_string(month);
	}
	filepath = this->outpath + filename + ".tif";

	//retrieve wrapper image from map
	ImageWrapper* image = this->getImageWrapper(filepath);

	//set the value of the particular image
	return image->setVal(index, val, hitNODATA);
}