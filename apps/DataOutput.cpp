#include "DataOutput.hpp"

DataOutput::ImageBuffer::ImageBuffer(std::string filepath, RefGridProperties refGrid) {
	//no synchronization required for the constructor
	this->image = std::make_unique<GDALRasterImage>(filepath, refGrid);

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

DataOutput::DataOutput(RefGridProperties& refGrid, std::string outpath, std::unordered_map<std::string, PPPG_OP_VAR> vars) {
	this->vars = vars;
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

int DataOutput::writeOutputGrids(const std::unordered_map<std::string, double>& opVarVals, long cellIndex) {
	//for each possible output variable
	for (auto& [pN, opV] : this->vars) {
		//determine value, name, and tell dataOutput class to write
		float val = (float)(opVarVals.at(pN));
		std::string name = opV.gridName;
		name = name.substr(0, name.find_last_of("."));
		this->setVal(-1, -1, name, cellIndex, val);
	}
	return EXIT_SUCCESS;
}

void DataOutput::writeMonthlyOutputGrids(const std::unordered_map<std::string, double>& opVarVals, int calYear, int calMonth, long cellIndex) {
	//for each possible output variable
	for (auto& [pN, opV] : this->vars) {

		//skip output variable if it is not marked for recurring output
		if (opV.recurYear == -1) {
			continue;
		}

		//skip output variable if it is not marked for recurring output
		if (!opV.recurStart) {
			continue;
		}

		// skip output variable if it is not at the recur interval
		if (((calYear - opV.recurStart) % opV.recurYear) != 0) {
			continue;
		}

		//skip output variable if we're not meant to be printing every month AND we're not on the month we're meant to be printing
		if (!opV.recurMonthly && opV.recurMonth != calMonth) {
			continue;
		}

		//determine value, name, and tell dataOutput class to write
		float val = (float)(opVarVals.at(pN));
		std::string name = opV.gridName;
		name = name.substr(0, name.find_last_of("."));
		this->setVal(calYear, calMonth, name, cellIndex, val);
	}
}