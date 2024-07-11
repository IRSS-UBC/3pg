#include "DataInput.hpp"

DataInput::DataInput(Logger& logger) {
	this->logger = &logger;
	this->refGrid = nullptr;
}

DataInput::~DataInput() {
	//check each parameter
	for (auto iterator = this->inputParams.begin(); iterator != this->inputParams.end(); iterator++) {
		PPPG_PARAM param = iterator->second;

		//if it's a grid parameter, we have to clean up it's corrosponding GDALRasterImage
		if (param.data.spType == pTif && param.data.g != nullptr) {
			delete param.data.g;
			param.data.g = nullptr;
		}
	}
}

bool DataInput::getScalar(std::vector<std::string> value, PPPG_PARAM& param) {
	try {
		//get the value from the array
		double val = std::stod(value.front());
		*(param.adr) = val;

		//mark the data as scalar
		param.data.spType = pScalar;

		//log input param acquired
		string output = "    " + param.id + "        constant: " + to_string(*(param.adr));
		std::cout << output << std::endl;
		this->logger->Log(output);

		return true;
	}
	catch (std::invalid_argument const& e) {
		//if we threw an error, it's not a double value so return false
		return false;
	}
}

bool DataInput::getGrid(std::vector<std::string> value, PPPG_PARAM& param) {
	try {
		//get the file path
		std::filesystem::path filePath = value.front();

		//ensure the file extension is .tif
		if (filePath.extension() != ".tif") {
			std::string errstr = filePath.filename().generic_string() + " is an invalid file type. File extension must be '.tif'";
			std::cout << errstr << std::endl;
			this->logger->Log(errstr);
			return false;
		}

		//ensure the filepath exists
		if (!std::filesystem::exists(filePath)) {
			std::string errstr = filePath.string() + " does not exist.";
			std::cout << errstr << std::endl;
			this->logger->Log(errstr);
			return false;
		}

		//set data type and file path
		param.data.spType = pTif;
		param.data.gridName = filePath.string();

		//try opening the grid and compare to the refgrid (or create the refgrid)
		if (!openCheckGrid(param.data)) {
			return false;
		}

		//log input param
		string output = "    " + param.id + "        raster: " + param.data.gridName;
		std::cout << output << std::endl;
		this->logger->Log(output);
		return true;
	}
	catch (std::filesystem::filesystem_error const& e) {
		//set and print/log error string
		std::string errstr = " " + value.front() + " could not be interpreted as a scalar or grid name";
		std::cout << errstr << std::endl;
		this->logger->Log(errstr);
		this->logger->Log(e.what());
		return false;
	}
}

bool DataInput::openCheckGrid(PPPG_VVAL& vval) {
	//ensure we can open and read from the file as a GDALRasterImage
	try {
		vval.g = new GDALRasterImage(vval.gridName);
	}
	catch (const std::exception& e) {
		std::string errstr = "failed to open " + vval.gridName + "\n" + e.what();
		std::cout << errstr << std::endl;
		this->logger->Log(errstr);
		return false;
	}

	//if the refgrid is null, this is grid becomes the refgrid
	//otherwise, check grid dimensions
	if (this->refGrid == nullptr) {
		this->refGrid = vval.g;
	}
	else {
		if (
			(fabs(this->refGrid->xMin - vval.g->xMin) > 0.0001) ||	//check xMin
			(fabs(this->refGrid->yMin - vval.g->yMin) > 0.0001) ||	//check yMin
			(fabs(this->refGrid->xMax - vval.g->xMax) > 0.0001) ||	//check xMax
			(fabs(this->refGrid->yMax - vval.g->yMax) > 0.0001) ||	//check yMax
			(this->refGrid->nRows != vval.g->nRows) ||				//check nRows
			(this->refGrid->nCols != vval.g->nCols)					//check nCols
			) {
			std::string errstr = "Grid dimensions of " + vval.gridName + " differs from " + this->refGrid->name;
			std::cout << errstr << std::endl;
			this->logger->Log(errstr);
			return false;
		}
	}

	return true;
}

bool DataInput::tryAddParam(std::string name, std::vector<std::string> value) {
	//TODO: actual 3pg doesn't care about case (I don't think), make sure that is reflected here
	PPPG_PARAM param;

	//if the string isn't the exact parameter name, see if it is in the parameter name map
	if (this->allParams.find(name) == this->allParams.end()) {
		auto search = this->paramNames.find(name);

		if (search == this->paramNames.end()) {
			return false;
		}

		param.id = search->second;
	}
	else {
		param.id = name;
	}

	//soilIndex special case
	if (param.id == "soilIndex") {
		if ("S" == value.front()) {
			value[0] = "1";
		}
		else if ("SL" == value.front()) {
			value[0] = "2";
		}
		else if ("CL" == value.front()) {
			value[0] = "3";
		}
		else if ("C" == value.front()) {
			value[0] = "4";
		}
	}

	//try to get the value as a scalar
	if (DataInput::getScalar(value, param)) {
		//if gotten, add to inputParams map
		this->inputParams.emplace(param.id, param);

		//remove from required params set
		this->requiredParams3PGS.erase(param.id);
		this->requiredParams3PG.erase(param.id);

		return true;
	}
	
	//try to get the value as a grid
	if (DataInput::getGrid(value, param)) {
		//if gotten, add to inputParams map and return
		this->inputParams.emplace(param.id, param);

		//remove from required params set
		this->requiredParams3PGS.erase(param.id);
		this->requiredParams3PG.erase(param.id);

		return true;
	}

	//if we've made it here, the user intended it as an input param but
	//we were unable to add it as one. Fail.
	exit(EXIT_FAILURE);
}

bool DataInput::inputFinished(bool modelMode3PGS) {
	bool haveSeedlingMass = this->inputParams.find("SeedlingMass") != this->inputParams.end();
	bool haveWFi = this->inputParams.find("WFi") != this->inputParams.end();
	bool haveWRi = this->inputParams.find("WRi") != this->inputParams.end();
	bool haveWSi = this->inputParams.find("WSi") != this->inputParams.end();

	//must have seedling mass and one of WFi, WRi, and WSi
	if (!haveSeedlingMass && (!haveWFi || !haveWRi || !haveWSi)) {
		std::string errstr;
		if (!haveSeedlingMass) {
			errstr = "Missing parameter: SeedlingMass.";
		}
		else {
			errstr = "Missing parameter: at least one of WFi, WRi, and WSi required.";
		}
		std::cout << errstr << std::endl;
		this->logger->Log(errstr);
		return false;
	}

	//note: we may not need the following if else statement. The previous code mentions it seems to cause
	//errors when not there. Something to look into (TODO)
	if (haveSeedlingMass) {
		//set WFi, WRi, and WSi if seedling mass is set
		*inputParams["WFi"].adr = 0;
		inputParams["WFi"].data.spType = pScalar;
		*inputParams["WRi"].adr = 0;
		inputParams["WRi"].data.spType = pScalar;
		*inputParams["WSi"].adr = 0;
		inputParams["WSi"].data.spType = pScalar;
	}
	else {
		//set SeedlingMass if it isn't already
		*inputParams["SeedlingMass"].adr = 0;
		inputParams["SeedlingMass"].data.spType = pScalar;
	}

	//if we're using 3PGS, ensure we have all required 3PGS parameters
	if (modelMode3PGS && this->requiredParams3PGS.size() != 0) {
		return false;
	}

	//if we're using 3PG (not 3PGS), ensure we have all required 3PG parameters
	if (!modelMode3PGS && this->requiredParams3PG.size() != 0) {
		return false;
	}

	finishedInput = true;
	return true;
}

InputParams DataInput::getInputParams(int row, int col) {
	//error checking
	if (!finishedInput) {
		throw std::exception("should NOT be able to call getInputParams() if input failed!");
	}

	//declare, copy values into, then return an instance of InputParams
	InputParams params;
	params.pFS2 = DataInput::getValFromParam("pFS2", row, col);
	params.pFS20 = DataInput::getValFromParam("pFS20", row, col);
	params.StemConst = DataInput::getValFromParam("StemConst", row, col);
	params.StemPower = DataInput::getValFromParam("StemPower", row, col);
	params.pRx = DataInput::getValFromParam("pRx", row, col);
	params.pRn = DataInput::getValFromParam("pRn", row, col);
	params.growthTmin = DataInput::getValFromParam("growthTmin", row, col);
	params.growthTopt = DataInput::getValFromParam("growthTopt", row, col);
	params.growthTmax = DataInput::getValFromParam("growthTmax", row, col);
	params.kF = DataInput::getValFromParam("kF", row, col);
	params.gammaFx = DataInput::getValFromParam("gammaFx", row, col);
	params.gammaF0 = DataInput::getValFromParam("gammaF0", row, col);
	params.tgammaF = DataInput::getValFromParam("tgammaF", row, col);
	params.Rttover = DataInput::getValFromParam("Rttover", row, col);
	params.MaxCond = DataInput::getValFromParam("MaxCond", row, col);
	params.CoeffCond = DataInput::getValFromParam("CoeffCond", row, col);
	params.BLcond = DataInput::getValFromParam("BLcond", row, col);
	params.m0 = DataInput::getValFromParam("m0", row, col);
	params.fN0 = DataInput::getValFromParam("fN0", row, col);
	params.fNn = DataInput::getValFromParam("fNn", row, col);
	params.thinPower = DataInput::getValFromParam("thinPower", row, col);
	params.mF = DataInput::getValFromParam("mF", row, col);
	params.mR = DataInput::getValFromParam("mR", row, col);
	params.mS = DataInput::getValFromParam("mS", row, col);
	params.SWconst0 = DataInput::getValFromParam("SWconst0", row, col);
	params.SWpower0 = DataInput::getValFromParam("SWpower0", row, col);
	params.wSx1000 = DataInput::getValFromParam("wSx1000", row, col);
	params.MaxAge = DataInput::getValFromParam("MaxAge", row, col);
	params.nAge = DataInput::getValFromParam("nAge", row, col);
	params.rAge = DataInput::getValFromParam("rAge", row, col);
	params.SLA0 = DataInput::getValFromParam("SLA0", row, col);
	params.SLA1 = DataInput::getValFromParam("SLA1", row, col);
	params.tSLA = DataInput::getValFromParam("tSLA", row, col);
	params.k = DataInput::getValFromParam("k", row, col);
	params.fullCanAge = DataInput::getValFromParam("fullCanAge", row, col);
	params.alpha = DataInput::getValFromParam("alpha", row, col);
	params.fracBB0 = DataInput::getValFromParam("fracBB0", row, col);
	params.fracBB1 = DataInput::getValFromParam("fracBB1", row, col);
	params.tBB = DataInput::getValFromParam("tBB", row, col);
	params.y = DataInput::getValFromParam("y", row, col);
	params.rhoMin = DataInput::getValFromParam("rhoMin", row, col);
	params.rhoMax = DataInput::getValFromParam("rhoMax", row, col);
	params.tRho = DataInput::getValFromParam("tRho", row, col);
	params.Qa = DataInput::getValFromParam("Qa", row, col);
	params.Qb = DataInput::getValFromParam("Qb", row, col);
	params.gDM_mol = DataInput::getValFromParam("gDM_mol", row, col);
	params.molPAR_MJ = DataInput::getValFromParam("molPAR_MJ", row, col);
	params.LAIgcx = DataInput::getValFromParam("LAIgcx", row, col);
	params.MaxIntcptn = DataInput::getValFromParam("MaxIntcptn", row, col);
	params.LAImaxIntcptn = DataInput::getValFromParam("LAImaxIntcptn", row, col);
	params.Lat = DataInput::getValFromParam("Lat", row, col);
	params.FRp = DataInput::getValFromParam("FRp", row, col);
	params.FRstart = DataInput::getValFromParam("FRstart", row, col);
	params.FRend = DataInput::getValFromParam("FRend", row, col);
	params.FRdec = DataInput::getValFromParam("FRdec", row, col);
	params.soilIndex = DataInput::getValFromParam("soilIndex", row, col);
	params.MaxASW = DataInput::getValFromParam("MaxASW", row, col);
	params.MinASWp = DataInput::getValFromParam("MinASWp", row, col);
	params.StartAge = DataInput::getValFromParam("StartAge", row, col);
	params.StartMonth = DataInput::getValFromParam("StartMonth", row, col);
	params.yearPlanted = DataInput::getValFromParam("yearPlanted", row, col);
	params.SeedlingMass = DataInput::getValFromParam("SeedlingMass", row, col);
	params.WFi = DataInput::getValFromParam("WFi", row, col);
	params.WRi = DataInput::getValFromParam("WRi", row, col);
	params.WSi = DataInput::getValFromParam("WSi", row, col);
	params.StemNoi = DataInput::getValFromParam("StemNoi", row, col);
	params.ASWi = DataInput::getValFromParam("ASWi", row, col);
	params.MinASWTG = DataInput::getValFromParam("MinASWTG", row, col);
	params.NDVI_FPAR_intercept = DataInput::getValFromParam("NDVI_FPAR_intercept", row, col);
	params.NDVI_FPAR_constant = DataInput::getValFromParam("NDVI_FPAR_constant", row, col);
	return params;
}

double DataInput::getValFromParam(std::string paramName, int row, int col) {
	auto search = this->inputParams.find(paramName);

	//if the param doesn't exist, set to -1
	if (search == this->inputParams.end()) {
		return -1.0;
	}

	//if the param is scalar, return it's value
	if (search->second.data.spType == pScalar) {
		return *search->second.adr;
	}

	//if the param is a grid, return the value at the row and column specified
	if (search->second.data.spType == pTif) {
		return search->second.data.g->GetVal(row, col);
	}

	throw std::exception("a parameter has been set incorrectly as neither a scalar or a grid.");
}

void DataInput::findRunPeriod(MYDate& minMY, MYDate& maxMY) {
	if (!finishedInput) {
		throw std::exception("should NOT be able to call findRunPeriod() if input failed!");
	}

	//get required parameters from parameter map
	PPPG_PARAM yearPlantedParam = this->inputParams.find("yearPlanted")->second;
	PPPG_PARAM startAgeParam = this->inputParams.find("StartAge")->second;
	PPPG_PARAM endYearParam = this->inputParams.find("EndYear")->second;
	PPPG_PARAM startMonthParam = this->inputParams.find("startMonth")->second;

	//get maxes and mins depending on whether they're scalar or grid parameters
	int yearPlantedMin = (yearPlantedParam.data.spType == pScalar) ? *yearPlantedParam.adr : yearPlantedParam.data.g->GetMin();
	int startAgeMin = (startAgeParam.data.spType == pScalar) ? *startAgeParam.adr : startAgeParam.data.g->GetMin();
	int endYearMax = (endYearParam.data.spType == pScalar) ? *endYearParam.adr : endYearParam.data.g->GetMax();
	int startMonthMax = (startMonthParam.data.spType == pScalar) ? *startMonthParam.adr : startMonthParam.data.g->GetMax();

	/* 
	determine minMY values 
	*/
	if (startAgeParam.data.spType != pScalar && yearPlantedParam.data.spType != pScalar) {
		//if both are raster, find the minimum sum of min year planted and min start age pixels
		//we do this because minimum year planted and minimum start age may never occur on the same pixels
		double overallMin = 0;
		for (int row = 0; row < yearPlantedParam.data.g->nRows; row++) {
			for (int col = 0; col < yearPlantedParam.data.g->nCols; col++) {
				//convert gotten values to double first so we don't have any overflow of floats as we add them
				double curMin = (double)yearPlantedParam.data.g->GetVal(row, col) + (double)startAgeParam.data.g->GetVal(row, col);
				
				//set the overall minimum accordingly
				overallMin = (curMin < overallMin) ? curMin : overallMin;
			}
		}
	}
	else {
		//otherwise, just add the mins together
		minMY.year = yearPlantedMin + startAgeMin;
	}
	//month isn't used so set to null
	minMY.mon = NULL;

	/* 
	determine maxMY values
	*/
	if (endYearParam.data.spType != pScalar && startMonthParam.data.spType != pScalar) {
		//find the maximum month considering only pixels that are the maximum year
		maxMY.mon = startMonthParam.data.g->maxFromIndices(endYearParam.data.g->getIndicesWhere(endYearMax));
	}
	else {
		//otherwise, just use the start month max
		maxMY.mon = startMonthMax;
	}
	//max year is just the max year
	maxMY.year = endYearMax;

	/*
	error check on years and months
	*/
	//check valid month
	if (maxMY.mon < 0 || maxMY.mon > 12) {
		//if month isn't within the range of 0 to 12, print and log error
		std::string errstr = "Invalid start month detected: " + to_string(minMY.mon);
		std::cout << errstr << std::endl;
		logger->Log(errstr);

		//then exit
		exit(EXIT_FAILURE);
	}

	//check valid years
	if (minMY.year > maxMY.year) {
		//if minimum year is larger than maximum year, print and log error
		std::string errstr = "min year (" + to_string(minMY.year) + ") is greater than max year (" + to_string(maxMY.year) + ")";
		std::cout << errstr << std::endl;
		logger->Log(errstr);

		//then exit
		exit(EXIT_FAILURE);
	}
}