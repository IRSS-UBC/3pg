#include "DataInput.hpp"
#include "util.hpp"

extern Logger logger;

DataInput::DataInput() {
	this->refGrid = nullptr;
}

DataInput::~DataInput() {
	//check each parameter
	for (auto iterator = this->inputParams.begin(); iterator != this->inputParams.end(); iterator++) {
		PPPG_PARAM param = iterator->second;

		//if it's a grid parameter, we have to clean up it's corrosponding GDALRasterImage
		if (param.spType == pTif && param.g != nullptr) {
			delete param.g;
			param.g = nullptr;
		}
	}
}

bool DataInput::getScalar(std::vector<std::string> value, PPPG_PARAM& param) {
	try {
		//get the value from the array
		double val = std::stod(value.front());
		param.val = val;

		//mark the data as scalar
		param.spType = pScalar;

		//log input param acquired
		std::string output = "    " + param.id + "        constant: " + std::to_string(param.val);
		std::cout << output << std::endl;
		logger.Log(output);

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
			logger.Log(errstr);
			return false;
		}

		//ensure the filepath exists
		if (!std::filesystem::exists(filePath)) {
			std::string errstr = filePath.string() + " does not exist.";
			std::cout << errstr << std::endl;
			logger.Log(errstr);
			return false;
		}

		//set data type
		param.spType = pTif;

		//try opening the grid and compare to the refgrid (or create the refgrid)
		if (!openCheckGrid(filePath.string(), param)) {
			return false;
		}

		//log input param
		std::string output = "    " + param.id + "        raster: " + param.g->name;
		std::cout << output << std::endl;
		logger.Log(output);
		return true;
	}
	catch (std::filesystem::filesystem_error const& e) {
		//set and print/log error string
		std::string errstr = " " + value.front() + " could not be interpreted as a scalar or grid name";
		std::cout << errstr << std::endl;
		logger.Log(errstr);
		logger.Log(e.what());
		return false;
	}
}

bool DataInput::openCheckGrid(std::string path, PPPG_PARAM& param) {
	//ensure we can open and read from the file as a GDALRasterImage
	try {
		param.g = new GDALRasterImage(path);
	}
	catch (const std::exception& e) {
		std::string errstr = "failed to open " + path + "\n" + e.what();
		std::cout << errstr << std::endl;
		logger.Log(errstr);
		return false;
	}

	//if the refgrid is null, this is grid becomes the refgrid
	//otherwise, check grid dimensions
	if (this->refGrid == nullptr) {
		this->refGrid = param.g;
	}
	else {
		if (
			(fabs(this->refGrid->xMin - param.g->xMin) > 0.0001) ||	//check xMin
			(fabs(this->refGrid->yMin - param.g->yMin) > 0.0001) ||	//check yMin
			(fabs(this->refGrid->xMax - param.g->xMax) > 0.0001) ||	//check xMax
			(fabs(this->refGrid->yMax - param.g->yMax) > 0.0001) ||	//check yMax
			(this->refGrid->nRows != param.g->nRows) ||				//check nRows
			(this->refGrid->nCols != param.g->nCols)					//check nCols
			) {
			std::string errstr = "Grid dimensions of " + path + " differs from " + this->refGrid->name;
			std::cout << errstr << std::endl;
			logger.Log(errstr);
			return false;
		}
	}

	return true;
}

bool DataInput::tryAddInputParam(std::string name, std::vector<std::string> value) {
	//TODO: actual 3pg doesn't care about case (I don't think), make sure that is reflected here
	PPPG_PARAM param;

	//if the string isn't the exact parameter name, see if it is in the parameter name map
	if (this->allInputParams.find(name) == this->allInputParams.end()) {
		//search for the actual (shortened) param name
		auto search = this->inputParamNames.find(name);

		//if the string passed isn't an input param return false
		if (search == this->inputParamNames.end()) {
			return false;
		}

		//the string passed is the long version of the input param: set the id accordingly
		param.id = search->second;
	}
	else {
		//if the string is the exact parameter name, set the id accordingly
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
		this->requiredInputParams3PGS.erase(param.id);
		this->requiredInputParams3PG.erase(param.id);

		return true;
	}
	
	//try to get the value as a grid
	if (DataInput::getGrid(value, param)) {
		//if gotten, add to inputParams map and return
		this->inputParams.emplace(param.id, param);

		//remove from required params set
		this->requiredInputParams3PGS.erase(param.id);
		this->requiredInputParams3PG.erase(param.id);

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

	//if we don't have seedling mass, need all of WFi, WRi, and WSi
	if (!haveSeedlingMass && (!haveWFi || !haveWRi || !haveWSi)) {
		std::string errstr;
		if (!haveSeedlingMass) {
			errstr = "Missing parameter: SeedlingMass.";
		}
		else {
			errstr = "Missing parameter: at least one of WFi, WRi, and WSi required.";
		}
		std::cout << errstr << std::endl;
		logger.Log(errstr);
		return false;
	}

	this->haveSeedlingMass = haveSeedlingMass;
	this->haveMinASWTG = this->inputParams.find("MinASWTG") != this->inputParams.end();
	this->haveAgeDepFert = (
		this->inputParams.find("FRstart") != this->inputParams.end() &&
		this->inputParams.find("FRend") != this->inputParams.end() &&
		this->inputParams.find("FRdec") != this->inputParams.end()
	);


	//if we're using 3PGS, ensure we have all required 3PGS parameters
	if (modelMode3PGS && this->requiredInputParams3PGS.size() != 0) {
		return false;
	}

	//if we're using 3PG (not 3PGS), ensure we have all required 3PG parameters
	if (!modelMode3PGS && this->requiredInputParams3PG.size() != 0) {
		return false;
	}

	finishedInput = true;
	return true;
}

bool DataInput::getInputParams(long cellIndex, InputParams& params) {
	//error checking
	if (!finishedInput) {
		throw std::exception("should NOT be able to call getInputParams() if input failed!");
	}

	//declare, copy values into, then return an instance of InputParams
	try {
		params.pFS2 = DataInput::getValFromParam("pFS2", cellIndex);
		params.pFS20 = DataInput::getValFromParam("pFS20", cellIndex);
		params.StemConst = DataInput::getValFromParam("StemConst", cellIndex);
		params.StemPower = DataInput::getValFromParam("StemPower", cellIndex);
		params.pRx = DataInput::getValFromParam("pRx", cellIndex);
		params.pRn = DataInput::getValFromParam("pRn", cellIndex);
		params.growthTmin = DataInput::getValFromParam("growthTmin", cellIndex);
		params.growthTopt = DataInput::getValFromParam("growthTopt", cellIndex);
		params.growthTmax = DataInput::getValFromParam("growthTmax", cellIndex);
		params.kF = DataInput::getValFromParam("kF", cellIndex);
		params.gammaFx = DataInput::getValFromParam("gammaFx", cellIndex);
		params.gammaF0 = DataInput::getValFromParam("gammaF0", cellIndex);
		params.tgammaF = DataInput::getValFromParam("tgammaF", cellIndex);
		params.Rttover = DataInput::getValFromParam("Rttover", cellIndex);
		params.MaxCond = DataInput::getValFromParam("MaxCond", cellIndex);
		params.CoeffCond = DataInput::getValFromParam("CoeffCond", cellIndex);
		params.BLcond = DataInput::getValFromParam("BLcond", cellIndex);
		params.m0 = DataInput::getValFromParam("m0", cellIndex);
		params.fN0 = DataInput::getValFromParam("fN0", cellIndex);
		params.fNn = DataInput::getValFromParam("fNn", cellIndex);
		params.thinPower = DataInput::getValFromParam("thinPower", cellIndex);
		params.mF = DataInput::getValFromParam("mF", cellIndex);
		params.mR = DataInput::getValFromParam("mR", cellIndex);
		params.mS = DataInput::getValFromParam("mS", cellIndex);
		params.SWconst0 = DataInput::getValFromParam("SWconst0", cellIndex);
		params.SWpower0 = DataInput::getValFromParam("SWpower0", cellIndex);
		params.wSx1000 = DataInput::getValFromParam("wSx1000", cellIndex);
		params.MaxAge = DataInput::getValFromParam("MaxAge", cellIndex);
		params.nAge = DataInput::getValFromParam("nAge", cellIndex);
		params.rAge = DataInput::getValFromParam("rAge", cellIndex);
		params.SLA0 = DataInput::getValFromParam("SLA0", cellIndex);
		params.SLA1 = DataInput::getValFromParam("SLA1", cellIndex);
		params.tSLA = DataInput::getValFromParam("tSLA", cellIndex);
		params.k = DataInput::getValFromParam("k", cellIndex);
		params.fullCanAge = DataInput::getValFromParam("fullCanAge", cellIndex);
		params.alpha = DataInput::getValFromParam("alpha", cellIndex);
		params.fracBB0 = DataInput::getValFromParam("fracBB0", cellIndex);
		params.fracBB1 = DataInput::getValFromParam("fracBB1", cellIndex);
		params.tBB = DataInput::getValFromParam("tBB", cellIndex);
		params.y = DataInput::getValFromParam("y", cellIndex);
		params.rhoMin = DataInput::getValFromParam("rhoMin", cellIndex);
		params.rhoMax = DataInput::getValFromParam("rhoMax", cellIndex);
		params.tRho = DataInput::getValFromParam("tRho", cellIndex);
		params.Qa = DataInput::getValFromParam("Qa", cellIndex);
		params.Qb = DataInput::getValFromParam("Qb", cellIndex);
		params.gDM_mol = DataInput::getValFromParam("gDM_mol", cellIndex);
		params.molPAR_MJ = DataInput::getValFromParam("molPAR_MJ", cellIndex);
		params.LAIgcx = DataInput::getValFromParam("LAIgcx", cellIndex);
		params.MaxIntcptn = DataInput::getValFromParam("MaxIntcptn", cellIndex);
		params.LAImaxIntcptn = DataInput::getValFromParam("LAImaxIntcptn", cellIndex);
		params.Lat = DataInput::getValFromParam("Lat", cellIndex);
		params.FRp = DataInput::getValFromParam("FRp", cellIndex);
		params.FRstart = DataInput::getValFromParam("FRstart", cellIndex);
		params.FRend = DataInput::getValFromParam("FRend", cellIndex);
		params.FRdec = DataInput::getValFromParam("FRdec", cellIndex);
		params.soilIndex = DataInput::getValFromParam("soilIndex", cellIndex);
		params.MaxASW = DataInput::getValFromParam("MaxASW", cellIndex);
		params.MinASWp = DataInput::getValFromParam("MinASWp", cellIndex);
		params.StartAge = DataInput::getValFromParam("StartAge", cellIndex);
		params.EndYear = DataInput::getValFromParam("EndYear", cellIndex);
		params.StartMonth = DataInput::getValFromParam("StartMonth", cellIndex);
		params.yearPlanted = DataInput::getValFromParam("yearPlanted", cellIndex);
		params.SeedlingMass = DataInput::getValFromParam("SeedlingMass", cellIndex);
		params.WFi = DataInput::getValFromParam("WFi", cellIndex);
		params.WRi = DataInput::getValFromParam("WRi", cellIndex);
		params.WSi = DataInput::getValFromParam("WSi", cellIndex);
		params.StemNoi = DataInput::getValFromParam("StemNoi", cellIndex);
		params.ASWi = DataInput::getValFromParam("ASWi", cellIndex);
		params.MinASWTG = DataInput::getValFromParam("MinASWTG", cellIndex);
		params.NDVI_FPAR_intercept = DataInput::getValFromParam("NDVI_FPAR_intercept", cellIndex);
		params.NDVI_FPAR_constant = DataInput::getValFromParam("NDVI_FPAR_constant", cellIndex);

		if (params.yearPlanted < 1 || isnan(params.yearPlanted) || params.StemNoi < 1) {
			return false;
		}

		if (params.StartAge < 1) {
			params.StartAge = 1;
		}

		if (params.EndYear < 1) {
			params.EndYear = 2;
		}

		if (params.StartMonth < 1) {
			params.StartMonth = 1;
		}

		return true;
	}
	catch (std::runtime_error e) {
		return false;
	}
}

double DataInput::getValFromParam(std::string paramName, long cellIndex) {
	auto search = this->inputParams.find(paramName);

	//if the param doesn't exist, set to -1
	if (search == this->inputParams.end()) {
		//use predefined default parameters for Lat, rhoMax, rhoMin, and tRho if not set by user
		if (paramName == "Lat")			{ return 1000; }
		else if (paramName == "rhoMax") { return 0.5; }
		else if (paramName == "rhoMin") { return 0.5; }
		else if (paramName == "tRho")	{ return 4; }
		else {
			//otherwise, set to 0
			return 0;
		}
	}

	//if the param is scalar, return it's value
	if (search->second.spType == pScalar) {
		return search->second.val;
	}

	//if the param is a grid, return the value at the row and column specified
	if (search->second.spType == pTif) {
		double val = search->second.g->GetVal(cellIndex);
		if (isnan(val)) {
			throw std::runtime_error("nan");
		}
		else {
			return val;
		}
	}

	throw std::exception("a parameter has been set incorrectly as neither a scalar or a grid.");
}

void DataInput::findRunPeriod(MYDate& minMY, MYDate& maxMY) {
	//error checking
	if (!finishedInput) {
		throw std::exception("should NOT be able to call findRunPeriod() if input failed!");
	}

	//get required parameters from parameter map
	PPPG_PARAM yearPlantedParam = this->inputParams.find("yearPlanted")->second;
	PPPG_PARAM startAgeParam = this->inputParams.find("StartAge")->second;
	PPPG_PARAM endYearParam = this->inputParams.find("EndYear")->second;
	PPPG_PARAM startMonthParam = this->inputParams.find("StartMonth")->second;

	//get maxes and mins depending on whether they're scalar or grid parameters
	int yearPlantedMin = (yearPlantedParam.spType == pScalar) ? yearPlantedParam.val : yearPlantedParam.g->GetMin();
	int startAgeMin = (startAgeParam.spType == pScalar) ? startAgeParam.val : startAgeParam.g->GetMin();
	int endYearMax = (endYearParam.spType == pScalar) ? endYearParam.val : endYearParam.g->GetMax();
	int startMonthMax = (startMonthParam.spType == pScalar) ? startMonthParam.val : startMonthParam.g->GetMax();

	/* 
	determine minMY values 
	*/
	if (startAgeParam.spType != pScalar && yearPlantedParam.spType != pScalar) {
		//if both are raster, find smallest sum of yearPlanted and startAge pixels.
		//we do this because the minimum sum (which is the year we should start on)
		//does not have to be at any of the pixels where yearPlanted or startAge are smallest
		double overallMin = std::numeric_limits<double>::max();
		for (int row = 0; row < yearPlantedParam.g->nRows; row++) {
			for (int col = 0; col < yearPlantedParam.g->nCols; col++) {
				//convert gotten values to double first so we don't have any overflow of floats as we add them
				double curMin = (double)yearPlantedParam.g->GetVal(row, col) + (double)startAgeParam.g->GetVal(row, col);
				
				//set the overall minimum accordingly
				overallMin = (curMin < overallMin) ? curMin : overallMin;
			}
		}

		minMY.year = overallMin;
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
	if (endYearParam.spType != pScalar && startMonthParam.spType != pScalar) {
		//find the maximum month considering only pixels that are the maximum year
		maxMY.mon = startMonthParam.g->maxFromIndices(endYearParam.g->getIndicesWhere(endYearMax));
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
		std::string errstr = "Invalid start month detected: " + std::to_string(minMY.mon);
		std::cout << errstr << std::endl;
		logger.Log(errstr);

		//then exit
		exit(EXIT_FAILURE);
	}

	//check valid years
	if (minMY.year > maxMY.year) {
		//if minimum year is larger than maximum year, print and log error
		std::string errstr = "min year (" + std::to_string(minMY.year) + ") is greater than max year (" + std::to_string(maxMY.year) + ")";
		std::cout << errstr << std::endl;
		logger.Log(errstr);

		//then exit
		exit(EXIT_FAILURE);
	}

	//valid run period successfully determined
	string runPeriodStr = "first run year = " + to_string(minMY.year) + ", last run mon/year = " + to_string(maxMY.mon) + "/" + to_string(maxMY.year);
	std::cout << runPeriodStr << std::endl;
	logger.Log(runPeriodStr);
}