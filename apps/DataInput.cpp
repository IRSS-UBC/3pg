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

bool DataInput::getScalar(std::string value, PPPG_PARAM& param) {
	try {
		//get the value from the array
		double val = std::stod(value);
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

bool DataInput::getGrid(std::string value, PPPG_PARAM& param) {
	try {
		//get the file path
		std::filesystem::path filePath = value;

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
		if (!openCheckGrid(filePath.string(), param.g)) {
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
		std::string errstr = " " + value + " could not be interpreted as a scalar or grid name";
		std::cout << errstr << std::endl;
		logger.Log(errstr);
		logger.Log(e.what());
		return false;
	}
}

bool DataInput::openCheckGrid(std::string path, GDALRasterImage*& grid) {
	//ensure we can open and read from the file as a GDALRasterImage
	try {
		grid = new GDALRasterImage(path);
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
		this->refGrid = grid;
	}
	else {
		if (
			(fabs(this->refGrid->xMin - grid->xMin) > 0.0001) ||	//check xMin
			(fabs(this->refGrid->yMin - grid->yMin) > 0.0001) ||	//check yMin
			(fabs(this->refGrid->xMax - grid->xMax) > 0.0001) ||	//check xMax
			(fabs(this->refGrid->yMax - grid->yMax) > 0.0001) ||	//check yMax
			(this->refGrid->nRows != grid->nRows) ||				//check nRows
			(this->refGrid->nCols != grid->nCols)					//check nCols
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
	if (DataInput::getScalar(value.front(), param)) {
		//if gotten, add to inputParams map
		this->inputParams.emplace(param.id, param);

		//remove from required params set
		this->requiredInputParams3PGS.erase(param.id);
		this->requiredInputParams3PG.erase(param.id);

		return true;
	}
	
	//try to get the value as a grid
	if (DataInput::getGrid(value.front(), param)) {
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

bool DataInput::tryAddSeriesParam(std::string name, std::vector<std::string> values, std::ifstream& paramFp, int& lineNo) {
	//if the name does not have a corrosponding series param, return
	auto search = this->seriesParamNameMap.find(name);
	if (search == this->seriesParamNameMap.end()) {
		return false;
	}
	
	//get the array index from seriesParamNameMap to get the according PPPG_SERIES_PARAM
	PPPG_SERIES_PARAM* param = &this->seriesParams[search->second];

	if (!values.empty()) {
		/*
		parameter is in 'one year' style, where a single year of monthly values is given
		*/

		//ensure all 12 monthly values are given
		if (values.size() != 12) {
			std::string errstr = "there must be 12 monthly values given for parameter " + name + ". There were only " + std::to_string(values.size()) + ".";
			std::cout << errstr << std::endl;
			logger.Log(errstr);
			exit(EXIT_FAILURE);
		}

		//parameters in the 'one year' style are stored in the seriesParamMap with key -1
		param->firstYear = -1;
		param->lastYear = -1;

		//use only the memory we need (12 indices);
		param->monthlyParams.resize(12);
		param->monthlyParams.shrink_to_fit();

		for (int i = 0; i < 12; i++) {
			param->monthlyParams[i].id = name + " month " + std::to_string(i);

			//try to add as a scalar param
			if (DataInput::getScalar(values[i], param->monthlyParams[i])) {
				continue;
			}

			//try to add as a grid param
			if (DataInput::getGrid(values[i], param->monthlyParams[i])) {
				continue;
			}

			//if we've made it here it's neither a usable scalar nor grid parameter
			exit(EXIT_FAILURE);
		}
	}
	else {
		//parameter is in 'multi year' style, where the lines following the current one
		//contain a year followed by monthly values
		std::string line;
		while (std::getline(paramFp, line)) {
			lineNo++;
			int year;
			std::vector<std::string> sTokens;

			//if empty we've iterated through every year, stop iterating
			if (line.size() == 0) {
				break;
			}

			//get the year and monthly tokens from the current line
			sTokens = boost::split(sTokens, line, boost::is_any_of(", \n\t"), boost::token_compress_on);

			//ensure all 12 monthly values are given
			if (sTokens.size() != 13) {
				std::string errstr = "there must be 12 monthly values given for parameter " + name + " at year " + sTokens.front() + ". There were only " + std::to_string(sTokens.size() - 1) + ".";
				std::cout << errstr << std::endl;
				logger.Log(errstr);
				exit(EXIT_FAILURE);
			}

			//try to convert first token to a year
			try {
				year = std::stoi(sTokens.front());
			}
			catch (const std::out_of_range& oor) {
				std::string errstr = "Year could not be converted to integer on line " + std::to_string(lineNo);
				std::cout << errstr << std::endl;
				logger.Log(errstr);
				exit(EXIT_FAILURE);
			}

			//set the first year if it hasn't been set
			if (param->firstYear == 0) {
				param->firstYear = year;
			}

			//set the last year and check to ensure years are given in consecutive order
			if (param->lastYear == 0) {
				param->lastYear = year;
			}
			else if (year == param->lastYear + 1) {
				param->lastYear = year;
			}
			else {
				std::string errstr = "yearly inputs for parameter " + name + " are not in consecutive order on line " + std::to_string(lineNo);
				std::cout << errstr << std::endl;
				logger.Log(errstr);
				exit(EXIT_FAILURE);
			}

			//resize the params vector to fit 12 more values
			param->monthlyParams.resize(param->monthlyParams.size() + 12);

			for (int i = 0; i < 12; i++) {
				int paramIndex = (param->lastYear - param->firstYear) * 12 + i;
				param->monthlyParams[paramIndex].id = name + " year " + sTokens.front() + " month " + std::to_string(i);

				//try to add as a scalar param
				if (DataInput::getScalar(sTokens[i + 1], param->monthlyParams[paramIndex])) {
					continue;
				}

				//try to add as a grid param
				if (DataInput::getGrid(sTokens[i + 1], param->monthlyParams[paramIndex])) {
					continue;
				}

				//if we've made it here it's neither a usable scalar nor grid parameter
				exit(EXIT_FAILURE);
			}
		}
	}

	//add series param name to set of acquired params
	acquiredSeriesParams.insert(name);
	return true;
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

	//series parameters
	bool haveTmax = this->acquiredSeriesParams.find("Tmax") != this->acquiredSeriesParams.end();
	bool haveTmin = this->acquiredSeriesParams.find("Tmin") != this->acquiredSeriesParams.end();
	bool haveTavg = this->acquiredSeriesParams.find("Tavg") != this->acquiredSeriesParams.end();
	bool haveVPD = this->acquiredSeriesParams.find("VPD") != this->acquiredSeriesParams.end();
	bool haveRain = this->acquiredSeriesParams.find("Rain") != this->acquiredSeriesParams.end();
	bool haveSolarRad = this->acquiredSeriesParams.find("Solar Radtn") != this->acquiredSeriesParams.end();
	bool haveNetRad = this->acquiredSeriesParams.find("Net radtn") != this->acquiredSeriesParams.end();
	bool haveFrost = this->acquiredSeriesParams.find("Frost days") != this->acquiredSeriesParams.end();
	bool haveNDVI = this->acquiredSeriesParams.find("NDVI_AVH") != this->acquiredSeriesParams.end();
	
	//check tmax/tavg
	if (!haveTmax && !haveTavg) {
		std::cout << "must have Tmax or Tavg" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//check tmin/tavg
	if (!haveTmin && !haveTavg) {
		std::cout << "must have Tmin or Tavg" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//check tavg/vpd
	if (haveTavg && !haveVPD) {
		std::cout << "must have VPD if using Tavg" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//check rain
	if (!haveRain) {
		std::cout << "must have Rain" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//check solar radation
	if (!haveSolarRad) {
		std::cout << "must have Solar Radiation" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	//check frost
	if (!haveFrost) {
		std::cout << "must have Frost" << std::endl;
		exit(EXIT_FAILURE);
	}

	//check model mode which requires NDVI series parameters
	if (modelMode3PGS && this->acquiredSeriesParams.find("NDVI_AVH") != this->acquiredSeriesParams.end()) {
		std::cout << "3PGS mode should have NDVI series parameters" << std::endl;
		exit(EXIT_FAILURE);
	}

	this->haveNDVI = haveNDVI;
	this->haveNetRad = haveNetRad;
	this->haveVPD = haveVPD;
	this->haveTavg = haveTavg;

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
		params.pFS2 = DataInput::getValFromInputParam("pFS2", cellIndex);
		params.pFS20 = DataInput::getValFromInputParam("pFS20", cellIndex);
		params.StemConst = DataInput::getValFromInputParam("StemConst", cellIndex);
		params.StemPower = DataInput::getValFromInputParam("StemPower", cellIndex);
		params.pRx = DataInput::getValFromInputParam("pRx", cellIndex);
		params.pRn = DataInput::getValFromInputParam("pRn", cellIndex);
		params.growthTmin = DataInput::getValFromInputParam("growthTmin", cellIndex);
		params.growthTopt = DataInput::getValFromInputParam("growthTopt", cellIndex);
		params.growthTmax = DataInput::getValFromInputParam("growthTmax", cellIndex);
		params.kF = DataInput::getValFromInputParam("kF", cellIndex);
		params.gammaFx = DataInput::getValFromInputParam("gammaFx", cellIndex);
		params.gammaF0 = DataInput::getValFromInputParam("gammaF0", cellIndex);
		params.tgammaF = DataInput::getValFromInputParam("tgammaF", cellIndex);
		params.Rttover = DataInput::getValFromInputParam("Rttover", cellIndex);
		params.MaxCond = DataInput::getValFromInputParam("MaxCond", cellIndex);
		params.CoeffCond = DataInput::getValFromInputParam("CoeffCond", cellIndex);
		params.BLcond = DataInput::getValFromInputParam("BLcond", cellIndex);
		params.m0 = DataInput::getValFromInputParam("m0", cellIndex);
		params.fN0 = DataInput::getValFromInputParam("fN0", cellIndex);
		params.fNn = DataInput::getValFromInputParam("fNn", cellIndex);
		params.thinPower = DataInput::getValFromInputParam("thinPower", cellIndex);
		params.mF = DataInput::getValFromInputParam("mF", cellIndex);
		params.mR = DataInput::getValFromInputParam("mR", cellIndex);
		params.mS = DataInput::getValFromInputParam("mS", cellIndex);
		params.SWconst0 = DataInput::getValFromInputParam("SWconst0", cellIndex);
		params.SWpower0 = DataInput::getValFromInputParam("SWpower0", cellIndex);
		params.wSx1000 = DataInput::getValFromInputParam("wSx1000", cellIndex);
		params.MaxAge = DataInput::getValFromInputParam("MaxAge", cellIndex);
		params.nAge = DataInput::getValFromInputParam("nAge", cellIndex);
		params.rAge = DataInput::getValFromInputParam("rAge", cellIndex);
		params.SLA0 = DataInput::getValFromInputParam("SLA0", cellIndex);
		params.SLA1 = DataInput::getValFromInputParam("SLA1", cellIndex);
		params.tSLA = DataInput::getValFromInputParam("tSLA", cellIndex);
		params.k = DataInput::getValFromInputParam("k", cellIndex);
		params.fullCanAge = DataInput::getValFromInputParam("fullCanAge", cellIndex);
		params.alpha = DataInput::getValFromInputParam("alpha", cellIndex);
		params.fracBB0 = DataInput::getValFromInputParam("fracBB0", cellIndex);
		params.fracBB1 = DataInput::getValFromInputParam("fracBB1", cellIndex);
		params.tBB = DataInput::getValFromInputParam("tBB", cellIndex);
		params.y = DataInput::getValFromInputParam("y", cellIndex);
		params.rhoMin = DataInput::getValFromInputParam("rhoMin", cellIndex);
		params.rhoMax = DataInput::getValFromInputParam("rhoMax", cellIndex);
		params.tRho = DataInput::getValFromInputParam("tRho", cellIndex);
		params.Qa = DataInput::getValFromInputParam("Qa", cellIndex);
		params.Qb = DataInput::getValFromInputParam("Qb", cellIndex);
		params.gDM_mol = DataInput::getValFromInputParam("gDM_mol", cellIndex);
		params.molPAR_MJ = DataInput::getValFromInputParam("molPAR_MJ", cellIndex);
		params.LAIgcx = DataInput::getValFromInputParam("LAIgcx", cellIndex);
		params.MaxIntcptn = DataInput::getValFromInputParam("MaxIntcptn", cellIndex);
		params.LAImaxIntcptn = DataInput::getValFromInputParam("LAImaxIntcptn", cellIndex);
		params.Lat = DataInput::getValFromInputParam("Lat", cellIndex);
		params.FRp = DataInput::getValFromInputParam("FRp", cellIndex);
		params.FRstart = DataInput::getValFromInputParam("FRstart", cellIndex);
		params.FRend = DataInput::getValFromInputParam("FRend", cellIndex);
		params.FRdec = DataInput::getValFromInputParam("FRdec", cellIndex);
		params.soilIndex = DataInput::getValFromInputParam("soilIndex", cellIndex);
		params.MaxASW = DataInput::getValFromInputParam("MaxASW", cellIndex);
		params.MinASWp = DataInput::getValFromInputParam("MinASWp", cellIndex);
		params.StartAge = DataInput::getValFromInputParam("StartAge", cellIndex);
		params.EndYear = DataInput::getValFromInputParam("EndYear", cellIndex);
		params.StartMonth = DataInput::getValFromInputParam("StartMonth", cellIndex);
		params.yearPlanted = DataInput::getValFromInputParam("yearPlanted", cellIndex);
		params.SeedlingMass = DataInput::getValFromInputParam("SeedlingMass", cellIndex);
		params.WFi = DataInput::getValFromInputParam("WFi", cellIndex);
		params.WRi = DataInput::getValFromInputParam("WRi", cellIndex);
		params.WSi = DataInput::getValFromInputParam("WSi", cellIndex);
		params.StemNoi = DataInput::getValFromInputParam("StemNoi", cellIndex);
		params.ASWi = DataInput::getValFromInputParam("ASWi", cellIndex);
		params.MinASWTG = DataInput::getValFromInputParam("MinASWTG", cellIndex);
		params.NDVI_FPAR_intercept = DataInput::getValFromInputParam("NDVI_FPAR_intercept", cellIndex);
		params.NDVI_FPAR_constant = DataInput::getValFromInputParam("NDVI_FPAR_constant", cellIndex);

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

double DataInput::getValFromInputParam(std::string paramName, long cellIndex) {
	auto search = this->inputParams.find(paramName);

	//if the param doesn't exist, set to -1
	if (search == this->inputParams.end()) {
		//use predefined default parameters for Lat, rhoMax, rhoMin, and tRho if not set by user
		if (paramName == "Lat") { return 1000; }
		else if (paramName == "rhoMax") { return 0.5; }
		else if (paramName == "rhoMin") { return 0.5; }
		else if (paramName == "tRho") { return 4; }
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

bool DataInput::getSeriesParams(long cellIndex, int year, int month, SeriesParams& params) {
	//error checking
	if (!finishedInput) {
		throw std::exception("should NOT be able to call getInputParams() if input failed!");
	}

	//if one of the calls to getValFromSeriesParams() throws an std::runtime_error
	//it means we're at a nodata pixel
	try {
		double Tmax;
		double Tmin;

		//get parameters we will always have
		params.Rain = DataInput::getValFromSeriesParam(SeriesIndex::RAIN, year, month, cellIndex);
		params.SolarRad = DataInput::getValFromSeriesParam(SeriesIndex::SOLAR_RAD, year, month, cellIndex);
		params.FrostDays = DataInput::getValFromSeriesParam(SeriesIndex::FROST_DAYS, year, month, cellIndex);

		//get optional parameters if we have them
		params.NDVI_AVH = this->haveNDVI ? DataInput::getValFromSeriesParam(SeriesIndex::NDVI_AVH, year, month, cellIndex) : 0;
		params.NetRad = this->haveNetRad ? DataInput::getValFromSeriesParam(SeriesIndex::NET_RAD, year, month, cellIndex) : 0;

		//if we need Tmax and Tmin, get their series params
		if (!this->haveTavg || !this->haveVPD) {
			Tmax = DataInput::getValFromSeriesParam(SeriesIndex::TMAX, year, month, cellIndex);
			Tmin = DataInput::getValFromSeriesParam(SeriesIndex::TMIN, year, month, cellIndex);
		}

		//calculate Tavg if we don't have it directly
		if (this->haveTavg) {
			params.Tavg = DataInput::getValFromSeriesParam(SeriesIndex::TAVG, year, month, cellIndex);
		}
		else {
			params.Tavg = (Tmax + Tmin) / 2;
		}

		//calculate VPD if we don't have it directly
		if (this->haveVPD) {
			params.VPD = DataInput::getValFromSeriesParam(SeriesIndex::VPD, year, month, cellIndex);
		}
		else {
			double VPDmax = 6.1078 * exp(17.269 * Tmax / (237.3 + Tmax));
			double VPDmin = 6.1078 * exp(17.269 * Tmin / (237.3 + Tmin));
			params.VPD = (VPDmax - VPDmin) / 2;
		}

		return true;
	}
	catch (std::runtime_error e) {
		return false;
	}
}

double DataInput::getValFromSeriesParam(int paramIndex, int year, int month, long cellIndex) {
	//get parameter from the seriesParam array using the seriesParamNameMap
	PPPG_SERIES_PARAM* param = &this->seriesParams[paramIndex];
	
	//get the location of the current year/month param based on whether we have multiple years in the vector
	int monthIndex = (param->firstYear == -1) ? (month - 1) : (year - param->firstYear) * 12 + (month - 1);
	
	//get the monthly parameter
	PPPG_PARAM monthParam = param->monthlyParams[monthIndex];

	if (monthParam.spType == pScalar) {
		return monthParam.val;
	}

	//get parameter value depending on whether it is scalar or grid
	if (monthParam.spType == pTif) {
		double val = monthParam.g->GetVal(cellIndex);
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

GDALRasterImage* DataInput::getRefGrid() {
	return this->refGrid;
}

bool DataInput::haveNetRadParam() {
	return this->haveNetRad;
}