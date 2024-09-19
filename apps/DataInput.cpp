#include "DataInput.hpp"

DataInput::DataInput(std::function<void(std::string)>& log) {
	this->log = log;
}

DataInput::DataInput() {
	this->log = [](std::string message) {
		//do nothing
	};
}

bool DataInput::getScalar(std::string value, PPPG_PARAM* param) {
	try {
		//get the value from the array
		double val = std::stod(value);
		param->val = val;

		//mark the data as scalar
		param->spType = pScalar;

		//log input param acquired
		std::string output = "    " + param->id + "        constant: " + std::to_string(param->val);
		this->log(output);

		return true;
	}
	catch (const std::invalid_argument&) {
		//if we threw an error, it's not a double value so return false
		return false;
	}
}

bool DataInput::getGrid(std::string value, PPPG_PARAM* param) {
	try {
		//get the file path
		std::filesystem::path filePath = value;

		//ensure the file extension is .tif
		if (filePath.extension() != ".tif") {
			std::string errstr = filePath.filename().generic_string() + " is an invalid file type. File extension must be '.tif'";
			std::cout << errstr << std::endl;
			this->log(errstr);
			return false;
		}

		//ensure the filepath exists
		if (!std::filesystem::exists(filePath)) {
			std::string errstr = filePath.string() + " does not exist.";
			std::cout << errstr << std::endl;
			this->log(errstr);
			return false;
		}

		//set data type
		param->spType = pTif;

		//try opening the grid and compare to the refgrid (or create the refgrid)
		if (!openCheckGrid(filePath.string(), param->g)) {
			return false;
		}

		//log input param
		std::string output = "    " + param->id + "        raster: " + param->g->name;
		this->log(output);
		return true;
	}
	catch (std::filesystem::filesystem_error& e) {
		//set and print/log error string
		std::string errstr = " " + value + " could not be interpreted as a scalar or grid name";
		std::cout << errstr << std::endl;
		this->log(errstr);
		this->log(e.what());
		return false;
	}
}

bool DataInput::openCheckGrid(std::string path, std::unique_ptr<GDALRasterImage>& grid) {
	//ensure we can open and read from the file as a GDALRasterImage
	try {
		grid = std::make_unique<GDALRasterImage>(path);
	}
	catch (const std::exception& e) {
		std::string errstr = "failed to open " + path + "\n" + e.what();
		std::cout << errstr << std::endl;
		this->log(errstr);
		return false;
	}

	//if the refgrid is null, this is grid becomes the refgrid
	//otherwise, check grid dimensions
	if (this->refGrid.name == "empty") {
		this->refGrid = grid->getRefGrid();
	}
	else {
		if (
			(fabs(this->refGrid.xMin - grid->xMin) > 0.0001) ||	//check xMin
			(fabs(this->refGrid.yMin - grid->yMin) > 0.0001) ||	//check yMin
			(fabs(this->refGrid.xMax - grid->xMax) > 0.0001) ||	//check xMax
			(fabs(this->refGrid.yMax - grid->yMax) > 0.0001) ||	//check yMax
			(this->refGrid.nRows != grid->nRows) ||				//check nRows
			(this->refGrid.nCols != grid->nCols)					//check nCols
			) {
			std::string errstr = "Grid dimensions of " + path + " differs from " + this->refGrid.name;
			std::cout << errstr << std::endl;
			this->log(errstr);
			return false;
		}
	}

	return true;
}

bool DataInput::tryAddInputParam(std::string name, std::vector<std::string> value) {
	//convert to lowercase for comparision
	boost::algorithm::to_lower(name);

	//TODO: actual 3pg doesn't care about case (I don't think), make sure that is reflected here
	std::unique_ptr<PPPG_PARAM> param = std::make_unique<PPPG_PARAM>();

	//if the string isn't the exact parameter name, see if it is in the parameter name map
	if (!this->allInputParams.contains(name)) {
		//if the string passed isn't an input param return false
		if (!this->inputParamNames.contains(name)) {
			return false;
		}

		//the string passed is the long version of the input param: set the id accordingly
		param->id = this->inputParamNames.at(name);
	}
	else {
		//if the string is the exact parameter name, set the id accordingly
		param->id = name;
	}

	//soilIndex special case
	if (param->id == "soilindex") {
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

	//try to get the value as a scalar then grid
	if (DataInput::getScalar(value.front(), param.get()) || DataInput::getGrid(value.front(), param.get())) {
		//remove from required params set
		this->requiredInputParams3PGS.erase(param->id);
		this->requiredInputParams3PG.erase(param->id);

		//if gotten, add to inputParams map
		this->inputParams.emplace(param->id, std::move(param));
		return true;
	}

	//if we've made it here, the user intended it as an input param but
	//we were unable to add it as one. Fail.
	exit(EXIT_FAILURE);
}

bool DataInput::tryAddSeriesParam(std::string name, std::vector<std::string> values, std::ifstream& paramFp, int& lineNo) {
	//convert to lowercase for comparision
	boost::algorithm::to_lower(name);

	//if the name does not have a corrosponding series param, return
	if (!this->seriesParamNameMap.contains(name)) {
		return false;
	}

	//get the index of the series param
	SeriesIndex index = this->seriesParamNameMap.at(name);
	
	//get the array index from seriesParamNameMap to get the according PPPG_SERIES_PARAM
	PPPG_SERIES_PARAM* param = &this->seriesParams[index];

	if (!values.empty()) {
		/*
		parameter is in 'one year' style, where a single year of monthly values is given
		*/

		//ensure all 12 monthly values are given
		if (values.size() != 12) {
			std::string errstr = "there must be 12 monthly values given for parameter " + name + ". There were only " + std::to_string(values.size()) + ".";
			std::cout << errstr << std::endl;
			this->log(errstr);
			exit(EXIT_FAILURE);
		}

		//parameters in the 'one year' style are stored in the seriesParamMap with key -1
		param->firstYear = -1;
		param->lastYear = -1;

		//use only the memory we need (12 indices);
		param->monthlyParams.resize(12);
		param->monthlyParams.shrink_to_fit();

		for (int i = 0; i < 12; i++) {
			std::unique_ptr<PPPG_PARAM> monthlyParam = std::make_unique<PPPG_PARAM>();
			monthlyParam->id = name + " month " + std::to_string(i);

			//try to add as a scalar param then grid
			if (DataInput::getScalar(values[i], monthlyParam.get()) || DataInput::getGrid(values[i], monthlyParam.get())) {
				param->monthlyParams[i] = std::move(monthlyParam);
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
				this->log(errstr);
				exit(EXIT_FAILURE);
			}

			//try to convert first token to a year
			try {
				year = std::stoi(sTokens.front());
			}
			catch (const std::out_of_range&) {
				std::string errstr = "Year could not be converted to integer on line " + std::to_string(lineNo);
				std::cout << errstr << std::endl;
				this->log(errstr);
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
				this->log(errstr);
				exit(EXIT_FAILURE);
			}

			//resize the params vector to fit 12 more values
			param->monthlyParams.resize(param->monthlyParams.size() + 12);

			for (int i = 0; i < 12; i++) {
				int paramIndex = (param->lastYear - param->firstYear) * 12 + i;
				std::unique_ptr<PPPG_PARAM> monthlyParam = std::make_unique<PPPG_PARAM>();
				monthlyParam->id = name + " year " + sTokens.front() + " month " + std::to_string(i);

				//try to add as a scalar then grid param
				if (DataInput::getScalar(sTokens[i + 1], monthlyParam.get()) || DataInput::getGrid(sTokens[i + 1], monthlyParam.get())) {
					param->monthlyParams[paramIndex] = std::move(monthlyParam);
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

bool DataInput::tryAddOutputParam(std::string name, std::vector<std::string> value, int lineNo) {
	boost::algorithm::to_lower(name);

	PPPG_OP_VAR opVar;

	//fr is a possible output param, although it's also an input param. The fr output param is indicated by frout, not fr.
	if (name == "fr") {
		return false;
	}

	//get name if it is an output param, return false if it isn't
	if (!this->allOutputParams.contains(name)) {
		if (!this->outputParamNames.contains(name)) {
			return false;
		}

		opVar.id = this->outputParamNames.at(name);
	}
	else {
		opVar.id = name;
	}

	//ensure we have enough tokens
	if (value.empty()) {
		std::string outstr = "No grid name for param " + opVar.id + " on line: " + std::to_string(lineNo);
		std::cout << outstr << std::endl;
		this->log(outstr);
		exit(EXIT_FAILURE);
	}

	//ensure we don't have too many tokens
	if (value.size() > 5) {
		std::string outstr = "More than 5 value elements detected for param " + opVar.id + " on line: " + std::to_string(lineNo);
		std::cout << outstr << std::endl;
		this->log(outstr);
		exit(EXIT_FAILURE);
	}

	//ensure correct file extensions
	if (!value.front().ends_with(".tif")) {
		std::string outstr = value.front() + " is an invalid filename. Found " + value.front() + " but must be '.tif'";
		std::cout << outstr << std::endl;
		this->log(outstr);
		exit(EXIT_FAILURE);
	}

	//set the gridname in the parameter
	opVar.gridName = value.front();
	if (opVar.gridName.substr(opVar.gridName.find_last_of(".") + 1) != "tif") {
		std::string outstr = "output type must be of type tif.";
		std::cout << outstr << std::endl;
		this->log(outstr);
	}

	// Check for optional second, third, fourth and fifth tokens; these are used to specify recurring output pattern.
	// The following parsing rules apply:
	//    - If second token exists, then a third and fourth token must also exist. A fifth token is optional.
	//    - The third token must be an integer, representing the start year of the recurrence pattern.
	//    - The fourth token must be 'monthly' or 'month'
	//    - If fourth token is 'monthly', then fifth token must not exist (assumed to be 1).
	//    - If fourth token is 'month', then fifth token must be an integer between 1 and 12.
	if (value.size() > 1) {
		//ensure the start year is an integer
		try {
			opVar.recurStart = std::stoi(value[1]);
		}
		catch (std::invalid_argument) {
			std::string outstr = "Expected an integer start year in recuring output specification on line " + std::to_string(lineNo);
			std::cout << outstr << std::endl;
			this->log(outstr);
			exit(EXIT_FAILURE);
		}

		//ensure we have an interval
		if (value.size() < 3) {
			std::string outstr = "Expected an integer start year in recuring output specification on line " + std::to_string(lineNo);
			std::cout << outstr << std::endl;
			this->log(outstr);
			exit(EXIT_FAILURE);
		}

		//ensure the interval is an integer
		try {
			opVar.recurYear = std::stoi(value[2]);
		}
		catch (std::invalid_argument) {
			std::string outstr = "Expected an integer interval in recuring output specification on line " + std::to_string(lineNo);
			std::cout << outstr << std::endl;
			this->log(outstr);
			exit(EXIT_FAILURE);
		}

		//ensure interval isn't zero
		if (opVar.recurYear == 0) {
			std::string outstr = "Found interval of zero years in recuring output specification on line " + std::to_string(lineNo) + ". Expected non-zero";
			std::cout << outstr << std::endl;
			this->log(outstr);
			exit(EXIT_FAILURE);
		}

		//ensure we have month/monthly keyword
		if (value.size() < 4) {
			std::string outstr = "Expected an integer interval in recuring output specification on line " + std::to_string(lineNo);
			std::cout << outstr << std::endl;
			this->log(outstr);
			exit(EXIT_FAILURE);
		}

		//ensure the keyword is one we can use
		if (value[3] != "month" && value[3] != "monthly") {
			std::string outstr = "Unrecognised keyword '" + value[3] + "' on line " + std::to_string(lineNo) + ". expecting 'month' or 'monthly'.";
			std::cout << outstr << std::endl;
			this->log(outstr);
			exit(EXIT_FAILURE);
		}

		//set monthly recurrence
		opVar.recurMonthly = (value[3] == "monthly");

		if (opVar.recurMonthly) {
			//ensure the user didn't add too many inputs
			if (value.size() > 4) {
				std::string outstr = "too many inputs were given on line " + std::to_string(lineNo) + ". For monthly outputs, no month needs to be given.";
				std::cout << outstr << std::endl;
				this->log(outstr);
				exit(EXIT_FAILURE);
			}
		}
		else {
			//ensure a month was given
			if (value.size() < 5) {
				std::string outstr = "Found 'month' keyword but no month in recuring output specification on line " + std::to_string(lineNo);
				std::cout << outstr << std::endl;
				this->log(outstr);
				exit(EXIT_FAILURE);
			}

			//ensure the month is an integer
			try {
				opVar.recurMonth = std::stoi(value[4]);
			}
			catch (std::invalid_argument) {
				std::string outstr = "Expected an integer month in recuring output specification on line " + std::to_string(lineNo);
				std::cout << outstr << std::endl;
				this->log(outstr);
				exit(EXIT_FAILURE);
			}

			//ensure the month isn't 0
			if (opVar.recurMonth == 0) {
				std::string outstr = "Found month of zero in recuring output specification on line " + std::to_string(lineNo) + ". Expected non-zero";
				std::cout << outstr << std::endl;
				this->log(outstr);
				exit(EXIT_FAILURE);
			}
		}
	}

	//log parameter
	this->log("   variable: " + opVar.id + "   grid: " + opVar.gridName);
	if (opVar.recurStart) {
		std::string outputGridString = "      starting in " + std::to_string(opVar.recurStart) + ", writing every " + std::to_string(opVar.recurYear) + " years";
		if (opVar.recurMonthly) {
			outputGridString += ", with monthly values.";
		}
		else if (opVar.recurMonth != 0){
			outputGridString += ", on the " + std::to_string(opVar.recurMonth) + " month.";
		}
		this->log(outputGridString);
	}

	//continuous check for 3PG and 3PGS specific parameters
	this->allow3PG = this->allow3PG && !this->only3PGS.contains(opVar.id);
	this->allow3PGS = this->allow3PGS && !this->only3PG.contains(opVar.id);
	this->outputParams.emplace(opVar.id, opVar);

	return true;
}

bool DataInput::tryAddManagementParam(std::string name, std::ifstream& inFile, int& lineNo) {
	//convert to lowercase for comparision
	boost::algorithm::to_lower(name);

	int index;
	if (name.compare("management: fertility") == 0) {
		index = ManagementIndex::FERTILITY;
	}
	else if (name.compare("management: minasw") == 0) {
		index = ManagementIndex::MINASW;
	}
	else if (name.compare("management: irrigation") == 0) {
		index = ManagementIndex::IRRIGATION;
	}
	else {
		return false;
	}
	
	PPPG_MT_PARAM* table = &this->managementTables[index];
	
	std::string line;
	while (std::getline(inFile, line)) {
		//increment line number
		lineNo++;
	
		//blank line terminates table
		if (line.empty()) {
			break;
		}
	
		//tokenize line
		std::vector<std::string> tokens;
		boost::split(tokens, line, boost::is_any_of(", \n\t"), boost::token_compress_on);
	
		//ensure line has exactly 2 tokens
		if (tokens.size() != 2) {
			//print and log error
			std::string errstr = "could not read management table at line " + std::to_string(lineNo);
			std::cout << errstr << std::endl;
			this->log(errstr);
	
			//exit
			exit(EXIT_FAILURE);
		}
	
		int year;
		std::unique_ptr<PPPG_PARAM> param = std::make_unique<PPPG_PARAM>();
	
		//ensure the year is an integer
		try {
			year = std::stoi(tokens[0]);
		}
		catch (std::invalid_argument const&) {
			//print and log error
			std::string errstr = "expected aninteger year in management table at line " + std::to_string(lineNo);
			std::cout << errstr << std::endl;
			this->log(errstr);
	
			//exit
			exit(EXIT_FAILURE);
		}
	
		if (DataInput::getScalar(tokens[1], param.get()) || DataInput::getGrid(tokens[1], param.get())) {
			//if this is the first entry in the management table, set the initial year
			if (table->firstYear == -1) {
				table->firstYear = year;
			}
	
			//if years were skipped in the management table, fill in those years with the same index as the last element
			if (year > table->firstYear + table->yearToIndex.size()) {
				while (year > table->firstYear + table->yearToIndex.size()) {
					table->yearToIndex.push_back(table->yearToIndex.back());
				}
			}
	
			//add new param to the back of the yearlyParams table, and add the index to the index table
			table->yearlyParams.emplace_back();
			table->yearlyParams.back() = std::move(param);
			table->yearToIndex.push_back((int)table->yearlyParams.size() - 1);
			
			//keep iterating
			continue;
		}
	
		//if the value could not be interpreted as either a scalar or a grid, then exit
		exit(EXIT_FAILURE);
	}

	return true;
}

bool DataInput::inputFinished(bool modelMode3PGS) {
	if (modelMode3PGS && !this->allow3PGS || !modelMode3PGS && !this->allow3PG) {
		std::string mode = modelMode3PGS ? "3PGS" : "3PG";
		std::string errstr = "output parameters selected are not compatable with " + mode + " mode.";

		std::cout << errstr << std::endl;
		this->log(errstr);
		return false;
	}

	bool haveSeedlingMass = this->inputParams.contains("seedlingmass");
	bool haveWFi = this->inputParams.contains("wfi");
	bool haveWRi = this->inputParams.contains("wri");
	bool haveWSi = this->inputParams.contains("wsi");

	//we must have either seedling mass, or all of WFi, WRi, and WSi
	if (!haveSeedlingMass && (!haveWFi || !haveWRi || !haveWSi)) {
		std::string errstr;
		if (!haveSeedlingMass) {
			errstr = "Missing parameter: SeedlingMass.";
		}
		else {
			errstr = "Missing parameter: at least one of WFi, WRi, and WSi required.";
		}
		std::cout << errstr << std::endl;
		this->log(errstr);
		return false;
	}

	this->haveSeedlingMass = haveSeedlingMass;
	this->haveMinASWTG = this->inputParams.contains("minaswtg");
	this->haveAgeDepFert = (
		this->inputParams.contains("frstart") &&
		this->inputParams.contains("frend") &&
		this->inputParams.contains("frdec")
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
	bool haveTmax = this->acquiredSeriesParams.contains("tmax");
	bool haveTmin = this->acquiredSeriesParams.contains("tmin");
	bool haveTavg = this->acquiredSeriesParams.contains("tavg");
	bool haveVPD = this->acquiredSeriesParams.contains("vpd");
	bool haveRain = this->acquiredSeriesParams.contains("rain");
	bool haveSolarRad = this->acquiredSeriesParams.contains("solar radtn");
	bool haveNetRad = this->acquiredSeriesParams.contains("net radtn");
	bool haveFrost = this->acquiredSeriesParams.contains("frost") || this->acquiredSeriesParams.contains("frost days");
	bool haveNDVI = this->acquiredSeriesParams.contains("ndvi_avh");
	
	//check Tavg
	if (!haveTavg && (!haveTmax || !haveTmin)) {
		std::cout << "must have both Tmax and Tmin if lacking Tavg" << std::endl;
		return false;
	}
	
	//check VPD
	if (!haveVPD && (!haveTmax || !haveTmin)) {
		std::cout << "must have both Tmax and Tmin if lacking VPD" << std::endl;
		return false;
	}
	
	//check rain
	if (!haveRain) {
		std::cout << "must have Rain" << std::endl;
		return false;
	}
	
	//check solar radation
	if (!haveSolarRad) {
		std::cout << "must have Solar Radiation" << std::endl;
		return false;
	}
	
	//check frost
	if (!haveFrost) {
		std::cout << "must have Frost" << std::endl;
		return false;
	}

	//check model mode which requires NDVI series parameters
	if (modelMode3PGS && !haveNDVI) {
		std::cout << "3PGS mode should have NDVI series parameters" << std::endl;
		return false;
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
		params.pFS2 = DataInput::getValFromInputParam("pfs2", cellIndex);
		params.pFS20 = DataInput::getValFromInputParam("pfs20", cellIndex);
		params.StemConst = DataInput::getValFromInputParam("stemconst", cellIndex);
		params.StemPower = DataInput::getValFromInputParam("stempower", cellIndex);
		params.pRx = DataInput::getValFromInputParam("prx", cellIndex);
		params.pRn = DataInput::getValFromInputParam("prn", cellIndex);
		params.growthTmin = DataInput::getValFromInputParam("growthtmin", cellIndex);
		params.growthTopt = DataInput::getValFromInputParam("growthtopt", cellIndex);
		params.growthTmax = DataInput::getValFromInputParam("growthtmax", cellIndex);
		params.kF = DataInput::getValFromInputParam("kf", cellIndex);
		params.gammaFx = DataInput::getValFromInputParam("gammafx", cellIndex);
		params.gammaF0 = DataInput::getValFromInputParam("gammaf0", cellIndex);
		params.tgammaF = DataInput::getValFromInputParam("tgammaf", cellIndex);
		params.Rttover = DataInput::getValFromInputParam("rttover", cellIndex);
		params.MaxCond = DataInput::getValFromInputParam("maxcond", cellIndex);
		params.CoeffCond = DataInput::getValFromInputParam("coeffcond", cellIndex);
		params.BLcond = DataInput::getValFromInputParam("blcond", cellIndex);
		params.m0 = DataInput::getValFromInputParam("m0", cellIndex);
		params.fN0 = DataInput::getValFromInputParam("fn0", cellIndex);
		params.fNn = DataInput::getValFromInputParam("fnn", cellIndex);
		params.thinPower = DataInput::getValFromInputParam("thinpower", cellIndex);
		params.mF = DataInput::getValFromInputParam("mf", cellIndex);
		params.mR = DataInput::getValFromInputParam("mr", cellIndex);
		params.mS = DataInput::getValFromInputParam("ms", cellIndex);
		params.SWconst0 = DataInput::getValFromInputParam("swconst0", cellIndex);
		params.SWpower0 = DataInput::getValFromInputParam("swpower0", cellIndex);
		params.wSx1000 = DataInput::getValFromInputParam("wsx1000", cellIndex);
		params.MaxAge = DataInput::getValFromInputParam("maxage", cellIndex);
		params.nAge = DataInput::getValFromInputParam("nage", cellIndex);
		params.rAge = DataInput::getValFromInputParam("rage", cellIndex);
		params.SLA0 = DataInput::getValFromInputParam("sla0", cellIndex);
		params.SLA1 = DataInput::getValFromInputParam("sla1", cellIndex);
		params.tSLA = DataInput::getValFromInputParam("tsla", cellIndex);
		params.k = DataInput::getValFromInputParam("k", cellIndex);
		params.fullCanAge = DataInput::getValFromInputParam("fullcanage", cellIndex);
		params.alpha = DataInput::getValFromInputParam("alpha", cellIndex);
		params.fracBB0 = DataInput::getValFromInputParam("fracbb0", cellIndex);
		params.fracBB1 = DataInput::getValFromInputParam("fracbb1", cellIndex);
		params.tBB = DataInput::getValFromInputParam("tbb", cellIndex);
		params.y = DataInput::getValFromInputParam("y", cellIndex);
		params.rhoMin = DataInput::getValFromInputParam("rhomin", cellIndex);
		params.rhoMax = DataInput::getValFromInputParam("rhomax", cellIndex);
		params.tRho = DataInput::getValFromInputParam("trho", cellIndex);
		params.Qa = DataInput::getValFromInputParam("qa", cellIndex);
		params.Qb = DataInput::getValFromInputParam("qb", cellIndex);
		params.gDM_mol = DataInput::getValFromInputParam("gdm_mol", cellIndex);
		params.molPAR_MJ = DataInput::getValFromInputParam("molpar_mj", cellIndex);
		params.LAIgcx = DataInput::getValFromInputParam("laigcx", cellIndex);
		params.MaxIntcptn = DataInput::getValFromInputParam("maxintcptn", cellIndex);
		params.LAImaxIntcptn = DataInput::getValFromInputParam("laimaxintcptn", cellIndex);
		params.Lat = DataInput::getValFromInputParam("lat", cellIndex);
		params.FRp = DataInput::getValFromInputParam("frp", cellIndex);
		params.FRstart = DataInput::getValFromInputParam("frstart", cellIndex);
		params.FRend = DataInput::getValFromInputParam("frend", cellIndex);
		params.FRdec = DataInput::getValFromInputParam("frdec", cellIndex);
		params.soilIndex = DataInput::getValFromInputParam("soilindex", cellIndex);
		params.MaxASW = DataInput::getValFromInputParam("maxasw", cellIndex);
		params.MinASWp = DataInput::getValFromInputParam("minaswp", cellIndex);
		params.StartAge = DataInput::getValFromInputParam("startage", cellIndex);
		params.EndYear = DataInput::getValFromInputParam("endyear", cellIndex);
		params.StartMonth = DataInput::getValFromInputParam("startmonth", cellIndex);
		params.yearPlanted = DataInput::getValFromInputParam("yearplanted", cellIndex);
		params.SeedlingMass = DataInput::getValFromInputParam("seedlingmass", cellIndex);
		params.WFi = DataInput::getValFromInputParam("wfi", cellIndex);
		params.WRi = DataInput::getValFromInputParam("wri", cellIndex);
		params.WSi = DataInput::getValFromInputParam("wsi", cellIndex);
		params.StemNoi = DataInput::getValFromInputParam("stemnoi", cellIndex);
		params.ASWi = DataInput::getValFromInputParam("aswi", cellIndex);
		params.MinASWTG = DataInput::getValFromInputParam("minaswtg", cellIndex);
		params.NDVI_FPAR_intercept = DataInput::getValFromInputParam("ndvi_fpar_intercept", cellIndex);
		params.NDVI_FPAR_constant = DataInput::getValFromInputParam("ndvi_fpar_constant", cellIndex);
		params.CO2Start = DataInput::getValFromInputParam("co2start", cellIndex);
		params.CO2End = DataInput::getValFromInputParam("co2end", cellIndex);
		params.fCalpha700 = DataInput::getValFromInputParam("fcalpha700", cellIndex);
		params.fCg700 = DataInput::getValFromInputParam("fcg700", cellIndex);

		if (params.yearPlanted < 1 || params.StemNoi < 1) {
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
	//if the param doesn't exist, set to -1
	if (!this->inputParams.contains(paramName)) {
		//use predefined default parameters for Lat, rhoMax, rhoMin, and tRho if not set by user
		if (paramName == "Lat") { return 1000; }
		else if (paramName == "rhomax") { return 0.5; }
		else if (paramName == "rhomin") { return 0.5; }
		else if (paramName == "trho") { return 4; }
		else {
			//otherwise, set to 0
			return 0;
		}
	}

	PPPG_PARAM* param = this->inputParams.at(paramName).get();

	//if the param is scalar, return it's value
	if (param->spType == pScalar) {
		return param->val;
	}

	//if the param is a grid, return the value at the row and column specified
	if (param->spType == pTif) {
		double val = param->g->GetVal(cellIndex);
		if (param->g->IsNoData((float)val)) {
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
	PPPG_PARAM* monthParam = param->monthlyParams[monthIndex].get();

	if (monthParam->spType == pScalar) {
		return monthParam->val;
	}

	//get parameter value depending on whether it is scalar or grid
	if (monthParam->spType == pTif) {
		double val = monthParam->g->GetVal(cellIndex);
		if (monthParam->g->IsNoData((float)val)) {
			throw std::runtime_error("nan");
		}
		else {
			return val;
		}
	}

	throw std::exception("a parameter has been set incorrectly as neither a scalar or a grid.");
}

bool DataInput::getManagementParam(ManagementIndex index, long cellIndex, int year, double& val) {
	PPPG_MT_PARAM* table = &this->managementTables[index];
	
	//if management parameter doesn't exist, return false
	if (table->firstYear == -1 || year < table->firstYear) {
		return false;
	}

	//get index to management table vector for correct year
	int yearIndex;
	if (year - table->firstYear >= table->yearToIndex.size()) {
		yearIndex = table->yearToIndex.back();
	}
	else {
		yearIndex = table->yearToIndex[year - table->firstYear];
	}
	
	//get parameter
	PPPG_PARAM* param = table->yearlyParams[yearIndex].get();

	//get val if scalar and return true
	if (param->spType == pScalar) {
		val = param->val;
		return true;
	}
	
	//get val from grid, return true if not nan
	if (param->spType == pTif) {
		val = param->g->GetVal(cellIndex);
		return !param->g->IsNoData((float)val);
	}
	
	throw std::exception("a parameter has been set incorrectly as neither a scalar or a grid.");
}

std::unordered_map<std::string, PPPG_OP_VAR> DataInput::getOpVars() {
	return this->outputParams;
}

RefGridProperties DataInput::getRefGrid() {
	return this->refGrid;
}

bool DataInput::haveNetRadParam() {
	return this->haveNetRad;
}

RunPeriod DataInput::getRunPeriod() {
	RunPeriod runPeriod;
	// don't re-compute the run period if it's already been found
	if (RPStartYear.has_value() || RPEndYear.has_value()) {
		runPeriod.StartYear = RPStartYear;
		runPeriod.EndYear = RPEndYear;
		return runPeriod
	}
	else {
		PPPG_PARAM* yearPlantedParam = this->inputParams.at("yearplanted").get();
		PPPG_PARAM* startAgeParam = this->inputParams.at("startage").get();
		PPPG_PARAM* endYearParam = this->inputParams.at("endyear").get();
		PPPG_PARAM* startMonthParam = this->inputParams.at("startmonth").get();
		//get maxes and mins depending on whether they're scalar or grid parameters
		int yearPlantedMin = (yearPlantedParam->spType == pScalar) ? static_cast<int>(yearPlantedParam->val) : static_cast<int>(yearPlantedParam->g->GetMin());
		int startAgeMin = (startAgeParam->spType == pScalar) ? static_cast<int>(startAgeParam->val) : static_cast<int>(startAgeParam->g->GetMin());
		int endYearMax = (endYearParam->spType == pScalar) ? static_cast<int>(endYearParam->val) : static_cast<int>(endYearParam->g->GetMax());
		int startMonthMax = (startMonthParam->spType == pScalar) ? static_cast<int>(startMonthParam->val) : static_cast<int>(startMonthParam->g->GetMax());
		/* 
		determine start year of the run period
		*/
		if (startAgeParam->spType != pScalar && yearPlantedParam->spType != pScalar) {
			//if both are raster, find smallest sum of yearPlanted and startAge pixels.
			//we do this because the minimum sum (which is the year we should start on)
			//does not have to be at any of the pixels where yearPlanted or startAge are smallest
			double overallMin = std::numeric_limits<double>::max();
			for (int row = 0; row < yearPlantedParam->g->nRows; row++) {
				for (int col = 0; col < yearPlantedParam->g->nCols; col++) {
					//convert gotten values to double first so we don't have any overflow of floats as we add them
					double curMin = (double)yearPlantedParam->g->GetVal(row, col) + (double)startAgeParam->g->GetVal(row, col);
					
					//set the overall minimum accordingly
					overallMin = (curMin < overallMin) ? curMin : overallMin;
				}
			}
			runPeriod.StartYear = static_cast<int>(overallMin);
		}
		else {
			//otherwise, just add the mins together
			runPeriod.StartYear = yearPlantedMin + startAgeMin;
		}
		//max year is just the max year
		runPeriod.EndYear = endYearMax;
		/*
		error check on years
		*/
		if (runPeriod.StartYear > runPeriod.EndYear) {
			//if minimum year is larger than maximum year, print and log error
			std::string errstr = "min year (" + std::to_string(runPeriod.StartYear) + ") is greater than max year (" + std::to_string(runPeriod.EndYear) + ")";
			std::cout << errstr << std::endl;
			this->log(errstr);
			//then exit
			exit(EXIT_FAILURE);
		}
		//valid run period successfully determined
	}
	string runPeriodStr = "first run year = " + to_string(runPeriod.StartYear) + ", last run year = " + to_string(runPeriod.EndYear);
	this->log(runPeriodStr);
	return runPeriod;
}