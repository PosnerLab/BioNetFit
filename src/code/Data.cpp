/*
 * Exp.cpp
 *
 *  Created on: Jul 17, 2015
 *      Author: brandon
 */


//TODO: Add checks for sim length equality
//TODO: Test all these methods
//TODO: Is sos correct terminology? Maybe objFunc would be better
#include "Data.hh"

using namespace std;

Data::Data(std::string path, Swarm * swarm, bool isExp) {
	swarm_ = swarm;
	isExp_ = isExp;

	dataPath = convertToAbsPath(path);
	Data::parseData();

	if ( (swarm_->options.objFunc == 4 && isExp)) {
		//cout << "soscalc is 4" << endl;
		getColumnAverages();
	}

	if (swarm_->options.divideByInit && !isExp) {
		//cout << "div by init" << endl;
		divideByInit();
		dataCurrent = &dataDividedByInit_;
	}

	if (swarm_->options.logTransformSimData > 0 && !isExp) {
		//cout << "logsim" << endl;
		logTransformData();
		dataCurrent = &dataLogTransformed_;
	}

	if (swarm_->options.standardizeSimData && !isExp) {
		cout << "stand sim" << endl;
		if (swarm_->options.objFunc != 4) {
			//cout << "gca sim" << endl;
			getColumnAverages();
		}
		standardizeData();
		dataCurrent = &dataStandardized_;
	}

	if (swarm_->options.standardizeExpData && isExp) {
		cout << "stand exp" << endl;
		if (swarm_->options.objFunc != 4) {
			//cout << "gca exp" << endl;
			getColumnAverages();
		}
		standardizeData();
		dataCurrent = &dataStandardized_;
	}
}

Data::Data(map<string, map<double, double>> &dataSet) {
	dataOrig_ = dataSet;
	dataCurrent = &dataOrig_;
}

Data::Data() {}

std::string Data::getPath() {
	if (!dataPath.empty())
		return dataPath;
	else
		return "empty";
}

void Data::parseData(){
	vector<string> allLines;

	if (swarm_->options.usePipes){
		char buf[5120];
		string line;

		int fd = open(dataPath.c_str(), O_RDONLY);

		if (fd < 0) {
			cout << "couldn't open " << dataPath << endl;
		}

		int len;
		while( (len = read(fd, buf, 5120) )){
			buf[len] = 0;
			line = string(buf);
			//cout << line << "!";
			split(line, allLines, "\n");
			//printf("%s", buf);
		}

		if (len < 0) {
			perror("read error");
		}
		close(fd);
	}
	else {

		ifstream dataFile(dataPath);

		std::string line;
		std::string errMsg;

		if (dataFile.is_open()) {
			while (getline(dataFile, line)) {
				allLines.push_back(line);
				//cout << "line: " << line << endl;
			}
		}
		else {
			errMsg = "Error: Couldn't open data file " + dataPath + " for parsing.";
			outputError(errMsg);
		}
		dataFile.close();
	}

	string replacement = "";
	string basename = boost::regex_replace(getFilename(dataPath),boost::regex("_\\d+_\\d+$"),replacement);

	if(allLines[0].at(0) != '#') {
		string errMsg = "Error: Your data file (" + dataPath + ") doesn't contain (#) as the first value.";
		outputError(errMsg);
	}

	// TODO: Do we need to store columns that aren't in .exp?  Maybe not.
	vector<string> columns;
	split(allLines[0], columns, " \t");

	for (vector<string>::iterator l = allLines.begin()+1; l != allLines.end(); ++l) {
		vector<string> values;
		split(*l, values," \t");

		int i = 1;
		for(vector<string>::iterator col = columns.begin(); col != columns.end(); ++col) {

			// Skip column if we encounter a hash, a 'time' column, or a column corresponding to a scan parameter name
			if (*col == "#" || *col == "time" || *col == swarm_->options.model->actions.at(basename).scanParam){
				continue;
			}

			// If our column name ends in "_SD" we need to store the value in our stdev map
			if (col->size() > 3 && col->find("_SD") == col->size()-3) {
				string colWithoutSD = boost::regex_replace(*col,boost::regex("_SD$"),string(""));
				//cout << "inserting sd at col " << colWithoutSD << ": " << stof(values[0]) << " " << stof(values[i]) << endl;
				//standardDeviations[col].insert(make_pair(stof(values[0]),stof(values[i])));
				standardDeviations[colWithoutSD][stod(values[0])] = stod(values[i]);
				//cout << "test " << colWithoutSD << ": " << standardDeviations[colWithoutSD][stof(values[0])] << endl;
				// TODO: Make sure all SD's have values if we're using objfunc 2, lest we div_by_0
				// TODO: Make sure all normal columns also have SD columns if we're using objfunc 2
			}
			// Insert value into data map
			else {
				//cout << "inserting val to col " << *col << ": " << values[0] << " " << values[i] << endl;
				dataOrig_[*col][stod(values[0])] = stod(values[i]);
				//dataOrig_[col].insert(make_pair(stof(values[0]),stof(values[i])));
			}
			i++;
		}
	}

	// Always keep dataCurrent pointing to the most recent version of the data.
	dataCurrent = &dataOrig_;
}

void Data::standardizeData() {
	// This function standardizes our data around a mean of 0

	double mean;
	double sqtotal;
	int counter;
	double stdev;

	// Loop through data map
	for(map<string,map<double,double> >::iterator c = dataCurrent->begin(); c != dataCurrent->end(); ++c) {

		// Grab our mean, which has already been calculated
		mean = colAverages[c->first];
		sqtotal = 0;
		counter = 0;

		// Skip the column if mean is already 0
		if (mean == 0) {
			continue;
		}

		// Here we get the riemann sum squared for use in the stdev calculation
		for(map<double,double>::iterator v = c->second.begin(); v != c->second.end(); ++v) {
			//cout << "sd loop" << endl;
			if (v->second == NaN) {
				continue;
			}
			sqtotal += pow(mean - v->second, 2);
			counter++;
		}

		// Cacualate the standard deviation
		stdev = pow((sqtotal / (counter - 1)), 0.5);

		// Loop through and set values according to formula: val_new = (val_old - column_mean) / (column_stdev)
		for(map<double,double>::iterator v = c->second.begin(); v != c->second.end(); ++v) {
			//cout << "stand loop" << endl;
			if (v->second == NaN) {
				dataStandardized_[c->first][v->first] = NaN;
				continue;
			}
			//cout << "inserting: " << (v->second - mean) / stdev << endl;
			dataStandardized_[c->first][v->first] = (v->second - mean) / stdev;
		}
	}
}

void Data::divideByInit() {
	for(map<string,map<double,double> >::iterator c = dataCurrent->begin(); c != dataCurrent->end(); ++c) {
		for(map<double,double>::iterator v = c->second.begin(); v != c->second.end(); ++v) {
			if (v->second == 0) {
				string errMsg = "You chose to divide_by_init, but the first value in the column '" + c->first + "' is 0. We cannot divide by 0.";
				outputError(errMsg);
			}
			if (v->second == NaN) {
				dataDividedByInit_[c->first][v->first] = NaN;
				continue;
			}
			dataDividedByInit_[c->first][v->first] = v->second / c->second.begin()->first;
			//cout << "dividing t=" << v->first << " of " << v->second << "by " << c->first << " init of " << c->second.begin()->first << endl;
		}
	}
}

void Data::getColumnAverages() {
	double sum;
	int counter;
	for(map<string,map<double,double> >::iterator c = dataCurrent->begin(); c != dataCurrent->end(); ++c) {
		cout << "gca col loop" << endl;

		sum = 0;
		counter = 0;

		for(map<double,double>::iterator v = c->second.begin(); v != c->second.end(); ++v) {
			//cout << "gca tp loop" << endl;
			if (v->second == NaN) {
				continue;
			}
			sum += v->second;
			counter++;
		}
		colAverages.insert(pair<string,double>(c->first,sum/counter));
	}
}

void Data::logTransformData() {
	for(map<string,map<double,double> >::iterator c = dataCurrent->begin(); c != dataCurrent->end(); ++c) {
		for(map<double,double>::iterator v = c->second.begin(); v != c->second.end(); ++v) {
			if (v->second == 0) {
				string errMsg = "You chose to log transform simulation output, but the first value in the column '" + c->first + "' is 0. We cannot take the log of 0.";
				outputError(errMsg);
			}
			else if (v->second == NaN) {
				dataLogTransformed_[c->first][v->first] = NaN;
				continue;
			}
			dataLogTransformed_[c->first][v->first] = log10(v->second)/log10(swarm_->options.logTransformSimData);
		}
	}
}
