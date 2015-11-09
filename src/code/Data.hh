/*
 * Exp.hh
 *
 *  Created on: Jul 17, 2015
 *      Author: brandon
 */

#ifndef DATA_HH_
#define DATA_HH_

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>

#include "Utils.hh"
#include "Swarm.hh"

#define NAN 123456789

class Swarm;
class Data {

public:
	Data(std::string path, Swarm * swarm, bool isExp);
	void parseData();
	std::string getPath();

	void standardizeData();
	void divideByInit();
	void getColumnAverages();
	void logTransformData();

	std::map<std::string,std::map<double,double> > * dataCurrent_; // This points to the most recently modified version of the data
	std::map<std::string,std::map<double,double> > standardDeviations_;
	std::map<std::string,double> colAverages_;

private:
	std::string dataPath_;

	std::map<std::string,std::map<double,double> > dataOrig_;
	std::map<std::string,std::map<double,double> > dataStandardized_;
	std::map<std::string,std::map<double,double> > dataLogTransformed_;
	std::map<std::string,std::map<double,double> > dataDividedByInit_;

	Swarm * swarm_;
	bool isExp_;
};

#endif /* DATA_HH_ */
