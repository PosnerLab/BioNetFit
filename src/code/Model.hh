/*
 * Model.hh
 *
 *  Created on: Jul 13, 2015
 *      Author: brandon
 */

#ifndef MODEL_HH_
#define MODEL_HH_

#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <regex>
#include <cstdlib>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <fcntl.h>

#include "Utils.hh"
#include "Data.hh"

class Model {
friend class Config;
friend class Swarm;
friend class Particle;

public:
	Model(std::string path);

	std::string getLocation();
	int makeCopyofOrig(std::string newLocation);
	void outputModelWithParams(std::map<std::string,double> params, std::string path, std::string filename, std::string suffix, bool stopAtNetGen, bool onlyActions, bool netAndBngl, bool usePipe, bool isNetFile);
	void parseNet(std::string path);
	bool getHasGenerateNetwork() {return hasGenerateNetwork_;}

	struct action {
			std::string full;

			//double t_start;
			double t_end;
			//double o_steps;

			std::string type;
			std::string scanParam = "";

			Data *dataSet;
		};

	// Map key contains the action prefix, map value contains the action information
	std::unordered_map<std::string,action > actions;

private:
	void parseModel();

	std::string modelPath_;
	std::vector<std::string> fullContents_;
	std::vector<std::string> netContents_;

	bool hasGenerateNetwork_ = false;
	std::map<std::string,std::string> freeParams_;
};

#endif /* MODEL_HH_ */
