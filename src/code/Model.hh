/*
 * Model.hh
 *
 *  Created on: Jul 13, 2015
 *      Author: brandon
 */

#ifndef MODEL_HH_
#define MODEL_HH_

#include <iostream>
#include <vector>
#include <fstream>
#include <boost/regex.hpp>
#include <cstdlib>

#include "Utils.hh"
#include "Data.hh"
#include "FreeParam.hh"

class FreeParam;
class Data;

class Model {
	friend class Config;
	friend class Swarm;
	friend class Particle;

public:
	Model(std::string path);
	Model();

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
		std::string scanParam;

		Data *dataSet;

		template<typename Archive>
		void serialize(Archive& ar, const unsigned version) {

			ar & full;
			ar & t_end;
			ar & type;
			ar & scanParam;

			ar & dataSet;
		}
	};

	const std::map<std::string, FreeParam*>& getfreeParams_() const {
		return freeParams_;
	}

	// Map key contains the action prefix, map value contains the action information
	std::map<std::string, action> actions;

private:
	friend class boost::serialization::access;

	void parseModel();

	std::string modelPath_;
	std::vector<std::string> fullContents_;
	std::vector<std::string> netContents_;

	bool hasGenerateNetwork_;
	//std::map<std::string,std::string> freeParams_;
	std::map<std::string, FreeParam*> freeParams_;

	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		//std::cout << " serializing model" << std::endl;

		ar & freeParams_;

		ar & actions;
		ar & modelPath_;
		ar & fullContents_;
		ar & netContents_;

		ar & hasGenerateNetwork_;

	}
};

#endif /* MODEL_HH_ */
