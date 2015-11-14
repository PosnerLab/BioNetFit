/*
 * Particle.hh
 *
 *  Created on: Jul 19, 2015
 *      Author: brandon
 */

#ifndef PARTICLE_HH_
#define PARTICLE_HH_

#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <chrono>

#include <dirent.h>

#include <boost/filesystem.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>

namespace io = boost::iostreams;

#include "Model.hh"
#include "Swarm.hh"

class Swarm;

class Particle {
public:
	Particle(Swarm * swarm, int id);
	void setModel(Model * model);
	void setParam(std::pair<std::string,double> myParams);
	std::map<std::string,double> getParams() { return simParams_; }
	void setID(int id);
	int getID() { return id_; }
	void doParticle();
	void calculateFit();
	void setBasePath(std::string path) { basePath_ = path; }
	void generateParams();

	std::unordered_map<int,double> fitCalcs;

private:
	typedef double (Particle::*objFunc)(double exp,double sim ,double stdev);
	objFunc objFuncPtr;

	Model * model_;
	std::map<std::string,double> simParams_;
	int id_;
	std::unordered_map<std::string,Data*> dataFiles_;
	std::string basePath_;
	std::string state_; // Stopped, running, simulating, analyzing, results, breeding, waiting
	Swarm * swarm_;

	double objFunc_chiSquare(double sim, double exp, double stdev);
	double objFunc_sumOfSquares(double sim, double exp, double dummyvar);
	double objFunc_divByMeasured(double sim, double exp, double dummyvar);
	double objFunc_divByMean(double sim, double exp, double mean);

	void initBreedWithParticle(int pID, int swapID);
	void rcvBreedWithParticle(std::vector<std::string>& params, int reciprocate, int swapID, int pID);
};
#endif /* PARTICLE_HH_ */
