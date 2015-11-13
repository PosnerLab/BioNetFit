/*
 * Swarm.hh
 *
 *  Created on: Jul 17, 2015
 *      Author: brandon
 */

#ifndef SWARM_HH_
#define SWARM_HH_

#include "Swarm.hh"

#include <unordered_map>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <random>
#include <iomanip>
#include <boost/filesystem.hpp>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <cstdio>
#include <sys/types.h>
#include <stdlib.h>
#include <unistd.h>

//#include "Model.hh"
class Model;


#include "Utils.hh"

class Data;
#include "Particle.hh"
#include "Pheromones.hh"

class Pheromones;
class Particle;

class Swarm {
public:
	Swarm(bool isMaster);
	void addExp(std::string path);

	void setModel(std::string path);
	Model * getModel() { return options_.model; };

	void setSwarmSize(int size);
	int getSwarmSize() { return options_.swarmSize; }

	void setIsMaster(bool isMaster) { isMaster_ = isMaster; }
	bool getIsMaster() { return isMaster_; }

	void setExePath(std::string path) { exePath_ = path; }
	void setConfigPath(std::string path) { configPath_ = path; }

	// TODO: Make get/set methods inline
	void setSwarmType(std::string type);
	void setSwarmSynchronicity(int synchronicity);
	void setSwarmGenerations(int generations);
	void setSwarmMinFit(float minfit);
	bool checkSwarmConsistency();

	void setVerbosity(int verbosity) { options_.verbosity = verbosity; }
	int getVerbosity() { return options_.verbosity; }

	void setUsePipes(bool usePipes) { options_.usePipes = usePipes; }
	bool getUsePipes() { return options_.usePipes; }

	void setUseCluster(bool useCluster);
	bool getUseCluster() { return options_.useCluster; }

	void setOutputDir(std::string path) { options_.outputDir = path; }
	std::string getOutputDir() { return options_.outputDir; }

	void setDivideByInit(bool divideByInit) { options_.divideByInit = divideByInit; }
	bool getDivideByInit() { return options_.divideByInit; }

	void setLogTransformSimData(bool logTransformSimData) { options_.logTransformSimData = logTransformSimData; }
	bool getLogTransformSimData() { return options_.logTransformSimData; }

	void setStandardizeSimData(bool standardizeSimData) { options_.standardizeSimData = standardizeSimData; }
	bool getStandardizeSimData() { return options_.standardizeSimData; }

	void setStandardizeExpData(bool standardizeExpData) { options_.standardizeExpData = standardizeExpData; }
	bool getStandardizeExpData() { return options_.standardizeExpData; }

	void setSosCalc(int objFunc) { options_.objFunc = objFunc; }
	int getSosCalc() { return options_.objFunc; }

	void setSwapRate(float swapRate) { options_.swapRate = swapRate; }
	void setExtraWeight(int extraWeight) { options_.extraWeight = extraWeight; }
	void setMaxRetryDifferentParents(int maxRetryDifferentParents) { options_.maxRetryDifferentParents = maxRetryDifferentParents; }
	void setForceDifferentParents(bool forceDifferentParents) { options_.forceDifferentParents = forceDifferentParents; }
	void setDeleteOldFiles(bool deleteOldFiles) { options_.deleteOldFiles = deleteOldFiles; }
	void setParallelCount(int parallelCount) { options_.parallelCount = parallelCount; }

	void setCurrentGen(int gen);

	void setType(std::string type) { type_ = type; }

	std::string getParticleBasePath() { return particleBasePath_; }

	void doSwarm();
	void runGeneration();
	void breedGeneration();
	void doParticle(int pID);
	void launchParticle(Particle *p);
	Particle *createParticle(int pID);

	std::string recvFromParticle(Particle *p);
	std::unordered_map<int,Particle*> generateInitParticles(int numParticles = -1);

	// Initialize our random number engine
	std::mt19937 randNumEngine;

	Pheromones *swarmComm_;
	int currentGeneration_ = 0;
	std::string type_;
	//std::unordered_map<Particle*,std::string> runningParticles_;

	std::set<int> runningParticles_;
	std::set<int>::iterator runningParticlesIterator_;

	std::set<int> failedParticles_;
	std::set<int>::iterator failedParticlesIterator_;

	std::string exePath_;
	std::string configPath_;
	const char * particleBasePath_ = "";

	std::unordered_map<int,Particle*> allParticles_;
	std::vector<std::vector<std::string>> currGenFits;
	std::multimap<double,std::string> allGenFits;
	//std::vector<std::map<double,std::unordered_map<std::string,double>>> allGenFits; // Vector of map of double/string pairs. Vector used for sorting. Map contains fit value mapped to string of params used to generate that fit value

	bool isMaster_;

	struct SwarmOpts {
		std::string swarmType; // genetic or swarm
		std::string outputDir; // Directory to use for output
		bool synchronicity; // 1 for synchronous

		Model * model; // the model file

		int maxGenerations = 10;
		int swarmSize = 10;
		int minFit = 0;
		int boostrap = 0;
		int parallelCount = 2;

		bool usePipes = false;
		bool useCluster = false;

		bool divideByInit = false;
		int logTransformSimData = 0;
		bool standardizeSimData = false;
		bool standardizeExpData = false;

		bool deleteOldFiles = true;

		int objFunc = 1;
		int extraWeight = 0;
		float swapRate = 0.5;
		bool forceDifferentParents = true;
		int maxRetryDifferentParents = 100;

		int verbosity;

		std::unordered_map<std::string,Data*> expFiles;
		std::unordered_map<std::string,std::vector<float> > mutationRates;
		std::unordered_map<std::string,std::vector<std::string> > initParams;
	};
	SwarmOpts options_;

private:
	void cleanupFiles(const char * path);
	bool sortFits(Particle * a, Particle * b);
	void finishFit();
	void getAllParticleParams();
	void outputRunSummary(std::string outputDir);
};

#endif /* SWARM_HH_ */
