/*
 * Swarm.hh
 *
 *  Created on: Jul 17, 2015
 *      Author: brandon
 */

#ifndef SWARM_HH_
#define SWARM_HH_

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
#include <chrono>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <cstdio>
#include <sys/types.h>
#include <stdlib.h>
#include <unistd.h>

//#include "Model.hh"
class Model;
class FreeParam;


#include "Utils.hh"
#include "Timer.hh"
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
	void setJobName(std::string jobName) { options_.jobName = jobName; }
	void addMutate(std::string mutateString);

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

	void getClusterInformation();
	std::string generateSlurmCommand(std::string cmd);

	Pheromones *swarmComm_;
	int currentGeneration_ = 0;	std::string type_;
	//std::unordered_map<Particle*,std::string> runningParticles_;

	std::set<int> runningParticles_;
	std::set<int>::iterator runningParticlesIterator_;

	std::set<int> failedParticles_;
	std::set<int>::iterator failedParticlesIterator_;

	std::string exePath_;
	std::string configPath_;
	const char * particleBasePath_ = "";

	std::unordered_map<int,Particle*> allParticles_;
	//std::vector<std::vector<std::string>> currGenFits;
	std::multimap<double,std::string> allGenFits;
	//std::vector<std::map<double,std::unordered_map<std::string,double>>> allGenFits; // Vector of map of double/string pairs. Vector used for sorting. Map contains fit value mapped to string of params used to generate that fit value

	bool isMaster_;

	struct SwarmOpts {
		std::string jobName;	// name of the job
		std::string swarmType;	// genetic or swarm
		std::string outputDir;	// root directory to use for output
		std::string jobOutputDir;// outputDir + jobName

		bool synchronicity;		// 1 for synchronous

		Model * model; 			// the model file

		int maxGenerations = 10;// maximum number of generations
		int swarmSize = 10;		// how many particles in the swarm
		float minFit = 0;		// we won't accept any fits in breeding if they are over this value // TODO: Implement this
		float maxFit = 0;		// we stop fitting if we reach this value // TODO: Implement this
		int boostrap = 0;		// how many times to bootstrap
		int parallelCount = 2;	// how many particles to run in parallel

		bool usePipes = false;	// whether or not to use pipes to gather simulation output
		bool useCluster = false;// whether or not we are running on a cluster

		bool divideByInit = false;// whether or not to divide simulation outputs by the value at t=0
		int logTransformSimData = 0;// whether or not to log transform simulation data. this value acts as the base.
		bool standardizeSimData = false;// whether or not to standardize simulation data
		bool standardizeExpData = false;// whether or not to standardize experimental data

		bool deleteOldFiles = true; // whether or not to delete unneeded files during the fitting run

		int objFunc = 1;		// which objective function to use
		int extraWeight = 0;	// how much extra weight to add while breeding in genetic algorithm
		float swapRate = 0.5;	// the rate at which to swap parent parameters during breeding
		bool forceDifferentParents = true;// whether or not to force difference parents when breeding
		int maxRetryDifferentParents = 100;// how many times to attempt selection of different parents if forceDifferentParents is true

		int verbosity = 1;		// terminal output verbosity

		bool hasMutate = false; // whether or not we should be mutating parameters during breeding in genetic algorithm

		std::string clusterSoftware;// which cluster software to use
		std::string clusterAccount;	// user account to specify in cluster submission commands // TODO: Parse
		bool saveClusterOutput;		// whether or not to save output during a cluster fit // TODO: Parse
		std::string clusterQueue;	// The cluster queue to submit to // TODO: Parse

		std::unordered_map<std::string,Data*> expFiles; // experimental data file
	};
	SwarmOpts options_;

	std::mt19937 randNumEngine;

private:
	void cleanupFiles(const char * path);
	bool sortFits(Particle * a, Particle * b);
	void finishFit();
	void getAllParticleParams();
	void outputRunSummary(std::string outputDir);
	void killAllParticles(int tag);

	Timer tmr_;
};

#endif /* SWARM_HH_ */
