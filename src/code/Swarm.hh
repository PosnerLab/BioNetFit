/*
 * Swarm.hh
 *
 *  Created on: Jul 17, 2015
 *      Author: brandon
 */

#ifndef SWARM_HH_
#define SWARM_HH_


#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <tr1/random>
#include <iomanip>
#include <chrono>
#include <cstdio>
#include <set>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "Utils.hh"
#include "Timer.hh"
#include "Particle.hh"
#include "Pheromones.hh"

class Model;
class FreeParam;
class Data;
class Pheromones;
class Particle;

#define MAX_LONG 9223372036854775807

// Forward declaration of class boost::serialization::access
namespace boost {
namespace serialization {
class access;
}
}

class Swarm {
public:
	Swarm(bool isMaster);
	Swarm();
	void addExp(std::string path);

	void setModel(std::string path);
	Model * getModel() { return options.model; };

	void setSwarmSize(int size);
	int getSwarmSize() { return options.swarmSize; }

	void setIsMaster(bool master) { isMaster = master; }
	bool getIsMaster() { return isMaster; }

	void setIsClusterInit(bool clusterInit) { isClusterInit = clusterInit; }

	void setExePath(std::string path) { exePath_ = path; }
	void setConfigPath(std::string path) { configPath_ = path; }

	void setSimPath(std::string path) { options.simPath = path; }

	// TODO: Make get/set methods inline
	void setSwarmType(std::string type);
	void setSwarmSynchronicity(int synchronicity);
	void setSwarmGenerations(int generations);
	void setSwarmMinFit(float minfit);
	bool checkSwarmConsistency();

	void setVerbosity(int verbosity) { options.verbosity = verbosity; }
	int getVerbosity() { return options.verbosity; }

	void setUsePipes(bool usePipes) { options.usePipes = usePipes; }
	bool getUsePipes() { return options.usePipes; }

	void setUseCluster(bool useCluster) { options.useCluster = useCluster; }
	bool getUseCluster() { return options.useCluster; }

	void setSaveClusterOutput(bool saveClusterOutput) { options.saveClusterOutput = saveClusterOutput; }

	void setOutputDir(std::string path) { options.outputDir = path; }
	std::string getOutputDir() { return options.outputDir; }

	void setDivideByInit(bool divideByInit) { options.divideByInit = divideByInit; }
	bool getDivideByInit() { return options.divideByInit; }

	void setLogTransformSimData(bool logTransformSimData) { options.logTransformSimData = logTransformSimData; }
	bool getLogTransformSimData() { return options.logTransformSimData; }

	void setStandardizeSimData(bool standardizeSimData) { options.standardizeSimData = standardizeSimData; }
	bool getStandardizeSimData() { return options.standardizeSimData; }

	void setStandardizeExpData(bool standardizeExpData) { options.standardizeExpData = standardizeExpData; }
	bool getStandardizeExpData() { return options.standardizeExpData; }

	void setSosCalc(int objFunc) { options.objFunc = objFunc; }
	int getSosCalc() { return options.objFunc; }

	void setSwapRate(float swapRate) { options.swapRate = swapRate; }
	void setExtraWeight(int extraWeight) { options.extraWeight = extraWeight; }
	void setMaxRetryDifferentParents(int maxRetryDifferentParents) { options.maxRetryDifferentParents = maxRetryDifferentParents; }
	void setForceDifferentParents(bool forceDifferentParents) { options.forceDifferentParents = forceDifferentParents; }
	void setDeleteOldFiles(bool deleteOldFiles) { options.deleteOldFiles = deleteOldFiles; }
	void setParallelCount(int parallelCount) { options.parallelCount = parallelCount; }
	void setJobName(std::string jobName) { options.jobName = jobName; }
	void addMutate(std::string mutateString);

	void setCurrentGen(int gen);

	void setType(std::string t) { type = t; }

	void doSwarm();
	void runGeneration();
	void breedGeneration();
	void doParticle(int pID);
	void launchParticle(int pID);

	Particle *createParticle(int pID);
	std::string recvFromParticle(Particle *p);
	std::map<int,Particle*> generateInitParticles(int numParticles = -1);
	std::vector<int> checkMasterMessages();

	void getClusterInformation();
	std::string generateSlurmCommand(std::string cmd, bool multiProg = true);
	std::string generateSlurmMultiProgCmd(std::string runCmd, std::string serializedSwarmPath);

	void initComm();
	void initFit();

	Pheromones *swarmComm;

	std::multimap<double,std::string> allGenFits;
	bool isMaster;
	bool isClusterInit;
	std::tr1::mt19937 randNumEngine;
	int currentGeneration;
	std::string type;

	struct SwarmOpts {
		// TODO: Need to define defaults AND check for required variables
		std::string jobName;	// name of the job
		std::string swarmType;	// genetic or swarm
		std::string outputDir;	// root directory to use for output
		std::string jobOutputDir;// outputDir + jobName
		std::string simPath;	// Path to simulators

		bool synchronicity;		// 1 for synchronous

		Model * model; 			// the model file

		int maxGenerations;// maximum number of generations
		int swarmSize;		// how many particles in the swarm
		float minFit;		// we won't accept any fits in breeding if they are over this value // TODO: Implement this
		float maxFit;		// we stop fitting if we reach this value // TODO: Implement this
		int boostrap;		// how many times to bootstrap
		int parallelCount;	// how many particles to run in parallel

		bool usePipes;	// whether or not to use pipes to gather simulation output
		bool useCluster;// whether or not we are running on a cluster

		bool divideByInit;// whether or not to divide simulation outputs by the value at t=0
		int logTransformSimData;// whether or not to log transform simulation data. this value acts as the base.
		bool standardizeSimData;// whether or not to standardize simulation data
		bool standardizeExpData;// whether or not to standardize experimental data

		bool deleteOldFiles; // whether or not to delete unneeded files during the fitting run

		int objFunc;		// which objective function to use
		int extraWeight;	// how much extra weight to add while breeding in genetic algorithm
		float swapRate;	// the rate at which to swap parent parameters during breeding
		bool forceDifferentParents;// whether or not to force difference parents when breeding
		int maxRetryDifferentParents;// how many times to attempt selection of different parents if forceDifferentParents is true
		long maxFitTime;
		long maxNumSimulations;

		int verbosity;		// terminal output verbosity

		bool hasMutate; // whether or not we should be mutating parameters during breeding in genetic algorithm

		std::string clusterSoftware;// which cluster software to use
		std::string clusterAccount;	// user account to specify in cluster submission commands // TODO: Parse
		bool saveClusterOutput;		// whether or not to save output during a cluster fit // TODO: Parse
		std::string clusterQueue;	// The cluster queue to submit to // TODO: Parse

		std::map<std::string,Data*> expFiles; // experimental data file

		template<class Archive>
		void serialize(Archive &ar, const unsigned int version)
		{
			//std::cout << " serializing options" << std::endl;

			ar & jobName;	// name of the job

			ar & swarmType;	// genetic or swarm

			ar & outputDir;	// root directory to use for output
			ar & jobOutputDir;// outputDir + jobName
			ar & simPath;

			ar & synchronicity;		// 1 for synchronous

			ar & model; 			// the model file

			ar & maxGenerations;// maximum number of generations
			ar & swarmSize;		// how many particles in the swarm
			ar & minFit;		// we won't accept any fits in breeding if they are over this value // TODO: Implement this
			ar & maxFit;		// we stop fitting if we reach this value // TODO: Implement this
			ar & boostrap;		// how many times to bootstrap
			ar & parallelCount;	// how many particles to run in parallel

			ar & usePipes;	// whether or not to use pipes to gather simulation output
			ar & useCluster;// whether or not we are running on a cluster

			ar & divideByInit;// whether or not to divide simulation outputs by the value at t=0
			ar & logTransformSimData;// whether or not to log transform simulation data. this value acts as the base.
			ar & standardizeSimData;// whether or not to standardize simulation data
			ar & standardizeExpData;// whether or not to standardize experimental data

			ar & deleteOldFiles; // whether or not to delete unneeded files during the fitting run

			ar & objFunc;		// which objective function to use
			ar & extraWeight;	// how much extra weight to add while breeding in genetic algorithm
			ar & swapRate;	// the rate at which to swap parent parameters during breeding
			ar & forceDifferentParents;// whether or not to force difference parents when breeding
			ar & maxRetryDifferentParents;// how many times to attempt selection of different parents if forceDifferentParents is true

			ar & verbosity;		// terminal output verbosity

			ar & hasMutate; // whether or not we should be mutating parameters during breeding in genetic algorithm

			ar & clusterSoftware;// which cluster software to use
			ar & clusterAccount;	// user account to specify in cluster submission commands // TODO: Parse
			ar & saveClusterOutput;		// whether or not to save output during a cluster fit // TODO: Parse
			ar & clusterQueue;	// The cluster queue to submit to // TODO: Parse

			ar & expFiles; // experimental data file

		}
	};
	SwarmOpts options;

private:
	friend class boost::serialization::access;

	void cleanupFiles(const char * path);
	bool sortFits(Particle * a, Particle * b);
	void finishFit();
	void getAllParticleParams();
	void outputRunSummary(std::string outputDir);
	void killAllParticles(int tag);

	std::set<int> runningParticles_;
	std::set<int>::iterator runningParticlesIterator_;

	std::set<int> failedParticles_;
	std::set<int>::iterator failedParticlesIterator_;

	std::string exePath_;
	std::string configPath_;
	//const char * particleBasePath_ = "";

	std::map<int,Particle*> allParticles_;

	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		//std::cout << " serializing swarm" << std::endl;

		ar & allGenFits;
		ar & configPath_;
		ar & currentGeneration;
	    ar & exePath_;
	    ar & options;
	}

	Timer tmr_;
};

#endif /* SWARM_HH_ */
