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
#include <string>
#include <map>

#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
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
	Swarm();

	void addExp(std::string path);
	void setModel(std::string path);
	void setExePath(std::string path) { exePath_ = path; }
	void setConfigPath(std::string path) { configPath_ = path; }
	void setfitType(std::string type);
	void addMutate(std::string mutateString);
	void setsConf(std::string sConf) { sConf_ = sConf; }
	std::string getsConf() { return sConf_; }
	void setJobOutputDir(std::string dir);

	void initComm();
	void doSwarm();
	Particle *createParticle(unsigned int pID);

	void getClusterInformation();
	std::string generateSlurmCommand(std::string cmd, bool multiProg = true);
	std::string generateSlurmMultiProgCmd(std::string runCmd);
	std::string generateSlurmBatchFile(std::string runCmd);

	Pheromones *swarmComm;

	std::multimap<double, std::string> allGenFits;

	bool isMaster;
	boost::random::mt19937 randNumEngine;
	unsigned int currentGeneration;

	struct SwarmOpts {
		// TODO: Need to define defaults AND check for required variables
		std::string jobName;	// name of the job
		std::string fitType;	// genetic or swarm
		std::string outputDir;	// root directory to use for output
		std::string jobOutputDir;// outputDir + jobName
		std::string bngCommand;	// Path to simulators

		unsigned int outputEvery;

		Model * model; 			// the model file

		bool synchronicity;		// 1 for synchronous
		unsigned int maxGenerations;// maximum number of generations
		unsigned int swarmSize;		// how many particles in the swarm
		float minFit;		// we won't accept any fits in breeding if they are over this value // TODO: Implement this
		float maxFit;		// we stop fitting if we reach this value // TODO: Implement this
		unsigned int boostrap;		// how many times to bootstrap
		unsigned int parallelCount;	// how many particles to run in parallel

		bool usePipes;	// whether or not to use pipes to gather simulation output
		bool useCluster;// whether or not we are running on a cluster

		bool divideByInit;// whether or not to divide simulation outputs by the value at t=0
		int logTransformSimData;// whether or not to log transform simulation data. this value acts as the base.
		bool standardizeSimData;// whether or not to standardize simulation data
		bool standardizeExpData;// whether or not to standardize experimental data

		bool deleteOldFiles; // whether or not to delete unneeded files during the fitting run

		unsigned int objFunc;		// which objective function to use

		// Genetic algorithm options
		unsigned int extraWeight;	// how much extra weight to add while breeding in genetic algorithm
		float swapRate;	// the rate at which to swap parent parameters during breeding
		bool forceDifferentParents;// whether or not to force difference parents when breeding
		unsigned int maxRetryDifferentParents;// how many times to attempt selection of different parents if forceDifferentParents is true
		unsigned int smoothing;
		unsigned int keepParents;

		unsigned long maxFitTime;	// Maximum amount of time to let the fit run
		unsigned long maxNumSimulations; // Maximum number of simulations to run
		unsigned long maxNumIterations; // Maximum number of iterations a particle can run // TODO: Implement

		// PSO options
		float inertia; // 0.72
		float cognitive; // 1.49
		float social; // 1.49
		unsigned int nmax; // 20
		unsigned int nmin; // 80
		float inertiaInit; // 1
		float inertiaFinal; // 0.1
		float absTolerance; // 10E-4
		float relTolerance; // 10E-4

		std::string topology; // fullyconnected
		std::string psoType; // bbpso

		bool enhancedStop; // true
		bool enhancedInertia; // true

		unsigned int verbosity;		// terminal output verbosity

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

			ar & fitType;	// genetic or swarm

			ar & outputDir;	// root directory to use for output
			ar & jobOutputDir;// outputDir + jobName
			ar & bngCommand;

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
			ar & smoothing;		// How many simulations to average

			ar & maxNumSimulations;
			ar & maxNumIterations;

			// PSO options
			ar & inertia; // 0.72
			ar & cognitive; // 1.49
			ar & social; // 1.49
			ar & nmax; // 20
			ar & nmin; // 80
			ar & inertiaInit; // 1
			ar & inertiaFinal; // 0.1
			ar & absTolerance; // 10E-4
			ar & relTolerance; // 10E-4
			ar & topology;
			ar & psoType;
			ar & enhancedStop;
			ar & enhancedInertia;

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

	void initFit();
	std::vector<std::vector<unsigned int>> generateInitParticles();
	void launchParticle(unsigned int pID);
	void runGeneration();
	void breedGeneration();

	void cleanupFiles(const char * path);
	void finishFit();
	void getAllParticleParams();
	void outputRunSummary(std::string outputDir);
	void outputRunSummary();
	void killAllParticles(int tag);
	std::vector<unsigned int> checkMasterMessages();
	void checkExternalMessages();

	void processParticlesPSO(std::vector<unsigned int> particles, bool nextFlight = false);
	void updateEnhancedStop();
	double getEuclidianNorm(double y, unsigned int n);
	void updateParticleWeights();
	double calcParticleWeight(unsigned int particle);
	double calcWeightedAveragePosition();
	void processParamsPSO(std::vector<double> &params, unsigned int pID, double fit);
	bool checkStopCriteria();
	void updateInertia();

	std::vector<double> calcParticlePosPSO(unsigned int particle);
	std::vector<double> calcParticlePosBBPSO(unsigned int particle, bool exp = false);
	std::vector<double> getNeighborhoodBestPositions(unsigned int particle);

	unsigned int pickWeighted(double weightSum, std::multimap<double, unsigned int> &weights, unsigned int extraWeight);
	double mutateParam(FreeParam* fp, double paramValue);

	void insertKeyByValue(std::map<double, int> &theMap, double key, int value);

	std::set<int> runningParticles_;
	std::set<int>::iterator runningParticlesIterator_;

	std::set<int> failedParticles_;
	std::set<int>::iterator failedParticlesIterator_;

	std::string exePath_;
	std::string configPath_;
	std::string sConf_;

	std::vector<std::vector<unsigned int> > allParticles_;

	// TODO: These need to be initialized with 0s
	// Maybe we can change them to vectors, too
	std::map<unsigned int, double> particleBestFits_;
	std::multimap<double, unsigned int> particleBestFitsByFit_;
	std::map<unsigned int, std::vector<double>> particleParamVelocities_;
	std::map<unsigned int, std::vector<double>> particleBestParamSets_;
	std::map<unsigned int, std::vector<double>> particleCurrParamSets_;
	std::map<unsigned int, double> particleWeights_;
	std::map<unsigned int, unsigned int> particleIterationCounter_;

	unsigned int permanenceCounter_; // 0
	unsigned int flightCounter_; // 0
	double weightedAvgPos_; // 0
	double optimum_; // 0
	unsigned int inertiaUpdateCounter_; // 0;

	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		//std::cout << " serializing swarm" << std::endl;

		ar & allGenFits;
		ar & configPath_;
		ar & currentGeneration;
		ar & exePath_;
		ar & options;
		ar & particleBestFits_;
		ar & sConf_;

		ar & allParticles_;
		ar & particleBestFitsByFit_;
		ar & particleParamVelocities_;
		ar & particleBestParamSets_;
		ar & particleCurrParamSets_;
		ar & particleWeights_;
		ar & particleBestFits_;

		ar & permanenceCounter_;
		ar & flightCounter_;
		ar & weightedAvgPos_;
		ar & optimum_;
		ar & inertiaUpdateCounter_;
	}

	Timer tmr_;
};

#endif /* SWARM_HH_ */
