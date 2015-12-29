/*
 * Particle.hh
 *
 *  Created on: Jul 19, 2015
 *      Author: brandon
 */

#ifndef PARTICLE_HH_
#define PARTICLE_HH_

#include "Model.hh"
#include "Swarm.hh"

class Swarm;
class Data;
class FreeParam;
class Model;
class FreeParam;
class Data;
class Pheromones;
class Particle;

class Particle {
public:
	Particle(Swarm * swarm, int id);
	void setModel(Model * model);
	void setParam(std::pair<std::string,double> myParams);
	std::map<std::string,double> getParams() { return simParams_; }
	void setID(int id);
	int getID() { return id_; }

	void doParticle();
	void generateParams();

	std::map<int, double> fitCalcs;

private:
	friend class boost::serialization::access;

	void runModel(int id = 0);
	void checkMessagesGenetic();
	void checkMessagesPSO();

	void calculateFit();
	double objFunc_chiSquare(double sim, double exp, double stdev);
	double objFunc_sumOfSquares(double sim, double exp, double dummyvar);
	double objFunc_divByMeasured(double sim, double exp, double dummyvar);
	double objFunc_divByMean(double sim, double exp, double mean);

	void initBreedWithParticle(int pID, int swapID);
	//void rcvBreedWithParticle(std::vector<std::string>& params, int reciprocateTo, int swapID, int pID);
	//double mutateParam(FreeParam* fp, double paramValue);
	void smoothRuns();
	void finalizeSim();

	typedef double (Particle::*objFunc)(double exp,double sim ,double stdev);

	objFunc objFuncPtr;

	Model * model_;
	std::map<std::string, double> simParams_;
	int id_;
	std::map<std::string, std::map<int, Data*>> dataFiles_;
	Swarm * swarm_;
	std::map<std::string, double> bestParams_;
	unsigned int currentGeneration_;

	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {

		ar & model_;
		ar & simParams_;
		ar & id_;
		ar & dataFiles_;
		ar & swarm_;
		ar & currentGeneration_;
	}
};
#endif /* PARTICLE_HH_ */
