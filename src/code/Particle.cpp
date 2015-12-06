/*
 * Particle.cpp
 *
 *  Created on: Jul 19, 2015
 *      Author: brandon
 */

// TODO: Let's replace all the static_cast's with a custom function that converts ints/doubles to strings using stringstream

#include "Particle.hh"

using namespace std;
using namespace std::chrono;

Particle::Particle(Swarm * swarm, int id) {
	id_ = id;
	swarm_ = swarm;
}

void Particle::setModel(Model * model) {
	model_ = model;
}

void Particle::setParam(pair<string,double> myParams) {
	simParams_.insert(myParams);
}

void Particle::setID(int id) {
	id_ = id;
}

void Particle::generateParams() {
	// TODO: See here: http://www.johndcook.com/blog/cpp_TR1_random/
	// for a possibly better way to generate numbers
	// freeParams_ is a map with the first element being a parameter name, and the second element being a pointer to a FreeParam object
	for (map<string,FreeParam*>::iterator i = model_->freeParams_.begin(); i != model_->freeParams_.end(); ++i) {

		//vector<string> values;

		// Split our generation method/range into parts.
		//split(i->second, values);

		// Make sure our array contains the generation method
		/*if (values[0].empty()){
			string errMsg = "Problem parsing initial parameter generator for param: " + i->first + ".";
			outputError (errMsg);
		}*/

		// Store the free parameter name that we're currently generating
		string paramName;
		paramName = i->first;

		// Store the generation method
		//string genType = values[0];
		string genType = i->second->getGenerationMethod();

		// Generate a number on a linear scale
		if (genType == "random_var"){
			//TODO: Replace the random number generators here with c++11 generators
			// Store our min and max values
			double min = i->second->getGenMin();
			double max = i->second->getGenMax();

			float myrand = min + static_cast<double> (rand()) /(static_cast<double> (RAND_MAX/(max-min)));

			pair<string,double> paramPair = make_pair(paramName, myrand);
			setParam(paramPair);
		}
		else if (genType == "loguniform_var") { //TODO: Test distribution of this generator

			// Store our min and max values
			double min = i->second->getGenMin();
			double max = i->second->getGenMax();

			//cout << "generating with " << min << ":" << max << endl;

			// Generate a random double between 0 and rand_max
			double myrand = static_cast <double> (rand()) / static_cast<double> (RAND_MAX);
			//cout << id_ << " rand is " << myrand << endl;

			//
			double exp = log10(min) + myrand * (log10(max) - log10(min) );
			//cout << id_ << " base is " << exp << endl;
			//cout << id_ << " final " << pow(10, exp) << endl << endl;

			pair<string,double> paramPair = make_pair(paramName, pow(10, exp));
			setParam(paramPair);
		}
		//TODO: Add the rest of the init parm generation options
	}
}

void Particle::doParticle() {

	if (swarm_->options.swarmType == "genetic") {
		runGenetic();
	}
	else if (swarm_->options.swarmType == "pso") {
		runPSO();
	}
}

void Particle::runGenetic() {
	//cout << "Particle " << id_ << " waiting to begin" << endl;
	swarm_->swarmComm->recvMessage(0, id_, 18, true, swarm_->swarmComm->univMessageReceiver);
	//cout << "Particle " << id_ << " starting" << endl;

	while(1) {
		// First get our path and filename variables set up for use in model generation, sim command, etc
		string bnglFilename = to_string(static_cast<long long int>(id_)) + ".bngl";
		string path = swarm_->options.jobOutputDir + "" + to_string(static_cast<long long int>(swarm_->currentGeneration)) + "/";
		string bnglFullPath = path + bnglFilename;

		string pipePath;

		// Only need to generate files if we're in the first generation. In subsequent generations
		// the model generation is handled by beeding parents
		if (swarm_->currentGeneration == 1) {
			if (swarm_->options.model->getHasGenerateNetwork()){
				string netFilename = "base.net";
				string netFullPath = swarm_->options.jobOutputDir + netFilename;

				swarm_->options.model->parseNet(netFullPath);
				model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(static_cast<long long int>(id_)), false, false, true, false, false);
			}
			else {
				model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(static_cast<long long int>(id_)), false, false, false, false, false);
			}
		}

		// Generate .gdat pipes if we're using pipes
		if (swarm_->options.usePipes) {
			for (std::map<std::string,Model::action>::iterator i = model_->actions.begin(); i != model_->actions.end(); ++i) {
				pipePath = path + i->first + "_" + to_string(static_cast<long long int>(id_)) + ".gdat";
				createParticlePipe(pipePath.c_str());
			}
		}

		// Construct our simulation command
		string command = swarm_->options.simPath + "BNG2.pl --outdir " + path + " " + bnglFullPath + ">> " + path + "/" + to_string(static_cast<long long int>(id_)) + ".BNG_OUT 2>&1";
		if (swarm_->options.usePipes) {
			command += " &";
		}

		if (swarm_->options.verbosity >= 3) {
			cout << "Running model with command: " << command << endl;
		}

		// Run simulation command
		int ret = system(command.c_str());

		// Check for simulation command success
		if (ret == 0) { // TODO: Need to check for simulation status when using pipes. Going by return code doesn't work there because we're using the & operator
			string dataPath;

			// Save our simulation outputs to data objects
			for (std::map<std::string,Model::action>::iterator i = model_->actions.begin(); i != model_->actions.end(); ++i) {
				dataPath = path + "/" + i->first + "_" + to_string(static_cast<long long int>(id_)) + ".gdat";
				dataFiles_[i->first] = new Data(dataPath, swarm_, false);
			}

			// Calculate our fit
			calculateFit();

			// Put our fit calc into the message vector
			swarm_->swarmComm->univMessageSender.push_back(to_string(static_cast<long double>(fitCalcs.at(swarm_->currentGeneration))));

			// Put our simulation params into the message vector
			for (map<string,double>::iterator i = simParams_.begin(); i != simParams_.end(); ++i){
				swarm_->swarmComm->univMessageSender.push_back(to_string(static_cast<long double>(i->second)));
			}

			// Tell the swarm master that we're finished
			//cout << id_ << " telling swarm we're finished" << endl;
			swarm_->swarmComm->sendToSwarm(id_, 0, SIMULATION_END, true, swarm_->swarmComm->univMessageSender);
			//cout << id_ << " done telling swarm" << endl;
			// Reset the message vector
			swarm_->swarmComm->univMessageSender.clear();
		}
		else {
			// If our return code is not 0, tell the master that the simulation failed
			swarm_->swarmComm->sendToSwarm(int(id_), 0, SIMULATION_FAIL, true, swarm_->swarmComm->univMessageSender);
		}

		// TODO: Breeding is still too slow. Need to speed things up here.
		// Wait for message from master telling us who to breed with
		bool doContinue = false;

		// swapTracker holds swapIDs and pIDs to keep track of who is breeding with who,
		// and which swaps are completed
		unordered_map<int,int> swapTracker;

		// Holds iterator ranges when finding items in the message holder
		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		while (!doContinue) {
			// Retrieve any messages
			int numCheckedMessages = 0;
			int numMessages = swarm_->swarmComm->recvMessage(-1, id_, -1, true, swarm_->swarmComm->univMessageReceiver, true);

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(INIT_BREEDING);
			if (smhRange.first != smhRange.second) {
				for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
					Timer tmr;

					// Convert the swapID to a unique negative number less than 1000
					int swapID = stoi(sm->second.message[0]);
					swapID += 1000;
					//cout << id_ << " init loop " << swapID << endl;

					// Store the particle with which to breed
					int pID = stoi(sm->second.message[1]);

					// Store our swap id so we know which swap we're working within
					swapTracker[swapID] = pID;

					//cout << id_ << " init breeding with " << pID << ". SwapID: " << swapID << endl;

					// Initiate breeding with that particle
					initBreedWithParticle(pID, swapID);

					//double t = tmr.elapsed();
					//cout << id_ << " INIT_BREEDING took " << t << " seconds" << endl;

					++numCheckedMessages;
				}
			}

			if (numCheckedMessages >= numMessages) {
				swarm_->swarmComm->univMessageReceiver.clear();
				continue;
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
			if (smhRange.first != smhRange.second) {
				//cout << id_ << "found final " << endl;
				for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
					//Timer tmr;

					int messageIndex = 0;
					for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
						//cout << id_ << " updating parameter " << p->first << " to " << *o << endl;
						p->second = stod(sm->second.message[messageIndex]);
						++messageIndex;
					}

					// Construct our filenames
					bnglFilename = to_string(static_cast<long long int>(id_)) + ".bngl";
					path = swarm_->options.jobOutputDir + "/" + to_string(static_cast<long long int>(swarm_->currentGeneration + 1));
					bnglFullPath = path + "/" + bnglFilename;

					// And generate our models
					if (swarm_->options.model->getHasGenerateNetwork()){
						// If we're using ODE solver, output .net and .bngl
						model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(static_cast<long long int>(id_)), false, false, true, false, false);
					}
					else {
						// If we're using network free simulation, output .bngl
						model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(static_cast<long long int>(id_)), false, false, false, false, false);
					}

					// Tell the master we have our new params and are ready for the next generation
					//cout << id_ << " telling master we're finished " << endl;
					swarm_->swarmComm->sendToSwarm(id_, 0, DONE_BREEDING, false, swarm_->swarmComm->univMessageSender);

					//double t = tmr.elapsed();
					//cout << "SEND_FINAL_PARAMS took " << t << " seconds" << endl;
					++numCheckedMessages;
				}
			}

			if (numCheckedMessages >= numMessages) {
				swarm_->swarmComm->univMessageReceiver.clear();
				continue;
			}

			// TODO: maybe we should search for all tags within this loop, rather then using equal_ranges
			for (Pheromones::swarmMsgHolderIt sm = swarm_->swarmComm->univMessageReceiver.begin(); sm != swarm_->swarmComm->univMessageReceiver.end(); ++sm) {
				if (sm->first > 1000) {
					//Timer tmr;

					// Store our swapID
					int swapID = sm->first;

					// Store the sender
					int reciprocateTo = sm->second.sender;

					vector<string> params;

					// Store the parameters
					for (auto p = sm->second.message.begin(); p != sm->second.message.end(); ++p) {
						params.push_back(*p);
					}

					int pID;
					// If we find our swapID in the swap tracker, it must mean that
					// we are receiving params from parent #2

					if (swapTracker.find(swapID) != swapTracker.end() && swapTracker.at(swapID) != id_) {
						//cout << id_ << " found " << swapTracker[swapID] << " at " << swapID << endl;
						//cout <<  id_ << " being given swapped parameters in swapID " << swapID << ". Receiving from the reciprocator: "<< reciprocateTo << endl;
						// Convert swapID to pID. pID is the "child" who receives the final parameter set
						pID = swapID - 1000;

						// Parse the received particles and integrate them with our own
						rcvBreedWithParticle(params, 0, swapID, pID);

					}
					// If we don't find our swapID in the swap tracker, it must mean that
					// we are receiving params from parent #1
					else {
						// If we are here and find our swapID in the tracker, it means we are breeding
						// with ourself. We must leave an entry in the swapTracker but change it from
						// our pID to ensure that we can receive from ourself
						if (swapTracker.find(swapID) != swapTracker.end()) {
							//cout << id_ << " found " << swapTracker[swapID] << " at " << swapID << ", changing to 0" <<  endl;
							swapTracker[swapID] = 0;
						}

						//cout << id_  << " being given swapped parameters in swapID " << swapID << ". Receiving from the initiator: "<< reciprocateTo << endl;
						// Convert swapID to pID. pID is the "child" who receives the final parameter set
						pID = swapID - 999;

						// Parse the received particles and integrate them with our own
						rcvBreedWithParticle(params, reciprocateTo, swapID, pID);
					}
					//cout << id_ << " is done breeding" << endl;
					//double t = tmr.elapsed();
					//cout << "RECEIVE_BREED took " << t << " seconds" << endl;
					++numCheckedMessages;
				}
			}

			if (numCheckedMessages >= numMessages) {
				swarm_->swarmComm->univMessageReceiver.clear();
				continue;
			}

			if (swarm_->swarmComm->univMessageReceiver.find(FIT_FINISHED) != swarm_->swarmComm->univMessageReceiver.end()) {
				//cout << id_ << " exiting " << endl;
				return;
			}

			if (swarm_->swarmComm->univMessageReceiver.find(NEXT_GENERATION) != swarm_->swarmComm->univMessageReceiver.end()) {
				doContinue = true;
			}

			swarm_->swarmComm->univMessageReceiver.clear();
		}
		// Reset the swapTracker so it is ready for the next generation
		swapTracker.clear();

		// Next generation
		swarm_->currentGeneration++;
	}
}

void Particle::runPSO() {

}

void Particle::calculateFit() {
	bool usingSD = false;
	bool usingMean = false;

	// Construct pointers to the relevant fitting function.
	// This lets us avoid a switch/case or excess conditionals
	// within the calculation loop
	if (swarm_->options.objFunc == 1) {
		objFuncPtr = &Particle::objFunc_sumOfSquares;
	}
	else if (swarm_->options.objFunc == 2) {
		objFuncPtr = &Particle::objFunc_chiSquare;
		usingSD = true;
	}
	else if (swarm_->options.objFunc == 3) {
		objFuncPtr = &Particle::objFunc_divByMeasured;
	}
	else if (swarm_->options.objFunc == 4) {
		objFuncPtr = &Particle::objFunc_divByMean;
		usingMean = true;
	}

	double colSum;
	double setSum;
	double totalSum = 0;
	double divisor = 0;

	// Loop through .exp files. Iterator points to string/dataset pair
	for (std::map<std::string,Data*>::iterator e = swarm_->options.expFiles.begin(); e != swarm_->options.expFiles.end(); ++e) {
		// Loop through .exp columns. Iterator points to column/map pair
		//cout << "exp loop " << e->first << endl;
		setSum = 0;
		//cout << "first: " << e->second->dataCurrent->begin()->first << endl;
		//cout << "second " <<  e->second->dataCurrent->begin()->second.begin()->first << endl;
		//cout << " third " << e->second->dataCurrent->begin()->second.begin()->second << endl;
		for (std::map<std::string,std::map<double,double> >::iterator exp_col = e->second->dataCurrent->begin(); exp_col != e->second->dataCurrent->end(); ++exp_col) {
			// Loop through timepoints of column.  Iterator points to a timepoint/value pair

			if (usingMean) {
				divisor = e->second->colAverages.at(exp_col->first);
			}

			//cout << "col loop " << exp_col->first << endl;
			colSum = 0;
			for (std::map<double,double>::iterator timepoint = exp_col->second.begin(); timepoint != exp_col->second.end(); ++timepoint) {
				//cout << "tp loop " << timepoint->first << endl;
				// TODO: Need to handle missing points

				//float exp = timepoint->second;
				//cout << "exp: " << exp << endl;
				//float sim = dataFiles_.at(e->first)->dataCurrent->at(exp_col->first).at(timepoint->first);
				//cout << "sim: " << sim << endl;
				//float SD = e->second->standardDeviations.at(exp_col->first).at(timepoint->first);
				//cout << "SD: " << SD << endl;

				if (usingSD) {
					//cout << "using sd" << endl;
					divisor = e->second->standardDeviations.at(exp_col->first).at(timepoint->first);
				}

				colSum += (this->*objFuncPtr)(dataFiles_.at(e->first)->dataCurrent->at(exp_col->first).at(timepoint->first), timepoint->second, divisor);
			}
			setSum += colSum;
		}
		totalSum += setSum;
	}
	fitCalcs[swarm_->currentGeneration] = pow(totalSum,0.5);
	//cout << id_ << " fitcalc: " << totalSum << endl;
}

// #1
double Particle::objFunc_chiSquare(double sim, double exp, double stdev) {
	return pow(((sim - exp)/stdev),2);
}

// #2
double Particle::objFunc_sumOfSquares(double sim, double exp, double dummyvar) {
	return pow((sim - exp),2);
}

// #3
double Particle::objFunc_divByMeasured(double sim, double exp, double dummyvar) {
	return pow(((sim - exp)/sim),2);
}

// #4
double Particle::objFunc_divByMean(double sim, double exp, double mean) {
	return pow(((sim - exp)/mean),2);
}

void Particle::initBreedWithParticle(int pID, int swapID) {
	//Timer tmr;
	//uniform_int_distribution<int> unif(1, 100);
	tr1::uniform_int<int> unif(1, 100);

	//for (auto p : simParams_) {
	for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
		// Swap
		if (unif(swarm_->randNumEngine) < (swarm_->options.swapRate * 100) ) {
			swarm_->swarmComm->univMessageSender.push_back(to_string(static_cast<long double>(p->second)));
		}
		// Don't swap
		else {
			swarm_->swarmComm->univMessageSender.push_back("");
		}
	}

	// Send message to parent 2 containing our param list
	swarm_->swarmComm->sendToSwarm(id_, pID, swapID, false, swarm_->swarmComm->univMessageSender);
	swarm_->swarmComm->univMessageSender.clear();

	//double t = tmr.elapsed();
	//cout << "Init took " << t << " seconds" << endl;
}

void Particle::rcvBreedWithParticle(vector<string>& params, int reciprocateTo, int swapID, int pID) {
	//Timer tmr;

	// If we're going to reciprocate, we need to store our params to send back to other parent
	// AND send a complete parameter set to the child.
	if (reciprocateTo) {
		vector<string> messageToChild;

		int pi = 0;
		//for (auto p : simParams_){
		for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {

			// If the parameter value is empty, we are not swapping that value. We give
			// parent #1 back an empty string, and give the child our value.
			if (params[pi] == "") {
				swarm_->swarmComm->univMessageSender.push_back("");
				messageToChild.push_back(to_string(static_cast<long double>(p->second)));

				//string vm = "DIDN'T SWAP " + to_string(static_cast<long double>(p->second));
			}

			// If the parameter is NOT empty, we are swapping. We give parent #1 our value
			// and give the child the (possible mutated) value from parent #2.
			else {
				swarm_->swarmComm->univMessageSender.push_back(to_string(static_cast<long double>(p->second)));

				double mutParam = stod(params[pi]);

				if (swarm_->options.hasMutate && swarm_->options.model->getfreeParams_().at(p->first)->isHasMutation()) {
					mutParam = mutateParam(swarm_->options.model->freeParams_.at(p->first), mutParam);
				}
				messageToChild.push_back(to_string(static_cast<long double>(mutParam)));
				//string vm = "SWAPPED " + to_string(mutParam);
			}
			++pi;
		}

		swarm_->swarmComm->sendToSwarm(id_, reciprocateTo, swapID, false, swarm_->swarmComm->univMessageSender);
		swarm_->swarmComm->univMessageSender.clear();
		swarm_->swarmComm->sendToSwarm(id_, pID, SEND_FINAL_PARAMS_TO_PARTICLE, false, messageToChild);
	}
	// If we aren't reciprocating, that means we're parent #1 and are receiving
	// parameter values from the reciprocator (parent #2). Create a message to
	// the child containing our values and the values given to us by the reciprocator
	else {
		int pi = 0;

		for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {

			// If the parameter is empty, we are keeping our parameter value. We
			// give that value to the child.
			if (params[pi] == "") {
				swarm_->swarmComm->univMessageSender.push_back(to_string(static_cast<long double>(p->second)));
			}

			// If the parameter is NOT empty, we are give the child the (possibly
			// mutated) value from parent #1.
			else {
				double finalParam = stod(params[pi]);
				if (swarm_->options.hasMutate && swarm_->options.model->getfreeParams_().at(p->first)->isHasMutation()) {
					finalParam = mutateParam(swarm_->options.model->freeParams_.at(p->first), finalParam);
				}

				swarm_->swarmComm->univMessageSender.push_back(to_string(static_cast<long double>(finalParam)));
			}
			++pi;
		}

		// Send the new parameter set to the child indicated by swapID
		//cout << id_ << " sending final param set to particle " << pID << endl;
		swarm_->swarmComm->sendToSwarm(id_, pID, SEND_FINAL_PARAMS_TO_PARTICLE, false, swarm_->swarmComm->univMessageSender);
		swarm_->swarmComm->univMessageSender.clear();
	}
	//double t = tmr.elapsed();
	//cout << "Rcv took " << t << " seconds" << endl;
}

double Particle::mutateParam(FreeParam* fp, double paramValue) {
	//Timer tmr;
	//uniform_real_distribution<double> unif(0,1);
	tr1::uniform_real<double> unif(0,1);

	// Generate a random number and see if it's less than our mutation rate.  If is, we mutate.
	if (unif(swarm_->randNumEngine) < fp->getMutationRate()) {
		// Store our mutation factor
		float maxChange = paramValue * fp->getMutationFactor();

		// Generate a new distribution between 0 and double maxChange
		//using param_t = uniform_real_distribution<>::param_type;
		//param_t p{0.0, maxChange * 2};
		//unif.param(p);

		tr1::uniform_real<double> unif(0.0, maxChange * 2);

		// Ger our new random number between 0 and maxChange, subtract maxChange.
		double change = unif(swarm_->randNumEngine) - maxChange;

		// Add/subtract the value from the parameter
		paramValue+= change;
		//cout << "factor: " << fp->getMutationFactor() << " max change: " << maxChange << " rand: " << rand << endl << endl;
	}
	//double t = tmr.elapsed();
	//cout << "Mutate took " << t << " seconds" << endl;
	return paramValue;
}
