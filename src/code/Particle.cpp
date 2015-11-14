/*
 * Particle.cpp
 *
 *  Created on: Jul 19, 2015
 *      Author: brandon
 */

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

			cout << "generating with " << min << ":" << max << endl;

			// Generate a random double between 0 and rand_max
			double myrand = static_cast <double> (rand()) / static_cast<double> (RAND_MAX);
			//cout << id_ << " rand is " << myrand << endl;

			//
			double exp = log10(min) + myrand * (log10(max) - log10(min) );
			//cout << id_ << " base is " << exp << endl;
			cout << id_ << " final " << pow(10, exp) << endl << endl;

			pair<string,double> paramPair = make_pair(paramName, pow(10, exp));
			setParam(paramPair);
		}
		//TODO: Add the rest of the init parm generation options
	}
}

void Particle::doParticle() {

	while(1) {

		// First get our path and filename variables set up for use in model generation, sim command, etc
		string bnglFilename = to_string(id_) + ".bngl";
		string path = swarm_->options_.outputDir + "/" + to_string(swarm_->currentGeneration_);
		string bnglFullPath = path + "/" + bnglFilename;

		string pipePath;

		// Only need to generate files if we're in the first generation. In subsequent generations
		// the model generation is handled by beeding parents
		if (swarm_->currentGeneration_ == 1) {
			if (swarm_->options_.model->getHasGenerateNetwork()){
				string netFilename = "base.net";
				string netFullPath = swarm_->options_.outputDir + "/" + netFilename;

				swarm_->options_.model->parseNet(netFullPath);
				model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(id_), false, false, true, false, false);
			}
			else {
				model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(id_), false, false, false, false, false);
			}
		}

		// Generate .gdat pipes if we're using pipes
		if (swarm_->options_.usePipes) {
			for (std::unordered_map<std::string,Model::action>::iterator i = model_->actions.begin(); i != model_->actions.end(); ++i) {
				pipePath = path + "/" + i->first + "_" + to_string(id_) + ".gdat";
				createParticlePipe(pipePath.c_str());
			}
		}

		// Construct our simulation command
		string command = "/home/brandon/projects/GenFit/Simulators/BNG2.pl --outdir " + path + " " + bnglFullPath + ">> " + path + "/" + to_string(id_) + ".BNG_OUT 2>&1";
		if (swarm_->options_.usePipes) {
			command += " &";
		}

		cout << "Running model with command: " << command << endl;
		// Run simulation command
		int ret = system(command.c_str());

		// Check for simulation command success
		if (ret == 0) { // TODO: Need to check for simulation status when using pipes. Going by return code doesn't work there because we're using the & operator
			string dataPath;

			// Save our simulation outputs to data objects
			for (std::unordered_map<std::string,Model::action>::iterator i = model_->actions.begin(); i != model_->actions.end(); ++i) {
				dataPath = path + "/" + i->first + "_" + to_string(id_) + ".gdat";
				//dataFiles_.insert(make_pair(i->first, new Data(dataPath, swarm_, false)));
				dataFiles_[i->first] = new Data(dataPath, swarm_, false);
				//cout << "loaded data for " << dataPath << endl;
			}

			// Calculate our fit
			calculateFit();

			// Put our fit calc into the message vector
			swarm_->swarmComm_->univMessageSender.push_back(to_string(fitCalcs.at(swarm_->currentGeneration_)));
			//cout << id_ << " adding calc of: " << fitCalcs.at(swarm_->currentGeneration_) << endl;
			// Put our simulation params into the message vector
			for (map<string,double>::iterator i = simParams_.begin(); i != simParams_.end(); ++i){
				swarm_->swarmComm_->univMessageSender.push_back(to_string(i->second));
				cout << "adding " << i->first << " of " << i->second << endl;
			}

			// Tell the swarm master that we're finished
			swarm_->swarmComm_->sendToSwarm(int(id_), 0, SIMULATION_END, false, swarm_->swarmComm_->univMessageSender);

			// Reset the message vector
			swarm_->swarmComm_->univMessageSender.clear();
		}
		else {
			swarm_->swarmComm_->sendToSwarm(int(id_), 0, SIMULATION_FAIL, false, swarm_->swarmComm_->univMessageSender);
		}

		// Wait for message from master telling us who to breed with
		bool doContinue = false;

		unordered_map<int,int> swapTracker;
		std::vector<std::vector<std::string>> messageReceiver;
		//cout << id_ << " waiting for instructions from master" << endl;

		while (!doContinue) {
			// TODO: Find out why breeding takes so long and fix it!!
			swarm_->swarmComm_->recvMessage(-1, id_, -1, true, messageReceiver, true);

			// First loop through individual messages
			for (std::vector<std::vector<std::string>>::iterator m = messageReceiver.begin(); m != messageReceiver.end(); ++m) {

				std::vector<std::string>::iterator o = (*m).begin();
				int tag = stoi(*o);

				//cout << id_ << " found tag: " << tag << endl;
				if (tag == INIT_BREEDING) {
					//Timer tmr;

					// Jump to the swapID
					o+=4;

					// Convert the swapID to a unique negative number less than 1000
					int swapID = stoi(*o);
					swapID = (swapID * -1) - 1000;

					// Jump to particle with which to breed
					o+=1;

					// Store our swap id so we know which swap we're working within
					swapTracker[swapID] = id_;

					cout << id_ << " init breeding with " << *o << ". SwapID: " << swapID << endl;

					// Initiate breeding with that particle
					initBreedWithParticle(stoi(*o), swapID);

					//double t = tmr.elapsed();
					//cout << "INIT_BREEDING took " << t << " seconds" << endl;
				}

				// Receive new params from parents
				else if (tag == SEND_FINAL_PARAMS_TO_PARTICLE) {
					//Timer tmr;
					// Jump to the first parameter value
					o+=4;

					for (std::map<std::string,double>::iterator p = simParams_.begin(); p != simParams_.end(); ++p) {
						//cout << id_ << " updating parameter " << p->first << " to " << *o << endl;
						p->second = stod(*o);
						o++;
					}

					milliseconds ms = duration_cast< milliseconds >(
							system_clock::now().time_since_epoch()
					);
					cout << ms.count() << " " << id_ << " generating model" << endl;
					bnglFilename = to_string(id_) + ".bngl";
					path = swarm_->options_.outputDir + "/" + to_string(swarm_->currentGeneration_ + 1);
					bnglFullPath = path + "/" + bnglFilename;

					if (swarm_->options_.model->getHasGenerateNetwork()){
						// If we're using ODE solver, output .net and .bngl
						model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(id_), false, false, true, false, false);
					}
					else {
						// If we're using network free simulation, output .bngl
						model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(id_), false, false, false, false, false);
					}

					// Tell the master we have our new params and are ready for the next generation
					cout << id_ << " telling master I'm finished " << endl;
					swarm_->swarmComm_->sendToSwarm(id_, 0, DONE_BREEDING, false, swarm_->swarmComm_->univMessageSender);

					//double t = tmr.elapsed();
					//cout << "SEND_FINAL_PARAMS took " << t << " seconds" << endl;
				}
				// We need to breed. Tag is equal to swap id.
				else if (tag < -1000) {
					//Timer tmr;

					int swapID = tag;

					// Jump to the sender
					o+=2;
					int reciprocateTo = stoi(*o);

					// Jump to the parameter list
					o+=2;
					vector<string> params;

					// Store the parameters
					while (o != (*m).end()) {
						params.push_back(*o);
						o++;
					}

					int pID;
					// If we find our swapID in the swap tracker, it must mean that
					// we are receiving params from parent #2
					if (swapTracker.find(swapID) != swapTracker.end() && swapTracker.at(swapID) != id_) {
						cout <<  id_ << " being given swapped parameters in swapID " << swapID << ". Receiving from the reciprocator: "<< reciprocateTo << endl;
						// Convert swapID to pID. pID is the "child" who receives the final parameter set
						pID = (swapID * - 1) - 1000;

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
							swapTracker[swapID] = 0;
						}

						cout << id_  << " being given swapped parameters in swapID " << swapID << ". Receiving from the initiator: "<< reciprocateTo << endl;
						// Convert swapID to pID. pID is the "child" who receives the final parameter set
						pID = (swapID * - 1) - 999;

						// Parse the received particles and integrate them with our own
						rcvBreedWithParticle(params, reciprocateTo, swapID, pID);
					}
					//cout << id_ << " is done breeding" << endl;
					//double t = tmr.elapsed();
					//cout << "RECEIVE_BREED took " << t << " seconds" << endl;
					cout << id_ << " received from " << reciprocateTo << endl;
				}
				else if (tag == NEXT_GENERATION) {
					cout << id_ << " done breeding" << endl;
					// TODO: We're generating models at signal for next generation, but things would go faster
					// if we did it on the fly as particles breed
					doContinue = true;
				}
				else if (tag == FIT_FINISHED) {
					return;
				}
			}
			messageReceiver.clear();
		}
		swapTracker.clear();
		swarm_->currentGeneration_++;
	}
}

void Particle::calculateFit() {
	bool usingSD = false;
	bool usingMean = false;

	if (swarm_->options_.objFunc == 1) {
		objFuncPtr = &Particle::objFunc_sumOfSquares;
	}
	else if (swarm_->options_.objFunc == 2) {
		objFuncPtr = &Particle::objFunc_chiSquare;
		usingSD = true;
	}
	else if (swarm_->options_.objFunc == 3) {
		objFuncPtr = &Particle::objFunc_divByMeasured;
	}
	else if (swarm_->options_.objFunc == 4) {
		objFuncPtr = &Particle::objFunc_divByMean;
		usingMean = true;
	}

	double colSum;
	double setSum;
	double totalSum = 0;
	double divisor = 0;

	// Loop through .exp files. Iterator points to string/dataset pair
	for (std::unordered_map<std::string,Data*>::iterator e = swarm_->options_.expFiles.begin(); e != swarm_->options_.expFiles.end(); ++e) {
		// Loop through .exp columns. Iterator points to column/map pair
		//cout << "exp loop " << e->first << endl;
		setSum = 0;
		//cout << "first: " << e->second->dataCurrent_->begin()->first << endl;
		//cout << "second " <<  e->second->dataCurrent_->begin()->second.begin()->first << endl;
		//cout << " third " << e->second->dataCurrent_->begin()->second.begin()->second << endl;
		for (std::map<std::string,std::map<double,double> >::iterator exp_col = e->second->dataCurrent_->begin(); exp_col != e->second->dataCurrent_->end(); ++exp_col) {
			// Loop through timepoints of column.  Iterator points to a timepoint/value pair

			if (usingMean) {
				divisor = e->second->colAverages_.at(exp_col->first);
			}

			//cout << "col loop " << exp_col->first << endl;
			colSum = 0;
			for (std::map<double,double>::iterator timepoint = exp_col->second.begin(); timepoint != exp_col->second.end(); ++timepoint) {
				//cout << "tp loop " << timepoint->first << endl;
				// TODO: Need to handle missing points

				//float exp = timepoint->second;
				//cout << "exp: " << exp << endl;
				//float sim = dataFiles_.at(e->first)->dataCurrent_->at(exp_col->first).at(timepoint->first);
				//cout << "sim: " << sim << endl;
				//float SD = e->second->standardDeviations_.at(exp_col->first).at(timepoint->first);
				//cout << "SD: " << SD << endl;

				if (usingSD) {
					//cout << "using sd" << endl;
					divisor = e->second->standardDeviations_.at(exp_col->first).at(timepoint->first);
				}

				colSum += (this->*objFuncPtr)(dataFiles_.at(e->first)->dataCurrent_->at(exp_col->first).at(timepoint->first), timepoint->second, divisor);
			}
			setSum += colSum;
		}
		totalSum += setSum;
	}
	fitCalcs[swarm_->currentGeneration_] = totalSum;
	cout << id_ << " fitcalc: " << totalSum << endl;
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
	uniform_int_distribution<int> unif(1, 100);

	for (auto p : simParams_) {
		// Swap
		if (unif(swarm_->randNumEngine) < (swarm_->options_.swapRate * 100) ) {
			double parameter = p.second;

			/*
			if (swarm_->options_.hasMutate && swarm_->options_.model->getFreeParams().at(p.first)->isHasMutation()) {
				parameter = mutateParam(swarm_->options_.model->freeParams_.at(p.first), parameter);
			}
			 */
			swarm_->swarmComm_->univMessageSender.push_back(to_string(parameter));
		}
		// Don't swap
		else {
			swarm_->swarmComm_->univMessageSender.push_back("");
		}
	}

	// Send message to parent 2 containing our param list
	swarm_->swarmComm_->sendToSwarm(id_, pID, swapID, false, swarm_->swarmComm_->univMessageSender);
	swarm_->swarmComm_->univMessageSender.clear();
}

void Particle::rcvBreedWithParticle(vector<string>& params, int reciprocate, int swapID, int pID) {
	// If we're going to reciprocate, we need to store our params to send back to other parent
	if (reciprocate) {
		for (auto p : simParams_) {
			swarm_->swarmComm_->univMessageSender.push_back(to_string(p.second));
		}

		//cout << "telling " << reciprocate << " that we're reciprocating in swapID: " << swapID << endl;

		swarm_->swarmComm_->sendToSwarm(id_, reciprocate, swapID, false, swarm_->swarmComm_->univMessageSender);
		swarm_->swarmComm_->univMessageSender.clear();
	}

	// Get new params from first parent and store them
	int pi = 0;

	for (auto p : simParams_){
		if (params[pi] != "") {
			//cout << "adding " << params[pi] << endl;
			swarm_->swarmComm_->univMessageSender.push_back(params[pi]);
		}
		else {
			//cout << "adding " << p.second << endl;
			swarm_->swarmComm_->univMessageSender.push_back(to_string(p.second));
		}
		++pi;
	}

	// Send the new parameter set to the particle indicated by swapID
	cout << id_ << " sending final param set to particle " << pID << endl;
	swarm_->swarmComm_->sendToSwarm(id_, pID, SEND_FINAL_PARAMS_TO_PARTICLE, false, swarm_->swarmComm_->univMessageSender);
	swarm_->swarmComm_->univMessageSender.clear();

}
