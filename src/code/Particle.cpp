/*
 * Particle.cpp
 *
 *  Created on: Jul 19, 2015
 *      Author: brandon
 */

#include "Particle.hh"

using namespace std;

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

	for (map<string,string>::iterator i = model_->freeParams_.begin(); i != model_->freeParams_.end(); ++i) {

		vector<string> values;

		split(i->second, values);

		if (values[0].empty()){
			string errMsg = "Problem parsing initial parameter generator for param: " + i->first + ".";
			outputError (errMsg);
		}
		string paramName;
		paramName = i->first;

		string genType = values[0];

		if (genType == "random_var"){
			if (values[1].empty() || values[2].empty()){
				string errMsg = "Problem parsing initial parameter generator for param: " + paramName + ".";
				outputError (errMsg);
			}

			double min = atof(values[1].c_str());
			double max = atof(values[2].c_str());

			float myrand = min + static_cast<double> (rand()) /(static_cast<double> (RAND_MAX/(max-min)));

			pair<string,double> paramPair = make_pair(paramName, myrand);
			setParam(paramPair);
		}
		else if (genType == "loguniform_var") { //TODO: Test distribution of this generator
			if (values[1].empty() || values[2].empty()){
				string errMsg = "Problem parsing initial parameter generator for param: " + paramName + ".";
				outputError (errMsg);
			}
			double min = atof(values[1].c_str());
			double max = atof(values[2].c_str());

			double myrand = static_cast <double> (rand()) / static_cast<double> (RAND_MAX);
			double base = ( log(min) / log(10) ) + myrand * ( ( log(max) / log(10) ) - ( log(min) / log(10) ) );
			//TODO: The above generator doesn't work. It outputs numbers outside the specified range.

			pair<string,double> paramPair = make_pair(paramName, pow(base, 10));
			setParam(paramPair);
		}
		//TODO: Add the rest of the init parm generation options
	}
}

void Particle::doParticle() {

	//while(1) {

	// First get our path and filename variables set up for use in model generation, sim command, etc
	string bnglFilename = to_string(id_) + ".bngl";
	string path = swarm_->options_.outputDir + "/" + to_string(swarm_->currentGeneration_);
	string bnglFullPath = path + "/" + bnglFilename;

	string pipePath;

	if (swarm_->options_.model->getHasGenerateNetwork()){
		string netFilename = "base.net";
		string netFullPath = swarm_->options_.outputDir + "/" + netFilename;

		swarm_->options_.model->parseNet(netFullPath);
		model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(id_), false, false, true, false, false);
	}
	else {
		model_->outputModelWithParams(simParams_, path, bnglFilename, to_string(id_), false, false, false, false, false);
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

	// Run simulation command
	int ret = system(command.c_str());

	// Check for simulation command success
	if (ret == 0) { // TODO: Need to check for simulation status when using pipes. Going by return code doesn't work there because we're using the & operator
		string dataPath;

		// Save our simulation outputs to data objects
		for (std::unordered_map<std::string,Model::action>::iterator i = model_->actions.begin(); i != model_->actions.end(); ++i) {
			dataPath = path + "/" + i->first + "_" + to_string(id_) + ".gdat";
			dataFiles_.insert(make_pair(i->first, new Data(dataPath, swarm_, false)));
		}

		// Calculate our fit
		calculateFit();

		// Put our fit calc into the message vector
		swarm_->swarmComm_->univMessageSender.push_back(to_string(fitCalcs.at(swarm_->currentGeneration_)));

		// Put our simulation params into the message vector
		for (map<string,double>::iterator i = simParams_.begin(); i != simParams_.end(); ++i){
			swarm_->swarmComm_->univMessageSender.push_back(to_string(i->second));
			//cout << "adding " << i->first << " of " << i->second << endl;
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
	cout << id_ << " waiting for instructions from master" << endl;

	while (!doContinue) {

		swarm_->swarmComm_->recvMessage(-1, id_, -1, true, messageReceiver, true);

		// First loop through individual messages
		for (std::vector<std::vector<std::string>>::iterator m = messageReceiver.begin(); m != messageReceiver.end(); ++m) {

			std::vector<std::string>::iterator o = (*m).begin();
			int tag = stoi(*o);

			//cout << id_ << " found tag: " << tag << endl;
			if (tag == INIT_BREEDING) {
				// Jump particle with which to breed
				o+=4;

				// Generate and store our swap id so we know which swap we're working within
				uniform_int_distribution<int> unif(-32767, -1001);
				int swapID = unif(swarm_->randNumEngine);
				swapTracker[swapID] = id_;

				cout << "Particle " << id_ << " init breeding with " << *o << ". SwapID: " << swapID << endl;

				// Initiate breeding with that particle
				initBreedWithParticle(stoi(*o), swapID);
			}
			// We need to breed. Tag is equal to swap id.
			else if (tag < -1000) {
				int swapID = tag;

				// Jump to the sender
				o+=2;
				int reciprocateTo = stoi(*o);

				// Jump to the parameter list
				o+=2;
				vector<string> params;

				while (o != (*m).end()) {
					params.push_back(*o);
					o++;
				}

				// If we find our swapID in the swap tracker, it must mean that
				// we are reciprocating
				if (swapTracker.find(swapID) != swapTracker.end()) {
					cout << "Particle " << id_ << " being given swapped parameters in swapID " << swapID << ". Receiving from the reciprocator: "<< reciprocateTo << endl;
					rcvBreedWithParticle(params, 0, swapID);
				}
				else {
					cout << "Particle " << id_  << " being given swapped parameters in swapID " << swapID << ". Receiving from the initiator: "<< reciprocateTo << endl;
					rcvBreedWithParticle(params, reciprocateTo, swapID);
				}
			}
			else if (tag == NEXT_GENERATION) {
				cout << "done breeding" << endl;
				doContinue = true;
			}
		}
		messageReceiver.clear();
	}
	swapTracker.clear();
	//}
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
	//cout << "fitcalc: " << totalSum << endl;
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

	bool hasSwaps = false;
	for (auto p : simParams_) {
		// Swap
		if (unif(swarm_->randNumEngine) < (swarm_->options_.swapRate * 100) ) {
			swarm_->swarmComm_->univMessageSender.push_back(to_string(p.second));
			hasSwaps = true;
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

void Particle::rcvBreedWithParticle(vector<string>& params, int reciprocate, int swapID = 0) {
	// If we're going to reciprocate, we need to store our params to send back to other parent
	if (reciprocate) {
		for (auto p : simParams_) {
			swarm_->swarmComm_->univMessageSender.push_back(to_string(p.second));
		}

		cout << "telling " << reciprocate << " that we're reciprocating in swapID: " << swapID << endl;

		swarm_->swarmComm_->sendToSwarm(id_, reciprocate, swapID, false, swarm_->swarmComm_->univMessageSender);
		swarm_->swarmComm_->univMessageSender.clear();
	}
	else {
		cout << "telling master we're done breeding" << endl;
		swarm_->swarmComm_->sendToSwarm(id_, 0, DONE_BREEDING, false, swarm_->swarmComm_->univMessageSender);
	}

	// Get new params from first parent and store them
	int pi = 0;
	for (auto p : simParams_){
		if (params[pi] != "") {
			p.second = stod(params[pi]);
		}
		++pi;
	}
}
