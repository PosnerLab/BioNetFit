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
	model_ = 0;
	objFuncPtr = 0;
	currentGeneration_ = 1;
	island_ = 1;
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
		else if (genType == "loguniform_var") {

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

	if (swarm_->options.verbosity >= 3) {
		cout << "In doParticle(), waiting for message from master" << endl;
	}

	if (swarm_->resumingSavedSwarm) {
		cout << "Waiting for flight count from master" << endl;
		swarm_->swarmComm->recvMessage(0, id_, SEND_NUMFLIGHTS_TO_PARTICLE, true, swarm_->swarmComm->univMessageReceiver);
		currentGeneration_ = stoi(*(swarm_->swarmComm->univMessageReceiver.find(SEND_NUMFLIGHTS_TO_PARTICLE)->second.message.begin()));

		swarm_->swarmComm->univMessageReceiver.clear();
		cout << "Received flight count from master" << endl;

		cout << "Waiting for param set from master" << endl;
		swarm_->swarmComm->recvMessage(0, id_, SEND_FINAL_PARAMS_TO_PARTICLE, true, swarm_->swarmComm->univMessageReceiver);

		auto freeParam = swarm_->options.model->getFreeParams_().begin();
		for (auto paramVal = swarm_->swarmComm->univMessageReceiver.find(SEND_FINAL_PARAMS_TO_PARTICLE)->second.message.begin(); paramVal != swarm_->swarmComm->univMessageReceiver.find(SEND_FINAL_PARAMS_TO_PARTICLE)->second.message.end(); ++paramVal) {
			cout << "updating " << freeParam->first << " to " << *paramVal << endl;
			simParams_.insert(pair<string, double>(freeParam->first, stod(*paramVal)));
			++freeParam;
		}
		swarm_->swarmComm->univMessageReceiver.clear();

	}

	if (swarm_->options.fitType == "de") {
		// Need to figure out which island I'm on
		for (unsigned int island = 1; island <= swarm_->options.numIslands; ++island) {
			for (unsigned int particle = 1; particle <= (swarm_->options.swarmSize / swarm_->options.numIslands) * island; ++particle) {
				if (particle == id_) {
					island_ = island;
					goto theEnd;
				}
			}
		}
		theEnd:;
	}

	//cout << "Particle " << id_ << " waiting to begin" << endl;
	swarm_->swarmComm->recvMessage(0, id_, NEXT_GENERATION, true, swarm_->swarmComm->univMessageReceiver);
	swarm_->swarmComm->univMessageReceiver.clear();
	//cout << "Particle " << id_ << " starting" << endl;

	if (swarm_->options.verbosity >= 3) {
		cout << "In doParticle(), entering main run loop" << endl;
	}

	bool doContinue = true;
	bool doRunModel = true;
	while(doContinue) {
		if (swarm_->bootstrapCounter > 0 && currentGeneration_ == 1) {
			swarm_->swarmComm->recvMessage(0, id_, NEXT_GENERATION, true, swarm_->swarmComm->univMessageReceiver);
			swarm_->swarmComm->univMessageReceiver.clear();
		}

		if (doRunModel) {
			for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
				runModel(i);
			}

			if (swarm_->options.smoothing > 1) {
				smoothRuns();
			}

			finalizeSim();
		}
		else {
			doRunModel = true;
		}

		if (swarm_->options.fitType == "ga") {
			checkMessagesGenetic();
		}
		else if (swarm_->options.fitType == "pso") {
			checkMessagesPSO();
		}
		else if (swarm_->options.fitType == "de") {
			checkMessagesDE();
		}
		else if (swarm_->options.fitType == "sa") {
			doRunModel = checkMessagesDE();
		}

		if (doRunModel == true) {
			// Next generation
			++currentGeneration_;
		}
	}
}

void Particle::runModel(int iteration, bool localSearch) {

	// First get our path and filename variables set up for use in model generation, sim command, etc
	string bnglFilename = to_string(id_) + "_" + to_string(iteration) + ".bngl";
	string path = swarm_->options.jobOutputDir + to_string(currentGeneration_) + "/";

	string bnglFullPath = path + bnglFilename;
	cout << "path is " << path << endl;

	string suffix = to_string(id_) + "_" + to_string(iteration);

	string pipePath;

	// Only need to generate files if we're in the first generation. In subsequent generations
	// the model generation is handled by breeding parents
	if (currentGeneration_ == 1 && swarm_->bootstrapCounter == 0 && !localSearch) {
		if (swarm_->options.model->getHasGenerateNetwork()){
			string netFilename = "base.net";
			string netFullPath = swarm_->options.jobOutputDir + netFilename;

			if (iteration == 1) {
				swarm_->options.model->parseNet(netFullPath);
			}
			model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, true, false, false);
		}
		else {
			model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
		}
	}
	else if (localSearch) {
		cout << "regenerating models.." << endl;
		// And generate our models
		if (swarm_->options.model->getHasGenerateNetwork()){
			// If we're using ODE solver, output .net file
			bnglFilename = boost::regex_replace(bnglFilename, boost::regex("bngl$"), string("net"));
			model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, true);
		}
		else {
			// If we're using network free simulation, output .bngl
			model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
		}
	}

	// Generate .gdat pipes if we're using pipes
	string outputSuffix;
	if (swarm_->options.usePipes) {
		for (std::map<std::string,Model::action>::iterator i = model_->actions.begin(); i != model_->actions.end(); ++i) {
			if (i->second.scanParam.size() > 0) {
				outputSuffix = ".scan";
			}
			else {
				outputSuffix = ".gdat";
			}

			pipePath = path + i->first + "_" + to_string(id_) + "_" + to_string(iteration) + outputSuffix;
			cout << "pp is: " << pipePath << endl;
			createParticlePipe(pipePath.c_str());
		}
	}

	// Construct our simulation command
	string command = swarm_->options.bngCommand + " --outdir " + path + " " + bnglFullPath + " >> " + path + to_string(id_) + ".BNG_OUT 2>&1";
	if (swarm_->options.usePipes) {
		command += " &";
	}

	if (swarm_->options.verbosity >= 3) {
		cout << "Running model with command: " << command << endl;
	}

	// Run simulation command
	//int ret = system(command.c_str());
	int ret = runCommand(command);

	// Check for simulation command success
	if (ret == 0) { // TODO: Need to check for simulation status when using pipes. Going by return code doesn't work there because we're using the & operator
		// really not sure how we can do this easily
		string outputSuffix;
		// Save our simulation outputs to data objects
		for (std::map<std::string,Model::action>::iterator action = model_->actions.begin(); action != model_->actions.end(); ++action) {
			if (action->second.scanParam.size() > 0) {
				outputSuffix = ".scan";
			}
			else {
				outputSuffix = ".gdat";
			}
			string dataPath = path + action->first + "_" + to_string(id_) + "_" + to_string(iteration) + outputSuffix;
			dataFiles_[action->first].insert(pair<int, Data*>(iteration, new Data(dataPath, swarm_, false)));
		}
	}
	else {
		// If our return code is not 0, tell the master that the simulation failed
		//cout << "I failed" << endl;
		swarm_->swarmComm->sendToSwarm(int(id_), 0, SIMULATION_FAIL, true, swarm_->swarmComm->univMessageSender);
	}
}

void Particle::checkMessagesGenetic() {

	// Holds iterator ranges when finding items in the message holder
	pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

	while (1) {
		// Retrieve any messages
		int numCheckedMessages = 0;
		int numMessages = swarm_->swarmComm->recvMessage(-1, id_, -1, true, swarm_->swarmComm->univMessageReceiver, true);

		smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
		if (smhRange.first != smhRange.second) {
			//cout << id_ << " found final " << endl;
			for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
				//Timer tmr;

				int messageIndex = 0;
				for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
					//cout << id_ << " updating parameter " << p->first << " to " << sm->second.message[messageIndex] << endl;
					p->second = stod(sm->second.message[messageIndex]);
					//cout << "updated param " << p->first << ": " << p->second << endl;
					++messageIndex;
				}

				string path = swarm_->options.jobOutputDir + to_string(currentGeneration_ + 1) + "/";

				if (!checkIfFileExists(path)) {
					runCommand("mkdir " + path);
				}

				for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
					// Construct our filenames
					string bnglFilename = to_string(id_) + "_" + to_string(i) + ".bngl";
					string bnglFullPath = path + bnglFilename;
					string suffix = to_string(id_) + "_" + to_string(i);

					// And generate our models
					if (swarm_->options.model->getHasGenerateNetwork()){
						// If we're using ODE solver, output .net and .bngl
						model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, true, false, false);
					}
					else {
						// If we're using network free simulation, output .bngl
						model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
					}
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


		if (swarm_->swarmComm->univMessageReceiver.find(FIT_FINISHED) != swarm_->swarmComm->univMessageReceiver.end()) {
			swarm_->swarmComm->~Pheromones();
			exit(0);
		}

		if (swarm_->swarmComm->univMessageReceiver.find(NEW_BOOTSTRAP) != swarm_->swarmComm->univMessageReceiver.end()) {
			if (swarm_->options.verbosity >= 3) {
				cout << id_ << ": Starting a new bootstrapping run" << endl;
			}

			++swarm_->bootstrapCounter;
			currentGeneration_ = 0;
			simParams_.clear();
			generateParams();

			swarm_->swarmComm->univMessageReceiver.clear();
			return;
		}

		if (swarm_->swarmComm->univMessageReceiver.find(NEXT_GENERATION) != swarm_->swarmComm->univMessageReceiver.end()) {
			return;
		}

		swarm_->swarmComm->univMessageReceiver.clear();
	}
}

void Particle::checkMessagesPSO() {

	while(1) {

		if (swarm_->options.verbosity >= 3) {
			cout << "Checking for messages from master" << endl;
		}

		int numCheckedMessages = 0;
		int numMessages = swarm_->swarmComm->recvMessage(-1, id_, -1, true, swarm_->swarmComm->univMessageReceiver, true);


		if (swarm_->options.verbosity >= 3) {
			//cout << "Found " << numMessages << " messages" << endl;

			/*
			for (auto sm = swarm_->swarmComm->univMessageReceiver.begin(); sm != swarm_->swarmComm->univMessageReceiver.end(); ++sm) {
				cout << "tag: " << sm->second.tag << endl;
			}
			 */
		}

		// Holds iterator ranges when finding items in the message holder
		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		while (numCheckedMessages < numMessages) {
			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
			if (smhRange.first != smhRange.second) {
				//cout << "Receiving parameter list from master" << endl;
				for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {

					int messageIndex = 0;
					for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
						//cout << id_ << " updating parameter " << p->first << " to " << sm->second.message[messageIndex] << endl;
						p->second = stod(sm->second.message[messageIndex]);
						++messageIndex;
					}

					for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {

						// Construct our filenames
						string bnglFilename = to_string(id_) + "_" + to_string(i) + ".bngl";
						string path = swarm_->options.jobOutputDir + to_string(currentGeneration_ + 1) + "/";
						string bnglFullPath = path + bnglFilename;
						string suffix = to_string(id_) + "_" + to_string(i);

						if (!checkIfFileExists(path)) {
							runCommand("mkdir " + path);
						}

						// And generate our models
						if (swarm_->options.model->getHasGenerateNetwork()){
							// If we're using ODE solver, output .net and .bngl
							model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, true, false, false);
						}
						else {
							// If we're using network free simulation, output .bngl
							model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
						}
					}

					++numCheckedMessages;
				}
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(FIT_FINISHED);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Master told me to die" << endl;
				}

				exit(0);
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEW_BOOTSTRAP);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Starting a new bootstrapping run" << endl;
				}

				++swarm_->bootstrapCounter;
				currentGeneration_ = 0;

				swarm_->swarmComm->univMessageReceiver.clear();
				return;
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEXT_GENERATION);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Master told me to move to the next iteration" << endl;
				}

				swarm_->swarmComm->univMessageReceiver.clear();
				return;
			}
		}
		swarm_->swarmComm->univMessageReceiver.clear();
	}
}

bool Particle::checkMessagesDE() {

	while(1) {

		if (swarm_->options.verbosity >= 3) {
			cout << "Checking for messages from master" << endl;
		}

		int numCheckedMessages = 0;
		int numMessages = swarm_->swarmComm->recvMessage(-1, id_, -1, true, swarm_->swarmComm->univMessageReceiver, true);


		if (swarm_->options.verbosity >= 3) {
			//cout << "Found " << numMessages << " messages" << endl;

			/*
			for (auto sm = swarm_->swarmComm->univMessageReceiver.begin(); sm != swarm_->swarmComm->univMessageReceiver.end(); ++sm) {
				cout << "tag: " << sm->second.tag << endl;
			}
			 */
		}

		// Holds iterator ranges when finding items in the message holder
		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		while (numCheckedMessages < numMessages) {
			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
			if (smhRange.first != smhRange.second) {
				cout << "Receiving parameter list from master" << endl;
				for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {

					int messageIndex = 0;
					for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
						//cout << id_ << " updating parameter " << p->first << " to " << sm->second.message[messageIndex] << endl;
						p->second = stod(sm->second.message[messageIndex]);
						++messageIndex;
					}

					for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {

						// Construct our filenames
						string bnglFilename = to_string(id_) + "_" + to_string(i) + ".bngl";
						string path = swarm_->options.jobOutputDir + to_string(currentGeneration_ + 1) + "/";
						string bnglFullPath = path + bnglFilename;
						string suffix = to_string(id_) + "_" + to_string(i);

						if (!checkIfFileExists(path)) {
							runCommand("mkdir " + path);
						}

						// And generate our models
						if (swarm_->options.model->getHasGenerateNetwork()){
							// If we're using ODE solver, output .net and .bngl
							model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, true, false, false);
						}
						else {
							// If we're using network free simulation, output .bngl
							model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
						}
					}

					++numCheckedMessages;
				}
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(FIT_FINISHED);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Master told me to die" << endl;
				}

				exit(0);
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEW_BOOTSTRAP);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Starting a new bootstrapping run" << endl;
				}

				++swarm_->bootstrapCounter;
				currentGeneration_ = 0;

				swarm_->swarmComm->univMessageReceiver.clear();
				return true;
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(BEGIN_NELDER_MEAD);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << "Beginning Nelder-Mead local search" << endl;
				}
				cout << "deserializing simplex" << endl;
				std::stringstream iss;
				iss.str(smhRange.first->second.message[0]); // Simplex is serialized in the first indice of the message vector
				boost::archive::text_iarchive ia(iss);
				map<double, vector<double>> simplex;
				ia >> simplex;

				cout << "running search" << endl;
				runNelderMead(simplex);
				cout << "nelder mead finished. old calc: " << simplex.begin()->first << " new calc: " << fitCalcs[-1] << endl;
				vector<string> paramsStr;
				paramsStr.push_back(to_string(fitCalcs[-1]));
				for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
					paramsStr.push_back(to_string(p->second));
					cout << "new param: " << p->second << endl;
				}
				swarm_->swarmComm->sendToSwarm(id_, 0, SIMULATION_END, false, paramsStr);
				swarm_->swarmComm->univMessageReceiver.clear();

				return false;
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEXT_GENERATION);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Master told me to move to the next iteration" << endl;
				}

				swarm_->swarmComm->univMessageReceiver.clear();
				return true;
			}
		}
		swarm_->swarmComm->univMessageReceiver.clear();
	}
}

void Particle::calculateFit(bool local) {
	if (swarm_->options.verbosity >= 3) {
		cout << "Calculating fit" << endl;
	}

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
	for (auto e = swarm_->options.expFiles.begin(); e != swarm_->options.expFiles.end(); ++e) {
		// Loop through .exp columns. Iterator points to column/map pair
		//cout << "exp loop " << e->first << endl;
		setSum = 0;
		//cout << "first: " << e->second->dataCurrent->begin()->first << endl;
		//cout << "second " <<  e->second->dataCurrent->begin()->second.begin()->first << endl;
		//cout << " third " << e->second->dataCurrent->begin()->second.begin()->second << endl;
		for (std::map<std::string, std::map<double,double> >::iterator exp_col = e->second->dataCurrent->begin(); exp_col != e->second->dataCurrent->end(); ++exp_col) {
			// Loop through timepoints of column.  Iterator points to a timepoint/value pair

			if (usingMean) {
				//cout << "trying to set mean" << endl;
				divisor = e->second->colAverages.at(exp_col->first);
				//cout << "mean set" << endl;
			}

			//cout << "col loop " << exp_col->first << endl;
			colSum = 0;
			for (std::map<double,double>::iterator timepoint = exp_col->second.begin(); timepoint != exp_col->second.end(); ++timepoint) {
				if (std::isnan(timepoint->second)) {
					continue;
				}

				//double exp = timepoint->second;
				//cout << "exp: " << exp << endl;

				if (swarm_->options.smoothing == 1) {
					//cout << "sim: " << dataFiles_.at(e->first).at(swarm_->options.smoothing)->dataCurrent->at(exp_col->first).at(timepoint->first) << endl;
				}
				else {
					//cout << "sim: " << dataFiles_.at(e->first).at(swarm_->options.smoothing + 1)->dataCurrent->at(exp_col->first).at(timepoint->first) << endl;
				}

				if (usingSD) {
					divisor = e->second->standardDeviations.at(exp_col->first).at(timepoint->first);
					//cout << "divisor: " << divisor << endl;
				}

				double sim;
				// TODO: Introduce fudge tolerance to account for precision loss in simulation control column
				if (swarm_->options.smoothing == 1) {
					sim = dataFiles_.at(e->first).at(swarm_->options.smoothing)->dataCurrent->at(exp_col->first).at(timepoint->first);
				}
				else {
					sim = dataFiles_.at(e->first).at(swarm_->options.smoothing+1)->dataCurrent->at(exp_col->first).at(timepoint->first);
				}

				double sum = (this->*objFuncPtr)(sim, timepoint->second, divisor);

				if (swarm_->options.bootstrap) {
					//cout << "all sets: "<< swarm_->bootstrapMaps.size() << endl;
					//cout << "files: " << swarm_->bootstrapMaps[swarm_->bootstrapCounter].size() << endl;
					//cout << "cols: " << (swarm_->bootstrapMaps[swarm_->bootstrapCounter].begin())->first << endl;
					//cout << "tps: " << swarm_->bootstrapMaps[swarm_->bootstrapCounter][exp_col->first][timepoint->first] << endl;
					//cout << "multiplying " << colSum << " by " << swarm_->bootstrapMaps[swarm_->bootstrapCounter - 1][e->first][exp_col->first][timepoint->first] << endl;
					//cout << "m: " << colSum << " by " << swarm_->bootstrapMaps[swarm_->bootstrapCounter].at(e->first).at(exp_col->first).at(timepoint->first) << endl;
					sum *= swarm_->bootstrapMaps[swarm_->bootstrapCounter][e->first][exp_col->first][timepoint->first];
				}
				colSum += sum;
			}
			setSum += colSum;
		}
		totalSum += setSum;
	}

	// Erase our data sets
	dataFiles_.clear();

	// Store our fit calc
	if (!local) {
		fitCalcs[currentGeneration_] = pow(totalSum, 0.5);
	}
	else {
		fitCalcs[-1] = pow(totalSum, 0.5);
	}

	if (swarm_->options.verbosity >= 3) {
		cout << "Fit calculation: " << pow(totalSum, 0.5) << endl;
	}
}

// #1
double Particle::objFunc_sumOfSquares(double sim, double exp, double dummyvar) {
	return pow((abs(sim) - exp), 2);
}

// #2
double Particle::objFunc_chiSquare(double sim, double exp, double stdev) {
	return pow(((abs(sim) - exp) / stdev), 2);

	// TODO: Missing stdev?
}

// #3
double Particle::objFunc_divByMeasured(double sim, double exp, double dummyvar) {
	cout << "sim: " << sim << endl;
	cout << "exp: " << exp << endl;

	return pow(((abs(sim) - exp) / sim), 2);

	// TODO: DIV BY 0??
}

// #4
double Particle::objFunc_divByMean(double sim, double exp, double mean) {
	return pow(((abs(sim) - exp) / mean), 2);
}

void Particle::finalizeSim() {

	if (swarm_->options.verbosity >= 3) {
		cout << "Finalizing simulation" << endl;
	}

	// Calculate our fit
	calculateFit();

	// Put our fit calc into the message vector
	swarm_->swarmComm->univMessageSender.push_back(to_string(fitCalcs.at(currentGeneration_)));
	//cout << "stored fit calc of " << fitCalcs.at(currentGeneration_) << " as " << swarm_->swarmComm->univMessageSender[0] << endl;

	// Put our simulation params into the message vector
	for (map<string,double>::iterator i = simParams_.begin(); i != simParams_.end(); ++i){
		//cout << "stored param of " << i->second << endl;
		swarm_->swarmComm->univMessageSender.push_back(to_string(i->second));
	}

	// Tell the swarm master that we're finished
	if (swarm_->options.verbosity >= 3) {
		cout << "Telling master that my simulation is finished" << endl;
	}

	swarm_->swarmComm->sendToSwarm(id_, 0, SIMULATION_END, false, swarm_->swarmComm->univMessageSender);
	// Reset the message vector
	swarm_->swarmComm->univMessageSender.clear();
}

void Particle::smoothRuns() {

	if (swarm_->options.verbosity >= 3) {
		cout << "Smoothing simulation outputs" << endl;
	}

	map<string, map<double, double>> dataSet;
	// For each action/prefix/exp file
	for (auto action = dataFiles_.begin(); action != dataFiles_.end(); ++action) {
		// For each column
		// Insert a new iteration into the prefix set
		//cout << "action loop " << action->first << endl;
		for (auto col = action->second.at(1)->dataCurrent->begin(); col != action->second.at(1)->dataCurrent->end(); ++col) {
			// For each timepoint
			// This map holds time/param value pairs
			//cout << "col loop " << col->first << endl;
			map<double, double> timePairs;
			for (auto time = col->second.begin(); time != col->second.end(); ++time) {
				double sum = 0;
				int i = 0;
				// For each iteration
				for (unsigned int iteration = 1; iteration <= swarm_->options.smoothing; ++iteration) {
					//cout << "it loop " << iteration << endl;
					sum += dataFiles_.at(action->first).at(iteration)->dataCurrent->at(col->first).at(time->first);
					++i;
				}
				double average = sum / (double)i;
				//cout << "average: " << average << endl;
				pair<double, double> timePair;
				timePair = make_pair(time->first, average);
				timePairs.insert(timePair);
			}
			dataSet.insert(pair<string, map<double, double>>(col->first, timePairs));
		}
		action->second.insert(pair<int, Data*>(swarm_->options.smoothing + 1, new Data(dataSet)));
		dataSet.clear();
	}


	/*
	std::cout.precision(18);

	for (auto action = dataFiles_.begin(); action != dataFiles_.end(); ++action) {
		cout << action->first << endl;
		Data * data = action->second.at(swarm_->options.smoothing + 1);
		for (auto col = data->dataCurrent->begin(); col != data->dataCurrent->end(); ++col) {
			cout << col->first << endl;
			for (auto tp = col->second.begin(); tp != col->second.end(); ++tp) {
				cout << tp->second << endl;
			}
		}
	}
	 */
	//cout << id_ << " done smoothing" << endl;
}

void Particle::runNelderMead(map<double, vector<double>> simplex) {
	// The transformation coefficients
	float reflection = 1.0;
	float expansion = 2.0;
	float contraction = 0.5;
	float shrink = 0.5;
	unsigned int simulationCount = 0;

	while (simulationCount < 20) {
		// Get our important vertices
		auto sIt = simplex.end();
		advance(sIt, - 1); // Last element
		cout << "worst is: " << sIt->first << endl;
		auto worst = sIt;
		advance(sIt, - 1); // Second to last element
		cout << "second worst is: " << sIt->first << endl;
		auto good = sIt;
		auto best = simplex.begin();
		cout << "best is: " << best->first << endl;
		bool invalid = false;

		vector<vector<double>> centroidVectors;
		for (auto params = simplex.begin(); params != worst; ++params) {
			centroidVectors.push_back(swarm_->normalizeParams(params->second));
		}

		// Calculate the centroid
		vector<double> centroid = getCentroid(centroidVectors);

		// Reflect
		cout << "reflecting" << endl;
		vector<double> R; // The transformation vector
		for (unsigned int d = 0; d < centroid.size(); ++d) {
			double r = centroid[d] + reflection * (centroid[d] - worst->second[d]);
			if (r <= 0 || r > 1) {
				invalid = true;
			}
			cout << "creating new pt with eq " << centroid[d] << " + " << reflection << " * (" << centroid[d] << " - " << worst->second[d] << "): " << r << endl;
			R.push_back(r);
		}

		vector<double> deNormalizedReflection = swarm_->deNormalizeParams(R);

		auto tIt = deNormalizedReflection.begin();
		for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
			p->second = *tIt;
			++tIt;
		}

		double rCalc;
		if (!invalid) {
			for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
				runModel(i, true);
				++simulationCount;
			}
			if (swarm_->options.smoothing > 1) {
				smoothRuns();
			}
			calculateFit(true);
			rCalc = fitCalcs[-1];
		}
		else {
			invalid = false;
			rCalc = pow(10, 10);
		}
		cout << "reflection: " << rCalc << endl;

		// R Better than good, but worse than best
		if (rCalc > best->first && rCalc < good->first) {
			cout << "reflection was worse than best and better than good. erasing worst of " << worst->first << endl;
			simplex.erase(worst);
			cout << "sim size: " << simplex.size() << endl;
			simplex.insert(pair<double, vector<double>> (rCalc, R));
			cout << "inserting rCalc of " << rCalc << endl;

			continue;
		}
		// R Better than best
		else if (rCalc < best->first) {
			cout << "reflection was better than best. expanding" << endl;
			// Expand

			vector<double> E;
			for (unsigned int d = 0; d < centroid.size(); ++d) {
				double e = centroid[d] + expansion * (R[d] - centroid[d]);
				if (e <= 0 || e > 1) {
					invalid = true;
				}
				E.push_back(e);
				cout << "creating new pt with eq " << centroid[d] << " + " << expansion << " * (" << R[d] << " - " << centroid[d] << "): " << e << endl;
			}

			vector<double> deNormalizedExpansion = swarm_->deNormalizeParams(E);

			auto tIt = deNormalizedExpansion.begin();
			for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
				p->second = *tIt;
				++tIt;
			}

			double eCalc;
			if (!invalid) {
				for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
					runModel(i, true);
					++simulationCount;
				}
				if (swarm_->options.smoothing > 1) {
					smoothRuns();
				}
				calculateFit(true);
				eCalc = fitCalcs[-1];
				cout << "expansion: " << eCalc << endl;
			}
			else {
				eCalc = pow(10,10);
				invalid = false;
			}

			if (eCalc < rCalc) {
				cout << "expansion was better than reflection. looping." << endl;
				simplex.erase(worst);
				simplex.insert(pair<double, vector<double>> (eCalc, E));
				continue;
			}
			else {
				cout << "expansion was worse than reflection. looping." << endl;
				simplex.erase(worst);
				simplex.insert(pair<double, vector<double>> (rCalc, R));
				continue;
			}
		}
		// R worse than good
		else if (rCalc > good->first) {
			cout << "reflection was worse than good. contracting." << endl;
			// Contraction
			vector<double> C;
			for (unsigned int d = 0; d < centroid.size(); ++d) {
				double c = centroid[d] + contraction * (worst->second[d] - centroid[d]);
				if (c <= 0 || c > 1) {
					invalid = true;
				}
				C.push_back(c);
				cout << "creating new pt with eq " << centroid[d] << " + " << contraction << " * (" << worst->second[d] << " - " << centroid[d] << "): " << c << endl;
			}

			vector<double> deNormalizedContraction = swarm_->deNormalizeParams(C);

			auto tIt = deNormalizedContraction.begin();
			for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
				p->second = *tIt;
				++tIt;
			}

			double cCalc;
			if (!invalid) {
				for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
					runModel(i, true);
					++simulationCount;
				}
				if (swarm_->options.smoothing > 1) {
					smoothRuns();
				}
				calculateFit(true);
				cCalc = fitCalcs[-1];
				cout << "contraction: " << cCalc << endl;
			}
			else {
				invalid = false;
				cCalc = pow(10, 10);
			}

			if (cCalc < worst->first) {
				cout << "contraction was better than worst. looping." << endl;
				simplex.erase(worst);
				simplex.insert(pair<double, vector<double>> (cCalc, C));
				continue;
			}
			else {
				//cout << "contraction was worse than worst. shrinking" << endl;
				// Shrink
				for (auto sIt = ++simplex.begin(); sIt != simplex.end(); ++sIt) {
					for (unsigned int d = 0; d < centroid.size(); ++d) {
						double s = simplex.begin()->second[d] + (shrink * (sIt->second[d] - best->second[d]));
						sIt->second[d] = s;
					}
				}
			}
		}
	}

	fitCalcs[-1] = simplex.begin()->first;
	unsigned int d = 0;
	for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
		p->second = simplex.begin()->second[d++];
	}
}

vector<double> Particle::getCentroid(vector<vector<double>> centroidVectors) {
	vector<double> centroid;

	for (unsigned int d = 0; d < centroidVectors[0].size(); ++d) {
		double sum = 0;
		for (unsigned int set = 0; set < centroidVectors.size(); ++set) {
			sum += centroidVectors[set][d];
		}
		centroid.push_back(sum / centroidVectors[0].size());
	}

	return centroid;
}

/*
bool Particle::checkNelderMeadTerminationCriteria(vector<vector<double>> simplex, unsigned int numEvaluations, unsigned int numIterations) {
	if (swarm_->options.maxLocalEvaluations && numEvaluations >= swarm_->options.maxLocalEvaluations) {
		return true;
	}
	else if (swarm_->options.maxLocalIterations && numEvaluations >= swarm_->options.maxLocalIterations) {
		return true;
	}
	else {

	}

	return false;
}
 */
