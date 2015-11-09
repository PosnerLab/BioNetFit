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

	while(1) {

		// First get our path and filename variables set up for use in model generation, sim command, etc
		string bnglFilename = to_string(id_) + ".bngl";
		string netFilename = "base.net";

		string path = swarm_->options_.outputDir + "/" + to_string(swarm_->currentGeneration_);

		string bnglFullPath = path + "/" + bnglFilename;
		string netFullPath = swarm_->options_.outputDir + "/" + netFilename;

		string pipePath;

		if (swarm_->options_.model->getHasGenerateNetwork()){
			//cout << "has generate network" << endl;
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
		string command = "/home/brandon/projects/GenFit/Simulators/BNG2.pl --outdir " + path + " " + bnglFullPath + ">> bngout 2>&1";
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
		while (1) {

		}
	}
}

void Particle::calculateFit() {
	bool usingSD = false;
	bool usingMean = false;

	if (swarm_->options_.sosCalc == 2) {
		objFuncPtr = &Particle::objFunc_chiSquare;
		usingSD = true;
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

double Particle::objFunc_chiSquare(double sim, double exp, double stdev) {
	return pow(((exp - sim)/stdev),2);
}
