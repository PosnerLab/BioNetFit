/*
 * Config.cpp
 *
 *  Created on: Jul 13, 2015
 *      Author: brandon
 */

#include <unordered_map>

#include "Config.hh"

using namespace std;

Config::Config(string configFile) {
	configPath_ = convertToAbsPath(configFile);
}

Swarm * Config::createSwarmFromConfig () {
	Swarm *s = new Swarm();

	string line;
	ifstream confFile(configPath_);
	boost::smatch matches;
	string name, value;
	unordered_multimap<string,string> pairs;

	s->setConfigPath(configPath_);

	//TODO: Why are we using atoi with c_str when we could go straight to stoi?

	// TODO: We really should be assigning iterators inside the conditional
	// statements to avoid the extra map find() each time we find a match
	if (confFile.is_open()) {
		while (getline(confFile, line)) {
			name = "";
			value = "";

			if (line.length() > 0){
				if (line.at(0) == '#')
					continue;

				boost::regex keyVal("^\\s*([a-zA-Z_]+)\\s*=\\s*(.*)\\s*$");
				if (boost::regex_search(line,matches,keyVal)){
					name = matches[1];
					value = matches[2];
					pairs.insert(make_pair(name,value));
				}
				else {
					boost::regex keyVal("^\\s*([a-zA-Z_]+)\\s+(.*)\\s*$");
					if (boost::regex_search(line,matches,keyVal)){
						name = matches[1];
						value = matches[2];
						pairs.insert(make_pair(name,value));
					}
				}
			}
		}
		confFile.close();
	}
	else {
		string errMsg = "Error: Couldn't open config file " + configPath_ + " for parsing.";
		outputError(errMsg);
	}

	// Add our model file to the swarm
	if(pairs.find("model") != pairs.end()) {

		s->setModel(pairs.find("model")->second);
	}

	// Update the swarm type
	if(pairs.find("fit_type") != pairs.end()) {

		s->setfitType(pairs.find("fit_type")->second);
	}

	// Update the swarm size
	if(pairs.find("swarm_size") != pairs.end()) {

		int swarmSize = atoi(pairs.find("swarm_size")->second.c_str());

		// If the swarm size isn't even, make it even by adding a particle
		if (swarmSize % 2 != 0 && s->options.fitType == "genetic") {
			++swarmSize;
		}
		s->options.swarmSize = swarmSize;
	}

	// Update the sim path
	if(pairs.find("bng_command") != pairs.end()) {

		s->options.bngCommand = pairs.find("bng_command")->second;
	}

	// Update the synchronicity
	if(pairs.find("synchronicity") != pairs.end()) {

		s->options.synchronicity = (stoi(pairs.find("synchronicity")->second));
	}

	// Update the maximum number of generations
	if(pairs.find("max_generations") != pairs.end()) {

		s->options.maxGenerations = (stoi(pairs.find("max_generations")->second));
	}

	if(pairs.find("output_every") != pairs.end()) {

		s->options.outputEvery = stoi(pairs.find("output_every")->second.c_str());
	}

	// Tell the swarm if we're using pipes
	if(pairs.find("use_pipes") != pairs.end()) {

		s->options.usePipes = (stoi(pairs.find("use_pipes")->second) == 1) ? true : false;
	}

	// Tell the swarm if we should delete old files
	if(pairs.find("delete_old_files") != pairs.end()) {

		s->options.deleteOldFiles = (stoi(pairs.find("delete_old_files")->second) == 1) ? true : false;
	}

	// Tell the swarm if we're using a cluster
	if(pairs.find("use_cluster") != pairs.end()) {
		if (stoi(pairs.find("use_cluster")->second)) {

			// Set cluster platform if it was specified
			if(pairs.find("cluster_software") != pairs.end()) {
				s->options.clusterSoftware = pairs.find("cluster_software")->second;
			}

			s->options.useCluster = (stoi(pairs.find("use_cluster")->second) == 1) ? true : false;
			s->getClusterInformation();
			//cout << "setting useCluster to " << s->options.useCluster << endl;
			// TODO: Set parallel count accordingly
		}
	}

	// TODO: Need to ensure path gets made
	// Tell the swarm if we should save cluster output
	if(pairs.find("save_cluster_output") != pairs.end()) {
		s->options.saveClusterOutput = (stoi(pairs.find("save_cluster_output")->second) == 1) ? true : false;
	}

	// Update swap rate
	if(pairs.find("swap_rate") != pairs.end()) {
		s->options.swapRate = stof(pairs.find("swap_rate")->second);
	}

	// Update extra weight
	if(pairs.find("extra_weight") != pairs.end()) {
		s->options.extraWeight = stoi(pairs.find("extra_weight")->second);
	}

	// Whether or not to force difference parents
	if(pairs.find("force_different_parents") != pairs.end()) {
		s->options.forceDifferentParents = (stoi(pairs.find("force_different_parents")->second) == 1) ? true : false;
	}

	// How many retries when breeding
	if(pairs.find("max_breeding_retries") != pairs.end()) {
		s->options.maxRetryDifferentParents = stoi(pairs.find("max_breeding_retries")->second);
	}

	// Set fit value that will cause fit to end
	if(pairs.find("smoothing") != pairs.end()) {
		s->options.smoothing = stod(pairs.find("smoothing")->second);
	}

	// Set output directory
	if(pairs.find("output_dir") != pairs.end()) {
		s->options.outputDir = convertToAbsPath(pairs.find("output_dir")->second);
	}

	// Set job name
	if(pairs.find("job_name") != pairs.end()) {
		s->options.jobName = pairs.find("job_name")->second;
	}

	// Set the job output directory
	s->setJobOutputDir(s->options.outputDir + "/" + s->options.jobName + "/");
	// TODO: Move directory creation to a swarm function.  Also need to check cmd success.

	// Set verbosity
	if(pairs.find("verbosity") != pairs.end()) {
		s->options.verbosity = stoi(pairs.find("verbosity")->second);
	}

	// Set fit value that will cause fit to end
	if(pairs.find("min_fit") != pairs.end()) {
		s->options.minFit = stod(pairs.find("min_fit")->second);
	}

	// Set the maximum fit value to consider in breeding
	// TODO: Implement this
	if(pairs.find("max_fit") != pairs.end()) {
		s->options.maxFit = stod(pairs.find("max_fit")->second);
	}

	// Set maximum number of simulations in an asynchronous genetic fit
	if(pairs.find("max_num_simulations") != pairs.end()) {
		s->options.maxNumSimulations = stol(pairs.find("max_num_simulations")->second);
	}

	// Set maximum fitting time
	// TODO: Implement this in synchronous genetic fitting
	if(pairs.find("max_fit_time") != pairs.end()) {
		vector<string> timeElements;
		split(pairs.find("max_fit_time")->second, timeElements, ":");
		long timeInSeconds = (stol(timeElements[0]) * 3600) + (stol(timeElements[1]) * 60) + stol(timeElements[2]);
		s->options.maxFitTime = timeInSeconds;
	}

	// Update the maximum number of parallel threads (non-cluster only)
	if(pairs.find("parallel_count") != pairs.end()) {
		// TODO: Make sure PC isn't higher than swarm size
		if (s->options.useCluster) {
			s->options.parallelCount = s->options.swarmSize;
		}
		else {
			s->options.parallelCount = stoi(pairs.find("parallel_count")->second);
		}

	}
	// Whether or not to divide by value at t=0
	if(pairs.find("divide_by_init") != pairs.end()) {
		s->options.divideByInit = (stoi(pairs.find("divide_by_init")->second.c_str()) == 1) ? true : false;
	}

	// Whether or not to log transform simulation output
	if(pairs.find("log_transform_sim_data") != pairs.end()) {
		s->options.logTransformSimData = stoi(pairs.find("log_transform_sim_data")->second);
	}

	// Whether or not to standardize simulation output
	if(pairs.find("standardize_sim_data") != pairs.end()) {
		s->options.standardizeSimData = (stoi(pairs.find("standardize_sim_data")->second) == 1) ? true : false;
	}

	// Whether or not to standardize exp data
	if(pairs.find("standardize_exp_data") != pairs.end()) {
		s->options.standardizeExpData = (stoi(pairs.find("standardize_exp_data")->second) == 1) ? true : false;
	}

	// Update fit calculation method
	if(pairs.find("objfunc") != pairs.end()) {
		s->options.objFunc = stoi(pairs.find("objfunc")->second);
	}


	//////////////////////////////////////////////////////////////////
	// PSO STUFF
	//////////////////////////////////////////////////////////////////

	if (pairs.find("inertia") != pairs.end()) {
		s->options.inertia = stof(pairs.find("inertia")->second);
	}

	if (pairs.find("cognitive") != pairs.end()) {
		s->options.cognitive = stof(pairs.find("cognitive")->second);
	}

	if (pairs.find("social") != pairs.end()) {
		s->options.social = stof(pairs.find("social")->second);
	}

	if (pairs.find("inertia") != pairs.end()) {
		s->options.inertia = stof(pairs.find("inertia")->second);
	}

	if (pairs.find("inertia") != pairs.end()) {
		s->options.inertia = stof(pairs.find("inertia")->second);
	}

	if (pairs.find("inertia") != pairs.end()) {
		s->options.inertia = stof(pairs.find("inertia")->second);
	}

	if (pairs.find("nmin") != pairs.end()) {
		s->options.nmin = stoi(pairs.find("nmin")->second);
	}

	if (pairs.find("nmax") != pairs.end()) {
		s->options.nmax = stoi(pairs.find("nmax")->second);
	}

	if (pairs.find("inertia_init") != pairs.end()) {
		s->options.inertiaInit = stof(pairs.find("inertia_init")->second);
	}

	if (pairs.find("inertia_final") != pairs.end()) {
		s->options.inertiaFinal = stof(pairs.find("inertia_final")->second);
	}

	if (pairs.find("abs_tolerance") != pairs.end()) {
		s->options.absTolerance = stof(pairs.find("abs_tolerance")->second);
	}

	if (pairs.find("rel_tolerance") != pairs.end()) {
		s->options.relTolerance = stof(pairs.find("rel_tolerance")->second);
	}

	if (pairs.find("topology") != pairs.end()) {
		s->options.topology = pairs.find("topology")->second;
	}

	if (pairs.find("pso_type") != pairs.end()) {
		s->options.psoType = pairs.find("pso_type")->second;
	}

	if(pairs.find("enhanced_stop") != pairs.end()) {

		s->options.enhancedStop = (stoi(pairs.find("enhanced_stop")->second) == 1) ? true : false;
	}

	if(pairs.find("enhanced_inertia") != pairs.end()) {

		s->options.enhancedInertia = (stoi(pairs.find("enhanced_inertia")->second) == 1) ? true : false;
		cout << "enhanced inertia set to " << s->options.enhancedInertia << endl;
	}

	// Add any init param generation options
	//for (auto pair: pairs) {
	for (unordered_multimap<string, string>::iterator pair = pairs.begin(); pair != pairs.end(); ++pair) {
		// TODO: Implement a map.equalrange to speed this part up
		if (pair->first == "random_var" || pair->first == "lognormrandom_var" || pair->first == "loguniform_var") {
			vector<string> paramComponents;
			split(pair->second,paramComponents);
			//string genString = pair->first + " " + paramComponents[1] + " " + paramComponents[2];

			// Make sure we have three components to work with
			if (paramComponents.size() == 3) {
				// Make sure first parameter name exists as a free parameter
				if (s->options.model->freeParams_.count(paramComponents[0]) > 0) {
					// Make sure 2nd and 3rd components are numeric
					if (isFloat(paramComponents[1]) && isFloat(paramComponents[2])) {
						s->options.model->freeParams_.at(paramComponents[0])->setGenerationMethod(pair->first);
						s->options.model->freeParams_.at(paramComponents[0])->setParameterName(paramComponents[0]);
						s->options.model->freeParams_.at(paramComponents[0])->setGenMin(stof(paramComponents[1]));
						s->options.model->freeParams_.at(paramComponents[0])->setGenMax(stof(paramComponents[2]));

						//cout << "setting " << paramComponents[0] << " to " << pair->first << ":" << paramComponents[1] << ":" << paramComponents[2] << endl;
					}
					else {
						outputError("Error: Problem parsing your free parameter generation option in your .conf file. The min and/or max values were non-numeric.");
					}
				}
				else {
					cout << "Warning: We found a parameter '" << paramComponents[0] << "' in your .conf file, but don't see a matching free parameter specification in your model file. We will ignore this parameter." << endl;
				}
			}
			else {
				outputError("Error: Problem parsing your free parameter generation option in your .conf file. Each parameter generation option requires three components: the parameter name, minimum, and maximum.");
			}
		}
	}

	//for (auto i : s->options.model->freeParams_) {
	for (map<string, FreeParam*>::iterator i = s->options.model->freeParams_.begin(); i != s->options.model->freeParams_.end(); ++i) {
		if (!i->second) {
			string errMsg = "Error: We found a free parameter '" + i->first + "' specified in your model file but can't find a matching parameter generator in your .conf file.";
			outputError(errMsg);
		}
	}

	// Add any .exp files to the swarm
	std::pair <std::unordered_multimap<string,string>::iterator, std::unordered_multimap<string,string>::iterator> it;
	it = pairs.equal_range("exp");

	for (unordered_multimap<string,string>::iterator exp = it.first; exp != it.second; ++exp) {
		s->addExp(exp->second);
	}

	// Add any .mutation rates to the swarm
	it = pairs.equal_range("mutate");

	for (unordered_multimap<string,string>::iterator exp = it.first; exp != it.second; ++exp) {
		s->addMutate(exp->second);
		s->options.hasMutate = true;

		if (boost::regex_search(exp->second, boost::regex("^default\\s"))) {
			//cout << "found default mutation rate. breaking" << endl;
			break;
		}
	}

	// Link all model actions with their corresponding data sets (.exp files)
	// Remove any .exp files that aren't pointed to in the model file
	// And Remove any Model::actions without corresponding .exp files specified in the .conf file
	vector<string> prefixedActions;
	vector<Model::action> toDeleteActs;

	for (map<string, Model::action>::iterator i = s->options.model->actions.begin(); i != s->options.model->actions.end();)
		//for (vector<Model::action>::iterator i = s->options.model->actions.begin(); i != s->options.model->actions.end();)
	{
		prefixedActions.push_back(i->first);

		if (s->options.expFiles.count(i->first) == 1) {
			if (s->options.verbosity >=3 ) {
				cout << "Linking action " << i->first << " with exp file: " << s->options.expFiles[i->first]->getPath() << endl;
			}
			i->second.dataSet = s->options.expFiles[i->first];
			++i;
		}
		else { // Have a prefix but no .exp file
			cout << "Warning: The model file specifies an action with the prefix '" << i->first << "' but there isn't a matching .exp file specified in your .conf file. We will ignore this action command." << endl;
			//i = s->options.model->actions.erase(i);
			s->options.model->actions.erase(++i);
		}
	}

	//vector<Data*> toDeleteExp;
	//for (auto &i : s->options.expFiles){ // Have .exp but no prefix
	for (map<string, Data*>::iterator i = s->options.expFiles.begin(); i != s->options.expFiles.end(); ++i) {
		if(std::find(prefixedActions.begin(), prefixedActions.end(), i->first) == prefixedActions.end() ) {
			cout << "Warning: The .conf file specifies an .exp file '" << i->first << "' but there isn't a matching action command in your model file specified with the prefix=> argument." << endl;
			//toDeleteExp.push_back(s->options.expFiles[i->first]);
			//i = s->options.expFiles.erase(i);
			s->options.expFiles.erase(++i);
		}
	}

	/*
	// TODO: Make sure we have either min fit, max time, or max sims when doing an asynchronous fit
	// This whole checker needs to be much more complex
	if (!s->options.synchronicity) {
		if (s->options.minFit == -1 && s->options.maxFitTime == MAX_LONG && s->options.maxNumSimulations == MAX_LONG && !s->options.enhancedStop) {
			outputError("Error: You are running an asynchronous fit, but failed to specify either a minimum fitting value, maximum fit time, or maximum number of simulations in the .conf file.");
		}
	}
	 */

	/*
	for (auto i : s->options.expFiles){
		cout << "Exp: " << i.second->getPath() << endl;
	}

	for (auto i : s->options.modelFile->actions_){
		cout << "Prefix: " << i.prefix << endl;
	}

	for (auto i : s->options.modelFile->actions_){
		cout << "DataSets: " << i.dataSet->getPath() << endl;
	}
	 */

	return s;
}

string Config::getLocation () {
	return configPath_;
}

void Config::checkConsistency() {
	/*
	 *	All Requires:
	 *	1. fit_type
	 *	2. job_name
	 *	3. model
	 *	4. exp
	 *	5. output_dir
	 *	6. bng_command
	 *	7. swarm_size
	 *
	 */

}
