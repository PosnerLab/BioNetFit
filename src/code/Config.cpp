/*
 * Config.cpp
 *
 *  Created on: Jul 13, 2015
 *      Author: brandon
 */

#include "Config.hh"

using namespace std;

Config::Config(string configFile) {
	configPath_ = convertToAbsPath(configFile);
}

Swarm * Config::createSwarmFromConfig (bool isMaster) {
	Swarm *s = new Swarm(isMaster);

	string line;
	ifstream confFile(configPath_);
	smatch matches;
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

				regex keyVal("^\\s*([a-zA-Z_]+)\\s*=\\s*(.*)\\s*$");
				if (regex_search(line.cbegin(),line.cend(),matches,keyVal)){
					name = matches[1];
					value = matches[2];
					pairs.insert(make_pair(name,value));
				}
				else {
					regex keyVal("^\\s*([a-zA-Z_]+)\\s+(.*)\\s*$");
					if (regex_search(line.cbegin(),line.cend(),matches,keyVal)){
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
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setModel(pairs.find("model")->second);
	}

	// Update the swarm size
	if(pairs.find("swarm_size") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		int swarmSize = atoi(pairs.find("swarm_size")->second.c_str());

		// If the swarm size isn't even, make it even by adding a particle
		if (swarmSize % 2 != 0) {
			swarmSize++;
		}
		s->setSwarmSize(swarmSize);
	}

	// Update the swarm type
	if(pairs.find("swarm_type") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setSwarmType(pairs.find("swarm_type")->second);
	}

	// Update the synchronicity
	if(pairs.find("synchronicity") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setSwarmSynchronicity(atoi(pairs.find("synchronicity")->second.c_str()));
	}

	// Update the maximum number of generations
	if(pairs.find("generations") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setSwarmGenerations(atoi(pairs.find("generations")->second.c_str()));
	}

	// Tell the swarm if we're using pipes
	if(pairs.find("use_pipes") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setUsePipes((atoi(pairs.find("use_pipes")->second.c_str()) == 1) ? true : false );
	}

	// Tell the swarm if we should delete old files
	if(pairs.find("delete_old_files") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setDeleteOldFiles((atoi(pairs.find("delete_old_files")->second.c_str()) == 1) ? true : false );
	}

	// Tell the swarm if we're using a cluster
	if(pairs.find("use_cluster") != pairs.end()) {
		if (stoi(pairs.find("use_cluster")->second)) {
			if(pairs.find("cluster_software") != pairs.end()) {
				s->options.clusterSoftware = pairs.find("cluster_software")->second;
			}

			s->setUseCluster((atoi(pairs.find("use_cluster")->second.c_str()) == 1) ? true : false );
			//s->setParallelCount(s->getSwarmSize());
			s->getClusterInformation();
		}
	}

	// Update swap rate
	if(pairs.find("swap_rate") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setSwapRate(atof(pairs.find("swap_rate")->second.c_str()));
	}

	// Update extra weight
	if(pairs.find("extra_weight") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setExtraWeight(atoi(pairs.find("extra_weight")->second.c_str()));
	}

	// Whether or not to force difference parents
	if(pairs.find("force_different_parents") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setUseCluster((atoi(pairs.find("force_different_parents")->second.c_str()) == 1) ? true : false );
		s->setParallelCount(s->getSwarmSize());
	}

	// How many retries when breeding
	if(pairs.find("max_breeding_retries") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setMaxRetryDifferentParents(atoi(pairs.find("max_breeding_retries")->second.c_str()));
	}

	// Set output directory
	if(pairs.find("output_dir") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setOutputDir(convertToAbsPath(pairs.find("output_dir")->second));
	}

	// Set job name
	if(pairs.find("job_name") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setJobName(pairs.find("job_name")->second);
	}

	// Set the job output directory
	s->options.jobOutputDir = s->options.outputDir + "/" + s->options.jobName + "/";

	// Set verbosity
	if(pairs.find("verbosity") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setVerbosity(atoi(pairs.find("verbosity")->second.c_str()));
	}

	// Set fit value that will cause fit to end
	if(pairs.find("min_fit") != pairs.end()) {
			//cout << "Adding model file: " << pairs.find("model")->second << endl;
			s->options.minFit = stod(pairs.find("min_fit")->second);
	}

	// Set the maximum fit value to consider in breeding
	// TODO: Implement this
	if(pairs.find("max_fit") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->options.maxFit = (stod(pairs.find("max_fit")->second));
	}

	// Set maximum number of simulations in an asynchronous genetic fit
	if(pairs.find("maxNumSimulations") != pairs.end()) {
		s->options.maxNumSimulations = stol(pairs.find("max_num_simulations")->second);
	}

	// Set maximum fitting time
	// TODO: Implement this in synchronous genetic fitting
	if(pairs.find("max_fit_time") != pairs.end()) {
		vector<string> timeElements;
		split(pairs.find("max_fit_time")->second, timeElements, ":");
		long timeInSeconds = (stol(timeElements[0]) * 3600) + (stol(timeElements[1]) * 60) + stol(timeElements[2]);
		s->options.maxFit = timeInSeconds;
	}

	// Update the maximimum number of parallel threads (non-cluster only)
	if(pairs.find("parallel_count") != pairs.end() && !s->getUseCluster()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setParallelCount(atoi(pairs.find("parallel_count")->second.c_str()));
	}
	// Whether or not to divide by value at t=0
	if(pairs.find("divide_by_init") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setDivideByInit((atoi(pairs.find("divide_by_init")->second.c_str()) == 1) ? true : false );
	}

	// Whether or not to log transform simulation output
	if(pairs.find("log_transform_sim_data") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setLogTransformSimData((atoi(pairs.find("log_transform_sim_data")->second.c_str())));
	}

	// Whether or not to standardize simulation output
	if(pairs.find("standardize_sim_data") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setStandardizeSimData((atoi(pairs.find("standardize_sim_data")->second.c_str()) == 1) ? true : false );
	}

	// Whether or not to standardize exp data
	if(pairs.find("standardize_exp_data") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setStandardizeExpData((atoi(pairs.find("standardize_exp_data")->second.c_str()) == 1) ? true : false );
	}

	// Update fit calculation method
	if(pairs.find("objfunc") != pairs.end()) {
		//cout << "Adding model file: " << pairs.find("model")->second << endl;
		s->setSosCalc((atoi(pairs.find("objfunc")->second.c_str())));
	}

	// Add any init param generation options
	for (auto pair : pairs) {
		// TODO: Implement a map.equalrange to speed this part up
		if (pair.first == "random_var" || pair.first == "lognormrandom_var" || pair.first == "loguniform_var") {
			vector<string> paramComponents;
			split(pair.second,paramComponents);
			//string genString = pair.first + " " + paramComponents[1] + " " + paramComponents[2];

			// Make sure we have three components to work with
			if (paramComponents.size() == 3) {
				// Make sure first parameter name exists as a free parameter
				if (s->options.model->freeParams_.count(paramComponents[0]) > 0) {
					// Make sure 2nd and 3rd components are numeric
					if (isFloat(paramComponents[1]) && isFloat(paramComponents[2])) {
						s->options.model->freeParams_.at(paramComponents[0])->setGenerationMethod(pair.first);
						s->options.model->freeParams_.at(paramComponents[0])->setParameterName(paramComponents[0]);
						s->options.model->freeParams_.at(paramComponents[0])->setGenMin(stof(paramComponents[1]));
						s->options.model->freeParams_.at(paramComponents[0])->setGenMax(stof(paramComponents[2]));

						//cout << "setting " << paramComponents[0] << " to " << pair.first << ":" << paramComponents[1] << ":" << paramComponents[2] << endl;
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

	for (auto i : s->options.model->freeParams_) {
		if (!i.second) {
			string errMsg = "Error: We found a free parameter '" + i.first + "' specified in your model file but can't find a matching parameter generator in your .conf file.";
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

		if (regex_search(exp->second, regex("^default\\s"))) {
			//cout << "found default mutation rate. breaking" << endl;
			break;
		}
	}

	/*
	for (auto pair : pairs) {
		if (pair.first == "exp") {
			//cout << "Adding exp file: " << pair.second << endl;
			s->addExp(pair.second);
		}
	}
	 */

	// Link all model actions with their corresponding data sets (.exp files)
	// Remove any .exp files that aren't pointed to in the model file
	// And Remove any Model::actions without corresponding .exp files specified in the .conf file
	vector<string> prefixedActions;
	vector<Model::action> toDeleteActs;

	for (std::map<string,Model::action>::iterator i = s->options.model->actions.begin(); i != s->options.model->actions.end();)
		//for (vector<Model::action>::iterator i = s->options.model->actions.begin(); i != s->options.model->actions.end();)
	{
		//if (!i->prefix.empty()) {
		//	prefixedActions.push_back(i->prefix);
		//}
		prefixedActions.push_back(i->first);

		if (s->options.expFiles.count(i->first) == 1) {
			if (s->getVerbosity() >=3 ) {
				cout << "Linking action " << i->first << " with exp file: " << s->options.expFiles[i->first]->getPath() << endl;
			}
			i->second.dataSet = s->options.expFiles[i->first];
			++i;
		}
		else { // Have a prefix but no .exp file
			cout << "Warning: The model file specifies an action with the prefix '" << i->first << "' but there isn't a matching .exp file specified in your .conf file. We will ignore this action command." << endl;
			i = s->options.model->actions.erase(i);
		}
	}

	vector<Data*> toDeleteExp;
	for (auto &i : s->options.expFiles){ // Have .exp but no prefix
		if(std::find(prefixedActions.begin(), prefixedActions.end(), i.first) == prefixedActions.end() ) {
			cout << "Warning: The .conf file specifies an .exp file '" << i.first << "' but there isn't a matching action command in your model file specified with the prefix=> argument." << endl;
			toDeleteExp.push_back(s->options.expFiles[i.first]);
		}
	}

	// TODO: We could probably delete in the loop
	// Delete the unmatched Exp objects
	for (auto i : toDeleteExp) {
		//cout << "Deleting object" << endl;
		s->options.expFiles.erase(getFilename(i->getPath()));
		delete i;
	}

	// TODO: Make sure we have either min fit, max time, or max sims when doing an asynchronous fit
	if (s->options.synchronicity) {
		if (s->options.minFit == -1 && s->options.maxFitTime == LONG_MAX && s->options.maxNumSimulations == LONG_MAX) {
			outputError("Error: You are running an asynchronous fit, but failed to specify either a minimum fitting value, maximum fit time, or maximum number of simulations in the .conf file.");
		}
	}

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

void Config::createConfigFromSwarm () {

}

string Config::getLocation () {
	return configPath_;
}
