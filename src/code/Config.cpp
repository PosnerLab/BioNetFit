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
	swarm_ = 0;

	configPath_ = convertToAbsPath(configFile);
}

// TODO: Make "max_flocks" synonymous with max_generations

Swarm * Config::createSwarmFromConfig () {
	swarm_ = new Swarm();

	string line;
	ifstream confFile(configPath_);
	boost::smatch matches;
	string name, value;
	unordered_multimap<string,string> pairs;

	swarm_->setConfigPath(configPath_);

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
		swarm_->outputError("Error: Couldn't open config file " + configPath_ + " for parsing.");
	}

	// Add our model file to the swarm
	if(pairs.find("model") != pairs.end()) {
		//cout << "Processing model" << endl;
		swarm_->setModel(pairs.find("model")->second);
	}

	// Update the swarm type
	if(pairs.find("fit_type") != pairs.end()) {
		//cout << "Processing fit type" << endl;
		swarm_->setfitType(pairs.find("fit_type")->second);
	}

	// Update the swarm size
	if(pairs.find("swarm_size") != pairs.end()) {
		//cout << "Processing swarm size" << endl;
		int swarmSize = stoi(pairs.find("swarm_size")->second);

		// If the swarm size isn't even, make it even by adding a particle
		//if (swarmSize % 2 != 0 && swarm_->options.fitType == "genetic") {
		//	++swarmSize;
		//}
		swarm_->options.swarmSize = swarmSize;
	}

	// Update the sim path
	if(pairs.find("bng_command") != pairs.end()) {
		//cout << "Processing bng command" << endl;
		swarm_->options.bngCommand = pairs.find("bng_command")->second;
		if (!checkIfFileExists(swarm_->options.bngCommand)) {
			swarm_->outputError("Error: The specified bng command: " + swarm_->options.bngCommand + " does not exist. Quitting.");
		}
	}

	// Update the synchronicity
	if(pairs.find("synchronicity") != pairs.end()) {
		//cout << "Processing synchronicity" << endl;
		swarm_->options.synchronicity = (stoi(pairs.find("synchronicity")->second));
	}

	// Update the maximum number of generations
	if(pairs.find("max_generations") != pairs.end()) {
		//cout << "Processing max generations" << endl;
		swarm_->options.maxGenerations = (stoi(pairs.find("max_generations")->second));
	}

	if(pairs.find("output_every") != pairs.end()) {
		//cout << "Processing outputevery" << endl;
		swarm_->options.outputEvery = stoi(pairs.find("output_every")->second.c_str());
	}

	// Tell the swarm if we're using pipes
	if(pairs.find("use_pipes") != pairs.end()) {
		//cout << "Processing usepipes" << endl;
		swarm_->options.usePipes = (stoi(pairs.find("use_pipes")->second) == 1) ? true : false;
	}

	// Tell the swarm if we should delete old files
	if(pairs.find("delete_old_files") != pairs.end()) {
		//cout << "Processing delete old files" << endl;
		swarm_->options.deleteOldFiles = (stoi(pairs.find("delete_old_files")->second) == 1) ? true : false;
	}

	// Tell the swarm if we're using a cluster
	if(pairs.find("use_cluster") != pairs.end()) {
		//cout << "Processing use cluster" << endl;
		if (stoi(pairs.find("use_cluster")->second)) {

			// Set cluster platform if it was specified
			if(pairs.find("cluster_software") != pairs.end()) {
				//cout << "Processing cluster software" << endl;
				swarm_->options.clusterSoftware = pairs.find("cluster_software")->second;
			}

			swarm_->options.useCluster = (stoi(pairs.find("use_cluster")->second) == 1) ? true : false;
			swarm_->getClusterInformation();
			//cout << "setting useCluster to " << swarm_->options.useCluster << endl;
			// TODO: Set parallel count accordingly

			if (pairs.find("email_when_finished") != pairs.end() && stoi(pairs.find("email_when_finished")->second)) {
				swarm_->options.emailWhenFinished = true;
			}

			if (pairs.find("email_address") != pairs.end()) {
				swarm_->options.emailAddress = pairs.find("email_address")->second;
			}
		}
	}

	// TODO: Need to ensure path gets made
	// Tell the swarm if we should save cluster output
	if(pairs.find("save_cluster_output") != pairs.end()) {
		cout << "Processing save cluster output" << endl;
		swarm_->options.saveClusterOutput = (stoi(pairs.find("save_cluster_output")->second) == 1) ? true : false;
	}

	// Update seed
	if(pairs.find("seed") != pairs.end()) {
		cout << "Processing seed" << endl;
		swarm_->options.seed = stoi(pairs.find("seed")->second);
	}

	// Update swap rate
	if(pairs.find("swap_rate") != pairs.end()) {
		cout << "Processing swap rate" << endl;
		swarm_->options.swapRate = stof(pairs.find("swap_rate")->second);
	}

	// Update number of parents to keep unchanged in breeding
	if(pairs.find("keep_parents") != pairs.end()) {
		cout << "Processing keep parents" << endl;
		swarm_->options.keepParents = stoi(pairs.find("keep_parents")->second);
	}

	// Update extra weight
	if(pairs.find("extra_weight") != pairs.end()) {
		cout << "Processing extra weight" << endl;
		swarm_->options.extraWeight = stoi(pairs.find("extra_weight")->second);
	}

	// Whether or not to force difference parents
	if(pairs.find("force_different_parents") != pairs.end()) {
		cout << "Processing force different parents" << endl;
		swarm_->options.forceDifferentParents = (stoi(pairs.find("force_different_parents")->second) == 1) ? true : false;
	}

	// How many retries when breeding
	if(pairs.find("max_breeding_retries") != pairs.end()) {
		cout << "Processing max breeding retries" << endl;
		swarm_->options.maxRetryDifferentParents = stoi(pairs.find("max_breeding_retries")->second);
	}

	// Set fit value that will cause fit to end
	if(pairs.find("smoothing") != pairs.end()) {
		cout << "Processing smoothing" << endl;
		swarm_->options.smoothing = stoi(pairs.find("smoothing")->second);

		if (swarm_->options.smoothing == 0) {
			swarm_->options.smoothing = 1;
		}
	}

	// Set output directory
	if(pairs.find("output_dir") != pairs.end()) {
		cout << "Processing output dir" << endl;
		swarm_->options.outputDir = convertToAbsPath(pairs.find("output_dir")->second);
	}

	// Set job name
	if(pairs.find("job_name") != pairs.end()) {
		cout << "Processing job name" << endl;
		swarm_->options.jobName = pairs.find("job_name")->second;
	}

	// Set the job output directory
	swarm_->setJobOutputDir(swarm_->options.outputDir + "/" + swarm_->options.jobName + "/");

	// Set verbosity
	if(pairs.find("verbosity") != pairs.end()) {
		cout << "Processing verbosity" << endl;
		swarm_->options.verbosity = stoi(pairs.find("verbosity")->second);
	}

	// Set fit value that will cause fit to end
	if(pairs.find("min_fit") != pairs.end()) {
		cout << "Processing minfit" << endl;
		swarm_->options.minFit = stod(pairs.find("min_fit")->second);
	}

	// Set the maximum fit value to consider in breeding
	if(pairs.find("max_fit") != pairs.end()) {
		cout << "Processing maxfit" << endl;
		swarm_->options.maxFit = stod(pairs.find("max_fit")->second);
	}

	// Set maximum number of simulations in an asynchronous genetic fit
	if(pairs.find("max_num_simulations") != pairs.end()) {
		cout << "Processing max num sims" << endl;
		swarm_->options.maxNumSimulations = stol(pairs.find("max_num_simulations")->second);
	}

	// Set maximum fitting time
	// TODO: Implement this in synchronous genetic fitting
	if(pairs.find("max_fit_time") != pairs.end()) {
		cout << "Processing max fit time" << endl;
		vector<string> timeElements;
		split(pairs.find("max_fit_time")->second, timeElements, ":");
		long timeInSeconds = (stol(timeElements[0]) * 3600) + (stol(timeElements[1]) * 60) + stol(timeElements[2]);
		swarm_->options.maxFitTime = timeInSeconds;
	}

	// Update the maximum number of parallel threads (non-cluster only)
	if(pairs.find("parallel_count") != pairs.end()) {
		cout << "Processing parallel count" << endl;
		// TODO: Make sure PC isn't higher than swarm size
		if (swarm_->options.useCluster) {
			swarm_->options.parallelCount = swarm_->options.swarmSize;
		}
		else {
			swarm_->options.parallelCount = stoi(pairs.find("parallel_count")->second);
		}

	}
	// Whether or not to divide by value at t=0
	if(pairs.find("divide_by_init") != pairs.end()) {
		cout << "Processing divide by init" << endl;
		swarm_->options.divideByInit = (stoi(pairs.find("divide_by_init")->second.c_str()) == 1) ? true : false;
	}

	// Whether or not to log transform simulation output
	if(pairs.find("log_transform_sim_data") != pairs.end()) {
		swarm_->options.logTransformSimData = stoi(pairs.find("log_transform_sim_data")->second);
	}

	// Whether or not to standardize simulation output
	if(pairs.find("standardize_sim_data") != pairs.end()) {
		swarm_->options.standardizeSimData = (stoi(pairs.find("standardize_sim_data")->second) == 1) ? true : false;
	}

	// Whether or not to standardize exp data
	if(pairs.find("standardize_exp_data") != pairs.end()) {
		swarm_->options.standardizeExpData = (stoi(pairs.find("standardize_exp_data")->second) == 1) ? true : false;
	}

	// Update fit calculation method
	if(pairs.find("objfunc") != pairs.end()) {
		cout << "Processing objfunc" << endl;
		swarm_->options.objFunc = stoi(pairs.find("objfunc")->second);
	}


	//////////////////////////////////////////////////////////////////
	// PSO STUFF
	//////////////////////////////////////////////////////////////////

	if (pairs.find("inertia") != pairs.end()) {
		cout << "Processing inertia" << endl;
		swarm_->options.inertia = stof(pairs.find("inertia")->second);
	}

	if (pairs.find("cognitive") != pairs.end()) {
		cout << "Processing cognitive" << endl;
		swarm_->options.cognitive = stof(pairs.find("cognitive")->second);
	}

	if (pairs.find("social") != pairs.end()) {
		cout << "Processing social" << endl;
		swarm_->options.social = stof(pairs.find("social")->second);
	}

	if (pairs.find("inertia") != pairs.end()) {
		cout << "Processing inertia" << endl;
		swarm_->options.inertia = stof(pairs.find("inertia")->second);
	}

	if (pairs.find("nmin") != pairs.end()) {
		cout << "Processing nmin" << endl;
		swarm_->options.nmin = stoi(pairs.find("nmin")->second);
	}

	if (pairs.find("nmax") != pairs.end()) {
		cout << "Processing nmax" << endl;
		swarm_->options.nmax = stoi(pairs.find("nmax")->second);
	}

	if (pairs.find("inertia_init") != pairs.end()) {
		cout << "Processing inertia_init" << endl;
		swarm_->options.inertiaInit = stof(pairs.find("inertia_init")->second);
	}

	if (pairs.find("inertia_final") != pairs.end()) {
		cout << "Processing inertia_final" << endl;
		swarm_->options.inertiaFinal = stof(pairs.find("inertia_final")->second);
	}

	if (pairs.find("abs_tolerance") != pairs.end()) {
		swarm_->options.absTolerance = stof(pairs.find("abs_tolerance")->second);
	}

	if (pairs.find("rel_tolerance") != pairs.end()) {
		swarm_->options.relTolerance = stof(pairs.find("rel_tolerance")->second);
	}

	if (pairs.find("mutate_qpso") != pairs.end()) {
		swarm_->options.mutateQPSO = (stoi(pairs.find("mutate_qpso")->second) == 1) ? true : false;
	}

	if (pairs.find("beta_min") != pairs.end()) {
		swarm_->options.betaMin = stof(pairs.find("beta_min")->second);
	}

	if (pairs.find("beta_max") != pairs.end()) {
		swarm_->options.betaMax = stof(pairs.find("beta_max")->second);
	}

	if (pairs.find("topology") != pairs.end()) {
		swarm_->options.topology = pairs.find("topology")->second;
	}

	if (pairs.find("pso_type") != pairs.end()) {
		swarm_->options.psoType = pairs.find("pso_type")->second;
	}

	if(pairs.find("enhanced_stop") != pairs.end()) {
		swarm_->options.enhancedStop = (stoi(pairs.find("enhanced_stop")->second) == 1) ? true : false;
	}

	if(pairs.find("enhanced_inertia") != pairs.end()) {
		swarm_->options.enhancedInertia = (stoi(pairs.find("enhanced_inertia")->second) == 1) ? true : false;

	}

	if(pairs.find("num_islands") != pairs.end()) {
		swarm_->options.numIslands = stoi(pairs.find("num_islands")->second);

		if (swarm_->options.swarmSize % swarm_->options.numIslands != 0 && swarm_->options.fitType == "de") {
			outputError("Error: The number of islands must divide evenly into the population size. Quitting.");
		}
	}

	if(pairs.find("mutate_type") != pairs.end()) {
		swarm_->options.mutateType = stoi(pairs.find("mutate_type")->second);
	}

	if(pairs.find("crossover_rate") != pairs.end()) {
		swarm_->options.cr = stof(pairs.find("crossover_rate")->second);
	}

	if(pairs.find("migration_frequency") != pairs.end()) {
		swarm_->options.migrationFrequency = stoi(pairs.find("migration_frequency")->second);
	}

	if(pairs.find("num_to_migrate") != pairs.end()) {
		swarm_->options.numToMigrate = stoi(pairs.find("num_to_migrate")->second);
	}

	if(pairs.find("min_temp") != pairs.end()) {
		swarm_->options.minTemp = stof(pairs.find("min_temp")->second);
	}

	if(pairs.find("min_radius") != pairs.end()) {
		swarm_->options.minRadius = stof(pairs.find("min_radius")->second);
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
				if (swarm_->options.model->freeParams_.count(paramComponents[0]) > 0) {
					// Make sure 2nd and 3rd components are numeric
					if (isFloat(paramComponents[1]) && isFloat(paramComponents[2])) {
						swarm_->options.model->freeParams_.at(paramComponents[0])->setGenerationMethod(pair->first);
						swarm_->options.model->freeParams_.at(paramComponents[0])->setParameterName(paramComponents[0]);
						swarm_->options.model->freeParams_.at(paramComponents[0])->setGenMin(stof(paramComponents[1]));
						swarm_->options.model->freeParams_.at(paramComponents[0])->setGenMax(stof(paramComponents[2]));

						//cout << "setting " << paramComponents[0] << " to " << pair->first << ":" << paramComponents[1] << ":" << paramComponents[2] << endl;
					}
					else {
						swarm_->outputError("Error: Problem parsing your free parameter generation option in your .conf file. The min and/or max values were non-numeric.");
					}
				}
				else {
					cout << "Warning: We found a parameter '" << paramComponents[0] << "' in your .conf file, but don't see a matching free parameter specification in your model file. We will ignore this parameter." << endl;
				}
			}
			else {
				swarm_->outputError("Error: Problem parsing your free parameter generation option in your .conf file. Each parameter generation option requires three components: the parameter name, minimum, and maximum.");
			}
		}
	}

	//for (auto i : swarm_->options.model->freeParams_) {
	for (map<string, FreeParam*>::iterator i = swarm_->options.model->freeParams_.begin(); i != swarm_->options.model->freeParams_.end(); ++i) {
		if (!i->second) {
			outputError("Error: We found a free parameter '" + i->first + "' specified in your model file but can't find a matching parameter generator in your .conf file.");
		}
	}

	// Add any .exp files to the swarm
	std::pair <std::unordered_multimap<string,string>::iterator, std::unordered_multimap<string,string>::iterator> it;
	it = pairs.equal_range("exp_file");

	for (unordered_multimap<string,string>::iterator exp = it.first; exp != it.second; ++exp) {
		swarm_->addExp(exp->second);
	}

	if(pairs.find("bootstrap") != pairs.end()) {
		swarm_->options.bootstrap = stoi(pairs.find("bootstrap")->second);
		if (swarm_->options.bootstrap) {
			vector<map<string, map<string, map<double, unsigned int>>>> bootstrapMaps;
			swarm_->generateBootstrapMaps(bootstrapMaps);
			swarm_->bootstrapMaps = bootstrapMaps;
		}
	}

	// Add any .mutation rates to the swarm
	it = pairs.equal_range("mutate");

	for (unordered_multimap<string,string>::iterator exp = it.first; exp != it.second; ++exp) {
		swarm_->addMutate(exp->second);
		swarm_->hasMutate = true;

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

	for (map<string, Model::action>::iterator i = swarm_->options.model->actions.begin(); i != swarm_->options.model->actions.end();)
		//for (vector<Model::action>::iterator i = swarm_->options.model->actions.begin(); i != swarm_->options.model->actions.end();)
	{
		prefixedActions.push_back(i->first);

		if (swarm_->options.expFiles.count(i->first) == 1) {
			if (swarm_->options.verbosity >=3 ) {
				cout << "Linking action " << i->first << " with exp file: " << swarm_->options.expFiles[i->first]->getPath() << endl;
			}
			i->second.dataSet = swarm_->options.expFiles[i->first];
			++i;
		}
		else { // Have a prefix but no .exp file
			cout << "Warning: The model file specifies an action with the prefix '" << i->first << "' but there isn't a matching .exp file specified in your .conf file. We will ignore this action command." << endl;
			//i = swarm_->options.model->actions.erase(i);
			swarm_->options.model->actions.erase(++i);
		}
	}

	//vector<Data*> toDeleteExp;
	//for (auto &i : swarm_->options.expFiles){ // Have .exp but no prefix
	for (map<string, Data*>::iterator i = swarm_->options.expFiles.begin(); i != swarm_->options.expFiles.end(); ++i) {
		if(std::find(prefixedActions.begin(), prefixedActions.end(), i->first) == prefixedActions.end() ) {
			cout << "Warning: The .conf file specifies an .exp file '" << i->first << "' but there isn't a matching action command in your model file specified with the prefix=> argument." << endl;
			//toDeleteExp.push_back(swarm_->options.expFiles[i->first]);
			//i = swarm_->options.expFiles.erase(i);
			swarm_->options.expFiles.erase(++i);
		}
	}

	/*
	// TODO: Make sure we have either min fit, max time, or max sims when doing an asynchronous fit
	// This whole checker needs to be much more complex
	if (!swarm_->options.synchronicity) {
		if (swarm_->options.minFit == -1 && swarm_->options.maxFitTime == MAX_LONG && swarm_->options.maxNumSimulations == MAX_LONG && !swarm_->options.enhancedStop) {
			outputError("Error: You are running an asynchronous fit, but failed to specify either a minimum fitting value, maximum fit time, or maximum number of simulations in the .conf file.");
		}
	}
	 */

	/*
	for (auto i : swarm_->options.expFiles){
		cout << "Exp: " << i.second->getPath() << endl;
	}

	for (auto i : swarm_->options.modelFile->actions_){
		cout << "Prefix: " << i.prefix << endl;
	}

	for (auto i : swarm_->options.modelFile->actions_){
		cout << "DataSets: " << i.dataSet->getPath() << endl;
	}
	 */
	checkConsistency();

	return swarm_;
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
	 *	Enhanced stop requires:
	 *	1. nmax
	 *
	 *	Enhanced inertia requires:
	 *	1. Enhanced stop
	 *
	 *	In DE, different mutate types require different numbers of minimum particles in each island
	 *
	 * Make sure we have SOME sort of stop criteria
	 * QPSO requires maxNumSimulations to update beta_. the max num of iterations should be close to the expected. would be nice if the user had to guess. is there an adaptive beta_ algorithm like enhancedInertia?
	 */

	if (swarm_->options.fitType.empty()) {
		swarm_->outputError("Error: You didn't specify a fit_type in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	// TODO: This doesn't work!!
	if (swarm_->options.jobName.empty()) {
		swarm_->outputError("Error: You didn't specify a job_name in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.model == 0) {
		swarm_->outputError("Error: You didn't specify a model in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.expFiles.empty()) {
		swarm_->outputError("Error: You didn't specify any .exp files in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.outputDir.empty()) {
		swarm_->outputError("Error: You didn't specify an output_dir in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.bngCommand.empty()) {
		swarm_->outputError("Error: You didn't specify a bng_command in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.swarmSize == 0) {
		swarm_->outputError("Error: You didn't specify a swarm_size in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.enhancedInertia) {
		if (!swarm_->options.enhancedStop) {
			swarm_->outputError("Error: You set enhanced_inertia to true, but enhanced_stop to false. Enhanced inertia requires enhanced_stop to be set to true. Quitting.");
		}
	}

	if (swarm_->options.enhancedStop) {
		if (swarm_->options.nmax == 0) {
			swarm_->outputError("Error: You set enhanced_stop to true, but did not set a maximum permance value (nmax) in your .conf file. Quitting.");
		}
	}

	if (swarm_->options.useCluster && swarm_->options.emailWhenFinished && swarm_->options.emailAddress.empty()) {
		swarm_->outputError("Error: You set email_when_finished to true, but did not specify an email address. Quitting.");
	}
}
