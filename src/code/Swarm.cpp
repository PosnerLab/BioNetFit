/*
 * Swarm.cpp
 *
 *  Created on: Jul 17, 2015
 *      Author: brandon
 */

#include "Swarm.hh"

using namespace std;
using namespace std::chrono;

Swarm::Swarm(bool isMaster) {

	// Whether or not we are master
	isMaster_ = isMaster;

	// Create the comunication class
	Pheromones *ph = new Pheromones();

	// Initialize the communication class
	ph->init(this);

	// Set our communication class
	swarmComm_ = ph;

	// Seed our random number engine
	randNumEngine.seed(std::random_device{}());
	randNumEngine.discard(700000);
	srand (std::random_device{}());
}

void Swarm::addExp(string path) {
	Timer tmr;

	path = convertToAbsPath(path);
	string basename = getFilename(path);

	//cout << "inserting EXP?" << endl;
	this->options_.expFiles.insert(make_pair(basename, new Data(path, this, true)));
	//cout << "exp inserted" << endl;

	double t = tmr.elapsed();
	cout << "Adding .exp took " << t << " seconds" << endl;
}

void Swarm::setModel(string path) {
	Timer tmr;

	path = convertToAbsPath(path);
	//cout << "setting model" << endl;
	this->options_.model = new Model(path);
	//cout << "model set" << endl;

	double t = tmr.elapsed();
	cout << "Adding .bngl took " << t << " seconds" << endl;
}

void Swarm::setSwarmSize(int size) {
	this->options_.swarmSize = size;
}

void Swarm::setSwarmType(string type) {
	this->options_.swarmType = type;
}

void Swarm::setSwarmSynchronicity(int synchronicity) {
	this->options_.synchronicity = synchronicity;
}

void Swarm::setSwarmGenerations(int generations) {
	this->options_.maxGenerations = generations;
}

void Swarm::setSwarmMinFit(float minfit) {
	this->options_.minFit = minfit;
}

bool Swarm::checkSwarmConsistency() {
	return true;
}

void Swarm::addMutate(std::string mutateString) {
	vector<string> mutComponents;
	split(mutateString, mutComponents);

	// Make sure we have three components to work with
	if (mutComponents.size() == 3) {
		// Make sure first parameter name exists as a free parameter
		if (options_.model->freeParams_.count(mutComponents[0]) > 0 || mutComponents[0] == "default") {
			// Make sure 2nd and 3rd components are numeric
			if (isFloat(mutComponents[1]) && isFloat(mutComponents[2])) {
				if (mutComponents[0] != "default") {
					cout << "not default" << endl;
					options_.model->freeParams_.at(mutComponents[0])->setMutationRate(stof(mutComponents[1]));
					options_.model->freeParams_.at(mutComponents[0])->setMutationFactor(stof(mutComponents[2]));
					options_.model->freeParams_.at(mutComponents[0])->setHasMutation(true);

					//cout << "setting " << mutComponents[0] << " to " << mutComponents[1] << ":" << mutComponents[2] << endl;
				}
				else {
					for (map<string,FreeParam*>::iterator fp = options_.model->freeParams_.begin(); fp != options_.model->freeParams_.end(); ++fp) {
						//cout << "setting " << mutComponents[0] << " to " << mutComponents[1] << ":" << mutComponents[2] << endl;

						fp->second->setMutationRate(stof(mutComponents[1]));
						fp->second->setMutationFactor(stof(mutComponents[2]));
						fp->second->setHasMutation(true);
					}
				}
			}
			else {
				outputError("Error: Problem parsing the mutation option in your .conf file. The mutation rate and/or factor are non-numeric.");
			}
		}
		else {
			cout << "Warning: We found a mutation option '" << mutComponents[0] << "' in your .conf file, but don't see a matching free parameter specification in your model file. We will ignore this mutation factor." << endl;
		}
	}
	else {
		outputError("Error: Problem parsing a mutation option in your .conf file. Each mutation option requires three components: the parameter name, mutation rate, and mutation factor.");
	}
}

void Swarm::doSwarm() {
	system ("exec rm -r /home/brandon/projects/GenFit2/Debug/output/*");
	system ("exec rm /home/brandon/projects/GenFit2/Debug/pOUT");

	// Generate all particles that will be present in the swarm
	allParticles_ = generateInitParticles();

	// If we are running ODE, we want to generate a network and network reader.
	// The network file has parameters that are replaced in every generation,
	// and the network reader is used to read those network files
	if (options_.model->getHasGenerateNetwork()) {

		// First create a particle which will generate the network
		// and fill it with dummy params
		allParticles_.at(1)->setModel(options_.model);
		allParticles_.at(1)->generateParams();

		// Output model file with dummy parameters. This file will be used to generate
		// our initial network
		options_.model->outputModelWithParams(allParticles_.at(1)->getParams(), options_.outputDir, "base.bngl", "", false, false, false, false, false);

		// Construct our simulation command and run the network generator
		string modelPath = options_.outputDir + "/base.bngl";
		string command = "/home/brandon/projects/GenFit/Simulators/BNG2.pl --outdir " + options_.outputDir + " " + modelPath + " >> " + options_.outputDir + "/netgen_output 2>&1";
		int ret = system(command.c_str());

		// Now that we have a .net file, replace the full .bngl with a .bngl
		// containing ONLY action commands and a .net file loader. This is
		// used later to run our .net files.
		options_.model->outputModelWithParams(allParticles_.at(1)->getParams(), options_.outputDir, "base.bngl", "", false, true, false, false, false);

		// Now store our .net file for later use
		string netPath = options_.outputDir + "/base.net";
		options_.model->parseNet(netPath);
	}

	// Main swarming loops
	if (options_.synchronicity == 1) {

		string createDirCmd = "mkdir " + options_.outputDir + "/1";
		system(createDirCmd.c_str());

		while (currentGeneration_ < options_.maxGenerations){
			runGeneration();

			string currentDirectory = options_.outputDir + "/" + to_string(currentGeneration_);
			if (options_.deleteOldFiles) {
				cleanupFiles(currentDirectory.c_str());
			}

			if (currentGeneration_ < options_.maxGenerations) {
				breedGeneration();
			}
		}
		cout << "fit finished " << endl;
		finishFit();
	}
}

Particle * Swarm::createParticle(int pID) {
	Particle * p = new Particle(this, pID);
	p->setBasePath(particleBasePath_);

	return p;
}

unordered_map<int,Particle*> Swarm::generateInitParticles(int pID) {
	Timer tmr;

	unordered_map<int,Particle*> allParticles;

	// TODO: Do we really need to create particle objects if we're a master?
	// Why not just keep track of them as integers? Maybe as backups in case
	// any of them fail??
	if (pID == -1) {
		for (int i = 1; i <= options_.swarmSize; i++) {
			//Particle *p = new Particle(i);
			//p->setBasePath(particleBasePath_);
			allParticles[i] = createParticle(i);
		}
	}
	return allParticles;

	double t = tmr.elapsed();
	cout << "Particle creation took " << t << " seconds" << endl;
}

void Swarm::launchParticle(Particle *p) {

	if (options_.useCluster) {

	}
	else {
		//string rm = "rm " + to_string(p->getID());
		//system(rm.c_str());
		int i;
		string command = exePath_ + " particle " + to_string(p->getID()) + " run " + to_string(currentGeneration_) + " " + configPath_;
		//if (p->getID() == 10) {
		//command = command + ">> " + to_string(p->getID()) + " 2>&1";
		command = command + ">> pOUT 2>&1";
		//}
		command = command + " &";
		// Create the pipe which will be used to communicate with the particle
		//string path = particleBasePath_ + to_string(p->getID());
		//createParticlePipe(path.c_str());

		if (options_.verbosity >= 3) {
			cout << "Running Particle " << p->getID() << endl;
		}

		i = system(command.c_str());
		//runningParticles_[p] = "";
		runningParticles_.insert(p->getID());
	}
}

void Swarm::setUseCluster(bool useCluster) {
	options_.useCluster = useCluster;
}

void Swarm::setCurrentGen(int gen) {
	currentGeneration_ = gen;
}

void Swarm::runGeneration () {
	// TODO: Implement walltime
	if(options_.verbosity >= 1) {
		cout << "Running generation " << currentGeneration_ + 1 << " with " << allParticles_.size() << " particles..." << endl;
	}

	currentGeneration_ += 1;

	std::vector<std::vector<std::string>> messageHolder;
	int numMessages = 0;
	int numLaunchedParticles = 0;
	int numFinishedParticles = 0;
	unordered_map<int,Particle*>::iterator p = allParticles_.begin();

	while (numFinishedParticles < options_.swarmSize) {
		// Launch particles (staying within bounds of parallel_count)
		if (runningParticles_.size() < options_.parallelCount && numLaunchedParticles < options_.swarmSize) {
			if (currentGeneration_ == 1) {
				launchParticle(p->second);
			}
			else {
				swarmComm_->sendToSwarm(0, p->first, NEXT_GENERATION, false, swarmComm_->univMessageSender);
				runningParticles_.insert(p->first);
			}
			numLaunchedParticles += 1;
			++p;
		}

		// Check for any messages from particles
		usleep(10000);
		//cout << "checking messages" << endl;
		numMessages = swarmComm_->recvMessage(-1, 0, -1, false, messageHolder);

		if (numMessages >= 1) {
			//cout << "found message" << endl;
			// First loop through individual messages
			for (std::vector<std::vector<std::string>>::iterator i = messageHolder.begin(); i != messageHolder.end(); ++i) {
				/// Loop through message contents
				//for (std::vector<std::string>::iterator o = (*i).begin(); o != (*i).end(); ++o) {
				std::vector<std::string>::iterator o = (*i).begin();

				if (stoi(*o) == SIMULATION_END || stoi(*o) == SIMULATION_FAIL) {

					// Jump to the particle ID and store it
					o+=2;
					int pID = stoi(*o);
					o-=2;

					//cout << "found simulation end for" << pID << endl;

					// Get an iterator to the particle in our list of running particles
					runningParticlesIterator_ = runningParticles_.find(pID);
					// Then remove it
					runningParticles_.erase(runningParticlesIterator_);
					// Increment our finished counter
					numFinishedParticles += 1;
					//cout << "particle " << pID << " finished simulation" << endl;
					// Jump back to the message tag

					if (stoi(*o) == SIMULATION_FAIL) {
						o+=2;
						cout << "Particle " << pID << " failed in gen " << currentGeneration_ << endl;

						// Store particle ID in our list of failed particles
						failedParticles_.insert(stoi(*o));

						// Give the particle a very large fit value
						//allParticles_.at(pID)->fitCalcs[currentGeneration_] = 1E10;
						//currGenFits.push_back(allParticles_.at(pID));
						//allGenFits.insert(pair<double,string>(-1,""));
					}
					else {

						// [TAG]
						// [ID]
						// [SENDER]
						// [MESSAGESIZE]
						// [ARR0]
						// [ARR1]
						// ...

						string params = "gen" + to_string(currentGeneration_) + "perm" + to_string(pID) + " ";

						o+=4; // Jump to the fit value
						double fitCalc = stod(*o);

						// Store the parameters given to us by the particle
						while(1) {
							o++;
							if (o != (*i).end()) {
								params += *o + " ";
								//cout << "adding " << *o<< endl;
							}
							else {
								//cout << "message end" << endl;
								break;
							}
						}

						// Then store it
						//allParticles_.at(pID)->fitCalcs[currentGeneration_] = stod(*o);
						//currGenFits.push_back(allParticles_.at(pID));
						//cout << "saving fit for " << pID << " of " << fitCalc << endl;
						allGenFits.insert(pair<double,string>(fitCalc,params));
						//cout << "particle " << pID << " has fit calc of " << allGenFits.find(fitCalc)->first << endl;
					}
					//cout << "finished processing simulation end" << endl;
				}
				else if (stoi(*o) == GET_RUNNING_PARTICLES) {
					vector<string> intParts;
					for (auto i: runningParticles_) {
						intParts.push_back(to_string(i));
					}
					o+=2;
					swarmComm_->sendToSwarm(0,stoi(*o),SEND_RUNNING_PARTICLES,true,intParts);
				}
			}
			messageHolder.clear();
			numMessages = 0;
		}
	}
	if (failedParticles_.size() > (options_.swarmSize - 3) ) {
		outputError("Error: You had too many failed runs. Check simulation output or adjust walltime.");
		//TODO: Cleanup and exit
	}

	// Sort current generation particle fit vector
	//sort(currGenFits.begin(), currGenFits.end(), [this](Particle * a, Particle * b) {return sortFits(a, b);});
}

void Swarm::cleanupFiles(const char * path) {
	DIR* dp;
	dirent* de;
	errno = 0;

	vector<string> filesToDel;

	dp = opendir(path);

	if (dp)
	{
		while (true)
		{
			//cout << "about to readdir" << endl;
			errno = 0;
			de = readdir( dp );
			if (de == NULL) break;
			//cout << "read " << de->d_name << endl;
			if (regex_match(string( de->d_name ),regex(".+xml$|.+species$|.+cdat$|.+gdat$|.+net$|.+BNG_OUT"))) {
				//cout << "found a file";
				filesToDel.push_back( string(de->d_name) );
			}
		}
		closedir( dp );
	}
	else {
		cout << "Warning: Couldn't open " << path << " to delete un-needed simulation files." << endl;
	}

	for (auto i: filesToDel) {
		string fullPath = string(path) + "/" + i;
		remove(fullPath.c_str());
	}
}

bool Swarm::sortFits(Particle * a, Particle * b) {
	return a->fitCalcs.at(currentGeneration_) < b->fitCalcs.at(currentGeneration_);
}

void Swarm::breedGeneration() {
	//Timer tmr;

	string createDirCmd = "mkdir " + options_.outputDir + "/" + to_string(currentGeneration_ + 1);
	system(createDirCmd.c_str());

	// We let particles do the actual breeding.  The master's role is to generate a breeding pattern and
	// tell which particles to breed with which
	/*
	cout << "fit\tperm\t";

	for (auto i: options_.model->freeParams_) {
		cout << i.first << "\t";
	}

	cout << endl;

	for (auto i: allGenFits) {
		cout << i.first << ": " << i.second << endl;
	}*/

	multimap<double,double> weights; // First element is original fit, second element is subtracted from max

	// Fill in the weight map with fit values
	double weightSum;

	for (multimap<double,string>::iterator f = allGenFits.begin(); f != allGenFits.end(); ++f) {
		weights.insert(pair<double,double>(f->first,0));
		weightSum += f->first;
		//cout << "f: " << f->first << endl;
	}

	double maxWeight = weights.rbegin()->first; // TODO: Is this true??
	//cout << "max: " << maxWeight << endl;

	// Fill the second element of the weight map with difference between maxWeight and fit value
	for (map<double,double>::iterator w = weights.begin(); w != weights.end(); ++w) {
		w->second = maxWeight - w->first;
	}
	/*
	for (auto i: weights) {
		cout << "weight diff: " << i.second << endl;
	}*/
	int parentPairs = options_.swarmSize / 2;

	//cout << "we have " << parentPairs << " parent pairs" << endl;

	double p1;
	double p2;
	smatch match;
	vector<string> parentVec;
	int numCurrBreeding = 0;
	int numFinishedBreeding = 0;

	//cout << "sum: " << weightSum << endl;

	int swapID = -1;
	for (int i = 0; i < parentPairs; ++i) {
		swapID += 2;

		// Pick the fit values (particle parents) used in breeding
		p1 = pickWeighted(weightSum, weights, options_.extraWeight, randNumEngine);
		p2 = pickWeighted(weightSum, weights, options_.extraWeight, randNumEngine);

		// If we want different parents used in breeding, make sure that happens
		int retryCount = 0;
		while (p1 == p2 && options_.forceDifferentParents) {
			retryCount++;
			if (retryCount > options_.maxRetryDifferentParents) {
				if (options_.verbosity >= 3) {
					cout << "Tried to many time to select different parents for breeding. Selecting the first two." << endl;
				}

				// Get iterator to the weight map
				multimap<double,double>::iterator w = weights.begin();

				// The weight map is sorted, so the first element will be the best fit
				p1 = w->first;

				// Increment the map iterator until we find a fit value that isn't ours
				while (p1 == p2 && w != weights.end()) {
					++w;
					p2 = w->first;
					sleep(1);
				}

				break;
			}
			p2 = pickWeighted(weightSum, weights, options_.extraWeight, randNumEngine);
		}

		//cout << "selected " << p1 << " and " << p2 << endl;

		// We always send the breeding information to the first parent. That particle will then
		// initiate breeding with its mate

		// Get the string containing gen number, perm number, and param values
		// We're searching a multimap which may contain duplicate fit values,
		// but it shouldn't matter which we pick because if two runs have the
		// exact same fit, they should have the same parameters sets as each other
		string parentString1 = allGenFits.find(p1)->second;
		string parentString2 = allGenFits.find(p2)->second;

		// Split that string into array. We're only interested in the first element which contains
		// the particle ID/permutation number
		split(parentString1,parentVec);

		// Extract the particle ID
		regex_search(parentVec[0], match, regex("gen\\d+perm(\\d+)"));
		string pID1 = match[1];

		// Clear the parent vector for use with the next parent
		parentVec.clear();

		// Do same as above, but for the second parent
		split(parentString2,parentVec);
		regex_search(parentVec[0], match, regex("gen\\d+perm(\\d+)"));
		string pID2 = match[1];
		parentVec.clear();

		//cout << "p1 is " << pID1 << " p2 is " << pID2 << endl;

		// Add second parent and swapID to the message
		swarmComm_->univMessageSender.push_back(to_string(swapID));
		swarmComm_->univMessageSender.push_back(pID2);

		// Send the message to the first parent
		//cout << "sending message to " << stoi(pID1) << " with id of " << to_string(swapID) << endl;

		// Tell parent 1 to initiate breeding
		swarmComm_->sendToSwarm(0, stoi(pID1), INIT_BREEDING, false, swarmComm_->univMessageSender);
		swarmComm_->univMessageSender.clear();
		numCurrBreeding+=2;

		//double t = tmr.elapsed();
		//cout << "BREEDING took " << t << " seconds" << endl;
		//tmr.reset();

		// While the current number of breeding particles is greater than twice the parallel count
		while ( numCurrBreeding >= (options_.parallelCount * 15)) {
			//cout << "checking for finished" << endl;
			int numMessages = swarmComm_->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm_->univMessageReceiver, true);
			//cout << "found " << numMessages << " finished" << endl;
			numFinishedBreeding+=numMessages;
			numCurrBreeding-=numMessages;
		}
	}
	while (numFinishedBreeding < options_.swarmSize) {
		int numMessages = swarmComm_->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm_->univMessageReceiver, true);
		numFinishedBreeding+=numMessages;
	}
	//cout << "got them all!" << endl;
}

void Swarm::finishFit() {
	string outputDir = options_.outputDir + "/Results";
	string command = "mkdir " + outputDir;

	if (system(command.c_str())) {
		string errMsg = "Error: Couldn't create directory: " + outputDir + " to contain final fitting results.";
		outputError(errMsg);
	}

	outputRunSummary(outputDir);
	//generateBestFitModels(outputDir);
	//copyBestFitToResults(outputDir);
	killAllParticles(FIT_FINISHED);

	cout << "Finished fitting in " << tmr_.elapsed() << " seconds. Results can be found in " << options_.outputDir << "/Results" << endl;
}

void Swarm::outputRunSummary(string outputDir) {
	string summaryFilename = outputDir + "/all_fits.txt";
	ofstream outputFile;

	outputFile.open(summaryFilename, ofstream::out);

	if (outputFile.is_open()) {
		outputFile.precision(8);
		outputFile.setf(ios::scientific);

		// Output first two fields of header
		outputFile << left << setw(16) << "Fit" << left << setw(16) << "Permutation";

		// Output parameter names
		for (auto i: options_.model->freeParams_) {
			outputFile << left << setw(16) << i.first;
		}

		outputFile << endl;

		vector<string> paramVals;

		for (auto i: allGenFits) {
			split(i.second, paramVals);

			outputFile << left << setw(16) << i.first << left << setw(16) << paramVals[0];

			for (int i = 1; i < paramVals.size(); i++) {
				outputFile << left << setw(16) << stod(paramVals[i]);
			}

			paramVals.clear();

			outputFile << endl;
		}

		outputFile.close();
	}
	else {
		cout << "Warning: couldn't open: " << outputFile << " to write fit summary." << endl;
	}
}

void Swarm::killAllParticles(int tag) {
	for (auto p: allParticles_) {
		swarmComm_->sendToSwarm(0, p.first, tag, false, swarmComm_->univMessageSender);
	}
}

void Swarm::getClusterInformation() {

	// If user didn't specify cluster platform, let's figure it out ourself
	if (options_.clusterSoftware.size() == 0) {
		// Test for slurm
		if (getOutputFromCommand("which srun").length() > 0) {
			options_.clusterSoftware = "slurm";
		}
		// Test for PBS-type
		else if (getOutputFromCommand("which qsub").length() > 0) {
			// Test for Torque/PBS
			if(getOutputFromCommand("which maui").length() > 0 || getOutputFromCommand("which moab").length() > 0) {
				options_.clusterSoftware = "torque";
			}
			// Test for SGE
			else if (getOutputFromCommand("which sge_execd").length() > 0 || getOutputFromCommand("which qconf").length() > 0 || getOutputFromCommand("which qmon").length() > 0) {
				outputError("Error: BioNetFit doesn't support for GridEngine clusters. If you are not running on a GridEngine cluster, specify the cluster platform in the .conf file using the 'cluster_software' option.");
			}
		}
	}
	else {
		if (options_.clusterSoftware != "slurm" && options_.clusterSoftware != "torque") {
			outputError("You specified an unrecognized cluster type in your .conf file. BioNetFit only supports 'torque' or 'slurm' cluster types.");
		}
	}

	// If we still don't know the cluster type, let's ask the user.
	if (options_.clusterSoftware.size() == 0) {
		string input;
		string clusterType;

		while (1) {
			cout << "BioNetFit couldn't determine which type of cluster software you are using. Specify (T) for Torque, or (S) for Slurm" << endl;
			getline(cin, input);
			stringstream(input) >> clusterType;

			if (clusterType == "S" || clusterType == "s") {
				options_.clusterSoftware = "slurm";
			}
			else if (clusterType == "T" || clusterType == "t") {
				options_.clusterSoftware = "torque";
			}
		}
	}
}

string Swarm::generateSlurmCommand(string cmd) {
	string command;

	// srun submits the job to the cluster
	command += "srun";

	// Add the job name
	command += " -J " + options_.jobName;

	// Only need one CPU per particle
	command += " -c 1";

	// Specify the cluster account if needed
	if (!options_.clusterAccount.empty()) {
		command += " -A " + options_.clusterAccount;
	}

	if (!options_.clusterQueue.empty()) {
		command += " -p	" + options_.clusterQueue;
	}

	// Specify output directory if needed
	if (options_.saveClusterOutput) {
		command += " -o " + options_.jobOutputDir + options_.jobName + "_cluster_output";
	}
	else {
		command += " -o /dev/null";
	}

	command+= " " + cmd;
}
