/*
 * Swarm.cpp
 *
 *  Created on: Jul 17, 2015
 *      Author: brandon
 */

#include "Swarm.hh"

using namespace std;

Swarm::Swarm(bool isMaster) {

	// Whether or not we are master
	isMaster_ = isMaster;

	// Create the comunication class
	Pheromones *ph = new Pheromones();

	// Initialize the communication class
	ph->init(this);

	// Set our communication class
	swarmComm_ = ph;
}

void Swarm::addExp(string path) {
	path = convertToAbsPath(path);
	string basename = getFilename(path);

	//cout << "inserting EXP?" << endl;
	this->options_.expFiles.insert(make_pair(basename, new Data(path, this, true)));
}

void Swarm::setModel(string path) {
	path = convertToAbsPath(path);
	this->options_.model = new Model(path);
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

		//allParticles_.at(2)->setModel(options_.model);
		//allParticles_.at(2)->generateParams();

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
		while (currentGeneration_ < options_.maxGenerations){
			runGeneration();

			string currentDirectory = options_.outputDir + "/" + to_string(currentGeneration_);
			cleanupFiles(currentDirectory.c_str());
			breedGeneration();
		}
	}
}

Particle * Swarm::createParticle(int pID) {
	Particle * p = new Particle(this, pID);
	p->setBasePath(particleBasePath_);

	return p;
}

unordered_map<int,Particle*> Swarm::generateInitParticles(int pID) {

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
}

void Swarm::launchParticle(Particle *p) {

	if (options_.useCluster) {

	}
	else {
		int i;
		string command = exePath_ + " particle " + to_string(p->getID()) + " run " + to_string(currentGeneration_) + " " + configPath_;
		//if (p->getID() == 1) {
		command = command + ">> pOUT 2>&1";
		//}
		command = command + " &";
		// Create the pipe which will be used to communicate with the particle
		//string path = particleBasePath_ + to_string(p->getID());
		//createParticlePipe(path.c_str());

		cout << "Running particle with command: " << command << endl;
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
		cout << "Running generation " << currentGeneration_ << " with " << allParticles_.size() << " particles..." << endl;
	}

	currentGeneration_ += 1;

	string createDirCmd = "mkdir " + options_.outputDir + "/" + to_string(currentGeneration_);
	system(createDirCmd.c_str());

	std::vector<std::vector<std::string>> messageHolder;
	int numMessages = 0;
	int numLaunchedParticles = 0;
	int numFinishedParticles = 0;
	unordered_map<int,Particle*>::iterator p = allParticles_.begin();

	while (numFinishedParticles < options_.swarmSize) {

		// Launch particles (staying within bounds of parallel_count)
		if (runningParticles_.size() < options_.parallelCount && numLaunchedParticles < options_.swarmSize) {
			launchParticle(p->second);
			numLaunchedParticles += 1;
			++p;
		}

		// Check for any messages from particles
		usleep(250000);
		numMessages = swarmComm_->recvMessage(-1, 0, -1, false, messageHolder);

		if (numMessages >= 1) {
			//cout << "found message" << endl;
			// First loop through individual messages
			for (std::vector<std::vector<std::string>>::iterator i = messageHolder.begin(); i != messageHolder.end(); ++i) {
				//cout << "message loop" << endl;
				/// Loop through message contents
				for (std::vector<std::string>::iterator o = (*i).begin(); o != (*i).end(); ++o) {
					//cout << "content loop" << endl;
					if (stoi(*o) == SIMULATION_END || stoi(*o) == SIMULATION_FAIL) {
						//cout << "found simulation end" << endl;
						// Jump to the particle ID and store it
						o+=2;
						int pID = stoi(*o);

						// Get an iterator to the particle in our list of running particles
						runningParticlesIterator_ = runningParticles_.find(pID);
						// Then remove it
						runningParticles_.erase(runningParticlesIterator_);
						// Increment our finished counter
						numFinishedParticles += 1;
						cout << "particle " << pID << " finished simulation" << endl;
						// Jump back to the message tag
						o-=2;
						if (stoi(*o) == SIMULATION_FAIL) {
							o+=2;
							cout << "Particle " << *o << " failed" << endl;

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
							std::vector<std::string>::iterator messageStart = o;
							string params = "gen" + to_string(currentGeneration_) + "perm" + to_string(pID) + " ";

							o+=4; // Jump to the fit value
							//map<string,string>::iterator fp = options_.model->freeParams_.begin();
							while(1) {
								o++;
								if (o != (*i).end()) {
									params += *o + " ";
									//cout << "adding " << *o << ". should be " << fp->first << endl;
								}
								else {
									//cout << "message end" << endl;
									break;
								}
								//fp++;
							}

							o = messageStart;

							o+=4;

							// Then store it
							allParticles_.at(pID)->fitCalcs[currentGeneration_] = stod(*o);
							//currGenFits.push_back(allParticles_.at(pID));
							allGenFits.insert(pair<double,string>(stod(*o),params));
							//cout << "particle " << pID << " has fit calc of " << *o;
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
			}
			messageHolder.clear();
			numMessages = 0;
		}
	}
	if (failedParticles_.size() > (options_.swarmSize - 3) ) {
		outputError("Error: You had too many failed runs. Check simulation output or adjust walltime.");
		//TODO: Cleanup and exit
		exit(EXIT_FAILURE);
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
			if (regex_match(string( de->d_name ),regex(".+xml$|.+species$|.+cdat$|.+gdat$|.+net$"))) {
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
	// We let particles do the actual breeding.  The master's role is to generate a breeding pattern and
	// tell which particles to breed with which

	cout << "fit       perm     ";

	for (auto i: options_.model->freeParams_) {
		cout << i.first << " ";
	}

	cout << endl;

	for (auto i: allGenFits) {
		cout << i.first << ": " << i.second << endl;
	}

	map<double,double> weights; // First element is original fit, second element is subtracted from max

	// Fill in the weight map with fit values
	double weightSum;
	for (map<double,string>::iterator f = allGenFits.begin(); f != allGenFits.end(); ++f) {
		weights.insert(pair<double,double>(f->first,0));
		weightSum += f->first;
	}

	double maxWeight = weights.end()->first; // TODO: Is this true??

	// Fill the second element of the weight map with difference between maxWeight and fit value
	for (map<double,double>::iterator w = weights.begin(); w != weights.end(); ++w) {
		w->second = maxWeight - w->first;
	}

	int parentPairs = options_.swarmSize / 2;

	cout << "we have " << parentPairs << " parent pairs" << endl;

	double p1;
	double p2;
	smatch match;
	vector<string> parentVec;

	for (int i = 0; i <= parentPairs; ++i) {
		// Pick the fit values (particle parents) used in breeding
		p1 = pickWeighted(weightSum, weights, options_.extraWeight);
		p2 = pickWeighted(weightSum, weights, options_.extraWeight);

		// If we want different parents used in breeding, make sure that happens
		int retryCount = 0;
		while (p1 == p2 && options_.forceDifferentParents) {
			retryCount++;
			if (retryCount > options_.maxRetryDifferentParents) {
				cout << "Tried to many time to select different parents for breeding. Selecting the first two." << endl;
				map<double,double>::iterator w = weights.begin();
				p1 = w->first;
				++w;
				p2 = w->first;

				break;
			}
			p2 = pickWeighted(weightSum, weights, options_.extraWeight);
		}

		cout << "selected " << p1 << " and " << p2 << endl;

		// We always send the breeding information to the first parent. That particle will then
		// initiate breeding with its mate

		// Get the string containing gen number, perm number, and param values
		string parentString1 = allGenFits.at(p1);
		string parentString2 = allGenFits.at(p2);

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

		cout << "p1 is " << pID1 << " p2 is " << pID2 << endl;

		// Add second parent to the message
		swarmComm_->univMessageSender.push_back(pID2);
		// Send the message to the first parent
		swarmComm_->sendToSwarm(0, stoi(pID1), BREED_PARENT, false, swarmComm_->univMessageSender);
		swarmComm_->univMessageSender.clear();
	}
}
