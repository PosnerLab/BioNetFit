/*
 * Swarm.cpp
 *
 *  Created on: Jul 17, 2015
 *      Author: brandon
 */

#include "Swarm.hh"

using namespace std;
using namespace std::chrono;

Swarm::Swarm() {

	// Whether or not we are master
	isMaster = false;
	sConf_ = "";

	options.simPath = "";
	options.maxGenerations = 10;
	options.swarmSize = 10;
	options.minFit = -1;
	options.maxFit = 0;
	options.boostrap = 0;
	options.parallelCount = 0;
	options.useCluster = false;
	options.saveClusterOutput = false;
	options.usePipes = false;

	options.divideByInit = false;
	options.logTransformSimData = false;
	options.standardizeSimData = false;
	options.standardizeExpData = false;

	options.deleteOldFiles = true;
	options.objFunc = 1;
	options.extraWeight = 0;
	options.swapRate = 0.5;
	options.forceDifferentParents = true;
	options.maxRetryDifferentParents = 100;

	options.maxFitTime = MAX_LONG;
	options.maxNumSimulations = MAX_LONG;

	// PSO options
	options.inertia = 0.72; // 0.72
	options.cognitive = 1.49; // 1.49
	options.social = 1.49; // 1.49
	options.nmax = 20; // 20
	options.nmin = 80; // 80
	options.inertiaInit = 1; // 1
	options.inertiaFinal = 0.1; // 0.1
	options.absTolerance = 10e-4; // 10E-4
	options.relTolerance = 10e-4; // 10E-4

	options.topology = "fullyconnected"; // fullyconnected
	options.psoType = "pso"; // pso

	options.enhancedStop = true; // true
	options.enhancedInertia = true; // true

	options.verbosity = 1;
	options.hasMutate = false;

	/*
	std::map<int, double> particleBestFits_;
	std::map<int, std::vector<double>> particleParamVelocities_;
	std::map<int, std::vector<double>> particleBestParamSets_;
	std::map<int, std::vector<double>> particleCurrParamSets_;
	std::map<int, double> particleWeights_;
	std::map<double, int> particleBestFitsByFit_;
	 */

	int permanenceCounter_ = 0; // 0
	int flightCounter_ = 0; // 0
	double weightedAvgPos_ = 0; // 0
	double optimum_ = 0; // 0
	int inertiaUpdateCounter_ = 0; // 0;

	// TODO: Make sure everything is being seeded properly and in the proper place. Also let's do away with rand()
	// Seed our random number engine
	randNumEngine.seed(std::tr1::random_device{}());
	//randNumEngine.discard(700000);

	srand (std::tr1::random_device{}());
}

void Swarm::initComm() {
	// Create the comunication class
	Pheromones *ph = new Pheromones();

	// Initialize the communication class
	ph->init(this);

	// Set our communication class
	swarmComm = ph;
}

void Swarm::addExp(string path) {
	Timer tmr;

	path = convertToAbsPath(path);
	string basename = getFilename(path);

	//cout << "inserting EXP?" << endl;
	this->options.expFiles.insert(make_pair(basename, new Data(path, this, true)));
	//cout << "exp inserted" << endl;

	//double t = tmr.elapsed();
	//cout << "Adding .exp took " << t << " seconds" << endl;
}

void Swarm::setModel(string path) {
	Timer tmr;

	path = convertToAbsPath(path);
	//cout << "setting model" << endl;
	this->options.model = new Model(path);
	//cout << "model set" << endl;

	double t = tmr.elapsed();
	cout << "Adding .bngl took " << t << " seconds" << endl;
}

void Swarm::setSwarmType(string type) {
	this->options.swarmType = type;

	if (options.swarmType == "swarm") {
		setParallelCount(options.swarmSize);
		cout << "setting parallelCount to swarmSize" << endl;
	}
}

void Swarm::setJobOutputDir(string path) {
	options.jobOutputDir = path;

	if (!checkIfFileExists(options.outputDir)) {
		string cmd = "mkdir " + options.outputDir;
		//cout << "running: " << cmd << endl;
		int ret = runCommand(cmd);
	}

	if (checkIfFileExists(options.jobOutputDir)) {
		string input;
		string answer;
		while (1) {
			cout << "Warning: Your output directory " << options.jobOutputDir << " already exists. Overwrite? (Y or N) ";
			getline(cin, input);
			stringstream myInp(input);
			myInp >> answer;

			if (answer == "Y" || answer == "y") {
				string cmd = "rm -r " + options.jobOutputDir + "*";
				//cout << "running: " << cmd << endl;
				int ret = runCommand(cmd);
				break;
			}
			else if (answer == "N" || answer == "n") {
				outputError("Error: Output directory already exists. Quitting.");
			}
		}
	}
	else {
		string cmd = "mkdir " + options.jobOutputDir;
		//cout << "else running: " << cmd << endl;
		//int ret = system(cmd.c_str());
		int ret = runCommand(cmd);
	}
}

void Swarm::addMutate(std::string mutateString) {
	vector<string> mutComponents;
	split(mutateString, mutComponents);

	// TODO: Should the parsing belong in the Config class? Maybe.
	// Make sure we have three components to work with
	if (mutComponents.size() == 3) {
		// Make sure first parameter name exists as a free parameter
		if (options.model->freeParams_.count(mutComponents[0]) > 0 || mutComponents[0] == "default") {
			// Make sure 2nd and 3rd components are numeric
			if (isFloat(mutComponents[1]) && isFloat(mutComponents[2])) {
				if (mutComponents[0] != "default") {
					cout << "not default" << endl;
					options.model->freeParams_.at(mutComponents[0])->setMutationRate(stof(mutComponents[1]));
					options.model->freeParams_.at(mutComponents[0])->setMutationFactor(stof(mutComponents[2]));
					options.model->freeParams_.at(mutComponents[0])->setHasMutation(true);

					//cout << "setting " << mutComponents[0] << " to " << mutComponents[1] << ":" << mutComponents[2] << endl;
				}
				else {
					for (map<string,FreeParam*>::iterator fp = options.model->freeParams_.begin(); fp != options.model->freeParams_.end(); ++fp) {
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

	//system ("exec rm -r /home/brandon/projects/GenFit2/Debug/output/*");
	system ("exec rm /home/brandon/projects/GenFit2/Debug/pOUT");

	initFit();

	// Main fit loops
	if (options.synchronicity == 1) {
		// Synchronous genetic
		if (options.swarmType == "genetic") {
			// TODO: Error checking. Make particles check and create next gen dir
			string createDirCmd = "mkdir " + options.jobOutputDir + "1";
			runCommand(createDirCmd);

			while (currentGeneration <= options.maxGenerations){
				runGeneration();

				string currentDirectory = options.jobOutputDir + "/" + to_string(static_cast<long long int>(currentGeneration));
				if (options.deleteOldFiles) {
					cleanupFiles(currentDirectory.c_str());
				}

				if (currentGeneration <= options.maxGenerations) {
					breedGeneration();
				}
			}
			cout << "fit finished " << endl;
			finishFit();
		}
	}
	// Asynchronous fit loops
	else {
		// Genetic fit
		if (options.swarmType == "genetic") {
			vector<int> finishedSimulations;
			int numFinishedSimulations = options.swarmSize;
			double totalFitTime;

			runGeneration();
			while(1) {
				finishedSimulations = checkMasterMessages();
				numFinishedSimulations += finishedSimulations.size();

				if (options.maxNumSimulations && numFinishedSimulations > options.maxNumSimulations) {
					break;
				}
				if (options.maxFitTime && totalFitTime > options.maxFitTime) {
					break;
				}
				if (options.minFit && allGenFits.begin()->first < options.maxFitTime) {
					break;
				}

				//for (auto pID: finishedSimulations) {
				for (auto pID = finishedSimulations.begin(); pID != finishedSimulations.end(); ++pID) {
					launchParticle(*pID);
				}
			}
			finishFit();
		}
		// PSO fit
		else if (options.swarmType == "pso") {
			// Generate all particles that will be present in the swarm
			allParticles_ = generateInitParticles();
			string createDirCmd = "mkdir " + options.jobOutputDir;

			bool stopCriteria = false;
			int clusterCheckCounter = 0;
			vector<int> finishedParticles;
			int numFinishedParticles;

			// Launch all particles
			for (int p = 1; p <= allParticles_.size(); ++p) {
				launchParticle(p);
			}

			// Wait until they are all finished so we can initialize
			// stuff for subsequent iterations
			while(numFinishedParticles < options.swarmSize) {
				// Check for our finished particles
				finishedParticles = checkMasterMessages();
				numFinishedParticles += finishedParticles.size();

				// Sleep for a second
				usleep(1000000);
			}

			// Fill a vector with all pID's
			vector<int> allParticles;
			for (int p = 1; p < options.swarmSize; ++p) {
				allParticles.push_back(p);
			}

			// Send all pID's in for processing (update velocities and positions)
			processParticlesPSO(allParticles, false);

			// Set our optimum
			optimum_ = particleBestFitsByFit_.begin()->first;

			// If we're using enhanced stop criteria, update the enhanced stop variables
			if (options.enhancedStop) {
				// Initialize our particle weights and weighted average position
				updateEnhancedStop();
			}

			// Re-launch the particles in to the Swarm proper
			for (int p = 1; p <= options.swarmSize; ++p) {
				launchParticle(p);
			}

			while (!stopCriteria) {
				usleep(250000);
				finishedParticles = checkMasterMessages();

				if (finishedParticles.size()) {
					// Process finished particles
					processParticlesPSO(finishedParticles, true);

					// Check for stop criteria
					stopCriteria = checkStopCriteria();
				}

				// Only check cluster queue every minute
				++clusterCheckCounter;
				if (clusterCheckCounter >= 240) {
					//checkClusterQueue();
					clusterCheckCounter = 0;
				}
			}

			// If we're out of the loop, we've finished the fit
			finishFit();
		}
	}
}

bool Swarm::checkStopCriteria() {
	if (options.enhancedStop) {

		// Eq 4 from Moraes et al
		// Algorithm 1 from Moraes et al
		if (abs(particleBestFitsByFit_.begin()->first - optimum_) < options.absTolerance + (options.absTolerance * particleBestFitsByFit_.begin()->first) ) {
			++permanenceCounter_;
			++inertiaUpdateCounter_;

			double quotient = permanenceCounter_ / options.nmin;
			if (floor(quotient) == quotient) {
				double oldAvgPos = weightedAvgPos_;
				updateEnhancedStop();

				if ( getEuclidianNorm(weightedAvgPos_ - oldAvgPos, options.model->getNumFreeParams()) < options.absTolerance ) {
					// Fit finished!
					return true;
				}
			}
		}
		else {
			permanenceCounter_ = 0;
			if (particleBestFitsByFit_.begin()->first < optimum_) {
				optimum_ = particleBestFitsByFit_.begin()->first;
			}
		}
	}

	// If we've run more than the maximum number of simulations
	if (flightCounter_ >= options.maxNumSimulations) {
		return true;
	}

	if (permanenceCounter_ >= options.maxNumSimulations) {
		return true;
	}

	/*
	// TIME
	if () {
	}
	*/
}

void Swarm::updateEnhancedStop() {
	// First calculate particle weights
	vector<float> particleWeights;

	updateParticleWeights();
	weightedAvgPos_ = calcWeightedAveragePosition();
}

double Swarm::getEuclidianNorm(double y, int n) {

	// Eq 10 in Moraes at al
	double sum;
	for (int i = 1; i <= n; ++i) {
		sum += pow(y, 2);
	}

	double norm = sqrt( (1/n) * sum);

	return norm;
}

void Swarm::updateParticleWeights() {
	for (int p = 1; p <= allParticles_.size(); ++p) {
		particleWeights_.at(p) = calcParticleWeight(p);
	}
}

double Swarm::calcWeightedAveragePosition() {

	double sum;
	for (int p = 1; p <= allParticles_.size(); ++p) {

		if (p == particleBestFitsByFit_.begin()->second) {
			continue;
		}

		sum += particleWeights_.at(p) * particleBestFits_[p];
	}

	return sum;
}

double Swarm::calcParticleWeight(int particle) {
	// Get reciprocal of euclidian norm of the difference between swarm best fit and particle best fit
	double numerator = 1 / getEuclidianNorm( (particleBestFits_.at(particle) - particleBestFitsByFit_.begin()->first), options.model->getNumFreeParams());

	// Eq 8 in Moraes et al
	double sum;
	for (int i = 1; i <= options.swarmSize; ++i) {
		// Make sure we're not using the particle with the best fit -- it will
		// result in a div_by_0
		if (i == particleBestFitsByFit_.begin()->second) {
			continue;
		}

		// Add 1/euclidian
		sum += 1 / numerator;
	}

	double weight = numerator / sum;

	return weight;
}

Particle * Swarm::createParticle(int pID) {
	Particle * p = new Particle(this, pID);
	//p->setBasePath(particleBasePath_);

	return p;
}

vector<vector<int>> Swarm::generateInitParticles(int pID) {
	Timer tmr;

	vector<vector<int>> allParticles (options.swarmSize+1);

	if (options.swarmType == "pso") {
		if (options.topology == "fullyconnected") {
			for (int p = 1; p <= options.swarmSize; ++p) {
				vector<int> connections;
				for (int c = 1; c <= options.swarmSize; ++c) {
					if (c != p) {
						connections.push_back(c);
					}
				}
				allParticles[p] = connections;
			}
		}
		else if (options.topology == "ring") {
			// Connect the first particle manually
			int firstConnection[] = {2, options.swarmSize};
			allParticles[1] = vector<int> (firstConnection, firstConnection + sizeof(firstConnection) / sizeof(firstConnection[0]));

			vector<int> connections;
			for (int p = 2; p <= options.swarmSize - 1; ++p) {
				connections.clear();
				// Connect to particle before, and particle after
				connections.push_back(p-1);
				connections.push_back(p+1);
				allParticles[p] = connections;
			}

			// Connect the last particle manually
			int lastConnection[] = {options.swarmSize - 1, 1};
			allParticles[options.swarmSize + 1] = vector<int> (lastConnection, lastConnection + sizeof(lastConnection) / sizeof(lastConnection[0]));
		}
		else if (options.topology == "star") {
			vector<int> connections;

			// First connect the central particle to all others
			for (int c = c; c <= options.swarmSize; ++c) {
				connections.push_back(c);
			}
			allParticles[1] = connections;

			// Then connect all particles to the central particles
			for (int p = 2; p <= options.swarmSize; ++p) {
				connections.clear();
				connections.push_back(1);
				allParticles[p] = connections;
			}
		}
		else if (options.topology == "mesh") {
			int desiredArea = options.swarmSize;
			int length;
			int width;

			float root = sqrt(desiredArea);

			// We have an integer. We will make a square.
			if (floor(root) == root) {
				length = root;
				width = root;
			}
			// sqrt is not an integer. Let's find the rectangle with the
			// least perimeter that will satisfy the area
			else {
				length = floor(root);
				width = ceil(root);
				float area = 0;

				while(desiredArea != area) {
					++length;
					--width;
					if (width == 0 || length == 0) {
						outputError("Error: Length or width is 0 when calculating mesh size.");
					}
					area = length*width;
				}
			}

			// Construct a matrix of dimenstion length x width
			// and fill it with particles
			int p = 0;
			vector<vector<int>> matrix(length, vector<int>(width));
			for (int x = 0; x < length; ++x) {
				for (int y = 0; y < width; ++y) {
					++p;
					matrix[x][y] = p;
				}
			}

			// Make our connections
			vector<int> connections;
			for (int x = 0; x < length; ++x) {
				for (int y = 0; y < width; ++y) {
					connections.clear();
					p = matrix[x][y];

					// If we're io the left bound
					if (x == 0) {
						// Add the particle to the right of us
						connections.push_back(matrix[x+1][y]);
					}
					// If we're on the right bound
					else if (x == length) {
						// Add the particle to the left of us
						connections.push_back(matrix[x-1][y]);
					}
					// If we're in a center
					else {
						// Add the particles on either side of us
						connections.push_back(matrix[x-1][y]);
						connections.push_back(matrix[x+1][y]);
					}

					if (y == 0) {
						connections.push_back(matrix[x][y+1]);
					}
					else if (y == width) {
						connections.push_back(matrix[x][y-1]);
					}
					else {
						connections.push_back(matrix[x][y-1]);
						connections.push_back(matrix[x][y+1]);
					}
					allParticles[p] = connections;
				}
			}
		}
		else if (options.topology == "toroidal") {
			int desiredArea = options.swarmSize;
			int length;
			int width;

			float root = sqrt(desiredArea);

			// We have an integer. We will make a square.
			if (floor(root) == root) {
				length = root;
				width = root;
			}
			// sqrt is not an integer. Let's find the rectangle with the
			// least perimeter that will satisfy the area
			else {
				length = floor(root);
				width = ceil(root);
				float area = 0;

				while(desiredArea != area) {
					++length;
					--width;
					if (width == 0 || length == 0) {
						outputError("Error: Length or width is 0 when calculating mesh size.");
					}
					area = length*width;
				}
			}

			// Construct a matrix of dimensions length x width
			// and fill it with particles
			int p = 0;
			vector<vector<int>> matrix(length, vector<int>(width));
			for (int x = 0; x < length; ++x) {
				for (int y = 0; y < width; ++y) {
					++p;
					matrix[x][y] = p;
				}
			}

			// Make our connections
			vector<int> connections;
			for (int x = 0; x < length; ++x) {
				for (int y = 0; y < width; ++y) {
					connections.clear();
					p = matrix[x][y];

					if (x == 0) {
						connections.push_back(matrix[x+1][y]);
						connections.push_back(matrix[length][y]);
					}
					else if (x == length) {
						connections.push_back(matrix[x-1][y]);
						connections.push_back(matrix[0][y]);
					}
					else {
						connections.push_back(matrix[x-1][y]);
						connections.push_back(matrix[x+1][y]);
					}

					if (y == 0) {
						connections.push_back(matrix[x][y+1]);
						connections.push_back(matrix[x][width]);
					}
					else if (y == width) {
						connections.push_back(matrix[x][y-1]);
						connections.push_back(matrix[x][0]);
					}
					else {
						connections.push_back(matrix[x][y-1]);
						connections.push_back(matrix[x][y+1]);
					}
					allParticles[p] = connections;
				}
			}
		}
		else if (options.topology == "tree") {
			int usedParticles = 1;
			int numLevels = 1;
			int previousLevel = 1;
			int currentLevel;

			// Determine number of levels in tree
			while (usedParticles < options.swarmSize) {
				currentLevel = previousLevel*2;
				usedParticles += currentLevel;
				previousLevel = currentLevel;
				++numLevels;
			}

			vector<vector<int>> tree(numLevels, vector<int>());
			usedParticles = 1;
			previousLevel = 1;
			int currentParticle = 2;
			bool doneFilling = false;

			// Construct tree and fill it with particles
			tree[0].push_back(1);

			// For each level in tree
			for (int level = 1; level <= numLevels; ++level) {
				// Current level's particle count is double that of previous level
				currentLevel = previousLevel * 2;

				// For each slot in current level
				for (int i = 0; i < currentLevel; ++i) {
					// Fill slot with particle
					tree[level].push_back(currentParticle);
					++currentParticle;

					// Make sure we're not filling past our swarm size
					if (currentParticle == options.swarmSize) {
						doneFilling = true;
						break;
					}
				}
				if (doneFilling) {
					break;
				}
			}

			// For each level in tree
			for (int level = 2; level < numLevels; ++level) {
				int prevGroupCounter = 0;
				int particleCounter = 0;

				// For each particle in level
				for (int p = 0; p <= tree[level].size(); ++level) {

					// Connect current particle (p) with proper particle in last level
					allParticles[p].push_back(tree[level-1][prevGroupCounter]);

					// Connect particle in last level to current particle (p)
					allParticles[tree[level-1][prevGroupCounter]].push_back(p);

					// If our particle counter reaches two, we're in a new pair in (or new particle
					// in previous level)
					if (++particleCounter == 2) {
						particleCounter = 0;
						++prevGroupCounter;
					}
				}
			}
		}
	}

	return allParticles;

	double t = tmr.elapsed();
	cout << "Particle creation took " << t << " seconds" << endl;
}

void Swarm::processParticlesPSO(vector<int> particles, bool newFlight) {
	// For each particle in our particle set
	for (auto particle = particles.begin(); particle != particles.end(); ++particle) {
		// We need to already have particle best position updated by the time we get here

		// Will hold positions for next iteration
		vector<double> nextPositions;

		// Calculate the next iteration positions according
		// to user preference
		if (options.psoType == "pso") {
			if (options.enhancedInertia) {
				updateInertia();
			}
			nextPositions = calcParticlePosPSO(*particle);
		}
		else if (options.psoType == "bbpso") {
			nextPositions = calcParticlePosBBPSO(*particle);
		}
		else if (options.psoType == "bbpsoexp") {
			nextPositions = calcParticlePosBBPSO(*particle, true);
		}

		// Convert positions to string so they can be sent to the particle
		vector<string> nextPositionsStr;
		for (auto param = nextPositions.begin(); param != nextPositions.end(); ++param) {
			nextPositionsStr.push_back(to_string(static_cast<long double>(*param)));
		}

		// Increment our flight counter
		++flightCounter_;

		// Finally, send the parameters
		swarmComm->sendToSwarm(0, *particle, SEND_FINAL_PARAMS_TO_PARTICLE, false, nextPositionsStr);
	}
}

vector<double> Swarm::calcParticlePosPSO(int particle) {
	int i = 0;

	// This vector holds the new positions to be sent to the particle
	vector<double> nextPositions(particleCurrParamSets_.size());

	// Get the best positions for particle's neighborhood
	vector<double> neighborhoodBestPositions = getNeighborhoodBestPositions(particle);

	// For each parameter in the current parameter set
	for (auto param = particleCurrParamSets_.begin(); param != particleCurrParamSets_.end(); ++param) {

		// Set up formula variables
		double currVelocity = options.inertia * particleParamVelocities_.at(particle)[i];
		double r1 = ((double) rand() / (RAND_MAX)); // TODO: These need to be inclusive
		double r2 = ((double) rand() / (RAND_MAX));
		double personalBestPos = particleBestParamSets_.at(particle)[i];
		double currPos = particleCurrParamSets_.at(particle)[i];

		// Set velocity
		double nextVelocity = (options.inertia * currVelocity) + options.cognitive * r1 * (personalBestPos - currPos) + options.social * r2 * (neighborhoodBestPositions[i] - currPos);

		// Set position
		nextPositions[i] = currPos + nextVelocity;
		++i;
	}

	return nextPositions;
}

vector<double> Swarm::calcParticlePosBBPSO(int particle, bool exp) {
	// Get the best positions for particle's neighborhood
	vector<double> neighborhoodBestPositions = getNeighborhoodBestPositions(particle);

	// This vector holds the new positions to be sent to the particle
	vector<double> nextPositions(particleCurrParamSets_.size());

	bool usePersonalBest = false;

	// For each parameter in the current parameter set
	int i = 0;
	for (auto param = particleCurrParamSets_.begin(); param != particleCurrParamSets_.end(); ++param) {
		if (exp) {
			// TODO: Does this need to be inclusive?
			if ( ((float) rand() / (RAND_MAX)) < 0.5 ) {
				usePersonalBest = true;
			}
		}

		double personalBestPos = particleBestParamSets_.at(particle)[i];

		if (usePersonalBest) {
			nextPositions[i] = personalBestPos;
			usePersonalBest = false;
		}
		else {
			// Calculate our mean and std
			double mean = (personalBestPos + neighborhoodBestPositions[i]) / 2;
			double std = personalBestPos + neighborhoodBestPositions[i];

			// Create the gaussian distribution
			std::tr1::normal_distribution<double> dist(mean, std);

			// Pick our next position randomly from distribution
			nextPositions[i] = dist(randNumEngine);
		}
		++i;
	}

	return nextPositions;
}

void Swarm::updateInertia() {
	// Eq 3 from Moraes et al
	options.inertia = options.inertiaInit + (options.inertiaFinal - options.inertiaInit) * (inertiaUpdateCounter_ / options.nmax + inertiaUpdateCounter_);
}

vector<double> Swarm::getNeighborhoodBestPositions(int particle) {

	// Set the current best fit to our own best fit
	double currBestFit = particleBestFits_.at(particle);
	int currentBestNeighbor = particle;

	// For every neighbor in this particle's neighborhood
	for (auto neighbor = allParticles_[particle].begin(); neighbor != allParticles_[particle].end(); ++neighbor) {
		// Set best fit of neighbor being tested
		double neighborBestFit = particleBestFits_.at(*neighbor);

		// If Neighbor's best fit is better than ours, update the best
		// neighbor. Also, update the current best fit value
		if (neighborBestFit < currBestFit) {
			currentBestNeighbor = *neighbor;
			currBestFit = neighborBestFit;
		}
	}

	return particleBestParamSets_.at(currentBestNeighbor);
}

void Swarm::processParamsPSO(vector<double> &params, int pID, double fit) {

	// First update the particles current parameter set
	int i = 0;
	for (auto param = params.begin(); param != params.end(); ++param) {
		particleCurrParamSets_.at(pID)[i] = *param;
		++i;
	}

	// The the fit value of this param set is less than the particles best fit
	// we should update the particle's best fit, then store the best fit params
	if (fit < particleBestFits_.at(pID)) {
		particleBestFits_[pID] = fit;

		i = 0;
		for (auto param = params.begin(); param != params.end(); ++param) {
			particleBestParamSets_.at(pID)[i] = *param;
			++i;
		}
	}
}

void Swarm::launchParticle(int pID) {

	if (currentGeneration == 1 && !options.useCluster) {

		// Construct command needed to run the particle
		string command = exePath_ + " particle " + to_string(static_cast<long long int>(pID)) + " run " + to_string(static_cast<long long int>(currentGeneration)) + " " + sConf_;
		command = command + " >> pOUT 2>&1";
		command = command + " &";

		// TODO: Check system return value for success
		int ret = runCommand(command);

		cout << "Command: " << command << endl;
	}

	if (options.verbosity >= 3) {
		cout << "Running Particle " << pID << endl;
	}

	runningParticles_.insert(pID);
	swarmComm->sendToSwarm(0, pID, NEXT_GENERATION, false, swarmComm->univMessageSender);
}

void Swarm::setCurrentGen(int gen) {
	currentGeneration = gen;
}

void Swarm::runGeneration () {
	// TODO: Implement walltime
	if(options.verbosity >= 1) {
		cout << "Running generation " << currentGeneration << " with " << options.swarmSize << " particles..." << endl;
	}

	int numLaunchedParticles = 0;
	int numFinishedParticles = 0;

	vector<int> finishedParticles;
	int p = 1;

	while (numFinishedParticles < options.swarmSize) {

		// If we're running on a cluster and submitting all particles at once
		/*if (options.useCluster && options.parallelCount == options.swarmSize) {
			cout << "3" << endl;
			while (p != allParticles_.end()) {
				launchParticle(p->first);
				numLaunchedParticles += 1;
				++p;
			}
		}
		else {*/

		if (runningParticles_.size() < options.parallelCount && numLaunchedParticles < options.swarmSize) {
			launchParticle(p);
			numLaunchedParticles += 1;
			++p;
		}
		//}

		// Check for any messages from particles
		usleep(10000);
		finishedParticles = checkMasterMessages();
		numFinishedParticles += finishedParticles.size();
		finishedParticles.clear();
	}
	if (failedParticles_.size() > (options.swarmSize - 3) ) {
		outputError("Error: You had too many failed runs. Check simulation output or adjust walltime.");
		//TODO: Cleanup and exit
	}

	currentGeneration += 1;
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
			if (boost::regex_match(string( de->d_name ),boost::regex(".+xml$|.+species$|.+cdat$|.+gdat$|.+net$|.+BNG_OUT"))) {
				//cout << "found a file";
				filesToDel.push_back( string(de->d_name) );
			}
		}
		closedir( dp );
	}
	else {
		cout << "Warning: Couldn't open " << path << " to delete un-needed simulation files." << endl;
	}

	//for (auto i: filesToDel) {
	for (auto i = filesToDel.begin(); i != filesToDel.end(); ++i) {
		string fullPath = string(path) + "/" + *i;
		remove(fullPath.c_str());
	}
}

bool Swarm::sortFits(Particle * a, Particle * b) {
	return a->fitCalcs.at(currentGeneration) < b->fitCalcs.at(currentGeneration);
}

void Swarm::breedGeneration() {
	//Timer tmr;

	// TODO: Check ret
	string createDirCmd = "mkdir " + options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration));
	int ret = runCommand(createDirCmd);

	// We let particles do the actual breeding.  The master's role is to generate a breeding pattern and
	// tell which particles to breed with which
	/*
	cout << "fit\tperm\t";

	for (auto i: options.model->freeParams_) {
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
	int parentPairs = options.swarmSize / 2;

	//cout << "we have " << parentPairs << " parent pairs" << endl;

	double p1;
	double p2;
	boost::smatch match;
	vector<string> parentVec;
	int numCurrBreeding = 0;
	int numFinishedBreeding = 0;

	//cout << "sum: " << weightSum << endl;

	int swapID = -1;
	for (int i = 0; i < parentPairs; ++i) {
		swapID += 2;

		// Pick the fit values (particle parents) used in breeding
		p1 = pickWeighted(weightSum, weights, options.extraWeight, randNumEngine);
		p2 = pickWeighted(weightSum, weights, options.extraWeight, randNumEngine);

		// If we want different parents used in breeding, make sure that happens
		int retryCount = 0;
		while (p1 == p2 && options.forceDifferentParents) {
			retryCount++;
			if (retryCount > options.maxRetryDifferentParents) {
				if (options.verbosity >= 3) {
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
			p2 = pickWeighted(weightSum, weights, options.extraWeight, randNumEngine);
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
		boost::regex_search(parentVec[0], match, boost::regex("gen\\d+perm(\\d+)"));
		string pID1 = match[1];

		// Clear the parent vector for use with the next parent
		parentVec.clear();

		// Do same as above, but for the second parent
		split(parentString2,parentVec);
		boost::regex_search(parentVec[0], match, boost::regex("gen\\d+perm(\\d+)"));
		string pID2 = match[1];
		parentVec.clear();

		//cout << "p1 is " << pID1 << " p2 is " << pID2 << endl;

		// Add second parent and swapID to the message
		swarmComm->univMessageSender.push_back(to_string(static_cast<long long int>(swapID)));
		swarmComm->univMessageSender.push_back(pID2);

		// Send the message to the first parent
		//cout << "sending message to " << stoi(pID1) << " with id of " << to_string(static_cast<long long int>(swapID)) << endl;

		// Tell parent 1 to initiate breeding
		swarmComm->sendToSwarm(0, stoi(pID1), INIT_BREEDING, false, swarmComm->univMessageSender);
		swarmComm->univMessageSender.clear();
		numCurrBreeding+=2;

		//double t = tmr.elapsed();
		//cout << "BREEDING took " << t << " seconds" << endl;
		//tmr.reset();

		// While the current number of breeding particles is greater than twice the parallel count
		while ( numCurrBreeding >= (options.parallelCount * 15)) {

			int numMessages = swarmComm->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm->univMessageReceiver, true);

			numFinishedBreeding+=numMessages;
			numCurrBreeding-=numMessages;
		}
	}
	while (numFinishedBreeding < options.swarmSize) {

		int numMessages = swarmComm->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm->univMessageReceiver, true);

		numFinishedBreeding+=numMessages;
	}
	//cout << "got them all!" << endl;
}

void Swarm::finishFit() {
	string outputDir = options.jobOutputDir + "/Results";
	string command = "mkdir " + outputDir;

	if (runCommand(command)) {
		string errMsg = "Error: Couldn't create directory: " + outputDir + " to contain final fitting results.";
		outputError(errMsg);
	}

	outputRunSummary(outputDir);
	//generateBestFitModels(outputDir);
	//copyBestFitToResults(outputDir);

	killAllParticles(FIT_FINISHED);

	cout << "Finished fitting in " << tmr_.elapsed() << " seconds. Results can be found in " << options.jobOutputDir << "Results" << endl;
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
		//for (auto i: options.model->freeParams_) {
		for (auto i = options.model->freeParams_.begin(); i != options.model->freeParams_.end(); ++i) {
			outputFile << left << setw(16) << i->first;
		}

		outputFile << endl;

		vector<string> paramVals;

		//for (auto i: allGenFits) {
		for (auto i = allGenFits.begin(); i != allGenFits.end(); ++i) {
			split(i->second, paramVals);

			outputFile << left << setw(16) << i->first << left << setw(16) << paramVals[0];

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
	for (int p = 1; p <= options.swarmSize; ++p) {
		swarmComm->sendToSwarm(0, p, tag, false, swarmComm->univMessageSender);
	}
}

void Swarm::getClusterInformation() {

	// If user didn't specify cluster platform, let's figure it out ourself
	if (options.clusterSoftware.size() == 0) {
		// Test for slurm
		string output;
		runCommand("which srun", output);
		if (output.length() > 0) {
			options.clusterSoftware = "slurm";
		}
		// Test for PBS-type

		else{
			runCommand("which qsub", output);
			if(output.length() > 0) {
				// Test for Torque/PBS
				string output2;
				runCommand("which maui", output);
				runCommand("which moab", output2);
				if(output.length() > 0 || output2.length() > 0) {
					options.clusterSoftware = "torque";
				}
				// Test for SGE
				else {
					string output3;
					runCommand("which sge_execd", output);
					runCommand("which qconf", output);
					runCommand("which qmon", output);
					if (output.length() > 0 || output2.length() > 0 || output3.length() > 0) {
						outputError("Error: BioNetFit doesn't support for GridEngine clusters. If you are not running on a GridEngine cluster, specify the cluster platform in the .conf file using the 'cluster_software' option.");
					}
				}
			}
		}
	}
	else {
		if (options.clusterSoftware != "slurm" && options.clusterSoftware != "torque") {
			outputError("You specified an unrecognized cluster type in your .conf file. BioNetFit only supports 'torque' or 'slurm' cluster types.");
		}
	}

	// If we still don't know the cluster type, let's ask the user.
	if (options.clusterSoftware.size() == 0) {
		string input;
		string clusterType;

		while (1) {
			cout << "BioNetFit couldn't determine which type of cluster software you are using. Specify (T) for Torque, or (S) for Slurm" << endl;
			getline(cin, input);
			stringstream myInp(input);
			myInp >> clusterType;

			if (clusterType == "S" || clusterType == "s") {
				options.clusterSoftware = "slurm";
				break;
			}
			else if (clusterType == "T" || clusterType == "t") {
				options.clusterSoftware = "torque";
				break;
			}
		}
	}
}

string Swarm::generateSlurmCommand(string cmd, bool multiProg) {
	string command;

	//TODO: Need to generate job submission name
	//TODO: Need to display terminal output from cluster jobs

	// srun submits the job to the cluster
	command += "srun";

	// Add the job name
	command += " -J " + options.jobName;

	// Specify the cluster account if needed
	if (!options.clusterAccount.empty()) {
		command += " -A " + options.clusterAccount;
	}

	if (!options.clusterQueue.empty()) {
		command += " -p	" + options.clusterQueue;
	}

	// Specify output directory if needed
	if (options.saveClusterOutput) {
		command += " -o " + options.outputDir + "/" + options.jobName + "_cluster_output/" + options.jobName;
	}
	else {
		command += " -o /dev/null";
	}

	cout << "pc is: " << options.parallelCount;
	command += " -c" + to_string(static_cast<long long int>(options.parallelCount));
	command += " -l";


	if (multiProg) {
		command += " --multi-prog " + cmd;
	}
	else {
		command+= " " + cmd;
	}

	return command;
}

void Swarm::initFit () {
	// If we are running ODE, we want to generate a network and network reader.
	// The network file has parameters that are replaced in every generation,
	// and the network reader is used to read those network files
	if (options.model->getHasGenerateNetwork()) {

		// First create a particle which will generate the network
		// and fill it with dummy params
		Particle p = Particle(this, 1);
		p.setModel(options.model);
		p.generateParams();

		// Output model file with dummy parameters. This file will be used to generate
		// our initial network
		options.model->outputModelWithParams(p.getParams(), options.jobOutputDir, "base.bngl", "", false, false, false, false, false);

		// Construct our simulation command and run the network generator
		string modelPath = options.jobOutputDir + "/base.bngl";
		string command = options.simPath + "BNG2.pl --outdir " + options.jobOutputDir + " " + modelPath + " >> " + options.jobOutputDir + "netgen_output 2>&1";
		cout << "Generating initial .net file with command: " << command << endl;

		// TODO: Check this return
		int ret = runCommand(command);

		// Now that we have a .net file, replace the full .bngl with a .bngl
		// containing ONLY action commands and a .net file loader. This is
		// used later to run our .net files.
		options.model->outputModelWithParams(p.getParams(), options.jobOutputDir, "base.bngl", "", false, true, false, false, false);

		// Now store our .net file for later use
		string netPath = options.jobOutputDir + "/base.net";
		options.model->parseNet(netPath);
	}
}

vector<int> Swarm::checkMasterMessages() {
	vector<int> finishedParticles;
	//int numFinishedParticles = 0;
	//cout << "checking messages" << endl;
	int numMessages = swarmComm->recvMessage(-1, 0, -1, false, swarmComm->univMessageReceiver);

	if (numMessages >= 1) {
		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		smhRange = swarmComm->univMessageReceiver.equal_range(SIMULATION_END);
		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
			int pID = sm->second.sender;

			if (options.verbosity >= 3) {
				cout << "Particle " << pID << " finished simulation" << endl;
			}

			// Get an iterator to the particle in our list of running particles
			runningParticlesIterator_ = runningParticles_.find(pID);
			// Then remove it
			if (runningParticlesIterator_ == runningParticles_.end()) {
				outputError("Error: Couldn't remove particle from runningParticle list.");
			}
			runningParticles_.erase(runningParticlesIterator_);
			// Increment our finished counter
			//numFinishedParticles += 1;


			finishedParticles.push_back(pID);
			//cout << "pushed pid" << endl;

			string paramsString = "gen" + to_string(static_cast<long long int>(currentGeneration)) + "perm" + to_string(static_cast<long long int>(pID)) + " ";
			vector<double> paramsVec;

			double fitCalc = stod(sm->second.message[0]);
			//cout << "stored calc" << endl;

			// Store the parameters given to us by the particle
			for (vector<string>::iterator m = sm->second.message.begin()+1; m != sm->second.message.end(); ++m) {
				paramsString += *m + " ";
			}
			//cout << "stored params" << endl;

			// Then store it
			//cout << "saving fit for " << pID << " of " << fitCalc << endl;
			allGenFits.insert(pair<double,string>(fitCalc,paramsString));
			//cout << "saved fit" << endl;

			if (options.swarmType == "pso") {
				processParamsPSO(paramsVec, pID, fitCalc);
			}
		}

		// TODO: When sending NEXT_GENERATION, make sure failed particles have actually run again. If not, they need re-launched.
		smhRange = swarmComm->univMessageReceiver.equal_range(SIMULATION_FAIL);
		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
			int pID = sm->second.sender;

			// Get an iterator to the particle in our list of running particles
			runningParticlesIterator_ = runningParticles_.find(pID);
			// Then remove it
			runningParticles_.erase(runningParticlesIterator_);
			// Increment our finished counter
			//numFinishedParticles += 1;
			finishedParticles.push_back(pID);

			cout << "Particle " << pID << " failed in gen " << currentGeneration << endl;

			// Store particle ID in our list of failed particles
			failedParticles_.insert(pID);
		}

		smhRange = swarmComm->univMessageReceiver.equal_range(GET_RUNNING_PARTICLES);
		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
			vector<string> intParts;
			//for (auto i: runningParticles_) {
			for (auto i = runningParticles_.begin(); i != runningParticles_.end(); ++i) {
				intParts.push_back(to_string(static_cast<long long int>(*i)));
			}
			swarmComm->sendToSwarm(0, sm->second.sender, SEND_RUNNING_PARTICLES, true, intParts);
		}

		swarmComm->univMessageReceiver.clear();
		numMessages = 0;
	}

	return finishedParticles;
}

string Swarm::generateSlurmMultiProgCmd(string runCmd, string serializedSwarmPath) {
	string multiProgConfPath = options.jobOutputDir + "multiprog.conf";
	ofstream multiprog(multiProgConfPath, ios::out);

	if (multiprog.is_open()) {
		multiprog << "0 " << runCmd << " " << configPath_ << endl;
		for (int id = 1; id <= options.swarmSize; ++id) {
			multiprog << to_string(static_cast<long long int>(id)) + " " << runCmd << " particle " << to_string(static_cast<long long int>(id)) << " run " << serializedSwarmPath << endl;
		}
		multiprog.close();
	}


	return generateSlurmCommand(multiProgConfPath);
}

string Swarm::generateSlurmBatchFile(string runCmd) {
	string sbatchPath = options.jobOutputDir + "slurm_script.sh";
	ofstream sbatch(sbatchPath, ios::out);

	if (sbatch.is_open()) {
		//TODO: Need to generate job submission name
		//TODO: Need to display terminal output from cluster jobs

		sbatch << "#!/bin/sh" << endl << endl;

		// Add the job name
		sbatch << "#SBATCH -J " + options.jobName << endl;

		// Specify the cluster account if needed
		if (!options.clusterAccount.empty()) {
			sbatch << "#SBATCH -A " + options.clusterAccount << endl;
		}

		if (!options.clusterQueue.empty()) {
			sbatch << "#SBATCH -p " + options.clusterQueue << endl;
		}

		// Specify output directory if needed
		if (options.saveClusterOutput) {
			sbatch << "#SBATCH -o " + options.outputDir + "/" + options.jobName + "_cluster_output/" + options.jobName << endl;
		}
		else {
			sbatch << "#SBATCH -o /dev/null" << endl;
		}

		sbatch << "#SBATCH -n" + to_string(static_cast<long long int>(options.parallelCount)) << endl;

		sbatch << endl;

		sbatch << "module load openmpi" << endl;
		//sbatch << "module load intel << endl;"
		sbatch << "echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH" << endl;
		sbatch << "echo PATH: $PATH" << endl;
		sbatch << "echo which mpirun:" << endl;
		sbatch << "which mpirun" << endl;
		sbatch << "echo pwd:" << endl;
		sbatch << "pwd"	<< endl;
		sbatch << "echo ldd GenFit2:" << endl;
		sbatch << "ldd /home/bt285/BioNetFit/bin/GenFit2" << endl;;
		sbatch << "mpirun -mca mca_component_show_load_errors 10 -v --tag-output -np 1 " << runCmd << " load " << sConf_ << " : -np " << options.swarmSize << " " << runCmd << " particle 0 run " << sConf_ << endl;
		//sbatch << "mpirun -prepend-rank -np 1 " << runCmd << " load " << sConf_ << " : -np " << options.swarmSize << " " << runCmd << " particle 0 run " << sConf_ << endl;

		sbatch.close();
	}
	//TODO: Error checking if file couldn't be opened

	runCmd = "sbatch " + sbatchPath;

	return runCmd;
}
