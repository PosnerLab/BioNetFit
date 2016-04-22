/*
 * Swarm.cpp
 *
 *  Created on: Jul 17, 2015
 *      Author: brandon
 */
// TODO: Auto average best fit gdat for results?

#include "Swarm.hh"

using namespace std;
using namespace std::chrono;

Swarm::Swarm() {

	// Whether or not we are master
	isMaster = false;
	resumingSavedSwarm = false;
	sConf_ = "";
	swarmComm = 0;
	fitCompareTolerance = 1e-6;
	bootstrapCounter = 0;

	options.jobName = "";
	options.fitType = "";
	options.outputDir = "";
	options.bngCommand = "";
	options.outputEvery = 100;

	options.model = 0;

	options.synchronicity = 0;
	options.maxGenerations = 0;
	options.swarmSize = 0;
	options.minFit = 0;
	options.maxFit = 0;
	options.bootstrap = 0;
	options.parallelCount = 0;
	options.seed = 0;

	options.useCluster = false;
	options.saveClusterOutput = false;
	options.emailWhenFinished = false;
	options.emailAddress = "";

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
	options.smoothing = 1;
	options.keepParents = 0;

	options.maxFitTime = MAX_LONG;
	options.maxNumSimulations = MAX_LONG;

	// PSO options
	options.inertia = 0.72; // 0.72
	options.cognitive = 1.49; // 1.49
	options.social = 1.49; // 1.49
	options.nmax = 0; // 20
	options.nmin = 80; // 80
	options.inertiaInit = 1; // 1
	options.inertiaFinal = 0.1; // 0.1
	options.absTolerance = 10e-4; // 10E-4
	options.relTolerance = 10e-4; // 10E-4
	options.mutateQPSO = false;
	options.betaMin = 0.5;
	options.betaMax = 1.0;

	options.topology = "fullyconnected"; // fullyconnected
	options.psoType = "pso"; // pso

	options.enhancedStop = false; // true
	options.enhancedInertia = false; // true

	options.minTemp = pow(10, -10);
	options.minRadius = pow(10, -6);
	options.localSearchProbability = 0.01;
	options.randParamsProbability = 0.1;

	options.verbosity = 1;
	hasMutate = false;

	currentGeneration = 0;

	/*
	std::map<int, double> particleBestFits_;
	std::map<int, std::vector<double>> particleBestParamSets_;
	std::map<int, std::vector<double>> particleCurrParamSets_;
	std::map<int, double> particleWeights_;
	std::map<double, int> particleBestFitsByFit_;
	 */

	permanenceCounter_ = 0; // 0
	flightCounter_ = 0; // 0
	weightedAvgPos_ = 0; // 0
	optimum_ = 0; // 0
	inertiaUpdateCounter_ = 0; // 0;
	beta_ = 0.7;
	cauchyMutator_ = 0.2;

	auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	randNumEngine.seed(seed);
	randNumEngine.discard(700000);

	/*
	if (options.seed) {
		randNumEngine.seed(options.seed);
		randNumEngine.discard(700000);

		srand(options.seed);
	}
	else {
		// TODO: Make sure everything is being seeded properly and in the proper place. Also let's do away with rand()
		// Seed our random number engine
		auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		randNumEngine.seed(seed);
		randNumEngine.discard(700000);

		srand (std::tr1::random_device{}());
	}
	 */
}

void Swarm::initRNGS(int seed) {
	if (seed) {
		srand(seed);
		//cout << "seeding rng with: " << seed << endl;
	}
	else {
		// TODO: Make sure everything is being seeded properly and in the proper place. Also let's do away with rand()
		// Seed our random number engine
		srand (std::tr1::random_device{}());
	}
}

void Swarm::initComm() {
	// Create the communication class
	Pheromones *ph = new Pheromones();

	// Initialize the communication class
	ph->init(this);

	// Store our communication class
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
	this->options.model = new Model(this, path);
	//cout << "model set" << endl;

	double t = tmr.elapsed();
	cout << "Adding .bngl took " << t << " seconds" << endl;
}

void Swarm::setfitType(string type) {
	this->options.fitType = type;

	if (options.fitType == "swarm") {
		options.parallelCount = options.swarmSize;
	}
}

void Swarm::setJobOutputDir(string path) {
	options.jobOutputDir = path;

	if (!checkIfFileExists(options.outputDir)) {
		string cmd = "mkdir " + options.outputDir;
		//cout << "running: " << cmd << endl;
		if(runCommand(cmd) != 0) {
			outputError("Error: Couldn't create base output directory with command: " + cmd + ". Quitting.");
		}
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
				string cmd = "rm -r " + options.jobOutputDir;
				if (options.verbosity >= 1) {
					cout << "Deleting old output directory, this may take a few minutes..." << endl;
				}

				if(runCommand(cmd) != 0) {
					outputError("Error: Couldn't delete existing output directory with the command: " + cmd + ". Quitting.");
				}
				break;
			}
			else if (answer == "N" || answer == "n") {
				outputError("Error: Output directory already exists. Quitting.");
			}
		}
	}
	string cmd = "mkdir " + options.jobOutputDir;
	//cout << "else running: " << cmd << endl;
	//int ret = system(cmd.c_str());
	if(runCommand(cmd) != 0) {
		outputError("Error: Couldn't create output directory with command: " + cmd + ". Quitting.");
	}
}

void Swarm::addMutate(std::string mutateString) {
	vector<string> mutComponents;
	split(mutateString, mutComponents);

	//cout << "mut" << endl;
	// TODO: Should the parsing belong in the Config class? Maybe.
	// Make sure we have three components to work with
	if (mutComponents.size() == 3) {
		// Make sure first parameter name exists as a free parameter
		if (options.model->freeParams_.count(mutComponents[0]) > 0 || mutComponents[0] == "default") {
			// Make sure 2nd and 3rd components are numeric
			if (isFloat(mutComponents[1]) && isFloat(mutComponents[2])) {
				if (mutComponents[0] != "default") {
					//cout << "not default" << endl;
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

	system ("exec rm /home/brandon/projects/BioNetFit/Debug/pOUT");

	for (bootstrapCounter = 0; bootstrapCounter <= options.bootstrap; ++bootstrapCounter) {
		initFit();
		// TODO: Test all synchronous fits on PC
		// Main fit loops
		if (options.synchronicity) {
			// Synchronous genetic
			if (options.fitType == "genetic") {
				runSGA();
			}
			else if (options.fitType == "pso") {
				runSPSO();
			}
			else if (options.fitType == "de") {
				runSDE();
			}
			else if (options.fitType == "sa") {
				runSSA();
			}
		}

		// TODO: Need to make asynchronous fits work in PCs
		// Asynchronous fit loops
		else {
			// Genetic fit
			if (options.fitType == "genetic") {
				runAGA();
			}
			// PSO fit
			else if (options.fitType == "pso") {
				runAPSO();
			}
			else if (options.fitType == "de") {
				runADE();
			}
			else if (options.fitType == "sa") {
				runASA();
			}
		}
		finishFit();
	}
}

bool Swarm::checkStopCriteria() {
	if (options.verbosity >= 3) {
		//cout << "Checking stop criteria. Flight count is " << flightCounter_ << " and max is " << options.maxNumSimulations << endl;
		//cout << "Fittype is " << options.fitType << endl;
	}

	if (options.fitType == "pso") {
		if (options.enhancedStop) {
			// Eq 4 from Moraes et al
			// Algorithm 1 from Moraes et al
			if (abs(particleBestFitsByFit_.begin()->first - optimum_) < options.absTolerance + (options.relTolerance * particleBestFitsByFit_.begin()->first) ) {
				if (options.verbosity >= 3) {
					cout << "Incrementing permanence counter" << endl;
				}

				++permanenceCounter_;
				++inertiaUpdateCounter_;

				double quotient = ((double)permanenceCounter_ / options.nmin);
				// Check to see if quotient is a real number
				if (floor(quotient) == quotient) {
					double oldAvgPos = weightedAvgPos_;
					updateEnhancedStop();

					if ( getEuclidianNorm(weightedAvgPos_ - oldAvgPos, options.model->getNumFreeParams()) < options.absTolerance ) {
						// Fit finished!

						if (options.verbosity >= 3) {
							cout << "Stopped according to enhanced stop criteria" << endl;
						}
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
		else if (options.nmax) {
			if (abs(particleBestFitsByFit_.begin()->first - optimum_) < options.absTolerance + (options.relTolerance * particleBestFitsByFit_.begin()->first) ) {
				if (options.verbosity >= 3) {
					cout << "Incrementing permanence counter. Optimum is " << optimum_ << " and f(g) is " << particleBestFitsByFit_.begin()->first << ". Difference is " << abs(particleBestFitsByFit_.begin()->first - optimum_) << " and conditional is " << options.absTolerance + (options.relTolerance * particleBestFitsByFit_.begin()->first) << endl;
				}

				++permanenceCounter_;
				++inertiaUpdateCounter_;

				if (permanenceCounter_ >= options.nmax) {

					if (options.verbosity >= 3) {
						cout << "Stopped according to permanence > n_max" << endl;
					}
					return true;
				}
			}
			else {
				permanenceCounter_ = 0;
				if (particleBestFitsByFit_.begin()->first < optimum_) {
					optimum_ = particleBestFitsByFit_.begin()->first;
				}
			}
		}
	}
	else if ((options.fitType == "genetic" || options.fitType == "de") && options.synchronicity) {
		if (options.verbosity >= 3) {
			cout << "Checking if we've reached max generation. Current is " << currentGeneration << " and max is " << options.maxGenerations << endl;
		}
		if (options.maxGenerations && currentGeneration > options.maxGenerations) {
			return true;
		}
	}

	// If we've run more than the maximum number of simulations
	if (flightCounter_ >= options.maxNumSimulations) {
		if (options.verbosity >= 3) {
			cout << "Stopped according to flightCounter (" << flightCounter_ << ") >= maxNumSimulations (" << options.maxNumSimulations << ")" << endl;
		}
		return true;
	}

	//cout << "Fitcompare" << endl;
	//cout << particleBestFitsByFit_.begin()->first << endl;
	//cout << options.minFit << endl;
	if (swarmBestFits_.begin()->first <= options.minFit) {
		cout << "Stopped according to swarmBestFit (" << particleBestFitsByFit_.begin()->first << ") <= options.minFit (" << options.minFit << ")" << endl;
		return true;
	}

	// TODO: Check if all particles have converged to same solution

	/*
			// TIME
			if () {
			}
	 */

	return false;
}

void Swarm::updateEnhancedStop() {
	// First calculate particle weights

	if (options.verbosity >= 3) {
		cout << "Updating enhanced stop" << endl;
	}

	updateParticleWeights();
	weightedAvgPos_ = calcWeightedAveragePosition();
}

double Swarm::getEuclidianNorm(double y, unsigned int n) {

	// Eq 10 in Moraes at al
	double sum;
	for (unsigned int i = 1; i <= n; ++i) {
		sum += pow(y, 2);
	}

	double norm = sqrt( (1/n) * sum);

	return norm;
}

void Swarm::updateParticleWeights() {

	if (options.verbosity >= 3) {
		cout << "Updating particle weights" << endl;
	}

	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		particleWeights_[p] = calcParticleWeight(p);
	}

	cout << "done here" << endl;
}

double Swarm::calcWeightedAveragePosition() {

	if (options.verbosity >= 3) {
		cout << "Calculating weighted average position" << endl;
	}

	double sum = 0;
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {

		if (p == particleBestFitsByFit_.begin()->second) {
			continue;
		}

		sum += particleWeights_.at(p) * particleBestFits_[p];
	}

	return sum;
}

vector<double> Swarm::calcQPSOmBests() {
	// Eq 2b Liu et al
	vector<double> mBests;

	for (unsigned int param = 0; param < options.model->getFreeParams_().size(); ++param) {
		double sum = 0;
		for (unsigned int particle = 1; particle <= options.swarmSize; ++particle) {
			sum += abs(particleBestParamSets_.at(particle)[param]);
		}

		double mean = sum / (double)options.swarmSize;

		// TODO: Adaptive mutation ala Liu 2005
		if (options.mutateQPSO) {
			boost::random::cauchy_distribution<double> dist(0, cauchyMutator_);
			double mutator = dist(randNumEngine);
			mean += mutator;
			//cout << mean << " mutated by: " << mutator << endl;
		}
		mBests.push_back(sum / (double)options.swarmSize);
	}

	return mBests;
}

double Swarm::calcParticleWeight(unsigned int particle) {
	// Get reciprocal of euclidian norm of the difference between swarm best fit and particle best fit
	double numerator = 1 / getEuclidianNorm( (particleBestFits_[particle] - particleBestFitsByFit_.begin()->first), options.model->getNumFreeParams());

	// Eq 8 in Moraes et al
	double sum;
	for (unsigned int i = 1; i <= options.swarmSize; ++i) {
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

Particle * Swarm::createParticle(unsigned int pID) {
	Particle * p = new Particle(this, pID);
	//p->setBasePath(particleBasePath_);

	return p;
}

vector<vector<unsigned int>> Swarm::generateTopology(unsigned int populationSize) {
	Timer tmr;

	vector<vector<unsigned int>> allParticles (populationSize + 1);

	if (options.fitType == "pso" || options.fitType == "de") {

		if (options.verbosity >= 3) {
			cout << "Generating initial particles with a " << options.topology << " topology" << endl;
		}

		if (options.topology == "fullyconnected") {
			for (unsigned int p = 1; p <= populationSize; ++p) {
				vector<unsigned int> connections;
				for (unsigned int c = 1; c <= populationSize; ++c) {
					if (c != p) {
						connections.push_back(c);
					}
				}
				allParticles[p] = connections;
			}
		}
		else if (options.topology == "ring") {

			// Connect the first particle manually
			unsigned int firstConnection[] = {2, populationSize};
			allParticles[1] = vector<unsigned int> (firstConnection, firstConnection + sizeof(firstConnection) / sizeof(firstConnection[0]));

			vector<unsigned int> connections;
			for (unsigned int p = 2; p <= populationSize - 1; ++p) {

				connections.clear();
				// Connect to particle before, and particle after
				connections.push_back(p-1);
				connections.push_back(p+1);
				allParticles[p] = connections;

			}

			// Connect the last particle manually
			unsigned int lastConnection[] = {populationSize - 1, 1};
			allParticles[populationSize] = vector<unsigned int> (lastConnection, lastConnection + sizeof(lastConnection) / sizeof(lastConnection[0]));
		}
		else if (options.topology == "star") {
			vector<unsigned int> connections;

			// First connect the central particle to all others
			for (unsigned int c = 2; c <= populationSize; ++c) {
				connections.push_back(c);
			}
			allParticles[1] = connections;

			// Then connect all particles to the central particles
			for (unsigned int p = 2; p <= populationSize; ++p) {
				connections.clear();
				connections.push_back(1);
				allParticles[p] = connections;
			}
		}
		else if (options.topology == "mesh") {

			unsigned int desiredArea = populationSize;
			unsigned int divisor = ceil(sqrt(desiredArea));
			while(desiredArea % divisor != 0) {
				++divisor;
			}
			unsigned int length = divisor;
			unsigned int width = desiredArea / divisor;

			// Construct a matrix of dimension length x width
			// and fill it with particles
			int p = 0;
			vector<vector<unsigned int>> matrix(length, vector<unsigned int>(width));
			for (unsigned int x = 0; x < length; ++x) {
				for (unsigned int y = 0; y < width; ++y) {
					matrix[x][y] = ++p;
				}
			}

			// Make our connections
			vector<unsigned int> connections;
			for (unsigned int x = 0; x < length; ++x) {
				for (unsigned int y = 0; y < width; ++y) {
					connections.clear();
					p = matrix[x][y];

					// If we're on the left bound
					if (x == 0) {
						// Add the particle to the right of us
						connections.push_back(matrix[x+1][y]);
					}
					// If we're on the right bound
					else if (x == length - 1) {
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
					else if (y == width - 1) {
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
			unsigned int desiredArea = populationSize;
			unsigned int divisor = ceil(sqrt(desiredArea));
			while(desiredArea % divisor != 0) {
				++divisor;
			}
			unsigned int length = divisor;
			unsigned int width = desiredArea / divisor;

			// Construct a matrix of dimensions length x width
			// and fill it with particles
			unsigned int p = 0;
			vector<vector<unsigned int>> matrix(length, vector<unsigned int>(width));
			for (unsigned int x = 0; x < length; ++x) {
				for (unsigned int y = 0; y < width; ++y) {
					matrix[x][y] = ++p;
				}
			}

			// Make our connections
			vector<unsigned int> connections;
			for (unsigned int x = 0; x < length; ++x) {
				for (unsigned int y = 0; y < width; ++y) {
					connections.clear();
					p = matrix[x][y];

					if (x == 0) {
						connections.push_back(matrix[x+1][y]);
						connections.push_back(matrix[length-1][y]);
					}
					else if (x == (length - 1)) {
						connections.push_back(matrix[x-1][y]);
						connections.push_back(matrix[0][y]);
					}
					else {
						connections.push_back(matrix[x-1][y]);
						connections.push_back(matrix[x+1][y]);
					}

					if (y == 0) {
						connections.push_back(matrix[x][y+1]);
						connections.push_back(matrix[x][width-1]);
					}
					else if (y == (width - 1)) {
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
			unsigned int usedParticles = 1;
			unsigned int numLevels = 1;
			unsigned int previousLevel = 1;
			unsigned int currentLevel;

			// Determine number of levels in tree
			while (usedParticles < populationSize) {
				currentLevel = previousLevel*2;
				usedParticles += currentLevel;
				previousLevel = currentLevel;
				++numLevels;
			}

			vector<vector<unsigned int>> tree(numLevels, vector<unsigned int>());
			usedParticles = 1;
			currentLevel = 1;
			unsigned int currentParticle = 2;
			bool doneFilling = false;

			// Construct tree and fill it with particles
			tree[0].push_back(1);

			// For each level in tree
			for (unsigned int level = 1; level <= numLevels; ++level) {
				// Current level's particle count is double that of previous level
				currentLevel = currentLevel * 2;

				// For each slot in current level
				for (unsigned int i = 0; i < currentLevel; ++i) {
					// Fill slot with particle
					tree[level].push_back(currentParticle);
					++currentParticle;

					// Make sure we're not filling past our swarm size
					if (currentParticle == populationSize + 1) {
						doneFilling = true;
						break;
					}
				}
				if (doneFilling) {
					break;
				}
			}

			// For each level in tree
			for (unsigned int level = 1; level < numLevels; ++level) {
				int prevGroupCounter = 0;
				int particleCounter = 0;

				// For each particle in level
				for (unsigned int p = 0; p < tree[level].size(); ++p) {

					// Connect current particle (p) with proper particle in last level
					allParticles[tree[level][p]].push_back(tree[level-1][prevGroupCounter]);

					// Connect particle in last level to current particle (p)
					allParticles[tree[level-1][prevGroupCounter]].push_back(tree[level][p]);

					// If our particle counter reaches two, we're in a new pair in (or new particle
					// in previous level)
					if (++particleCounter == 2) {
						particleCounter = 0;
						++prevGroupCounter;
					}
				}
			}
		}
		else {
			outputError("Error: BioNetFit did not recognize your specified topology of: " + options.topology);
		}
	}

	/*
	int p = 0;
	for (auto o = allParticles.begin(); o != allParticles.end(); ++o) {
		cout << p << " is connected to:" << endl;
		for (auto i = (*o).begin(); i != (*o).end(); ++i) {
			cout << *i << endl;
		}
		++p;
	}
	 */

	return allParticles;

	double t = tmr.elapsed();
	cout << "Particle creation took " << t << " seconds" << endl;
}

void Swarm::processParticlesPSO(vector<unsigned int> particles, bool newFlight) {

	if (options.verbosity >= 3) {
		cout << "Processing " << particles.size() << " particles" << endl;
	}

	vector<double> mBests;
	if (options.psoType == "qpso") {
		mBests = calcQPSOmBests();
	}

	// For each particle in our particle set
	for (auto particle = particles.begin(); particle != particles.end(); ++particle) {
		if (options.verbosity >= 3) {
			cout << "Processing particle " << *particle << endl;
		}

		// Calculate the next iteration positions according
		// to user preference
		if (options.psoType == "pso") {
			if (options.enhancedInertia) {
				updateInertia();
			}
			particleCurrParamSets_[*particle] = calcParticlePosPSO(*particle);
		}
		else if (options.psoType == "bbpso") {
			particleCurrParamSets_[*particle] = calcParticlePosBBPSO(*particle);
		}
		else if (options.psoType == "bbpsoexp") {
			particleCurrParamSets_[*particle] = calcParticlePosBBPSO(*particle, true);
		}
		else if (options.psoType == "qpso") {
			updateContractionExpansionCoefficient();
			particleCurrParamSets_[*particle] = calcParticlePosQPSO(*particle, mBests);
		}

		if (options.verbosity >= 3) {
			cout << "Got next positions for particle " << *particle << endl;
		}
		// Convert positions to string so they can be sent to the particle
		vector<string> nextPositionsStr;
		for (auto param = particleCurrParamSets_.at(*particle).begin(); param != particleCurrParamSets_.at(*particle).end(); ++param) {
			nextPositionsStr.push_back(to_string(static_cast<long double>(*param)));
			//cout << *param << endl;
		}

		if (options.verbosity >= 3) {
			cout << "Sending new positions to " << *particle << endl;
		}

		// Finally, send the parameters
		swarmComm->sendToSwarm(0, *particle, SEND_FINAL_PARAMS_TO_PARTICLE, false, nextPositionsStr);
	}
}

vector<double> Swarm::calcParticlePosPSO(unsigned int particle) {

	if (options.verbosity >= 3) {
		cout << "Calculating velocity and position of particle " << particle << endl;
	}

	// This vector holds the new positions to be sent to the particle
	vector<double> nextPositions(particleCurrParamSets_.at(particle).size());
	vector<double> nextVelocities(particleCurrParamSets_.at(particle).size());

	// Get the best positions for particle's neighborhood
	vector<double> neighborhoodBestPositions = getNeighborhoodBestPositions(particle);

	int i = 0;

	//cout << "inertia: " << options.inertia << endl;
	//cout << "cognitive: " << options.cognitive << endl;
	//cout << "social: " << options.social << endl << endl;

	// For each parameter in the current parameter set
	for (auto param = particleCurrParamSets_.at(particle).begin(); param != particleCurrParamSets_.at(particle).end(); ++param) {
		//cout << "before " << *param << endl;
		// Set up formula variables
		double currVelocity = options.inertia * particleParamVelocities_.at(particle)[i];
		//cout << "cv: " << currVelocity << endl;
		double r1 = ((double) rand() / (RAND_MAX)); // TODO: These need to be inclusive
		//cout << "r1: " << r1 << endl;
		double r2 = ((double) rand() / (RAND_MAX));
		//cout << "r2: " << r2 << endl;
		double personalBestPos = particleBestParamSets_.at(particle)[i];
		//cout << "pb: " << personalBestPos << endl;
		double currPos = particleCurrParamSets_.at(particle)[i];
		//cout << "cp: " << currPos << endl;
		//cout << "nbp: " << neighborhoodBestPositions[i] << endl;
		// Set velocity

		// TODO: Look into constriction factor - Clerc
		// Also, according to Eberhart and Shi, constriction factor + velocity clamping
		// results in fastest convergence
		double nextVelocity = (options.inertia * currVelocity) + options.cognitive * r1 * (personalBestPos - currPos) + options.social * r2 * (neighborhoodBestPositions[i] - currPos);
		//cout << "nv: " << nextVelocity << endl;

		// Set velocity
		particleParamVelocities_.at(particle)[i] = nextVelocity;
		// Set position
		nextPositions[i] = currPos + nextVelocity;

		/*
		if (nextPosition <= 0) {
			cout << "NP less than 0. Setting to Very Small Number" << endl;
			nextPosition = 0.00000001;
			particleParamVelocities_.at(particle)[i] = 0;
		}
		 */

		//cout << "after " << nextPositions[i] << endl << endl;
		++i;
	}

	return nextPositions;
}

vector<double> Swarm::calcParticlePosBBPSO(unsigned int particle, bool exp) {

	if (options.verbosity >= 3) {
		cout << "Calculating BBPSO for particle " << particle << endl;
	}

	// Get the best positions for particle's neighborhood
	vector<double> neighborhoodBestPositions = getNeighborhoodBestPositions(particle);

	/*
	cout << "nbps:" << endl;
	for (auto p = neighborhoodBestPositions.begin(); p != neighborhoodBestPositions.end(); ++p) {
		cout << "nbp: " << *p << endl;
	}
	 */

	// This vector holds the new positions to be sent to the particle
	vector<double> nextPositions(particleCurrParamSets_.size());

	bool usePersonalBest = false;

	// For each parameter in the current parameter set
	int i = 0;
	for (auto param = particleBestParamSets_.at(particle).begin(); param != particleBestParamSets_.at(particle).end(); ++param) {
		if (exp) {
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
			cout << "personalbest: " << personalBestPos << " neighborbest: " << neighborhoodBestPositions[i] << endl;
			double mean = (abs(personalBestPos) + abs(neighborhoodBestPositions[i])) / 2;
			cout << "mean: " << mean << endl;
			double std = abs(abs(personalBestPos) - abs(neighborhoodBestPositions[i]));
			cout << "std: " << std << endl;

			// Create the gaussian distribution
			boost::random::normal_distribution<double> dist(mean, std);

			// Pick our next position randomly from distribution
			nextPositions[i] = dist(randNumEngine);
			cout << "picked: " << nextPositions[i] << endl;
		}
		++i;
	}

	return nextPositions;
}

vector<double> Swarm::calcParticlePosQPSO(unsigned int particle, vector<double> mBests) {

	vector<double> nextPositions;
	vector<double> neighborhoodBests = getNeighborhoodBestPositions(particle);
	for (unsigned int d = 0; d < mBests.size(); ++d) {
		//cout << particle << " before: " << particleCurrParamSets_[particle][d] << endl;
		cout << particle << " mbest: " << mBests[d] << endl;
		//cout << particle << " best: " << particleBestParamSets_[particle][d] << endl;
		//cout << particle << " swarmbest: " << getNeighborhoodBestPositions(particle)[d] << endl;
		double fi1 = ((double) rand() / (RAND_MAX));
		//cout << "f1: " << fi1 << endl;
		//double fi2 = ((double) rand() / (RAND_MAX));
		//double p = ((fi1 * particleBestFits_.at(particle)) + (fi2 * swarmBestFits_.begin()->first)) / (fi1 + fi2);
		// TODO: Work out whether or not we should be using abs() for this and below
		double p = fi1 * abs(particleBestParamSets_[particle][d]) + (1 - fi1) * neighborhoodBests[d];
		double u = ((double) rand() / (RAND_MAX));

		// TODO: Linearly decrease beta from 1.0 to 0.5? See Liu et al.

		// Liu et all, 2005, eq 2a
		if (u > 0.5) {
			nextPositions.push_back(p - beta_ * abs(mBests[d] - abs(particleCurrParamSets_.at(particle)[d])) * log(1/u));
			//cout << particle << " p: " << p << endl;
			//cout << "mbest - curr: " << mBests[d] - particleCurrParamSets_.at(particle)[d] << endl;
			//cout << "log: " << log(1/u) << endl;
			//cout << particle << " after: " << p - beta_ * abs(mBests[d] - particleCurrParamSets_.at(particle)[d]) * log(1/u) << endl;
		}
		else {
			nextPositions.push_back(p + beta_ * abs(mBests[d] - abs(particleCurrParamSets_.at(particle)[d])) * log(1/u));
			//cout << particle << " after: " << p + beta_ * abs(mBests[d] - particleCurrParamSets_.at(particle)[d]) * log(1/u) << endl;
		}
	}

	return nextPositions;
}

void Swarm::updateInertia() {
	// Eq 3 from Moraes et al
	options.inertia = options.inertiaInit + (options.inertiaFinal - options.inertiaInit) * (inertiaUpdateCounter_ / (float)(options.nmax + inertiaUpdateCounter_));
	cout << "Setting inertia to " << options.inertia << endl;
}

void Swarm::updateContractionExpansionCoefficient() {
	beta_ = options.betaMax - (((float)flightCounter_ / options.maxNumSimulations) * (options.betaMax - options.betaMin));
	cout << "updating beta_ to " << beta_ << endl;
}

vector<double> Swarm::getNeighborhoodBestPositions(unsigned int particle) {

	if (options.verbosity >= 3) {
		cout << "Getting neighborhood best for particle " << particle << endl;
	}

	// Set the current best fit to our own best fit
	double currBestFit = particleBestFits_.at(particle);
	//cout << "cbf: " << currBestFit << endl;

	int currentBestNeighbor = particle;
	//cout << "set initial nb to " << currentBestNeighbor << endl;

	//cout << "allParticles size: " << populationTopology_.size() << endl;
	//cout << "allParticles neighbor size: " << populationTopology_[particle].size() << endl;

	// For every neighbor in this particle's neighborhood
	for (auto neighbor = populationTopology_[particle].begin(); neighbor != populationTopology_[particle].end(); ++neighbor) {

		auto it = particleBestFits_.find(*neighbor);

		// Skip this neighbor if it doesn't contain a fit value
		if (it->second == 0) {
			continue;
		}

		//cout << "checking if neighbor " << *neighbor << " has a better fit of " << it->second << " than " << currBestFit << endl;

		// If Neighbor's best fit is better than ours, update the best
		// neighbor. Also, update the current best fit value
		if (it->second < currBestFit) {
			//cout << "it does! setting current best fit of particle " << *neighbor << " of " << it->second << endl;
			currentBestNeighbor = *neighbor;
			currBestFit = it->second;
		}
	}

	//cout << "cbn: " << currentBestNeighbor << endl;
	return particleBestParamSets_.at(currentBestNeighbor);
}

void Swarm::processParamsPSO(vector<double> &params, unsigned int pID, double fit) {

	if (options.verbosity >= 3) {
		cout << "Processing finished params for particle " << pID << " with fit of " << fit << endl;
	}

	unsigned int i = 0;
	// If this is the particle's first iteration, we need to store its params
	if (particleIterationCounter_.at(pID) == 1) {
		//cout << "Storing " << params.size() << " params for particle " << pID << endl;
		for (auto param = params.begin(); param != params.end(); ++param) {
			//cout << *param << endl;
			particleCurrParamSets_[pID][i] = *param;
			++i;
		}
	}

	// The the fit value of this param set is less than the particles best fit
	// we should update the particle's best fit, then store the best fit params
	if (particleBestFits_.at(pID) == 0 || fit < particleBestFits_.at(pID)) {
		if (options.verbosity >= 3) {
			cout << "Updating best fit and params for particle " << pID << endl;
		}

		// Insert into best fit lists, erasing in the case of the map with fits as keys
		particleBestFits_[pID] = fit;

		map<double, unsigned int>::iterator toDelIt = particleBestFitsByFit_.end();
		for (auto it = particleBestFitsByFit_.begin(); it != particleBestFitsByFit_.end(); ++it) {
			//cout << "loop: " << it->second << endl;
			if (it->second == pID) {
				//cout << "Erasing old best fit for particle " << pID << endl;
				toDelIt = it;
			}
		}
		if (toDelIt != particleBestFitsByFit_.end()) {
			particleBestFitsByFit_.erase(toDelIt);
		}

		particleBestFitsByFit_.insert(pair<double, unsigned int>(fit, pID));

		unsigned int i = 0;
		for (auto param = params.begin(); param != params.end(); ++param) {
			//cout << "Updating best param for particle " << pID << ": " << *param << endl;
			particleBestParamSets_[pID][i] = *param;
			++i;
		}
	}
}

void Swarm::launchParticle(unsigned int pID, bool nextGen) {

	if (currentGeneration == 1 && !options.useCluster && !nextGen && bootstrapCounter == 0) {

		// Construct command needed to run the particle
		string command = exePath_ + " particle " + to_string(static_cast<long long int>(pID)) + " run " + to_string(static_cast<long long int>(currentGeneration)) + " " + sConf_;
		command = command + " >> pOUT 2>&1";
		command = command + " &";

		if (runCommand(command) != 0) {
			cout << "Warning: Couldn't launch particle " << pID << " with command: " << command <<  endl;
			failedParticles_.insert(pID);
			return;
		}
	}

	if (options.verbosity >= 3) {
		cout << "Running Particle " << pID << endl;
	}

	runningParticles_.insert(pID);
	swarmComm->sendToSwarm(0, pID, NEXT_GENERATION, false, swarmComm->univMessageSender);
}

void Swarm::runGeneration () {
	// TODO: Implement walltime
	if(options.verbosity >= 1) {
		cout << "Running generation " << currentGeneration << " with " << options.swarmSize << " particles..." << endl;
	}

	unsigned int numLaunchedParticles = 0;
	unsigned int numFinishedParticles = 0;

	vector<unsigned int> finishedParticles;
	unsigned int p = 1;

	while (numFinishedParticles < options.swarmSize) {

		if (runningParticles_.size() < options.parallelCount && numLaunchedParticles < options.swarmSize) {
			launchParticle(p);
			numLaunchedParticles += 1;
			++p;
		}

		// Check for any messages from particles
		usleep(10000);
		finishedParticles = checkMasterMessages();
		numFinishedParticles += finishedParticles.size();
		finishedParticles.clear();
	}
	if (failedParticles_.size() > (options.swarmSize - 3) ) {
		finishFit();
		outputError("Error: You had too many failed runs. Check simulation output (.BNG_OUT files) or adjust walltime.");
	}
	currentGeneration += 1;
}

void Swarm::cleanupFiles(const char * path) {
	// TODO: Should we fork to do this? ...yes.
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

void Swarm::breedGenerationGA(vector<unsigned int> children) {
	// TODO: Need to take into consideration failed particles? Definitely need to in first generation.
	if (options.verbosity >= 3) {
		cout << "Breeding generation" << endl;
	}

	if (children.size() == 0) {
		for (unsigned int i = 1; i <= options.swarmSize; ++i) {
			children.push_back(i);
		}
	}

	std::map<unsigned int, std::vector<double>> particleNewParamSets;

	// Create the output directory for the next generation
	if (!checkIfFileExists(options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration)))) {
		string createDirCmd = "mkdir " + options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration));
		int retryCounter = 0;
		while (runCommand(createDirCmd) != 0 && !checkIfFileExists(options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration)))) {
			if(++retryCounter >= 100) {
				outputError("Error: Couldn't create " + options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration)) + " to hold next generation's output.");
			}
			sleep(1);
			cout << "Trying again to create dir" << endl;
		}
	}

	unsigned int parentPoolSize = options.swarmSize;
	parentPoolSize = parentPoolSize - options.keepParents;

	// Create an iterator to our fit list, use it to get our maximum fit value, then reset the iterator to the beginning of the list
	multimap<double, unsigned int>::iterator w = swarmBestFits_.begin();
	advance(w, options.swarmSize - 1);
	double maxWeight = w->first;
	w = swarmBestFits_.begin();

	// Fill the second element of the weight map with difference between maxWeight and fit value
	multimap<double, unsigned int> weightDiffs;
	double diff;
	double weightSum = 0;
	for (unsigned int i = 0; i < options.swarmSize; ++i) {
		diff = maxWeight - w->first;
		weightSum += diff;
		weightDiffs.insert(pair<double, unsigned int>(diff, w->second));
		++w;
	}

	if (weightSum == 0) {
		finishFit();
		outputError("Your population has converged. Quitting.");
	}

	/*
	cout << "max: " << maxWeight << endl;
	cout << "weight sum: " << weightSum << endl;
	for (auto i = weightDiffs.begin(); i != weightDiffs.end(); ++i) {
		cout << "weight diff: " << i->first << " " << i->second << endl;
	}

	for (auto i = swarmBestFits_.begin(); i != swarmBestFits_.end(); ++i) {
		cout << "sbf: " << i->first << " " << i->second << endl;
	}

	for (auto i = allGenFits.begin(); i != allGenFits.end(); ++i) {
		cout << "agf: " << i->first << " " << i->second << endl;
	}
	 */


	// If we want to keep any parents unchanged, send unchanged param sets to children.
	// We start with the global best fit params and iterate through the fit list from there.
	unsigned int childCounter = 1;

	// Only keep parents if we're doing an entire generation at once
	if (children.size() == options.swarmSize) {
		auto parent = allGenFits.begin();
		for (unsigned int p = 1; p <= options.keepParents; ++p) {

			vector<string> params;
			split(parent->second, params);
			params.erase(params.begin());

			for (auto param = params.begin(); param != params.end(); ++param) {
				particleNewParamSets[childCounter].push_back(stod(*param));
			}

			//cout << "Sending unchanged params to " << childCounter << endl;
			swarmComm->sendToSwarm(0, childCounter, SEND_FINAL_PARAMS_TO_PARTICLE, false, params);
			++parent;
			++childCounter;
		}
	}

	float parentPairs;
	//cout << "children.size(): " << children.size() << endl;
	if (children.size() == options.swarmSize) {
		parentPairs = (float)parentPoolSize / 2;
	}
	else {
		parentPairs = (float)children.size() / 2;
	}

	//cout << "we have " << parentPairs << " parent pairs" << endl;

	boost::random::uniform_int_distribution<int> unif(1, 100);
	for (unsigned int i = 0; i < parentPairs; ++i) {

		unsigned int p1;
		unsigned int p2;
		//unsigned int maxFitCounter = 0;
		// TODO: maxfit needs to be implemented. to do this we need to be able to lookup the chosen particles fit value...
		// Pick the fit values (particle parents) used in breeding
		//do {
		p1 = pickWeighted(weightSum, weightDiffs, options.extraWeight);
		p2 = pickWeighted(weightSum, weightDiffs, options.extraWeight);

		// Quit if we try to many times to select suitable parents
		//if (++maxFitCounter >= 10000) {
		//	outputError("Error: Tried too many times to select parents that didn't exceed the specified max_fit value of " + to_string(static_cast<long double>(options.maxFit)) + ". Quitting.");
		//}
		// Make sure we don't exceed max_fit
		//} while (options.maxFit != 0 && particleBestFits_.at(p1) >= options.maxFit && particleBestFits_.at(p2) >= options.maxFit);

		// If we want different parents used in breeding, make sure that happens
		unsigned int retryCount = 0;
		while (p1 == p2 && options.forceDifferentParents) {
			retryCount++;
			if (retryCount > options.maxRetryDifferentParents) {
				if (options.verbosity >= 1) {
					cout << "Tried too many time to select different parents for breeding. Selecting the first two." << endl;
				}

				// Get reverse iterator to the weight map (best parents are at the end)
				multimap<double, unsigned int>::reverse_iterator w = weightDiffs.rbegin();

				// The weight map is sorted, so the first element will be the best fit
				p1 = w->second;
				//cout << "1: " << w->second;

				// Increment the map iterator until we find a fit value that isn't ours
				while (p1 == p2 && w != weightDiffs.rend()) {
					++w;
					p2 = w->second;
					//cout << "tried to set p2 as " << w->second;
				}

				break;
			}
			p2 = pickWeighted(weightSum, weightDiffs, options.extraWeight);
			//cout << "retry2: " << p2 << endl;
		}

		//cout << "chose " << p1 << " and " << p2 << endl;

		auto p1It = allGenFits.begin();
		advance(p1It, p1);
		vector<string> p1Vec;
		split(p1It->second, p1Vec);
		p1Vec.erase(p1Vec.begin());

		auto p2It = allGenFits.begin();
		advance(p2It, p2);
		vector<string> p2Vec;
		split(p2It->second, p2Vec);
		p2Vec.erase(p2Vec.begin());

		vector<string> c1Vec, c2Vec;
		unsigned int pi = 0;
		for (auto p = options.model->getFreeParams_().begin(); p != options.model->getFreeParams_().end(); ++p) {
			string p1Param, p2Param;
			//cout << "param: " << p->first << endl;
			p1Param = p1Vec[pi];
			p2Param = p2Vec[pi];

			//cout << "p1: " << p1Param << endl;
			//cout << "p2: " << p2Param << endl;
			if (unif(randNumEngine) < (options.swapRate * 100) ) {

				// TODO: Make sure individual mutation rates work
				if (hasMutate && p->second->isHasMutation()) {
					//cout << "about to mutate" << endl;
					p1Param = mutateParamGA(p->second, stod(p1Param));
					p2Param = mutateParamGA(p->second, stod(p2Param));
					//cout << "mutated" << endl;
				}

				//cout << "swapping p1: " << p1Param << " with p2: " << p2Param << endl;
				c1Vec.push_back(p2Param);
				particleNewParamSets[children[childCounter - 1]].push_back(stod(p2Param));
				//cout << "c1 final: " << p2Param << endl;
				c2Vec.push_back(p1Param);
				particleNewParamSets[children[childCounter]].push_back(stod(p1Param));
			}
			else {
				c1Vec.push_back(p1Vec[pi]);
				particleNewParamSets[children[childCounter - 1]].push_back(stod(p1Vec[pi]));
				//cout << "c1 final: " << p1Vec[pi] << endl;
				c2Vec.push_back(p2Vec[pi]);
				particleNewParamSets[children[childCounter]].push_back(stod(p2Vec[pi]));
			}
			++pi;
		}

		//cout << "sending to " << children[childCounter - 1] << endl;
		swarmComm->sendToSwarm(0, children[childCounter - 1], SEND_FINAL_PARAMS_TO_PARTICLE, false, c1Vec);
		++childCounter;

		// Make sure we don't send to too many parents (only relevant in last breeding with odd number of parents)
		if ( !( (fmod(parentPairs * 2, 2)) == 1 && parentPairs - i == 0.5 ) ) {
			//cout << "sending to " << children[childCounter - 1] << endl;
			swarmComm->sendToSwarm(0, children[childCounter - 1], SEND_FINAL_PARAMS_TO_PARTICLE, false, c2Vec);
			++childCounter;
		}
	}

	//particleCurrParamSets_ = particleNewParamSets;

	// Replace any changed param sets in the master set
	for (auto child = particleNewParamSets.begin(); child != particleNewParamSets.end(); ++child) {
		particleCurrParamSets_[child->first] = child->second;
	}

	unsigned int numFinishedBreeding = 0;
	//cout << "waiting for " << childCounter - 1 << endl;
	while (numFinishedBreeding < (childCounter - 1)) {
		//cout << "checking for DONEBREED" << endl;
		unsigned int numMessages = swarmComm->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm->univMessageReceiver, true);
		numFinishedBreeding += numMessages;
		//cout << numFinishedBreeding << endl;
	}
	//cout << "done?" << endl;
}

void Swarm::finishFit() {

	string outputDir;
	if (options.bootstrap) {
		outputDir = options.outputDir + "/" + options.jobName + "_bootstrap/";

		if (!checkIfFileExists(outputDir)) {
			string command = "mkdir " + outputDir + " && cp " + configPath_ + " " + outputDir;

			int errorCounter = 0;
			while (runCommand(command) != 0 && errorCounter < 10) {
				sleep(1);
				++errorCounter;
			}

			if (errorCounter > 10) {
				cout << "Error: Couldn't create directory: " + outputDir + " to contain final bootstrap results.";
				//outputBootstrapSummary();
			}
		}

		// Output to bootstrap file
		outputBootstrapSummary();

		// Output fit summary
		string outputFilePath = outputDir + to_string(static_cast<long long int>(bootstrapCounter + 1)) + "_all_fits.txt";
		outputRunSummary(outputFilePath);

		// Reset variables
		resetVariables();

		cout << "counter: " << bootstrapCounter << " bootstrap " << options.bootstrap << endl;
		if ((bootstrapCounter + 1) < options.bootstrap) {
			cout << "Deleting last fitting run and beginning next fit.." << endl;

			string command = "cd " + options.jobOutputDir + " && find -mindepth 1 -maxdepth 1 -type d -exec rm -r {} \\";
			runCommand(command);

			killAllParticles(NEW_BOOTSTRAP);
		}
		else {
			killAllParticles(FIT_FINISHED);
		}
	}
	else {
		outputDir = options.jobOutputDir + "Results/";
		string command = "mkdir " + outputDir + " && cp " + configPath_ + " " + outputDir;

		int errorCounter = 0;
		while (runCommand(command) != 0 && errorCounter < 10) {
			sleep(1);
			++errorCounter;
		}

		if (errorCounter > 10) {
			cout << "Error: Couldn't create directory: " + outputDir + " to contain final fitting results. Outputting results to screen.";
			outputRunSummary();
		}

		string outputFilePath = outputDir + "all_fits.txt";
		outputRunSummary(outputFilePath);

		killAllParticles(FIT_FINISHED);

		cout << "Finished fitting in " << tmr_.elapsed() << " seconds. Results can be found in " << options.jobOutputDir << "Results" << endl;
	}


	//generateBestFitModels(outputDir);
	//copyBestFitToResults(outputDir);


}

void Swarm::resetVariables() {
	allGenFits.clear();
	currentGeneration = 1;
	flightCounter_ = 0;
}

void Swarm::outputRunSummary(string outputPath) {
	ofstream outputFile;

	outputFile.open(outputPath, ofstream::out);

	if (outputFile.is_open()) {
		outputFile.precision(8);
		outputFile.setf(ios::scientific);

		// Output first two fields of header
		outputFile << left << setw(16) << "Fit" << left << setw(16) << "Iteration";

		// Output parameter names
		for (auto i = options.model->freeParams_.begin(); i != options.model->freeParams_.end(); ++i) {
			outputFile << left << setw(16) << i->first;
		}

		outputFile << endl;

		vector<string> paramVals;

		for (auto i = allGenFits.begin(); i != allGenFits.end(); ++i) {
			//cout << "about to split" << endl;
			split(i->second, paramVals);
			//cout << "split done" << endl;
			outputFile << left << setw(16) << i->first << left << setw(16) << paramVals[0];
			for (unsigned int i = 1; i < paramVals.size(); i++) {
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

void Swarm::outputRunSummary() {
	// Output first two fields of header
	cout << left << setw(16) << "Fit" << left << setw(16) << "Iteration";

	// Output parameter names
	//for (autoDE i: options.model->freeParams_) {
	for (auto i = options.model->freeParams_.begin(); i != options.model->freeParams_.end(); ++i) {
		cout << left << setw(16) << i->first;
	}

	cout << endl;

	vector<string> paramVals;

	//for (auto i: allGenFits) {
	for (auto i = allGenFits.begin(); i != allGenFits.end(); ++i) {
		split(i->second, paramVals);

		cout << left << setw(16) << i->first << left << setw(16) << paramVals[0];

		for (unsigned int i = 1; i < paramVals.size(); i++) {
			cout << left << setw(16) << stod(paramVals[i]);
		}

		paramVals.clear();

		cout << endl;
	}
}

void Swarm::outputBootstrapSummary() {
	string outputPath = options.outputDir + "/" + options.jobName + "_bootstrap/bootstrap_results.txt";
	ofstream outputFile;
	outputFile.open(outputPath, ofstream::out | std::ofstream::app);

	if (outputFile.is_open()) {
		outputFile.precision(8);
		outputFile.setf(ios::scientific);

		if (bootstrapCounter == 0) {
			// Output first two fields of header
			outputFile << left << setw(16) << "Fit";

			// Output parameter names
			for (auto i = options.model->freeParams_.begin(); i != options.model->freeParams_.end(); ++i) {
				outputFile << left << setw(16) << i->first;
			}
			outputFile << endl;
		}

		vector<string> paramVals;
		split(allGenFits.begin()->second, paramVals);
		outputFile << left << setw(16) << allGenFits.begin()->first;
		for (unsigned int i = 1; i < paramVals.size(); i++) {
			outputFile << left << setw(16) << stod(paramVals[i]);
		}
		outputFile << endl;
		outputFile.close();
	}
	else {
		cout << "Warning: couldn't open: " << outputFile << " to write fit summary." << endl;
	}
}

void Swarm::killAllParticles(int tag) {
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		cout << "killing " << p << " with tag: " << tag << endl;
		swarmComm->sendToSwarm(0, p, tag, false, swarmComm->univMessageSender);
	}
	swarmComm->univMessageSender.clear();
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

void Swarm::initFit () {

	if (resumingSavedSwarm) {
		if (options.verbosity >= 3) {
			cout << "Resuming fitting run" << endl;
		}

		for (auto particle = particleIterationCounter_.begin(); particle != particleIterationCounter_.end(); ++particle) {
			swarmComm->univMessageSender.push_back(to_string(static_cast<long long int>(currentGeneration - 1)));
			swarmComm->sendToSwarm(0, particle->first, SEND_NUMFLIGHTS_TO_PARTICLE, false, swarmComm->univMessageSender);
			swarmComm->univMessageSender.clear();

			for (auto param = particleCurrParamSets_.at(particle->first).begin(); param != particleCurrParamSets_.at(particle->first).end(); ++param) {
				//cout << "Sending " << *param << " to " << particle->first << endl;
				swarmComm->univMessageSender.push_back(to_string(static_cast<long double>(*param)));
			}
			swarmComm->sendToSwarm(0, particle->first, SEND_FINAL_PARAMS_TO_PARTICLE, false, swarmComm->univMessageSender);
			swarmComm->univMessageSender.clear();
		}

		// We need to set current generation to the iteration counter because the saved swarm
		// currentGeneration may not be accurate due to it changed between runGeneration() and
		// breedGenerationGA()
		if (options.fitType == "genetic") {
			currentGeneration = currentGeneration - 1;
		}
	}
	else {
		if (options.verbosity >= 3) {
			cout << "Initializing fitting run" << endl;
		}

		// If we are running ODE, we want to generate a network and network reader.
		// The network file has parameters that are replaced in every generation,
		// and the network reader is used to read those network files
		if (options.model->getHasGenerateNetwork() && bootstrapCounter == 0) {

			if (options.verbosity >= 3) {
				cout << "Creating dummy .bngl, running to get a .net file" << endl;
			}

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
			string command = options.bngCommand + " --outdir " + options.jobOutputDir + " " + modelPath + " >> " + options.jobOutputDir + "netgen_output 2>&1";
			cout << "Generating initial .net file with command: " << command << endl;

			if (options.useCluster) {
				command = generateSlurmCommand(command, false, 1);
			}

			if (runCommand(command) != 0) {
				outputError("Error: Couldn't generate initial .net file with command: " + command + ". Quitting.");
			}

			// Now that we have a .net file, replace the full .bngl with a .bngl
			// containing ONLY action commands and a .net file loader. This is
			// used later to run our .net files.
			options.model->outputModelWithParams(p.getParams(), options.jobOutputDir, "base.bngl", "", false, true, false, false, false);

			// Now store our .net file for later use
			string netPath = options.jobOutputDir + "/base.net";
			options.model->parseNet(netPath);
		}

		vector<double> params;
		for (unsigned int param = 0; param < options.model->getNumFreeParams(); ++param) {
			params.push_back(0);
		}
		for (unsigned int particle = 1; particle <= options.swarmSize; ++particle) {
			particleCurrParamSets_.insert(pair<int, vector<double>>(particle, params));
		}

		if (options.fitType == "pso") {
			for (unsigned int particle = 1; particle <= options.swarmSize; ++particle) {
				particleParamVelocities_[particle] = params;
				particleBestParamSets_[particle] = params;
				particleWeights_[particle] = 0;
				particleIterationCounter_[particle] = 0;
			}
		}
		if (options.fitType == "pso" || "de") {
			for (unsigned int particle = 1; particle <= options.swarmSize; ++particle) {
				particleBestFits_[particle] = 0;
			}
		}
		if(options.fitType == "de") {
			for (unsigned int island = 0; island <= options.numIslands; ++island) {
				islandToParticle_.push_back(vector<unsigned int>(options.swarmSize / options.numIslands,0));
			}
			for (unsigned int particle = 1; particle <= options.swarmSize; ++particle) {
				particleToIsland_.push_back(0);
				particleBestFitsByFit_.insert(pair<double, unsigned int>(0, particle));
			}
		}
		if (options.fitType == "sa") {
			for (unsigned int particle = 1; particle <= options.swarmSize; ++particle) {
				particleBestFitsByFit_.insert(pair<double, unsigned int>(0, particle));
			}
		}
	}
}

vector<unsigned int> Swarm::checkMasterMessages() {
	vector<unsigned int> finishedParticles;

	if (options.verbosity >= 3) {
		//cout << "Checking messages" << endl;
	}

	// First let's check interswarm communication
	int numMessages = swarmComm->recvMessage(-1, 0, -1, false, swarmComm->univMessageReceiver);

	if (numMessages >= 1) {
		if (options.verbosity >= 3) {
			//cout << "Found " << numMessages << " messages" << endl;
		}

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
				string errMsg = "Error: Couldn't remove particle " + to_string(static_cast<long long int>(pID)) + " from runningParticle list.";
				outputError(errMsg);
			}
			runningParticles_.erase(runningParticlesIterator_);

			++particleIterationCounter_[pID];
			++flightCounter_;
			finishedParticles.push_back(pID);

			unsigned int gen = currentGeneration;
			string paramsString;
			if (options.synchronicity == 1) {
				if (options.fitType == "genetic") {
					paramsString = "gen" + to_string(static_cast<long long int>(gen)) + "perm" + to_string(static_cast<long long int>(pID)) + " ";
				}
				else if (options.fitType == "pso") {
					paramsString = "flock" + to_string(static_cast<long long int>(gen)) + "particle" + to_string(static_cast<long long int>(pID)) + " ";
				}

			}

			else if (options.synchronicity == 0) {
				// Increment our flight counter

				// TODO: This run summary contains 1 less particle than it should.
				// Output a run summary every outputEvery flights
				if (flightCounter_ % options.outputEvery == 0) {
					string outputPath = options.jobOutputDir + to_string(static_cast<long long int>(flightCounter_)) + "_summary.txt";
					outputRunSummary(outputPath);
				}

				paramsString = to_string(static_cast<long long int>(flightCounter_)) + " ";
			}

			// Store the parameters given to us by the particle
			vector<double> paramsVec;
			int i = 0;
			for (vector<string>::iterator m = sm->second.message.begin() + 1; m != sm->second.message.end(); ++m) {
				paramsString += *m + " ";
				paramsVec.push_back(stod(*m));
				if (options.fitType == "genetic") {
					particleCurrParamSets_[pID][i] = stod(*m);
				}
				//cout << "pushed back: " << *m << endl;
				++i;
			}
			//cout << "stored params" << endl;

			double fitCalc = stod(sm->second.message[0]);
			allGenFits.insert(pair<double, string>(fitCalc, paramsString));
			swarmBestFits_.insert(pair<double, unsigned int>(fitCalc, pID));

			if (options.fitType == "pso") {
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
			if (runningParticlesIterator_ == runningParticles_.end()) {
				string errMsg = "Error: Couldn't remove particle " + to_string(static_cast<long long int>(pID)) + " from runningParticle list.";
				outputError(errMsg);
			}
			// Increment our finished counter
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

	}

	// Now let's check for any external messages
	//checkExternalMessages();

	return finishedParticles;
}

unordered_map<unsigned int, vector<double>> Swarm::checkMasterMessagesDE() {

	if (options.verbosity >= 3) {
		//cout << "Checking messages" << endl;
	}

	unordered_map<unsigned int, vector<double>> particleParams;

	// First let's check interswarm communication
	int numMessages = swarmComm->recvMessage(-1, 0, -1, false, swarmComm->univMessageReceiver);

	if (numMessages >= 1) {
		if (options.verbosity >= 3) {
			//cout << "Found " << numMessages << " messages" << endl;
		}

		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		smhRange = swarmComm->univMessageReceiver.equal_range(SIMULATION_END);
		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
			int pID = sm->second.sender;

			if (options.verbosity >= 3) {
				//cout << "Particle " << pID << " finished simulation" << endl;
			}

			// Get an iterator to the particle in our list of running particles
			runningParticlesIterator_ = runningParticles_.find(pID);
			// Then remove it
			if (runningParticlesIterator_ == runningParticles_.end()) {
				string errMsg = "Error: Couldn't remove particle " + to_string(static_cast<long long int>(pID)) + " from runningParticle list.";
				outputError(errMsg);
			}
			runningParticles_.erase(runningParticlesIterator_);

			++particleIterationCounter_[pID];

			/*
			// Only increment flight counter if we're finishing a trial set
			if (trial == true) {
				++flightCounter_;
			}
			 */

			/*
			unsigned int gen = currentGeneration;
			string paramsString;
			// Only output a run summary after a trial set has finished
			if (options.synchronicity == 1 && trial == true) {
				paramsString = "gen" + to_string(static_cast<long long int>(gen)) + "perm" + to_string(static_cast<long long int>(pID)) + " ";
			}
			else if (options.synchronicity == 0 && trial == true) {
				// TODO: This run summary contains 1 less particle than it should.
				// Output a run summary every outputEvery flights
				if (flightCounter_ % options.outputEvery == 0) {
					string outputPath = options.jobOutputDir + to_string(static_cast<long long int>(flightCounter_)) + "_summary.txt";
					outputRunSummary(outputPath);
				}

				paramsString = to_string(static_cast<long long int>(flightCounter_)) + " ";
			}
			 */

			double fitCalc = stod(sm->second.message[0]);

			/*
			bool replaceParams = false;
			if (trial == false) {
				particleBestFits_[pID] = fitCalc;
				insertKeyByValue(particleBestFitsByFit_, fitCalc, pID);
				replaceParams = true;
			}
			else if (trial == true && fitCalc < particleBestFits_.at(pID)) {
				cout << "replacing" << endl;
				particleBestFits_[pID] = fitCalc;
				insertKeyByValue(particleBestFitsByFit_, fitCalc, pID);
				replaceParams = true;
			}
			 */

			// Store the parameters given to us by the particle
			vector<double> paramsVec;
			paramsVec.push_back(fitCalc);
			//int i = 0;
			for (vector<string>::iterator m = sm->second.message.begin() + 1; m != sm->second.message.end(); ++m) {

				paramsVec.push_back(stod(*m));

				/*
				if (replaceParams) {
					//cout << "stored " << stod(*m) << " for particle " << pID << " as " << particleCurrParamSets_[pID][i] << endl;
					particleCurrParamSets_[pID][i] = stod(*m);
					paramsString += *m + " ";
				}
				else if (trial) {
					paramsString += to_string(static_cast<long double>(particleCurrParamSets_[pID][i])) + " ";
				}
				 */

				//++i;
			}

			/*
			if (trial == true) {
				allGenFits.insert(pair<double, string>(particleBestFits_[pID], paramsString));
				swarmBestFits_.insert(pair<double, unsigned int>(fitCalc, pID));
			}
			 */

			particleParams.insert(pair<unsigned int, vector<double>>(pID, paramsVec));
			//cout << "done processing particle " << pID << endl;
		}

		// TODO: When sending NEXT_GENERATION, make sure failed particles have actually run again. If not, they need re-launched.
		smhRange = swarmComm->univMessageReceiver.equal_range(SIMULATION_FAIL);
		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {

			int pID = sm->second.sender;

			// Get an iterator to the particle in our list of running particles
			runningParticlesIterator_ = runningParticles_.find(pID);

			// Then remove it
			if (runningParticlesIterator_ == runningParticles_.end()) {
				string errMsg = "Error: Couldn't remove particle " + to_string(static_cast<long long int>(pID)) + " from runningParticle list.";
				outputError(errMsg);
			}
			// Increment our finished counter
			//finishedParticles.push_back(pID);
			particleParams.insert(pair<unsigned int, vector<double>>(pID, vector<double>()));

			cout << "Particle " << pID << " failed in gen " << currentGeneration << endl;

			// Store particle ID in our list of failed particles
			failedParticles_.insert(pID);
		}

		swarmComm->univMessageReceiver.clear();

	}

	return particleParams;
}

void Swarm::checkExternalMessages() {
	string messagePath = options.jobOutputDir + ".req";
	if (checkIfFileExists(messagePath)) {
		ifstream messageFile(messagePath);

		std::string line;

		if (messageFile.is_open()) {
			while (getline(messageFile, line)) {
				if (line == "output results") {
					string outputDir = options.jobOutputDir + "Results";
					string command = "mkdir " + outputDir + " && cp " + configPath_ + " " + outputDir;

					int errorCounter = 0;
					while (!runCommand(command) && errorCounter < 10) {
						sleep(1);
						++errorCounter;
					}

					if (errorCounter > 10) {
						cout << "Error: Couldn't create directory: " + outputDir + " to contain final fitting results. Outputting results to screen.";
						outputRunSummary();
					}

					outputRunSummary(outputDir);
				}
			}
		}
		else {
			string errMsg = "Error: Couldn't open data file " + messagePath + " for parsing.";
			outputError(errMsg);
		}
		messageFile.close();
		string delCmd = "rm " + messagePath;
		runCommand(delCmd);
	}
}

string Swarm::generateSlurmMultiProgCmd(string runCmd) {
	string multiProgConfPath = options.jobOutputDir + "multiprog.conf";
	ofstream multiprog(multiProgConfPath, ios::out);

	if (multiprog.is_open()) {
		multiprog << "0 " << runCmd << " load " << sConf_ << endl;
		for (unsigned int id = 1; id <= options.swarmSize; ++id) {
			multiprog << to_string(static_cast<long long int>(id)) + " " << runCmd << " particle " << to_string(static_cast<long long int>(id)) << " run " << sConf_ << endl;
		}
		multiprog.close();
	}

	return generateSlurmCommand(multiProgConfPath);
}

string Swarm::generateTorqueBatchScript(string cmd) {
	string batchScriptPath = options.jobOutputDir + "batch_script.sh";
	ofstream batchScript(batchScriptPath, ios::out);

	if (batchScript.is_open()) {
		batchScript << "#!/bin/bash" << endl << endl;

		batchScript << "#PBS -N " << options.jobName << endl;

		if (options.saveClusterOutput) {
			batchScript << "#PBS -o " + options.outputDir + "/" + options.jobName + "_cluster_output/" + options.jobName << endl;
			batchScript << "#PBS -j oe" << endl;
		}
		else {
			batchScript << "#PBS -o /dev/null" << endl;
		}

		if (!options.clusterQueue.empty()) {
			batchScript << "#PBS -p	" + options.clusterQueue << endl;
		}

		// Specify the cluster account if needed
		if (!options.clusterAccount.empty()) {
			batchScript << "#PBS -A " + options.clusterAccount << endl;
		}

		batchScript << "#PBS -l procs=" << options.swarmSize + 1 << ",";

		if (options.maxFitTime != MAX_LONG) {
			batchScript << " walltime=00:00:" << to_string(static_cast<long long int>(options.maxFitTime)) << ",";
		}

		batchScript << endl << endl;

		batchScript << "mpirun -np1" << cmd << " load " << sConf_ << " : cmd ";

		batchScript.close();
	}

	return cmd;
}

string Swarm::generateSlurmCommand(string cmd, bool multiProg, int nCPU) {
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

	if (options.emailWhenFinished) {
		command += " --mail-type=END,FAIL --mail-user=" + options.emailAddress;
	}

	if (options.maxFitTime != MAX_LONG) {
		command += " -t 00:00:" + to_string(static_cast<long long int>(options.maxFitTime));
	}

	// Specify output directory if needed
	if (options.saveClusterOutput) {
		command += " -o " + options.outputDir + "/" + options.jobName + "_cluster_output/" + options.jobName;
	}
	else {
		command += " -o /dev/null";
	}

	if (nCPU == 0) {
		command += " -n" + to_string(static_cast<long long int>(options.parallelCount) + 1);
	}
	else {
		command += " -n" + to_string(static_cast<long long int>(nCPU));
	}

	command += " -l";


	if (multiProg) {
		command += " --multi-prog " + cmd;
	}
	else {
		command+= " " + cmd;
	}

	return command;
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

		sbatch << "#SBATCH -n" + to_string(static_cast<long long int>(options.parallelCount) + 1) << endl;

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
	else {
		outputError("Error: Couldn't generate slurm batch file: " + sbatchPath + ". Quitting.");
	}

	runCmd = "sbatch " + sbatchPath;

	return runCmd;
}

unsigned int Swarm::pickWeighted(double weightSum, multimap<double, unsigned int> &weights, unsigned int extraWeight) {
	double lowerBound = 0;
	double upperBound = weightSum;

	if (upperBound <= 0) {
		return 0;
	}

	boost::random::uniform_real_distribution<double> unif(lowerBound, upperBound);

	double random = unif(randNumEngine);
	double chosen = random * ( 1 - (extraWeight / 10 ));

	//cout << "chosen: " << chosen << endl;

	double currentSum = 0;
	unsigned int indexCounter = 0;
	for (multimap<double, unsigned int>::reverse_iterator w = weights.rbegin(); w != weights.rend(); ++w) {
		currentSum += w->first;
		//cout << "adding " << w->first << endl;

		if (currentSum >= chosen) {
			return indexCounter;
		}
		++indexCounter;
	}
	//cout << "fell off the end" << endl;
	return 0;
}

void Swarm::insertKeyByValue(multimap<double, unsigned int> &theMap, double key, unsigned int value) {
	//cout << "looking for " << value << endl;
	for (auto thePair = theMap.begin(); thePair != theMap.end(); ++thePair) {
		if (thePair->second == value) {
			//cout << "found " << value << " inserting " << key << endl;
			theMap.erase(thePair);
			theMap.insert(pair<double, unsigned int>(key, value));
			return;
		}
	}
}

string Swarm::mutateParamGA(FreeParam* fp, double paramValue) {
	//Timer tmr;
	//uniform_real_distribution<double> unif(0,1);
	boost::random::uniform_real_distribution<double> unif(0,1);

	// Generate a random number and see if it's less than our mutation rate.  If is, we mutate.
	if (unif(randNumEngine) < fp->getMutationRate()) {
		// Store our mutation factor
		float maxChange = paramValue * fp->getMutationFactor();

		// Generate a new distribution between 0 and double maxChange
		//using param_t = uniform_real_distribution<>::param_type;
		//param_t p{0.0, maxChange * 2};
		//unif.param(p);

		boost::random::uniform_real_distribution<double> unif(0.0, maxChange * 2);

		// Ger our new random number between 0 and maxChange, subtract maxChange.
		double change = unif(randNumEngine) - maxChange;

		// Add/subtract the value from the parameter
		paramValue+= change;
		//cout << "factor: " << fp->getMutationFactor() << " max change: " << maxChange << " rand: " << rand << endl << endl;
	}
	//double t = tmr.elapsed();
	//cout << "Mutate took " << t << " seconds" << endl;
	return to_string(static_cast<long double>(paramValue));
}

void Swarm::saveSwarmState() {
	string serializedSwarmPath = options.jobOutputDir + "swarmState.sconf";

	std::ofstream ofs(serializedSwarmPath);
	if (ofs.is_open()) {

		boost::archive::binary_oarchive ar(ofs);
		ar & this;
		ofs.close();
	}
	else {
		cout << "Warning: Couldn't save swarm state to " << serializedSwarmPath << endl;
	}
}

void Swarm::initPSOswarm(bool resumeFit) {
	vector<unsigned int> finishedParticles;
	unsigned int numFinishedParticles = 0;

	if (options.verbosity >= 3) {
		cout << "Launching swarm" << endl;
	}

	// Launch all particles
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		launchParticle(p);
	}

	if (options.verbosity >= 3) {
		cout << "Waiting for particles to finish" << endl;
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

	// Fill a vector with all pID's to send for processing
	vector<unsigned int> allParticles;
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
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
}

void Swarm::outputError(string errorMessage) {
	cout << errorMessage << endl;

	killAllParticles(FIT_FINISHED);

	exit (1);
}

vector<double> Swarm::mutateParticleDE(unsigned int particle, float mutateFactor) {

	if (options.verbosity >= 3) {
		cout << "Mutating particle " << particle << endl;
	}

	vector<double> mutatedParams;
	unsigned int currIsland;;
	unsigned int particlesPerIsland;

	if (options.fitType == "de") {
		currIsland = particleToIsland_.at(particle);
		particlesPerIsland = options.swarmSize / options.numIslands;
	}

	// f between 0.1 and 0.9
	//float f = 0.1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.9-0.1)));
	float f;
	if (mutateFactor) {
		f = mutateFactor;
	}
	else {
		f = 0.9;
	}

	boost::random::uniform_int_distribution<int> unif(0, particlesPerIsland - 1);

	if (options.mutateType == 1) {
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0 && particleToIsland_.at(fitVal->second) == currIsland) {
				bestParticle = fitVal->second;
				//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
				break;
			}
		}

		vector<double> bestParamSet = particleCurrParamSets_.at(bestParticle);
		unsigned int pi = 0;
		for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;

			while (p1 == p2) {
				p1 = unif(randNumEngine);
				p2 = unif(randNumEngine);
			}

			double p1Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p1])[pi];
			double p2Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p2])[pi];

			mutatedParams.push_back(bestParamSet[pi] + f * (p1Param - p2Param));
			++pi;
		}
	}
	else if (options.mutateType == 2) {
		unsigned int pi = 0;
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0 && particleToIsland_.at(fitVal->second) == currIsland) {
				bestParticle = fitVal->second;
			}
		}
		vector<double> bestParamSet = particleCurrParamSets_.at(bestParticle);
		for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			int p4 = 0;
			while (p1 == p2 || p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4 || p3 == p4) {
				p1 = unif(randNumEngine);
				p2 = unif(randNumEngine);
				p3 = unif(randNumEngine);
				p4 = unif(randNumEngine);
			}

			double p1Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p1])[pi];
			double p2Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p2])[pi];
			double p3Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p3])[pi];
			double p4Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p4])[pi];

			mutatedParams.push_back(bestParamSet[pi] + f * (p1Param - p2Param) + f * (p3Param - p4Param));
			++pi;
		}
	}
	else if (options.mutateType == 3) {
		unsigned int pi = 0;
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0 && particleToIsland_.at(fitVal->second) == currIsland) {
				bestParticle = fitVal->second;
			}
		}
		vector<double> bestParamSet = particleCurrParamSets_.at(bestParticle);
		for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			while (p1 == p2 || p1 == p3 || p2 == p3) {
				p1 = unif(randNumEngine);
				p2 = unif(randNumEngine);
				p3 = unif(randNumEngine);
			}

			double p1Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p1])[pi];
			double p2Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p2])[pi];
			double p3Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p3])[pi];

			mutatedParams.push_back(p1Param + f * (p2Param - p3Param));
			++pi;
		}
	}
	else if (options.mutateType == 4) {
		unsigned int pi = 0;
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0 && particleToIsland_.at(fitVal->second) == currIsland) {
				bestParticle = fitVal->second;
			}
		}
		vector<double> bestParamSet = particleCurrParamSets_.at(bestParticle);
		for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			int p4 = 0;
			int p5 = 0;
			// This (and above instances) should be optimized
			while (p1 == p2 || p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4 || p3 == p4 || p1 == p5 || p2 == p5 || p3 == p5 || p4 == p5) {
				p1 = unif(randNumEngine);
				p2 = unif(randNumEngine);
				p3 = unif(randNumEngine);
				p4 = unif(randNumEngine);
				p5 = unif(randNumEngine);
			}

			double p1Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p1])[pi];
			double p2Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p2])[pi];
			double p3Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p3])[pi];
			double p4Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p4])[pi];
			double p5Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p5])[pi];

			mutatedParams.push_back(p1Param + f * (p2Param - p3Param) + f * (p4Param - p5Param));
			++pi;
		}
	}

	return mutatedParams;
}

vector<double> Swarm::mutateParticleSA(unsigned int particle, float mutateFactor) {

	if (options.verbosity >= 3) {
		cout << "Mutating particle " << particle << endl;
	}

	vector<double> mutatedParams;

	// f between 0.1 and 0.9
	//float f = 0.1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.9-0.1)));
	float f;
	if (mutateFactor) {
		f = mutateFactor;
	}
	else {
		f = 0.9;
	}

	boost::random::uniform_int_distribution<int> unif(1, options.swarmSize);

	if (options.mutateType == 1) {
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0 ) {
				bestParticle = fitVal->second;
				//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
				break;
			}
		}

		vector<double> bestParamSet = particleCurrParamSets_.at(bestParticle);
		unsigned int pi = 0;
		for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;

			while (p1 == p2) {
				p1 = unif(randNumEngine);
				p2 = unif(randNumEngine);
			}

			//cout << "p1: " << p1 << " p2: " << p2 << " pi: " << pi << endl;
			double p1Param = particleCurrParamSets_.at(p1)[pi];
			double p2Param = particleCurrParamSets_.at(p2)[pi];

			mutatedParams.push_back(bestParamSet[pi] + f * (p1Param - p2Param));
			++pi;
		}
	}
	else if (options.mutateType == 2) {
		unsigned int pi = 0;
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0) {
				bestParticle = fitVal->second;
			}
		}
		vector<double> bestParamSet = particleCurrParamSets_.at(bestParticle);
		for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			int p4 = 0;
			while (p1 == p2 || p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4 || p3 == p4) {
				p1 = unif(randNumEngine);
				p2 = unif(randNumEngine);
				p3 = unif(randNumEngine);
				p4 = unif(randNumEngine);
			}

			double p1Param = particleCurrParamSets_.at(p1)[pi];
			double p2Param = particleCurrParamSets_.at(p2)[pi];
			double p3Param = particleCurrParamSets_.at(p3)[pi];
			double p4Param = particleCurrParamSets_.at(p4)[pi];

			mutatedParams.push_back(bestParamSet[pi] + f * (p1Param - p2Param) + f * (p3Param - p4Param));
			++pi;
		}
	}
	else if (options.mutateType == 3) {
		unsigned int pi = 0;
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0) {
				bestParticle = fitVal->second;
			}
		}
		vector<double> bestParamSet = particleCurrParamSets_.at(bestParticle);
		for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			while (p1 == p2 || p1 == p3 || p2 == p3) {
				p1 = unif(randNumEngine);
				p2 = unif(randNumEngine);
				p3 = unif(randNumEngine);
			}

			double p1Param = particleCurrParamSets_.at(p1)[pi];
			double p2Param = particleCurrParamSets_.at(p2)[pi];
			double p3Param = particleCurrParamSets_.at(p3)[pi];

			mutatedParams.push_back(p1Param + f * (p2Param - p3Param));
			++pi;
		}
	}
	else if (options.mutateType == 4) {
		unsigned int pi = 0;
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0) {
				bestParticle = fitVal->second;
			}
		}
		vector<double> bestParamSet = particleCurrParamSets_.at(bestParticle);
		for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			int p4 = 0;
			int p5 = 0;
			// This (and above instances) should be optimized
			while (p1 == p2 || p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4 || p3 == p4 || p1 == p5 || p2 == p5 || p3 == p5 || p4 == p5) {
				p1 = unif(randNumEngine);
				p2 = unif(randNumEngine);
				p3 = unif(randNumEngine);
				p4 = unif(randNumEngine);
				p5 = unif(randNumEngine);
			}

			double p1Param = particleCurrParamSets_.at(p1)[pi];
			double p2Param = particleCurrParamSets_.at(p2)[pi];
			double p3Param = particleCurrParamSets_.at(p3)[pi];
			double p4Param = particleCurrParamSets_.at(p4)[pi];
			double p5Param = particleCurrParamSets_.at(p5)[pi];

			mutatedParams.push_back(p1Param + f * (p2Param - p3Param) + f * (p4Param - p5Param));
			++pi;
		}
	}

	return mutatedParams;
}

vector<double> Swarm::crossoverParticleDE(unsigned int particle, vector<double> mutationSet, float crossoverRate) {
	if (options.verbosity >= 3) {
		cout << "Crossing over particle " << particle << endl;
	}

	// Uses binomial crossover
	boost::random::uniform_int_distribution<int> intDist(0, options.model->getNumFreeParams() - 1);
	unsigned int randParam = intDist(randNumEngine);
	boost::random::uniform_real_distribution<float> floatDist(0, 1);

	float cr;
	if (crossoverRate) {
		cr = crossoverRate;
	}
	else {
		cr = options.cr;
	}

	//cout << "crossing over particle " << particle << ". randParam is " << randParam << endl;
	vector<double> newParamSet;

	for (unsigned int p = 0; p < options.model->getNumFreeParams(); ++p) {
		float rand = floatDist(randNumEngine);
		//cout << "p: " << p << ". rand: " << rand << endl;

		if (rand < cr || p == randParam) {
			newParamSet.push_back(mutationSet[p]);
			//cout << "pushing back " << mutationSet[p] << " from mutation set" << endl;
		}
		else {
			newParamSet.push_back(particleCurrParamSets_.at(particle)[p]);
			//cout << "pushing back " << particleCurrParamSets_.at(particle)[p] << " from current set" << endl;
		}
	}

	return newParamSet;
}

void Swarm::sendMigrationSetDE(unsigned int island, vector<vector<unsigned int>> islandTopology, map<unsigned int, vector<vector<double>>> &migrationSets) {
	if (options.verbosity >= 3) {
		cout << "Sending migration set from island " << island << endl;
	}

	// Generate a random number between 1 and the number of neighbors
	// in the island topology
	boost::random::uniform_int_distribution<int> unif(0, islandTopology[island].size() - 1);

	// Fill a vector with particles from this island, starting with best fits and ending with worst
	vector<unsigned int> particlesToSend;
	for (auto fitIt = particleBestFitsByFit_.begin(); fitIt != particleBestFitsByFit_.end(); ++fitIt) {
		if (particleToIsland_.at(fitIt->second) == island) {
			particlesToSend.push_back(fitIt->second);
			if (particlesToSend.size() == options.numToMigrate) {
				break;
			}
		}
	}

	// Choose the index of our receiver
	unsigned int receivingIslandIndex = unif(randNumEngine);

	// Loop through particles to send
	vector<double> migrationSet;
	//for (auto particle = particlesToSend.begin(); particle != particlesToSend.end(); ++particle) {
	for (unsigned int i = 0; i < options.numToMigrate; ++i) {
		//cout << "sending set " << i + 1 << " from particle " << particlesToSend[i] << endl;
		// Fill migrationSet with this particles current parameters
		for (auto param = particleCurrParamSets_.at(particlesToSend[i]).begin(); param != particleCurrParamSets_.at(particlesToSend[i]).end(); ++param) {
			migrationSet.push_back(*param);
		}
		// Add that migration set to the list of the receiving island's
		// migration sets
		migrationSets[islandTopology[island][receivingIslandIndex]].push_back(migrationSet);
		migrationSet.clear();
		//cout << island << " sending migration set to island " << islandTopology[island][receivingIslandIndex] << endl;
	}
}

void Swarm::recvMigrationSetDE(unsigned int island, map<unsigned int, vector<vector<double>>> &migrationSets) {
	// If we have any sets to receive..
	if (migrationSets.find(island) != migrationSets.end() && migrationSets.at(island).size()) {

		if (options.verbosity >= 3) {
			cout << "Receiving " << migrationSets.at(island).size() << " migration sets for island " << island << endl;
		}

		// Create a list of particles from the receiving island that will
		// receive the migration set
		vector<unsigned int> particlesToRecv;
		for (map<double, unsigned int>::reverse_iterator fitIt = particleBestFitsByFit_.rbegin(); fitIt != particleBestFitsByFit_.rend(); ++fitIt) {
			cout << fitIt->second << endl;
			if (particleToIsland_.at(fitIt->second) == island) {
				particlesToRecv.push_back(fitIt->second);
			}
		}

		// Iterates through the list of receiving particles
		auto recvIt = particlesToRecv.begin();
		cout << 999 << endl;
		// For each migration set destined for this island
		unsigned int replacementCounter = 0;
		for (auto migrationSet = migrationSets[island].begin(); migrationSet != migrationSets[island].end();) {
			cout << "Replacing " << migrationSet->size() << " params for particle " << *recvIt << " with fit value of " << particleBestFits_.at(*recvIt) << endl;
			vector<string> paramStr;
			unsigned int i = 0;
			// For each parameter in the migration set
			for (auto param = migrationSet->begin(); param != migrationSet->end(); ++param) {
				// Insert parameter to the currentParamSets tracker
				cout << "param: " << *param << endl;
				particleCurrParamSets_.at(*recvIt)[i] = *param;
				paramStr.push_back(to_string(static_cast<long double>(*param)));
				++i;
			}
			// Send the parameters to the particle
			swarmComm->sendToSwarm(0, *recvIt, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramStr);
			// Next migration set, choose a new receiver from the island
			++recvIt;
			++replacementCounter;
			migrationSet = migrationSets[island].erase(migrationSet);

			if (replacementCounter >= options.swarmSize / options.numIslands) {
				break;
			}
		}
	}
}

void Swarm::runSGA() {
	if (options.verbosity >= 3) {
		cout << "Running a synchronous genetic fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + "1")) {
		string createDirCmd = "mkdir " + options.jobOutputDir + "1";
		if (runCommand(createDirCmd) != 0) {
			outputError("Error: Couldn't create first generation output directory with command: " + createDirCmd + ". Quitting.");
		}
	}

	bool stopCriteria = false;
	while (!stopCriteria){
		runGeneration();
		stopCriteria = checkStopCriteria();
		saveSwarmState();

		string currentDirectory = options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration));
		if (options.deleteOldFiles) {
			cleanupFiles(currentDirectory.c_str());
		}

		if (!stopCriteria) {
			string outputPath = options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration - 1 )) + "_summary.txt";
			outputRunSummary(outputPath);
			breedGenerationGA();
		}
	}
}

void Swarm::runSPSO() {
	if (options.verbosity >= 3) {
		cout << "Running a synchronous PSO fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + "1")) {
		string createDirCmd = "mkdir " + options.jobOutputDir + "1";
		runCommand(createDirCmd);
	}

	if (!resumingSavedSwarm) {
		// Generate all particles that will be present in the swarm
		populationTopology_ = generateTopology(options.swarmSize);
	}

	saveSwarmState();

	vector<unsigned int> finishedParticles;
	bool stopCriteria = false;
	while (!stopCriteria) {

		if (options.verbosity >= 1) {
			cout << "Launching flock " << currentGeneration << endl;
		}

		unsigned int p = 1;
		while (finishedParticles.size() < options.swarmSize && p <= options.swarmSize) {
			if (runningParticles_.size() < options.parallelCount) {
				launchParticle(p++);
			}

			usleep(250000);

			vector<unsigned int> currFinishedParticles;
			currFinishedParticles = checkMasterMessages();

			if (currFinishedParticles.size()) {
				finishedParticles.insert(finishedParticles.end(), currFinishedParticles.begin(), currFinishedParticles.end());
			}
		}

		// TODO: With the current structure and ordering of checkStopCriteria()
		// and processParticlesPSO(), we wind up with an extra directory containing
		// models after fitting is finished.

		// Process particles
		processParticlesPSO(finishedParticles, true);

		// Check for stop criteria
		stopCriteria = checkStopCriteria();

		if (!stopCriteria) {
			string outputPath = options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration)) + "_summary.txt";
			outputRunSummary(outputPath);
			++currentGeneration;
		}

		// Save swarm state
		saveSwarmState();

		// Clear our finished particles for the next round
		finishedParticles.clear();
	}
}

void Swarm::runSDE() {
	if (options.verbosity >= 3) {
		cout << "Running a synchronous DE fit" << endl;
	}

	// Create an output directory for each island
	if (!checkIfFileExists(options.jobOutputDir + to_string(static_cast<long long int>(1)))) {
		string createDirCmd = "mkdir " + options.jobOutputDir + to_string(static_cast<long long int>(1));
		runCommand(createDirCmd);
	}

	if (options.verbosity >= 3) {
		cout << "Generating island topology" << endl;
	}
	vector<vector<unsigned int>> islandTopology;
	if (!resumingSavedSwarm) {
		// Generate all particles that will be present in the swarm
		islandTopology = generateTopology(options.numIslands);
		// Map all particles and islands
		unsigned int particlesPerIsland = options.swarmSize / options.numIslands;
		unsigned int currParticle = 1;
		for (unsigned int i = 1; i <= options.numIslands; ++i) {
			for (unsigned int p = 0; p < particlesPerIsland; ++p) {
				particleToIsland_[currParticle] = i;
				islandToParticle_[i][p] = currParticle;
				cout << i << " " << p << " " << islandToParticle_[i][p] << endl;
				++currParticle;
			}
		}
	}

	saveSwarmState();

	if (options.verbosity >= 3) {
		cout << "Launching first generation" << endl;
	}

	unordered_map<unsigned int, vector<double>> finishedParticles;
	bool stopCriteria = false;
	unsigned int numFinishedParticles = 0;
	vector<unsigned int> islandFinishedParticles(options.numIslands + 1, 0);
	map<unsigned int, vector<vector<double>>> migrationSets;

	bool trialLoop = false;
	while(!stopCriteria) {

		// Generation loop
		unsigned int p = 1;
		while (numFinishedParticles < options.swarmSize) {
			// Launch next particle
			if (runningParticles_.size() < options.parallelCount && p <= options.swarmSize) {
				launchParticle(p++, trialLoop);
			}

			unordered_map<unsigned int, vector<double>> finished = checkMasterMessagesDE();

			if (finished.size()) {
				finishedParticles.insert(finished.begin(), finished.end());
				numFinishedParticles += finished.size();
			}
		}

		// For each finished particle
		for (auto particle = finishedParticles.begin(); particle != finishedParticles.end(); ++particle) {
			// Increment the counter that tracks the number of particles finished
			// for a given island
			islandFinishedParticles[particleToIsland_.at(particle->first)] += 1;

			if (options.verbosity >= 3) {
				cout << "Particle " << particle->first << " in island " << particleToIsland_.at(particle->first) << " finished. total finished is " << numFinishedParticles << endl;
			}

			// If we're in a trial loop..
			if (trialLoop && particle->second.size()) {
				// Create the beginning of our param string for use in summary and result outputs
				string paramsString = "gen" + to_string(static_cast<long long int>(currentGeneration)) + "perm" + to_string(static_cast<long long int>(particle->first)) + " ";

				// Increment the flight counter only at the end of a trial
				++flightCounter_;

				bool replaceParams = false;

				// If the trail sim fit is better than our current best fit, we need
				// to replace the old fit value with the new one
				if (particle->second[0] < particleBestFits_.at(particle->first)) {
					particleBestFits_[particle->first] = particle->second[0];
					insertKeyByValue(particleBestFitsByFit_, particle->second[0], particle->first);

					replaceParams = true;
				}

				// Loop through parameters sent to us by the particle
				unsigned int i = 0;
				for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {

					// If we need to replace old params with new ones (because fit calc
					// was better in the trial)...
					if (replaceParams) {
						// Store the parameter in the global param set and
						// concatenate the param string for output
						particleCurrParamSets_[particle->first][i] = *param;
						paramsString += to_string(static_cast<long double>(*param)) + " ";
						//cout << "added " << paramsString << endl;
					}
					// Otherwise, do nothing except use the current/old param
					// set in the output
					else {
						paramsString += to_string(static_cast<long double>(particleCurrParamSets_[particle->first][i])) + " ";
						//cout << "added: " << paramsString << endl;
					}
					++i;
				}

				// Insert fit values and parameters into our global trackers
				allGenFits.insert(pair<double, string>(particleBestFits_[particle->first], paramsString));
				swarmBestFits_.insert(pair<double, unsigned int>(particle->second[0], particle->first));
			}
			// If we're not in a trial loop...
			else if (particle->second.size()){
				// Replace old best fit value with the new one
				particleBestFits_[particle->first] = particle->second[0];
				insertKeyByValue(particleBestFitsByFit_, particle->second[0], particle->first);

				// Loop through the params sent to us by the particle
				unsigned int i = 0;
				for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
					// Replace current params with the new ones and add to
					// our output string
					particleCurrParamSets_[particle->first][i] = *param;
					++i;
				}
			}
		}

		finishedParticles.clear();
		numFinishedParticles = 0;

		// Done processing params. Now we need to...
		// Check each island to see if it has completed its generation
		for (unsigned int island = 1; island <= options.numIslands; ++island) {
			// If the number finished in this island is equal to the total size of the island
			if (islandFinishedParticles.at(island) == (options.swarmSize / options.numIslands)) {
				cout << "Island " << island << " finished" << endl;
				// Loop through the particles in the island
				for (auto particle = islandToParticle_.at(island).begin(); particle != islandToParticle_.at(island).end(); ++particle) {

					// Will hold mutated and crossed over parameter sets
					vector<double> newParamSet;
					// If we're not in a trial loop, we should mutate and crossover
					if (trialLoop == false) {
						// Create a mutation set for the particle
						newParamSet = mutateParticleDE(*particle);

						// Run crossover for the particle
						newParamSet = crossoverParticleDE(*particle, newParamSet);
					}
					// If we're in a main loop, we just send the current parameter
					// set back to the particle
					else {
						newParamSet = particleCurrParamSets_.at(*particle);
					}

					// Convert our param set to string for sending
					vector<string> paramVecStr;
					for (auto param = newParamSet.begin(); param != newParamSet.end(); ++param) {
						paramVecStr.push_back(to_string(static_cast<long double>(*param)));
					}

					// Send new param sets to particles for next generation
					swarmComm->sendToSwarm(0, *particle, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramVecStr);

					// Empty our finished particle counter for the next loop
					islandFinishedParticles[island] = 0;
				}
			}
		}

		// Switch from trial/main loops
		if (trialLoop == false) {
			trialLoop = true;

			cout << "Switching to trial loop" << endl;
		}
		else {
			// Only want to check stop criteria at end of trial set
			stopCriteria = checkStopCriteria();
			if (!stopCriteria) {
				// Send/receive migration sets
				if (currentGeneration % options.migrationFrequency == 0) {

					for (unsigned int island = 1; island <= options.numIslands; ++island) {
						sendMigrationSetDE(island, islandTopology, migrationSets);
					}

					for (unsigned int island = 1; island <= options.numIslands; ++island) {
						recvMigrationSetDE(island, migrationSets);
					}
				}

				string outputPath = options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration)) + "_summary.txt";
				outputRunSummary(outputPath);
				++currentGeneration;
				cout << "Switching to main loop" << endl;
				trialLoop = false;
			}
		}
	}
}

void Swarm::runAGA() {
	if (options.verbosity >= 3) {
		cout << "Running an asynchronous genetic fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + "1")) {
		string createDirCmd = "mkdir " + options.jobOutputDir + "1";
		if (runCommand(createDirCmd) != 0) {
			outputError("Error: Couldn't create first generation output directory with command: " + createDirCmd + ". Quitting.");
		}
	}

	// Holds finished particles
	vector<unsigned int> finishedParticles;

	// Run the first generation
	runGeneration();

	// Save swarm state
	saveSwarmState();

	// Breed generation
	breedGenerationGA();

	// Re-launch the particles in to the Swarm proper
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		launchParticle(p);
	}

	bool stopCriteria = false;
	while (!stopCriteria){
		usleep(250000);

		// Check for any messages from particles and store finished particles
		finishedParticles = checkMasterMessages();

		// Check stop criteria
		stopCriteria = checkStopCriteria();

		// Cleanup old files if needed
		string currentDirectory = options.jobOutputDir + to_string(static_cast<long long int>(currentGeneration));
		if (options.deleteOldFiles) {
			cleanupFiles(currentDirectory.c_str());
		}

		// If we haven't reached stop criteria, breed and re-launch finished particles
		if (!stopCriteria && finishedParticles.size()) {
			// Get new params for any finished particles
			breedGenerationGA(finishedParticles);

			// Re-launch any finished particles
			for (auto pID = finishedParticles.begin(); pID != finishedParticles.end(); ++pID) {
				launchParticle(*pID);
			}
		}
	}
}

void Swarm::runAPSO() {
	if (options.verbosity >= 3) {
		cout << "Running an asynchronous PSO fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + "1")) {
		string createDirCmd = "mkdir " + options.jobOutputDir + "1";
		runCommand(createDirCmd);
	}

	if (!resumingSavedSwarm) {
		// Generate all particles that will be present in the swarm
		populationTopology_ = generateTopology(options.swarmSize);

		// Run first flight
		initPSOswarm();
	}

	saveSwarmState();

	if (options.verbosity >= 3) {
		cout << "Re-launching swarm" << endl;
	}

	// Re-launch the particles in to the Swarm proper
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		launchParticle(p);
	}

	if (options.verbosity >= 3) {
		cout << "Entering main swarm loop" << endl;
	}

	vector<unsigned int> finishedParticles;
	bool stopCriteria = false;
	unsigned int clusterCheckCounter = 0;

	while (!stopCriteria) {
		usleep(250000);
		finishedParticles = checkMasterMessages();

		if (finishedParticles.size()) {
			// Process finished particles
			processParticlesPSO(finishedParticles, true);

			// Check for stop criteria
			stopCriteria = checkStopCriteria();

			// Save swarm state
			saveSwarmState();
		}

		// Only check cluster queue every minute
		++clusterCheckCounter;
		if (clusterCheckCounter >= 240) {
			//checkClusterQueue();
			clusterCheckCounter = 0;
		}

		if (!stopCriteria) {
			// Re-launch the particles in to the Swarm
			for (auto particle = finishedParticles.begin(); particle != finishedParticles.end(); ++particle) {
				launchParticle(*particle);
			}
		}
	}
}

void Swarm::runADE() {
	if (options.verbosity >= 3) {
		cout << "Running an asynchronous DE fit" << endl;
	}

	// Create an output directory for each island
	if (!checkIfFileExists(options.jobOutputDir + to_string(static_cast<long long int>(1)))) {
		string createDirCmd = "mkdir " + options.jobOutputDir + to_string(static_cast<long long int>(1));
		runCommand(createDirCmd);
	}

	if (options.verbosity >= 3) {
		cout << "Generating island topology" << endl;
	}

	vector<vector<unsigned int>> islandTopology;
	vector<bool> islandIsTrial = vector<bool>(options.numIslands + 1, false);

	if (!resumingSavedSwarm) {
		// Generate all particles that will be present in the swarm
		islandTopology = generateTopology(options.numIslands);
		// Map all particles and islands
		unsigned int particlesPerIsland = options.swarmSize / options.numIslands;
		unsigned int currParticle = 1;
		for (unsigned int i = 1; i <= options.numIslands; ++i) {
			for (unsigned int p = 0; p < particlesPerIsland; ++p) {
				particleToIsland_[currParticle] = i;
				islandToParticle_[i][p] = currParticle;
				cout << i << " " << p << " " << islandToParticle_[i][p] << endl;
				++currParticle;
			}
		}
	}

	saveSwarmState();

	if (options.verbosity >= 3) {
		cout << "Launching first generation" << endl;
	}

	unordered_map<unsigned int, vector<double>> finishedParticles;
	bool stopCriteria = false;
	vector<unsigned int> islandFinishedParticles(options.numIslands + 1, 0);
	vector<unsigned int> islandGenerationCounter(options.numIslands + 1, 1);
	map<unsigned int, vector<vector<double>>> migrationSets;

	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		launchParticle(p);
	}

	while(!stopCriteria) {

		finishedParticles = checkMasterMessagesDE();

		if (finishedParticles.size()) {
			// Increment the counter that tracks the number of particles finished
			// for a given island
			for (auto particle = finishedParticles.begin(); particle != finishedParticles.end(); ++particle) {
				islandFinishedParticles[particleToIsland_.at(particle->first)] += 1;
				cout << particle->first << " in island " << particleToIsland_.at(particle->first) << " finished." << endl;

				// TODO: Need to put this (and synch version) in a function
				// Let's process params and update any fit values
				string paramsString;
				paramsString = "gen" + to_string(static_cast<long long int>(currentGeneration)) + "perm" + to_string(static_cast<long long int>(particle->first)) + " ";
				if (islandIsTrial[particleToIsland_[particle->first]]) {

					if (flightCounter_ && flightCounter_ % options.outputEvery == 0) {
						string outputPath = options.jobOutputDir + to_string(static_cast<long long int>(flightCounter_)) + "_summary.txt";
						outputRunSummary(outputPath);
					}

					++flightCounter_;

					paramsString = to_string(static_cast<long long int>(flightCounter_)) + " ";

					bool replaceParams = false;

					if (particle->second[0] < particleBestFits_.at(particle->first) || particleBestFits_.at(particle->first) == 0) {
						particleBestFits_[particle->first] = particle->second[0];
						insertKeyByValue(particleBestFitsByFit_, particle->second[0], particle->first);

						replaceParams = true;
					}

					unsigned int i = 0;
					for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
						if (replaceParams) {
							//cout << "stored " << stod(*m) << " for particle " << pID << " as " << particleCurrParamSets_[pID][i] << endl;
							particleCurrParamSets_[particle->first][i] = *param;
							paramsString += to_string(static_cast<long double>(*param)) + " ";
						}
						else {
							paramsString += to_string(static_cast<long double>(particleCurrParamSets_[particle->first][i])) + " ";
						}
						++i;
					}

					allGenFits.insert(pair<double, string>(particleBestFits_[particle->first], paramsString));
					swarmBestFits_.insert(pair<double, unsigned int>(particle->second[0], particle->first));
				}
				else {
					particleBestFits_[particle->first] = particle->second[0];
					insertKeyByValue(particleBestFitsByFit_, particle->second[0], particle->first);

					unsigned int i = 0;
					for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
						particleCurrParamSets_[particle->first][i] = *param;
						++i;
					}
				}
			}

			// Check each island to see if it has completed its generation
			for (unsigned int island = 1; island <= options.numIslands; ++island) {
				// If the number finished in this island is equal to the total size of the island
				if (islandFinishedParticles.at(island) == (options.swarmSize / options.numIslands)) {
					cout << "Island " << island << " finished" << endl;

					for (auto particle = islandToParticle_.at(island).begin(); particle != islandToParticle_.at(island).end(); ++particle) {
						vector<double> newParamSet;

						if (islandIsTrial[island] == false) {
							// Create a mutation set for the particle
							newParamSet = mutateParticleDE(*particle);

							// Run crossover for the particle
							newParamSet = crossoverParticleDE(*particle, newParamSet);
						}
						else { // Finishing the trial set
							newParamSet = particleCurrParamSets_.at(*particle);
						}

						// Convert our param set to string for sending
						vector<string> paramVecStr;
						//for (auto param : particleCurrParamSets_.at(*particle)) {
						for (auto param = newParamSet.begin(); param != newParamSet.end(); ++param) {
							paramVecStr.push_back(to_string(static_cast<long double>(*param)));
						}

						// Send new param sets to particles for next generation
						swarmComm->sendToSwarm(0, *particle, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramVecStr);
						launchParticle(*particle);

						// Empty our finished particle counter for the next loop
						islandFinishedParticles[island] = 0;
					}

					if (islandIsTrial[island]) {
						stopCriteria = checkStopCriteria();

						if (!stopCriteria) {
							islandIsTrial[island] = false;

							if (islandGenerationCounter[island] % options.migrationFrequency == 0) {
								sendMigrationSetDE(island, islandTopology, migrationSets);
							}

							recvMigrationSetDE(island, migrationSets);
							++islandGenerationCounter[island];
						}
					}
					else {
						islandIsTrial[island] = true;
					}
				}
			}
		}
	}
}

void Swarm::runASA() {
	if (options.verbosity >= 3) {
		cout << "Running an asynchronous PSADE fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + "1")) {
		string createDirCmd = "mkdir " + options.jobOutputDir + "1";
		if (runCommand(createDirCmd) != 0) {
			outputError("Error: Couldn't create first generation output directory with command: " + createDirCmd + ". Quitting.");
		}
	}

	// Initialize and/or fill our parameter vectors
	vector<unsigned int> isLocal = vector<unsigned int>(options.swarmSize, 0);
	vector<double> particleTemps = vector<double>(options.swarmSize + 1, 0);
	vector<double> particleRadii = vector<double>(options.swarmSize + 1, 0);
	vector<float> particleFs = generateParticleFs();
	vector<float> particleCRs = generateParticleCRs();
	vector<float> cpuToParticle = vector<float>(options.swarmSize + 1, 0);
	vector<vector<float>> trialParams = vector<vector<float>>(options.swarmSize + 1, vector<float>(2, 0));

	// Launch the initialization population and map
	// CPUs to particles
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		launchParticle(p);
		cpuToParticle[p] = p;
	}

	unordered_map<unsigned int, vector<double>> finishedParticles;
	unordered_map<unsigned int, vector<double>> initFinishedParticles;
	unsigned int numFinishedParticles = 0;

	// Wait for initialization population to finish simulations
	while (numFinishedParticles < options.swarmSize) {
		initFinishedParticles = checkMasterMessagesDE();
		// Save finished particles in another map, we still need to process them in the main loop later on
		if (initFinishedParticles.size()) {
			finishedParticles.insert(initFinishedParticles.begin(), initFinishedParticles.end());
			numFinishedParticles += initFinishedParticles.size();
		}

		// Process each finished particle
		for (auto particle = initFinishedParticles.begin(); particle != initFinishedParticles.end(); ++particle) {
			++flightCounter_;
			cout << "particle " << cpuToParticle[particle->first] << " on cpu " << particle->first << " finished" << endl;
			string paramString = to_string(static_cast<long long int>(flightCounter_)) + " ";

			// Update parameter values
			unsigned int i = 0;
			for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
				particleCurrParamSets_[cpuToParticle[particle->first]][i++] = *param;
				cout << "updating particle " << cpuToParticle[particle->first] << " to " << *param << " at " << i << endl;
				paramString += to_string(static_cast<long long int>(*param)) + " ";
			}

			// Update fit calcs
			particleBestFits_[cpuToParticle[particle->first]] = particle->second[0];
			insertKeyByValue(particleBestFitsByFit_, particle->second[0], cpuToParticle[particle->first]);
			allGenFits.insert(pair<double, string>(particle->second[0], paramString));
			cout << "updating particle " << cpuToParticle[particle->first] << " calc to " << particle->second[0] << endl;
			//cout << "count: " << particleBestFitsByFit_.count(particle->second[0]) << endl;

			// Update F and CR
			trialParams[cpuToParticle[particle->first]][0] = particleFs[cpuToParticle[particle->first]];
			trialParams[cpuToParticle[particle->first]][1] = particleCRs[cpuToParticle[particle->first]];
		}
	}

	// Generate our initial temps and radii
	particleTemps = generateParticleTemps();
	particleRadii = generateParticleRadii();

	// Main loop
	cout << "entering main loop" << endl;
	bool stopCriteria = false;
	while (!stopCriteria) {

		// Check to make sure we aren't processing init swarm before checking for new
		if (!finishedParticles.size()) {
			finishedParticles = checkMasterMessagesDE();
		}

		// Process any finished particles
		for (auto particle = finishedParticles.begin(); particle != finishedParticles.end(); ++particle) {
			cout << "particle " << cpuToParticle[particle->first] << " on cpu " << particle->first << " finished" << endl;

			string paramString = to_string(static_cast<long long int>(flightCounter_ + 1)) + " ";

			// If the particle finished their local search
			if (isLocal[particle->first]) {
				cout << "particle just finished local search" << endl;

				cout << "greedy acceptance of fit of " << particle->second[0] << endl;
				// Greedy acceptance
				particleBestFits_[cpuToParticle[particle->first]] = particle->second[0];
				insertKeyByValue(particleBestFitsByFit_, particle->second[0], cpuToParticle[particle->first]);

				unsigned int i = 0;
				for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
					particleCurrParamSets_[cpuToParticle[particle->first]][i++] = *param;
					cout << "updating particle " << cpuToParticle[particle->first] << " to " << *param << " at " << i << endl;
					paramString += to_string(static_cast<long long int>(*param)) + " ";
				}

				allGenFits.insert(pair<double, string>(particle->second[0], paramString));

				isLocal[cpuToParticle[particle->first]] = false;
			}
			else {
				// If we pass metropolis selection, store the params for assignment to receiver
				if (metropolisSelection(particle->first, particle->second[0], particleTemps[particle->first])) {
					cout << "particle was accepted through metropolis selection" << endl;
					unsigned int i = 0;

					for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
						particleCurrParamSets_[cpuToParticle[particle->first]][i++] = *param;
						cout << "updating particle " << cpuToParticle[particle->first] << " to " << *param << " at " << i << endl;
						paramString += to_string(static_cast<long long int>(*param)) + " ";
					}

					// Update fitness
					particleBestFits_[cpuToParticle[particle->first]] = particle->second[0];
					insertKeyByValue(particleBestFitsByFit_, particle->second[0], cpuToParticle[particle->first]);
					allGenFits.insert(pair<double, string>(particle->second[0], paramString));

					cout << "accepting params and fit of " << particle->second[0] << endl;
					// Update CR and F

					cout << "updating F to " << trialParams[cpuToParticle[particle->first]][0] << " and CR to " << trialParams[cpuToParticle[particle->first]][1] << endl;
					particleFs[cpuToParticle[particle->first]] = trialParams[cpuToParticle[particle->first]][0];
					particleCRs[cpuToParticle[particle->first]] = trialParams[cpuToParticle[particle->first]][1];
				}

				// Do a local search at some probability
				// Or if particle is best in swarm
				if (options.localSearchProbability > ((float)rand() / (float)RAND_MAX) || cpuToParticle[particle->first] == particleBestFitsByFit_.begin()->second) {
					cout << "going to do a local search" << endl;
					isLocal[cpuToParticle[particle->first]] = true;
				}
			}

			if (flightCounter_ && flightCounter_ % options.outputEvery == 0) {
				string outputPath = options.jobOutputDir + to_string(static_cast<long long int>(flightCounter_)) + "_summary.txt";
				outputRunSummary(outputPath);
				//cout << "fc is " << flightCounter_ << ", outputting" << endl;
			}

			++flightCounter_;

			unsigned int receiver = rand() % options.swarmSize + 1;
			cout << "chose receiver of " << receiver << endl;
			vector<double> newParams;

			// If we aren't going to do a local search, generate a new trial point
			if (!isLocal[cpuToParticle[particle->first]] && particleBestFitsByFit_.begin()->second != receiver) {
				// Pick the controller
				unsigned int controller = pickWeightedSA();
				cout << "not doing local search. chose controller of " << controller << endl;
				// Generate a new trial vector
				newParams = generateTrialPointSA(controller, receiver, particleRadii, particleCRs, particleFs, trialParams);
				cout << "generated new trial points" << endl;

				// Make sure we know this CPU will be running this Receiver
				cpuToParticle[particle->first] = receiver;
				cout << "setting cpu " << particle->first << " as receiver " << receiver << endl;

				vector<string> newParamsStr;
				for (auto param = newParams.begin(); param != newParams.end(); ++param) {
					newParamsStr.push_back(to_string(static_cast<long double>(*param)));
				}

				cout << "running " << cpuToParticle[particle->first] << " with new params on cpu " << particle->first << endl;
				// Send the params to the CPU
				swarmComm->sendToSwarm(0, particle->first, SEND_FINAL_PARAMS_TO_PARTICLE, false, newParamsStr);
				launchParticle(particle->first);
			}
			else {
				cout << "running local search on cpu " << particle->first << endl;
				// Do local search
				runNelderMead(receiver, particle->first);
			}
		}
		finishedParticles.clear();
	}
}

void Swarm::runSSA() {
	if (options.verbosity >= 3) {
		cout << "Running an asynchronous PSADE fit" << endl;
	}

	// Initialize and/or fill our parameter vectors
	vector<unsigned int> isLocal = vector<unsigned int>(options.swarmSize, 0);
	vector<float> particleTemps = vector<float>(options.swarmSize + 1, 0);
	vector<float> particleRadii = vector<float>(options.swarmSize + 1, 0);
	vector<float> particleFs = generateParticleFs();
	vector<float> particleCRs = generateParticleCRs();
	vector<float> cpuToParticle = vector<float>(options.swarmSize + 1, 0);
	vector<vector<float>> trialParams = vector<vector<float>>(options.swarmSize + 1, vector<float>(2, 0));

	// Launch the initialization population and map
	// CPUs to particles
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		launchParticle(p);
		cpuToParticle[p] = p;
	}

	unordered_map<unsigned int, vector<double>> finishedParticles;
	unsigned int numFinishedParticles = 0;

	// Wait for initialization population to finish simulations
	while (numFinishedParticles < options.swarmSize) {


	}
}

vector<double> Swarm::generateParticleTemps() {
	// Eq 4 in Olensek et al

	// maxTemp is the difference between highest and lowest fits
	double maxTemp = particleBestFitsByFit_.rbegin()->first - particleBestFitsByFit_.begin()->first;
	//cout << "maxtemp is " << maxTemp << endl;
	//cout << "minTemp is " << options.minTemp << endl;
	//cout << "ss is " << options.swarmSize - 1 << endl;
	double ct = (1 / ((double)options.swarmSize - 1)) * log(maxTemp / options.minTemp);


	//cout << fixed << setprecision(10) << "ct is " << ct << endl;

	vector<double> temps = vector<double>(options.swarmSize + 1, 0);
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		temps[p] = exp(0 - (ct * (p - 1)));
		//cout << fixed << setprecision(10) << "temp of " << p << " is " << temps[p] << endl;
	}

	return temps;
}

vector<double> Swarm::generateParticleRadii() {
	// Eq 4 in Olensek et al

	float rMax = 1;
	float cr = (1 / ((double)options.swarmSize - 1)) * log(rMax / options.minRadius);

	vector<double> radii = vector<double>(options.swarmSize + 1, 0);
	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		radii[p] = exp(0 - (cr * (p - 1)));
		//cout << "radius of " << p << " is " << radii[p] << endl;
	}

	return radii;
}

vector<float> Swarm::generateParticleFs() {
	float min = 0.5;
	float max = 1.5;

	vector<float> Fs = vector<float>(options.swarmSize + 1, 0);

	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		float r = (float)rand() / (float)RAND_MAX;
		float diff = max - min;

		Fs[p] = (r * diff) + min;
	}

	return Fs;
}

vector<float> Swarm::generateParticleCRs() {
	float min = 0.1;
	float max = 0.9;

	vector<float> CRs = vector<float>(options.swarmSize + 1, 0);

	for (unsigned int p = 1; p <= options.swarmSize; ++p) {
		float r = (float)rand() / (float)RAND_MAX;
		float diff = max - min;

		CRs[p] = (r * diff) + min;
	}

	return CRs;
}

unsigned int Swarm::pickWeightedSA() {

	//Eq 6 Olensek et al

	// Calculate the reimann sum
	float sum = 0;
	for (unsigned int r = 1; r <= options.swarmSize; ++r) {
		sum += exp(0 - r);
	}

	// Fill vector with probabilities of particles being selected
	map<float, unsigned int> probabilities;
	unsigned int r = 0;
	for (auto particle = particleBestFitsByFit_.begin(); particle != particleBestFitsByFit_.end(); ++particle) {
		probabilities[exp(0 - ++r) / sum] = particle->second;
	}

	float rnd = ((float) rand() / (RAND_MAX));

	// Go through and choose the particle
	for (auto probability = probabilities.rbegin(); probability != probabilities.rend(); ++probability) {
		if (rnd < probability->first) {
			return probability->second;
		}

		rnd -= probability->first;
	}

	// Fell off the end..return the worst particle.
	return probabilities.rend()->second;
}

bool Swarm::metropolisSelection(unsigned int particle, double fit, float particleTemp) {
	float p = min <float>(1, exp( (0 - (fit - particleBestFits_.at(particle))) / particleTemp ) );
	float r = (float)rand() / (float)RAND_MAX;

	if (r < p) {
		return true;
	}
	else {
		return false;
	}
}

void Swarm::swapTR(vector<float> particleRadii, vector<float> particleTemps) {
	// Pick two random individuals
	unsigned int r1 = rand() % options.swarmSize + 1;
	unsigned int r2 = rand() % options.swarmSize + 1;

	// Make sure they're different
	while (r1 == r2) {
		r2 = rand() % options.swarmSize + 1;
	}

	// Choose probability and a random number [0,1]
	float p = min <float>(1, exp( (particleBestFits_.at(r1) - particleBestFits_.at(r2) ) * ((1/particleRadii[r1]) - (1/particleRadii[r2]))));
	float r = (float)rand() / (float)RAND_MAX;

	// Swap
	if (r < p) {
		float r1Temp = particleRadii[r1];
		float r1Radius = particleRadii[r1];

		particleTemps[r1] = particleTemps[r2];
		particleRadii[r1] = particleRadii[r2];

		particleTemps[r2] = r1Temp;
		particleRadii[r2] = r1Radius;
	}
}
vector<double> Swarm::generateTrialPointSA(unsigned int controller, unsigned int receiver, vector<double> particleRadii, vector<float>particleCRs, vector<float>particleFs, vector<vector<float>> &trialParams) {
	vector<double> currParams = normalizeParams(particleCurrParamSets_.at(controller));
	float cr;
	float f;

	// Choose either random F and CR or take from receiver
	if ( ((float)rand() / (float)RAND_MAX) < options.randParamsProbability) {
		cout << "choosing random F and CR" << endl;
		float min = 0.5;
		float max = 1.5;

		float r = (float)rand() / (float)RAND_MAX;
		float diff = max - min;
		f = (r * diff) + min;

		min = 0.1;
		max = 0.9;
		r = (float)rand() / (float)RAND_MAX;
		diff = max - min;
		cr = (r * diff) + min;

		//cout << "f is " << f << " and cr is " << cr << endl;
	}
	else {
		f = particleFs[receiver];
		cr = particleCRs[receiver];
		//cout << "using receiver f and cr of " << f << " and " << cr << endl;
	}

	trialParams[receiver][0] = f;
	trialParams[receiver][1] = cr;

	// Mutation functions pulls parameters from the global list
	currParams = mutateParticleSA(controller, f);
	// Crossover function uses params returned from mutation function
	currParams = crossoverParticleDE(controller, currParams, cr);

	//cout << "mutated and crossed over params" << endl;

	vector<double> newParams;
	for (auto param = currParams.begin(); param != currParams.end(); ++param) {
		float r = (float)rand() / (float)RAND_MAX;
		//cout << "before: " << *param << endl;
		float newParam = *param + particleRadii[controller] * tan(3.141592654 * (r - 0.5));
		//cout << "after: " << newParam << endl;

		if (newParam <= 0 || newParam > 1) {
			newParam = (double)rand() / (double)RAND_MAX;
			//cout << "new: " << newParam << endl;
		}

		newParams.push_back(newParam);
	}

	return deNormalizeParams(newParams);
}

void Swarm::generateBootstrapMaps(vector<map<string, map<string, map<double,unsigned int>>>> &bootStrapMaps) {
	for (unsigned int i = 0; i < options.bootstrap; ++i) {
		//cout << "i: " << i << endl;
		// For each .exp file
		map<string, map<string, map<double, unsigned int>>> maps;
		for (auto dataSet = options.expFiles.begin(); dataSet != options.expFiles.end(); ++dataSet) {
			//cout << "set: " << dataSet->first << endl;
			// For each column
			map<string, map<double, unsigned int>> bsMap;
			for (auto col = dataSet->second->dataCurrent->begin(); col != dataSet->second->dataCurrent->end(); ++col) {
				//cout << "col: " << col->first << endl;
				// Fill the map with 0's
				map<double, unsigned int> colVals;
				for (auto tp = col->second.begin(); tp != col->second.end(); ++tp) {
					//cout << "tp: " << tp->first << endl;
					colVals[tp->first] = 0;
				}

				// Select datapoints at random. If a datapoint is selected,
				// increment it's integer value in the vals map
				for (unsigned int o = 0; o < col->second.size(); ++o) {
					int i = rand() % (col->second.size() - 1);
					auto it = col->second.begin();
					advance(it, i);
					colVals[it->first] += 1;
					//cout << "chose " << i << " (" << it->first << ")" << endl;
				}

				// Insert this column into the map
				bsMap.insert(pair<string, map<double, unsigned int>>(col->first, colVals));
			}
			// Insert this dataset into the maps set
			maps.insert(pair<string, map<string, map<double, unsigned int>>>(dataSet->first, bsMap));
		}
		// Insert this set of maps into the master set
		bootStrapMaps.push_back(maps);
	}

	/*
	unsigned int setCounter = 1;
	for (auto bsSet = bootStrapMaps.begin(); bsSet != bootStrapMaps.end(); ++bsSet) {
		cout << "Set " << setCounter << endl;
		for (auto dataSet = bsSet->begin(); dataSet != bsSet->end(); ++dataSet) {
			cout << "Dataset " << dataSet->first << endl;
			for (auto col = dataSet->second.begin(); col != dataSet->second.end(); ++col) {
				cout << "Col " << col->first << endl;
				for (auto tp = col->second.begin(); tp != col->second.end(); ++tp) {
					cout << tp->first << " " << tp -> second << endl;
				}
			}
		}
	}
	 */
}

void Swarm::runNelderMead(unsigned int it, unsigned int cpu) {
	if (options.verbosity >= 3) {
		cout << "Running Nelder-Meade local search for particle " << it << " on cpu " << cpu << endl;
	}
	// On the master side, we just need to construct the simplex, serialize it,
	// then send it to the CPU..

	// calc -> params
	map<double, vector<double>> simplex {pair<double, vector<double>>(particleBestFits_.at(it), particleCurrParamSets_.at(it))};
	vector<unsigned int> usedParticles {it};

	cout << "constructing simplex" << endl;
	// First fill simplex with param sets (n+1 vertices)
	for (unsigned int i = 0; i < options.model->getNumFreeParams(); ++i) {
		unsigned int particle = rand() % options.swarmSize + 1;
		bool isDuplicate = false;

		// Make sure we aren't selecting duplicates
		do {
			isDuplicate = false;
			for (unsigned int p = 0; p < usedParticles.size(); ++p) {
				if (particle == p || simplex.find(particleBestFits_.at(particle)) != simplex.end()) {
					isDuplicate = true;
					break;
				}
			}
			particle = rand() % options.swarmSize + 1;
		} while (isDuplicate);

		usedParticles.push_back(particle);
		simplex[particleBestFits_.at(particle)] = particleCurrParamSets_.at(particle);
	}

	cout << "serializing simplex" << endl;
	std::stringstream oss;
	boost::archive::text_oarchive oa(oss);
	oa << simplex;
	std::string serializedSimplex(oss.str());

	vector<string> message;
	message.push_back(serializedSimplex);

	cout << "sending simplex to cpu " << cpu << " with starting calc of " << endl;
	runningParticles_.insert(cpu);
	swarmComm->sendToSwarm(0, cpu, BEGIN_NELDER_MEAD, false, message);
}

vector<double> Swarm::normalizeParams(vector<double> oldParams) {
	vector<double> newParams;

	unsigned int d = 0;
	for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
		//cout << "denormalized: " << oldParams[d] << endl;
		//cout << "max: " << param->second->getGenMax() << " min: " << param->second->getGenMin() << " diff: " << (param->second->getGenMax() - param->second->getGenMin()) << endl;
		newParams.push_back(oldParams[d++] / (param->second->getGenMax() - param->second->getGenMin()));
		//cout << "normalized: " << *(newParams.end() - 1) << endl;
	}

	return newParams;
}

vector<double> Swarm::deNormalizeParams(vector<double> oldParams) {
	vector<double> newParams;

	unsigned int d = 0;
	for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
		//cout << "normalized: " << oldParams[d] << endl;
		newParams.push_back( param->second->getGenMin() + oldParams[d++] * (param->second->getGenMax() - param->second->getGenMin()));
		//cout << "denormalized: " << *(newParams.end() - 1) << endl;
	}

	return newParams;
}
