//============================================================================
// Name        : GenFit2.cpp
// Author      : Brandon Thomas
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include "GenFit2.hh"

using namespace std;

int main(int argc, char *argv[]) {

	srand(clock());

	int generation = 0;
	string action;
	string configFile;
	string type;
	int pID;

	// GenFit2 [conf_file]
	if (argc == 2) {
		configFile = argv[1];
	}
	// GenFit2 [action] [conf_file]
	else if (argc == 3) {
		action = argv[1];
		configFile = argv[2];
	}
	// GenFit2 [type] [action] [conf_file]
	// Only used internally -- user should never need to specify type
	else if (argc == 4) {
		type = argv[1];
		action = argv[2];
		configFile = argv[3];
	}
	// GenFit2 [type] [pID] [action] [conf_file]
	else if (argc == 5) {
		type = argv[1];
		pID = atoi(argv[2]);
		action = argv[3];
		configFile = argv[4];
	}
	// GenFit2 [type] [pID] [action] [generation] [conf_file]
	else if (argc == 6) {
		type = argv[1];
		pID = atoi(argv[2]);
		action = argv[3];
		generation = atoi(argv[4]);
		configFile = argv[5];
	}
	// Output help if we have too too few or too many arguments
	else {
		outputHelp();
		return 0;
	}

	// Default action is to 'run'
	if (action.empty())
		action = "run";

	// Default type is 'master'
	if (type.empty())
		type = "master";

	// Be sure action and type are valid
	if (action != "run" && action != "cluster" && action != "load" && action != "results" && action != "resume") {
		outputError("Error: Couldn't find a valid 'action' in your arguments.");
	}
	if (type != "master" && type != "particle") {
		outputError("Error: Couldn't find a valid 'type' in your arguments.");
	}

	/*
	cout << "type: " << type << endl;
	cout << "action: " << action << endl;
	cout << "config: " << configFile << endl;
	cout << "pID: " << pID << endl;
	cout << "generation: " << generation << endl;
	*/

	// Regardless of type or action, we need to set up the Swarm
	//Timer tmr;

	//double t = tmr.elapsed();
	//cout << "Adding .conf took " << t << " seconds" << endl;

	//tmr.reset();

	Swarm *s;
	if (type == "master") {
		Config myconfig(configFile);
		if (action != "load" && action != "resume") {
			s = myconfig.createSwarmFromConfig();

			if (s->options.useCluster) {
				action = "cluster";
			}

			//t = tmr.elapsed();
			//cout << "Processing .conf took " << t << " seconds" << endl;

			s->currentGeneration = 1;
			s->setExePath(convertToAbsPath(argv[0]));
			s->isMaster = true;
			s->initComm();
			s->isMaster = false;
		}

		if (action == "cluster" || action == "run") {
			int randNum = rand();
			string serializedSwarmPath = to_string(static_cast<long long int>(randNum)) + ".sconf";

			std::ofstream ofs(serializedSwarmPath);
			if (ofs.is_open()) {
				s->setsConf(convertToAbsPath(serializedSwarmPath));

				boost::archive::binary_oarchive ar(ofs);
				ar & s;
				ofs.close();
			}

			s->isMaster = true;
		}

		if (action == "cluster") {
			string runCmd = s->getClusterCommand(string(convertToAbsPath(argv[0])));
			//string runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));

			if (s->options.saveClusterOutput) {
				string outputPath = s->options.outputDir + "/" + s->options.jobName + "_cluster_output";
				//cout << "string: " << outputPath << endl;
				if (!checkIfFileExists(outputPath)) {
					string makeClusterOutputDirCmd = "mkdir " + outputPath;
					if (runCommand(makeClusterOutputDirCmd) != 0) {
						cout << "Warning: Couldn't create cluster output directory with command: " << makeClusterOutputDirCmd << ". Turning off save_cluster_output" << endl;
						s->options.saveClusterOutput = false;
						runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));
					}
				}
			}

			cout << "Running BioNetFit on cluster with command: " << runCmd << endl;

			if(runCommand(runCmd) != 0) {
				outputError("Error: Couldn't launch BioNetFit on cluster with command: " + runCmd + ". Quitting.");
			}

			return 0;
		}
		else if (action == "load") {
			std::ifstream ifs(configFile);

			if (ifs.is_open()) {
				boost::archive::binary_iarchive ar(ifs);

				// Load data
				ar & s;
				ifs.close();
			}
			else {
				outputError("Error: Couldn't load config file: " + configFile + ".");
			}

			if (s->options.useCluster) {
				setenv("OMPI_MCA_mpi_warn_on_fork","0",1);
			}

			s->isMaster = true;
			s->setExePath(convertToAbsPath(argv[0]));
			s->initComm();
		}
		else if (action == "resume") {
			string swarmState = configFile + "/swarmState.sconf";
			std::ifstream ifs(swarmState);

			if (ifs.is_open()) {
				cout << "trying to load swarm..." << endl;
				boost::archive::binary_iarchive ar(ifs);

				// Load data
				ar & s;
				ifs.close();
			}
			else {
				outputError("Error: Couldn't load swarm state: " + swarmState + ".");
			}

			cout << "loaded" << endl;

			s->resumingSavedSwarm = true;

			if (s->options.useCluster) {
				setenv("OMPI_MCA_mpi_warn_on_fork","0", 0);
				int randNum = rand();
				string serializedSwarmPath = to_string(static_cast<long long int>(randNum)) + ".sconf";

				std::ofstream ofs(serializedSwarmPath);
				if (ofs.is_open()) {
					s->setsConf(convertToAbsPath(serializedSwarmPath));
					//cout << "Path is: " << s->getsConf() << endl;

					boost::archive::binary_oarchive ar(ofs);
					ar & s;
					ofs.close();
				}

				s->isMaster = true;
				string runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));

				if (s->options.saveClusterOutput) {
					string outputPath = s->options.outputDir + "/" + s->options.jobName + "_cluster_output";
					if (!checkIfFileExists(outputPath)) {
						string makeClusterOutputDirCmd = "mkdir " + outputPath;
						if (runCommand(makeClusterOutputDirCmd) != 0) {
							cout << "Warning: Couldn't create cluster output directory with command: " << makeClusterOutputDirCmd << ". Turning off save_cluster_output" << endl;
							s->options.saveClusterOutput = false;
							runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));
						}
					}
				}

				cout << "Running BioNetFit on cluster with command: " << runCmd << endl;

				if(runCommand(runCmd) != 0) {
					outputError("Error: Couldn't launch BioNetFit on cluster with command: " + runCmd + ". Quitting.");
				}
			}
			else {
				s->initComm();
				s->initRNGS(s->options.seed);
				s->isMaster = true;
				s->doSwarm();
			}

			return 0;
		}

		if (action != "results") {
			s->doSwarm();
		}
		else {
			string messageFilePath = s->options.jobOutputDir + ".req";

			ofstream outFile;
			outFile.open(messageFilePath);

			if (outFile.is_open()) {
				outFile << "output results";
				outFile.close();
			}
			else {
				string errMsg = "Error: Couldn't open file " + messageFilePath + " to request results from the swarm master.";
				outputError(errMsg);
			}
		}
	}

	// We are a particle
	else if (type == "particle"){
		// Try to open the serialized swarm
		while(1) {
			// Create and input archive
			std::ifstream ifs(configFile);

			if (ifs.is_open()) {
				boost::archive::binary_iarchive ar(ifs);

				// Load data
				ar & s;
				ifs.close();
				break;
			}
		}

		s->isMaster = false;
		s->setExePath(convertToAbsPath(argv[0]));

		s->initComm();

		if (pID == 0) {
			pID = s->swarmComm->getRank();
		}

		s->initRNGS(s->options.seed + pID);

		Particle *p = s->createParticle(pID);
		p->setModel(s->options.model);

		if (s->currentGeneration == 1) {
			p->generateParams();
		}

		if (s->options.useCluster) {
			setenv("OMPI_MCA_mpi_warn_on_fork","0",1);
		}

		p->doParticle();
	}

	s->swarmComm->~Pheromones();

	return 0;
}
