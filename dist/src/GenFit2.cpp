//============================================================================
// Name        : GenFit2.cpp
// Author      : Brandon Thomas
// Version     :
// Copyright   : 
// Description :
//============================================================================

#include <iostream>

//#include <mpi.h>
//#include <sys/shm.h>
//#include <sys/types.h>
//#include <fcntl.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <sys/stat.h>
//#include <sys/wait.h>
//#include <unistd.h>

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
	}

	// Default action is to 'run'
	if (action.empty())
		action = "run";

	// Default type is 'master'
	if (type.empty())
		type = "master";

	// Be sure action and type are valid
	if (action != "run") {
		outputError("Error: Couldn't find a valid 'action' in your arguments.");
	}
	if (type != "master" && type != "particle") {
		outputError("Error: Couldn't find a valid 'type' in your arguments.");
	}

	// Regardless of type or action, we need to set up the Swarm
	//Timer tmr;

	Config myconfig(configFile);

	//double t = tmr.elapsed();
	//cout << "Adding .conf took " << t << " seconds" << endl;

	//tmr.reset();

	Swarm *s;
	if (type == "master") {
		s = myconfig.createSwarmFromConfig((type=="master") ? true : false);

		//t = tmr.elapsed();
		//cout << "Processing .conf took " << t << " seconds" << endl;
		s->setType(type);
		s->setExePath(argv[0]);
		s->initComm();

		if (generation) {
			cout << "setting gen to " << generation << endl;
			s->currentGeneration = generation;
		}

		cout << "trying to serialize swarm class" << endl;

		std::ofstream ofs("swarm.txt");
		if (ofs.is_open()) {
			boost::archive::text_oarchive ar(ofs);
			ar & s;
			ofs.close();
		}

		s->setIsMaster(true);
		s->doSwarm();
	}

	// We are a particle
	else if (type == "particle"){

		// Create and input archive
		std::ifstream ifs("swarm.txt");
		if (ifs.is_open()) {
			boost::archive::text_iarchive ar(ifs);

			cout << "particle trying to load swarm";
			// Load data
			ar & s;
			ifs.close();
		}

		s->setIsMaster(false);
		s->setExePath(argv[0]);
		s->initComm();
		Particle *p = s->createParticle(pID);
		p->setModel(s->getModel());

		if (s->currentGeneration == 1) {
			p->generateParams();
		}
		p->doParticle();
	}

	return 0;
}
