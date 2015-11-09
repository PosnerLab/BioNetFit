//============================================================================
// Name        : GenFit2.cpp
// Author      : Brandon Thomas
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <mpi.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#include "GenFit2.hh"

#define MAX_BUF 1024

using namespace std;

int main(int argc, char *argv[]) {

	srand(clock());

	int generation;
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
	Config myconfig(configFile);
	Swarm *s = myconfig.createSwarmFromConfig((type=="master") ? true : false);
	s->setType(type);
	s->setExePath(argv[0]);
	if (generation) {
		s->currentGeneration_ = generation;
	}

	// We are the master
	if (type == "master") {
		s->setIsMaster(true);
		s->doSwarm();
	}
	// We are a particle
	else if (type == "particle"){
		s->setIsMaster(false);
		Particle *p = s->createParticle(pID);
		p->setModel(s->getModel());
		if (generation == 1) {
			p->generateParams();
		}
		p->doParticle();
	}

	/*int i;
	char * command = "/home/brandon/projects/nfsim_free/Debug/nfsim_free -xml /home/brandon/projects/nfsim_free/Debug/egfr.xml  >> output";
	pid_t pid;

	int fd;
	char * fifo = "/home/brandon/projects/nfsim_free/Debug/egfr_nf.gdat";
	char buf[MAX_BUF];

	int fifo_status;
	//if( remove( S ) != 0 )
	//	perror( "Error deleting file!" );

	cout << "Creating pipe with path: " << fifo << endl;

	unlink(fifo);
	fifo_status = mkfifo(fifo, 0666);

	if (fifo_status) {
		cout << "Error: Couldn't create output pipe: " <<endl;
	}

	pid = fork();

	if (pid < 0)
	{
		cout << "Fork Error";
		return -1;
	}
	else if (pid == 0) {
		cout << "Child " << getpid() << endl;
		i = system (command);
	}
	else
	{
		fd = open(fifo, O_RDONLY);

		if (fd < 0) {
			cout << "couldn't open " << fifo << endl;
		}

		int len;
		while( len = read(fd, buf, MAX_BUF)){
			buf[len] = 0;
			printf("Buffer recived: %s\n", buf);
		}

		if (len < 0) {
		    perror("read error");
		}

		else {
		    buf[len] = 0;
		    printf("Buffer recived: %s\n", buf);
		}

		close(fd);

		wait(NULL);
		cout << "Parent " << getpid() << endl;
		return 0;
	}

	if (i != 0){
		cout<<"Return code is "<<i<<endl;
	}

	unlink(fifo);
	 */
	return 0;
}
