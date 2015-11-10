/*
 * Utils.cpp
 *
 *  Created on: Jul 14, 2015
 *      Author: brandon
 */

#include "Utils.hh"

using namespace boost::filesystem;
using namespace std;

void outputError(string errorMessage) {
	cout << errorMessage << endl;

	exit (1);
}

string convertToAbsPath(string relPath) {
	path fullPath;

	try{
		fullPath = canonical(relPath);
	}
	catch(...){}

	if (!exists(fullPath)){
		string errMsg = "Error: Can't find the file: " + relPath;
		outputError(errMsg);
	}

	return fullPath.string();
}

int checkIfFileExists(string path) {
	if (exists(path))
		return 1;
	else
		return 0;
}

void split(const string& str, vector<string>& tokens, const string& delimiters) {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

string getFilename(string path) {
	return boost::filesystem::path(path).stem().string();
}

void outputHelp() {
	cout << "GenFit2 Usage:" << endl;
	cout << "GenFit2 [config_file]" << endl;
}

bool createParticlePipe(const char * path) {

	unlink(path);
	int fifo_status = mkfifo(path, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
	if (fifo_status) {
		cout << "Warning: Couldn't create pipe with path: " << path << endl;
		return false;
	}
	else {
		return true;
	}
}

double pickWeighted(double weightSum, map<double,double> &weights, int extraWeight, mt19937 &randNumEngine) {
	double lowerBound = 0;
	double upperBound = weightSum;

	uniform_real_distribution<double> unif(lowerBound,upperBound);

	double random = unif(randNumEngine);
	double chosen = random * ( 1 - (extraWeight / 10 ));

	cout << "random: " << random << " chosen: " << chosen << endl;
	double currentSum = 0;
	for (map<double,double>::iterator w = weights.begin(); w != weights.end(); ++w) {
		currentSum += w->second;

		if (currentSum >= chosen) {
			return w->first;
		}
	}
	//cout << "fell off the end" << endl;
	return weights.rbegin()->first;
}

/*
void saveParticle(const Particle &p, const char * filename) {
    // make an archive
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << p;
}
 */
