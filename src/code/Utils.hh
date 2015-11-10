/*
 * Utils.hh
 *
 *  Created on: Jul 14, 2015
 *      Author: brandon
 */

#ifndef UTILS_HH_
#define UTILS_HH_

#include <iostream>
#include <map>
#include <random>
#include <boost/filesystem.hpp>.

#include <sys/stat.h>

void outputError(std::string errorMessage);
std::string convertToAbsPath(std::string relPath);
std::string getFilename(std::string path);
int checkIfFileExists(std::string path);
void split(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ");
void outputHelp();
bool createParticlePipe(const char * path);
double pickWeighted(double weightSum, std::map<double,double>& weights, int extraWeight, std::mt19937 &randNumEngine);

#endif /* UTILS_HH_ */
