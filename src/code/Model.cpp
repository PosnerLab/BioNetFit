/*
 * Model.cpp
 *
 *  Created on: Jul 13, 2015
 *      Author: brandon
 */

#include "Model.hh"

using namespace std;

Model::Model(Swarm * swarm, string path) {

	swarm_ = swarm;
	hasGenerateNetwork_ = false;
	modelPath_ = convertToAbsPath(path);
	parseModel();
}

Model::Model() {
	swarm_ = 0;
	hasGenerateNetwork_ = false;
	modelPath_ = "";
}

void Model::parseNet(string path) {

	if (swarm_->options.verbosity >= 3) {
		cout << "Parsing .net file" << endl;
	}

	string line;
	ifstream modelFile(path);

	if (modelFile.is_open()) {
		while (getline(modelFile, line)) {
			if (line.at(0) == '#') {
				continue;
			}
			string line_with_newline = line + "\n";
			netContents_.push_back(line_with_newline);
			//cout << line_with_newline;
		}
	}
	else {
		swarm_->outputError("Error: Couldn't open model file " + path + " for parsing.");
	}
	modelFile.close();
}

void Model::parseModel() {
	string line;
	ifstream modelFile(modelPath_);
	vector<string> model_vector;

	if (modelFile.is_open()) {
		while (getline(modelFile, line)) {
			string line_with_newline = line + "\n";
			model_vector.push_back(line_with_newline);
		}
	}
	else {
		swarm_->outputError("Error: Couldn't open model file " + modelPath_ + " for parsing.");
	}
	modelFile.close();

	// Find any line-continuation backslashes and merge the lines together
	for (unsigned i = 0 ; i < model_vector.size(); i++) {
		int j = 0;
		while(1) {
			if (boost::regex_search(model_vector[i-j],boost::regex("\\\\"))) {
				model_vector[i-j] = boost::regex_replace(model_vector[i-j],boost::regex("\\\\\\W*"),string(""));
				model_vector[i-j] = model_vector[i-j] + model_vector[i+1];
				model_vector[i-j] = boost::regex_replace(model_vector[i-j],boost::regex("\\s*"),string(""));
				if (!boost::regex_search(model_vector[i-j],boost::regex("\\\\"))) {
					fullContents_.push_back(model_vector[i-j]);
				}
				i++;
				j++;
				if (!boost::regex_search(model_vector[i-j],boost::regex("\\\\"))) {
					break;
				}
			}
			else {
				break;
			}
		}
		if (j == 0){
			fullContents_.push_back(model_vector[i]);
		}
	}

	// Go through the model and store important things
	boost::smatch smatches;
	//for (string i : fullContents_) {
	for (vector<string>::iterator i = fullContents_.begin(); i != fullContents_.end(); ++i) {
		if (boost::regex_search(*i, boost::regex("^simulate|^simulate_nf|^simulate_ode|^simulate_ssa|^simulate_pla|^parameter_scan"))) {
			action newAction;
			string prefix;
			bool isParScan = 0;

			// If the action is a normal "simulate", find the method in the arguments and store it
			if (boost::regex_search(*i, smatches, boost::regex("^simulate\\(\\{.*method=>(\"|')(\\w{2,3})(\"|')"))) {
				newAction.type = smatches[2];
			}
			// If we the action command is simulate_xx, store the sim type
			else if (boost::regex_search(*i, smatches, boost::regex("^simulate_(\\w{2})\\(\\{"))) {
				newAction.type = smatches[1];
			}
			// If we found a parameter_scan command, store the type, scan parameter, and t_end
			else if (boost::regex_search(*i, smatches, boost::regex("^parameter_scan\\(\\{.*method=>(\"|')(\\w{2,3})(\"|')"))) {
				isParScan = 1;
				newAction.type = smatches[2];

				if (boost::regex_search(*i, smatches, boost::regex("parameter=>('|\")(\\w+)('|\")"))) {
					newAction.scanParam = smatches[2];
				}
				if (boost::regex_search(*i, smatches, boost::regex("par_max=>(\\w+)"))) {
					newAction.t_end = atof(smatches[1].str().c_str());
				}
				if (boost::regex_search(*i, smatches, boost::regex("par_scan_vals=>\\[(.+)\\]"))) {
					vector<string> values;
					split(smatches[1].str(), values, ",");
					newAction.t_end = atof(values.back().c_str());
				}
			}

			// Remove any suffixes from the command
			*i = boost::regex_replace(*i,boost::regex(",suffix=>('|\")\\w+('|\")"), string(""));

			// Find any prefixes and store them
			if (boost::regex_search(*i, smatches, boost::regex("prefix=>('|\")(\\w+)('|\")"))) {
				prefix = smatches[2];
			}

			// If we find a non parameter_scan action, store the t_end value
			if (boost::regex_search(*i, smatches, boost::regex("t_end=>(\\w+)")) && !isParScan) {
				newAction.t_end = atof(smatches[1].str().c_str());
			}

			newAction.full = *i;
			if (!prefix.empty()) {
				actions.insert(pair<string,action>(prefix,newAction));
			}
		}
		// Save any free parameters
		else if (boost::regex_search(*i, smatches, boost::regex("(\\s+|=\\s*)(\\w+)__FREE__"))) {
			//freeParams_.insert(pair<string,string>(smatches[2],""));
			FreeParam * fp = new FreeParam(smatches[2]);
			freeParams_.insert(pair<string,FreeParam*>(smatches[2],fp));
		}
		// Make sure we know if we need to do network generation
		else if (boost::regex_search(*i, smatches, boost::regex("^generate_network"))) {
			hasGenerateNetwork_ = true;
		}
	}

	if (actions.size() == 0) {
		swarm_->outputError("Error: We didn't find any action commands containing prefixes in your model file. Your model must contain an action command that uses a prefix which corresponds to the name of the .exp file to be fit to the data generated by that action command.");
	}
}

void Model::outputModelWithParams(map<string, double> params, string path, string filename, string suffix, bool stopAtNetGen=false, bool onlyActions=false, bool netAndBngl=false, bool usePipe=false, bool isNetFile=false) {

	if (netAndBngl) {
		// First output the .bngl file (containing only action commands)
		outputModelWithParams(params, path, filename, suffix, false, true, false, false, false);
		filename = boost::regex_replace(filename, boost::regex("bngl$"), string("net"));
		// Then output the .net file
		outputModelWithParams(params, path, filename, "", false, false, false, false, true);

		return;
	}

	if (!usePipe) {
		string fullPath = path + "/" + filename;
		unlink(fullPath.c_str());

		ofstream outFile;
		outFile.open(fullPath);
		boost::smatch matches;

		if (outFile.is_open()) {
			if (isNetFile) {
				bool inParameterBlock = true;
				int numParamsToReplace = params.size();
				int numReplacedParams = 0;

				//for (string line : netContents_){
				for (auto line = netContents_.begin(); line != netContents_.end(); ++line) {
					if (*line == "end parameters" || numReplacedParams == numParamsToReplace) {
						inParameterBlock = false;
					}
					else if (inParameterBlock) {
						// Replace free param with generated param
						//double tt = 0;
						//for (auto p : params) { // TODO: Is there a faster way to do this than loop through params over and over? It's still too slow.
						for (map<string, double>::iterator p = params.begin(); p != params.end(); ++p) {
							//cout << "p is " << p.first << endl;
							//Timer tmr;
							if(boost::regex_match(*line, matches, boost::regex("\\s+\\d+\\s+(\\w+)\\s+(.+)\\s+"))) {
								//double t = tmr.elapsed();
								//cout << "Match took " << t << " seconds" << endl;
								//tt += t;
								if (matches[1] == p->first) {
									string match = p->first + "\\s+.+";
									string replacement = p->first + " " + to_string(static_cast<long double>(abs(p->second))) + "\n";
									*line = boost::regex_replace(*line, boost::regex(match), replacement);
									numReplacedParams++;
								}
							}
						}
						//cout << "matches took " << tt << " seconds" << endl;
					}
					//cout << line;
					outFile << *line;
				}
			}
			else {
				if (onlyActions) {
					filename = boost::regex_replace(filename, boost::regex("bngl$"), string("net"));
					string line = "readFile({file=>\"" + path + filename + "\"})\n";
					outFile << line;
				}

				//for (string line : fullContents_){
				for (auto line = fullContents_.begin(); line != fullContents_.end(); ++line) {
					if (onlyActions) {
						string actionLine = *line;
						if (boost::regex_search(*line, boost::regex("^simulate|^simulate_nf|^simulate_ode|^simulate_ssa|^simulate_pla|^parameter_scan|^setConcentration|^addConcentration|^saveConcentration|^resetConcentrations|^setParameter|^saveParameters|^resetParameters|^quit|^substanceUnits|^version|^setOption"))) {
							if (!suffix.empty()) {
								string suffixLine = ",suffix=>\"" + suffix + "\"})";
								actionLine = boost::regex_replace(actionLine, boost::regex("\\}\\)"), suffixLine);
							}
						}
						else {
							continue;
						}
						outFile << actionLine;
					}
					else {
						string newLine = *line;

						bool inParameterBlock = true;
						int numReplacedParams = 0;
						int numParamsToReplace = params.size();
						// Skip line if it's empty or if it is a comment
						if (boost::regex_match(*line, boost::regex("^\\s*$")) || boost::regex_match(*line, boost::regex("^#"))) {
							continue;
						}
						// TODO: Next line is never working. Same as above for .net file?
						if (*line == "end parameters" || numReplacedParams == numParamsToReplace) {
							inParameterBlock = false;
							cout << "out of parameter block" << endl;
						}

						//cout << *line << endl;
						if (inParameterBlock) {
							// Replace free param with generated param
							if (boost::regex_search(*line, matches, boost::regex("(\\s+|=\\s*)(\\w+)__FREE__"))) {
								// Older version of C++ don't support double overload to to_string, so we have to cast to long double
								cout << "found a FP: " << matches[2] << endl;
								newLine = boost::regex_replace(*line, boost::regex("\\w+__FREE__"), to_string(static_cast<long double>(params[matches[2]])));
								cout << "replacing with: " << params[matches[2]] << endl;
								++numReplacedParams;
							}
						}
						// Add in the unique suffix
						if (!suffix.empty()) {
							if (boost::regex_search(*line, boost::regex("^simulate|^simulate_nf|^simulate_ode|^simulate_ssa|^simulate_pla|^parameter_scan"))) {
								string suffixLine = ",suffix=>\"" + suffix + "\"})";
								newLine = boost::regex_replace(*line, boost::regex("\\}\\)"), suffixLine);
								newLine = newLine + "\n";
							}
						}

						outFile << newLine;
						if (stopAtNetGen) {
							if (boost::regex_search(*line,matches,boost::regex("^generate_network"))) {
								break;
							}
						}
					}
				}
			}
			outFile.close();
		}
		else {
			swarm_->outputError("Couldn't open file '" + fullPath + "' for output");
		}
	}
	else {
		const char * fifo = path.c_str();

		int fifo_status = mkfifo(fifo, 0666);
		if (fifo_status) {
			swarm_->outputError("Couldn't create pipe '" + path + "' for output");
		}

		int fd = open(fifo, O_RDONLY);
		if (fd < 0) {
			cout << "Couldn't open " << fifo << endl;
		}

		//for (auto line : fullContents_) {
		for (vector<string>::iterator line = fullContents_.begin(); line != fullContents_.end(); ++line) {
			write(fd,line->c_str(),line->size());
		}
		close(fd);
	}
}
