/*
 * FreeParam.cpp
 *
 *  Created on: Nov 13, 2015
 *      Author: brandon
 */

#include "FreeParam.hh"

using namespace std;

FreeParam::FreeParam(string parameterName) {
	parameterName_ = parameterName;
	generationMethod_ = "";
	hasMutation_ = false;
	mutationRate_ = 0;
	mutationFactor_ = 0;
	genMin_ = 0;
	genMax_ = 0;
}

FreeParam::FreeParam() {
	parameterName_ = "";
	generationMethod_ = "";
	hasMutation_ = false;
	mutationRate_ = 0;
	mutationFactor_ = 0;
	genMin_ = 0;
	genMax_ = 0;
}
