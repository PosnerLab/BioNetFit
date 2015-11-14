/*
 * FreeParam.hh
 *
 *  Created on: Nov 13, 2015
 *      Author: brandon
 */

#ifndef CODE_FREEPARAM_HH_
#define CODE_FREEPARAM_HH_

#include <string>

class FreeParam {
public:
	FreeParam(std::string parameterName);

	const std::string& getGenerationMethod() const {
		return generationMethod_;
	}

	void setGenerationMethod(const std::string& generationMethod) {
		generationMethod_ = generationMethod;
	}

	float getMutationFactor() const {
		return mutationFactor_;
	}

	void setMutationFactor(float mutationFactor) {
		mutationFactor_ = mutationFactor;
	}

	float getMutationRate() const {
		return mutationRate_;
	}

	void setMutationRate(float mutationRate) {
		mutationRate_ = mutationRate;
	}

	const std::string& getParameterName() const {
		return parameterName_;
	}

	void setParameterName(const std::string& parameterName) {
		parameterName_ = parameterName;
	}

	float getGenMax() const {
		return genMax_;
	}

	void setGenMax(float genMax) {
		genMax_ = genMax;
	}

	float getGenMin() const {
		return genMin_;
	}

	void setGenMin(float genMin) {
		genMin_ = genMin;
	}

	bool isHasMutation() const {
		return hasMutation_;
	}

	void setHasMutation(bool hasMutation = false) {
		hasMutation_ = hasMutation;
	}

private:
	std::string parameterName_;
	std::string generationMethod_;
	float mutationRate_;
	float mutationFactor_;
	float genMin_;
	float genMax_;
	bool hasMutation_ = false;
};

#endif /* CODE_FREEPARAM_HH_ */
