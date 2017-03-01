/*============================================================================
// Name        : FreeParam.hh
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Last Update : 2017-01-07
// Copyright   :
// Description :
//============================================================================*/

#ifndef CODE_FREEPARAM_HH_
#define CODE_FREEPARAM_HH_

#include "Swarm.hh"

class FreeParam {
public:
	FreeParam(std::string parameterName);
	FreeParam();

	const std::string& getGenerationMethod() const {
		return generationMethod_;
	}

	void setGenerationMethod(const std::string& generationMethod) {
		generationMethod_ = generationMethod;
	}

	std::string getGenerationMethod() {
		return generationMethod_;
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

	void setIsLog(bool log) {
		isLog_ = log;
	}
	bool getIsLog() { return isLog_; }

private:
	friend class boost::serialization::access;

	std::string parameterName_;
	std::string generationMethod_;
	float mutationRate_;
	float mutationFactor_;
	float genMin_;
	float genMax_;
	bool hasMutation_;
	bool isLog_;

	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		//std::cout << " serializing fp" << std::endl;

		ar & parameterName_;
		ar & generationMethod_;
		ar & mutationRate_;
		ar & mutationFactor_;
		ar & genMin_;
		ar & genMax_;
		ar & hasMutation_;
		ar & isLog_;
	}
};

#endif /* CODE_FREEPARAM_HH_ */
