/*============================================================================
// Name        : Model.hh
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Last Update : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/

#ifndef MODEL_HH_
#define MODEL_HH_

#include <iostream>
#include <vector>
#include <fstream>
#include <boost/regex.hpp>
#include <cstdlib>

#include "Utils.hh"
#include "Data.hh"
#include "FreeParam.hh"
#include "Setting.hh"

class FreeParam;
class Data;
class Swarm;

class Model {
	friend class Config;
	friend class Swarm;
	friend class Particle;

public:
	Model(Swarm *swarm, std::string path);
	Model();

	std::string getLocation();
	int makeCopyofOrig(std::string newLocation);
	void outputModelWithParams(std::map<std::string,double> params, std::string path, std::string filename, std::string suffix, bool stopAtNetGen, bool onlyActions, bool netAndBngl, bool usePipe, bool isNetFile);
	void parseNet(std::string path);
	bool getHasGenerateNetwork() {return hasGenerateNetwork_;}
	unsigned int getNumFreeParams() { return freeParams_.size(); }
	std::string getName(){return modelPath_;} //razi added

	struct action {
		std::string full;

		//double t_start;
		double t_end;
		//double o_steps;

		std::string type;
		std::string scanParam;

		Data *dataSet;

		template<typename Archive>
		void serialize(Archive& ar, const unsigned version) {

			ar & full;
			ar & t_end;
			ar & type;
			ar & scanParam;

			ar & dataSet;
		}
	};

	const std::map<std::string, FreeParam*>& getFreeParams_() const {
		return freeParams_;
	}

	// Map key contains the action prefix, map value contains the action information
	std::map<std::string, action> actions;

	std::map<std::string, FreeParam*> freeParams_;


private:
	friend class boost::serialization::access;

	void parseModel();

	Swarm *swarm_;

	std::string modelPath_;
	std::vector<std::string> fullContents_;
	std::vector<std::string> netContents_;

	bool hasGenerateNetwork_;

	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		//std::cout << " serializing model" << std::endl;
		ar & swarm_;

		ar & freeParams_;

		ar & actions;
		ar & modelPath_;
		ar & fullContents_;
		ar & netContents_;

		ar & hasGenerateNetwork_;

	}
};

bool check_model_consistency(std::vector<Model *> models);//razi added

#endif /* MODEL_HH_ */
