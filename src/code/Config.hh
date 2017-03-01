/*============================================================================
// Name        : Config.hh
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Last Update : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/


#ifndef CONFIG_HH_
#define CONFIG_HH_

#include "Parser.hh"
#include "Utils.hh"
#include "Data.hh"
#include "Swarm.hh"
#include "Model.hh"

class Swarm;

class Config {
public:
	Config(std::string configFile);

	std::string getLocation();
	int makeCopy(std::string newLocation);
	Swarm * createSwarmFromConfig();

	int printDetails();   //razi added for debuging
	void verbose_on(){verbose=true; if (swarm_!=NULL) swarm_->verbose_on(); }  //razi added
	void verbose_off(){verbose=false; if (swarm_!=NULL) swarm_->verbose_off(); } //razi added

	bool verbose;
private:
	friend class boost::serialization::access;

	void checkConsistency();

	Swarm * swarm_;
	std::string configPath_;

	/*
	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & configPath_;
	}
	*/
};

#endif /* CONFIG_HH_ */
