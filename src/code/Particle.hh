/*============================================================================
// Name        : Particle.hh
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Lat Update: : 2017-01-7
// Copyright   :
// Description :
//============================================================================*/

#ifndef PARTICLE_HH_
#define PARTICLE_HH_

#include "Model.hh"
#include "Swarm.hh"
#include "Setting.hh"
#include <random> //Fixing support for other versions of gcc/++ and mpicc/++

class Swarm;
class Data;
class FreeParam;
class Model;
class FreeParam;
class Data;
class Pheromones;
class Particle;
	class subParticle;  //razi added to support multiple files in parallel
//	unsigned int fcalcsubParID(unsigned int ParID, unsigned int mid, unsigned int nModels){return 1+mid+(ParID-1)* nModels;} //razi: moved to Utils.hh
//	unsigned int fcalcMID(unsigned int subParID, unsigned int nModels){return (unsigned int)((subParID-1) % nModels);} //razi: moved to Utils.hh
//	unsigned int fcalcParID(unsigned int subParID, unsigned int nModels){return (unsigned int)((subParID-1)/ nModels)+1;} //razi: moved to Utils.hh

class Particle {
public:
	Particle(Swarm * swarm, int id);

	void setID(int id);
	int getID() { return id_; }

	void generateParams();

	std::vector<std::map<std::string, double> > simParams_;
	void setParam(std::pair<std::string, double> myParams, unsigned int mid);
	std::map<std::string,double> getParams(unsigned int mid) { return simParams_[mid]; }
	void setModel(Model * model, unsigned int mid);
	std::vector<std::map<int, double> > fitCalcs;

	unsigned int calcsubParID(unsigned int mid);
	unsigned int calcParID(unsigned int subParID, unsigned int mid);
	unsigned int calcModelID(unsigned int subParID);





//protected:    //was private:

	friend class boost::serialization::access;
	void runModel(unsigned int id = 0, bool localSearch = false);

	void runNelderMead(std::map<double, std::vector<double> > simplex, unsigned int mid);
	void checkMessagesGenetic(unsigned int mid);
	void checkMessagesPSO(unsigned int mid);
	bool checkMessagesDE(unsigned int mid);
	void checkMessagesGenetic(){checkMessagesGenetic(0);}
	void checkMessagesPSO(){checkMessagesPSO(0);}
	bool checkMessagesDE(){ return checkMessagesDE(0);}
	void calculateFit(bool local, unsigned int mid);
	void finalizeSim(unsigned int mid);
	void smoothRuns(unsigned int mid);

	std::vector<double> getCentroid(std::vector<std::vector<double> >);

	double objFunc_chiSquare(double sim, double exp, double stdev);
	double objFunc_sumOfSquares(double sim, double exp, double dummyvar);
	double objFunc_divByMeasured(double sim, double exp, double dummyvar);
	double objFunc_divByMean(double sim, double exp, double mean);


	typedef double (Particle::*objFunc)(double exp,double sim ,double stdev);

	objFunc objFuncPtr;

	unsigned int id_;

	std::vector <subParticle *> subParticles;
	std::vector <unsigned int> subParIDs;

	void addSubParticle(subParticle* subparticle, unsigned int mid, bool overwrite);
	bool subParExist(unsigned int subParID);
	unsigned int getNumSubParticle();
	subParticle * getSubParticle(unsigned int subParIndex){return subParticles.at(subParIndex);} //razi: in master it is equAL TO MODEL ID
	void deleteSubParticle(){subParticles.clear();}
	void removeSubParticle(unsigned int subParIndex){subParticles.erase(subParticles.begin()+subParIndex);}
	std::vector<Model *> models;
	void doParticle(unsigned int mid);
	std::vector<std::map<std::string, std::map<int, Data*> > > dataFiles_;

	Swarm * swarm_;
	unsigned int currentGeneration_;
	unsigned int island_;

	bool ignore_ = false;

	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {

		ar & models;
		ar & subParIDs;
		ar & simParams_;
		ar & id_;
		ar & swarm_;
		ar & currentGeneration_;
	}
};



class subParticle{
public:
	subParticle(Particle * p, unsigned int modelId, unsigned int subParID);
	void setParent(Particle * p) {parParticle = p;}
	void runModel(unsigned int iteration, bool localSearch);
	void setModel(unsigned int mid);
	Particle * parParticle; //razi: model id startuing from 0
	unsigned int mid_;  //razi: global subPar Id
	Model * model;
	unsigned int subParID_; //global subParticleID  =1+ (parID-1)*(number of models) + model_id

	void runNelderMead(std::map<double, std::vector<double> > simplex);
	void doParticle();
	void setModel(Model * model);
	void setParam(std::pair<std::string, double> myParams);
	void calculateFit(bool local);
	void finalizeSim();
	void smoothRuns();

	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & parParticle;
		ar & mid_;
		ar & model;
		ar & subParID_;
	}

};

#endif /* PARTICLE_HH_ */
