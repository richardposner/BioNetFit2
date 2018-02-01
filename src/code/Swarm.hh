/*============================================================================
// Name        : swarm.hh
// Authors     : Brandon Thomas, Abolfazl Razi,Raquel Dias
// Version     : 2.0
// Lat Update: : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/

#ifndef SWARM_HH_
#define SWARM_HH_


#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
//razi changed was #include <tr1/random>
//#if(defined _WIN32 || defined  _WIN64)
	//#include <boost/tr1/tr1/random>
//#else
//	#include <tr1/random>
//#endif

//#include <boost/tr1/tr1/random>
#include <boost/random.hpp>//new Boost 1.65
#include <iomanip>
#include <chrono>
#include <cstdio>
#include <set>
#include <string>
#include <map>

#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/cauchy_distribution.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/utility.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "Setting.hh"    //razi added
#include "Utils.hh"
#include "Timer.hh"
#include "Particle.hh"
#include "Pheromones.hh"
#include "FreeParam.hh"  //razi added
#include "Evaluate.hh" //Raquel: added


class Model;
class FreeParam;
class Data;
class Pheromones;
class Particle;

#define MAX_LONG 9223372036854775807

// Forward declaration of class boost::serialization::access
namespace boost {
namespace serialization {
class access;
}
}

float evaluateResults(string inputFile1, string inputFile2, map<int,string> constraints,  string outname);

class Swarm {
public:
	Swarm();
	std::vector<unsigned int> checkMasterMessages();
	void fixRunningParticle(int subParID);
	void addRunningParticle(int subParID);

	void setExePath(std::string path) { exePath_ = path; }
	void setConfigPath(std::string path) { configPath_ = path; }

#ifdef VER2	//razi added to support multiple files
	void addExp(std::string path, int mid);
	void setExpPath(std::string prefixPath, std::string path, int mid);  //call addExp
	void setModels(std::string prefixPath, std::string path, bool overwrite);
	void setModel(std::string path, int mid, bool overwrite);


	std::string getModelName(unsigned int modelId, bool FullPath); //std::string getModelPath(unsigned int modelId) {return modelPaths_.at(modelId);} //std::string getModel(unsigned int modelId);

	std::string getExpPath(unsigned int modelId) {return getExpPath(modelId, 0);}//razi: full exp file name
	std::string getExpPath(unsigned int modelId, unsigned int expid) {return expPaths_[modelId][expid];}
	std::string getExp(unsigned int modelId); //razi: abs file name
	std::string getExp(unsigned int modelId, unsigned int expid);

	bool check_model_exp_consistency();
	unsigned int getNumModels(){return options.models.size();}
	unsigned int getNumExp(unsigned mid){return expPaths_[mid].size();}
	void setDefaultModelIndex(int mid); //set default model
	int getdefaultModelIndex(){return options.defaultModel;}

	unsigned int getNumFreeParams() { return options.freeParams_.size(); }
	const std::map<std::string, FreeParam*>& getFreeParams_() const {return options.freeParams_;}

	/* razi: old version,  does not support multiple files
	void setModelPath(std::string path) {modelPath_ = path;}
	void setExpPath(std::string path) {expPath_ = path;}
	std::string getExpPath() {return expPath_;}
	std::string getModelPath() {return modelPath_;}
	*/

	std::string getExePath() {return exePath_;}
	std::string getConfigPath() {return configPath_;}

	void consolidate_model_params();
	void update_cur_particle_params(unsigned int pID, unsigned int mid, bool overwrite);
	void update_fitCalcs();   //razi: later develop to calculate the particle fitting based on subparticle results
	std::string order_params(std::string paramsString, unsigned int mid);

	void setsConf(std::string sConf, int mid);
	std::string getsConf(unsigned int mid) { return sConf_.at(mid); }

	vector<pair<int,float>> resultChecking(); //Raquel: result checking function added
	//void processLateParticles(int subParID, bool breed, int currGen);
	void processLateParticles(std::map<std::string, double> simPar, int subParID, bool breed, int currGen);

#else
	void addExp(std::string path);
	void setModel(std::string path);

	void setsConf(std::string sConf) { sConf_ = sConf; }
	std::string getsConf() { return sConf_; }
#endif

	void setfitType(std::string type);
	void addMutate(std::string mutateString);
	void setJobOutputDir(std::string dir);
	void generateBootstrapMaps(std::vector<std::map<std::string, std::map<std::string, std::map<double,unsigned int>>>> &bootStrapMaps);

	void outputError(std::string errorMessage);

	void initComm();
	void initRNGS(int seed=0);
	void doSwarm();

	void runSGA();
	void runSPSO();
	void runSDE();
	void runAGA();
	void runAPSO();
	void runADE();
	void runSSA();
	void runASA();


	void verbose_on(){verbose=true; options.verbosity=4;}  //razi added
	void verbose_off(){verbose=false;} //razi added
	bool check_verbose(){return verbose;} //razi added


	int printDetails(); //razi added for debugging

	std::vector<double> normalizeParams(std::vector<double> params);
	std::vector<double> deNormalizeParams(std::vector<double> params);

	Particle *createParticle(unsigned int pID);

	void getClusterInformation();
	std::string getClusterCommand(std::string cmd);
	std::string generateTorqueBatchScript(std::string cmd);
	std::string generateSlurmCommand(std::string cmd, bool multiProg = true, unsigned int nCPU = 0);
	std::string generateSlurmMultiProgCmd(std::string runCmd);
	std::string generateSlurmBatchFile(std::string runCmd);
	std::string generateMPICommand(std::string runCmd);
	std::string generateBNF2MPICommand(std::string runCmd);

	Pheromones *swarmComm;

	std::multimap<double, std::string> allGenFits;
	std::vector<std::map<std::string, std::map<std::string, std::map<double,unsigned int>>>> bootstrapMaps;


	bool isMaster;
	boost::random::mt19937 generalRand;
	boost::random::mt19937 parameterGenRand;
	unsigned int currentGeneration;
	unsigned int bootstrapCounter;
	bool resumingSavedSwarm;
	bool hasMutate;
	float fitCompareTolerance;
	bool commInit;
	bool verbose;  //razi added for debug


	// Holds the current parameter set being used by each particle
	std::map<unsigned int, std::vector<double>> particleCurrParamSets_;
	std::map<unsigned int, std::vector<double>> particleCurrParamSets_copy;

	// Holds the running best parameter set for each particle
	//razi changed, so we keep currrent parameters for each model separately
	std::map<unsigned int, std::map<unsigned int, std::vector<double>>> subparticleCurrParamSets_;   //razi: PID, MID, PARAMA VALUESS
	std::map<unsigned int, std::map<unsigned int, std::vector<double>>> subparticleCurrParamSets_copy;   //razi: PID, MID, PARAMA VALUESS

public:

	struct SwarmOpts {
		std::string jobName;	// name of the job
		std::string fitType;	// genetic or swarm
		std::string outputDir;	// root directory to use for output
		std::string jobOutputDir;// outputDir + jobName
		std::string bngCommand;	// Path to simulators
		std::map<std::string,Data*> expFiles; // experimental data file
		std::map<std::string, FreeParam*> freeParams_;
#ifdef VER2
		std::vector<Model *> models; 			// the model files
		int defaultModel;
		std::map<int,string> constraints_; //Raquel: added constraint options support
		float constraintWeight = 0; //Raquel: added constraint weight option support
#else
		Model * model; 			// the model file
#endif

		// General options
		unsigned int verbosity;		// terminal output verbosity
		bool synchronicity;		// 1 for synchronous
		unsigned int maxGenerations;// maximum number of generations
		int swarmSize;		// how many particles in the swarm //Raquel removing Wsigned warnings
		float minFit;		// we stop fitting if we reach this value
		unsigned int parallelCount;	// how many particles to run in parallel
		unsigned int objFunc;		// which objective function to use
		bool usePipes;	// whether or not to use pipes to gather simulation output
		bool useCluster;// whether or not we are running on a cluster
		int seed; // seed for the random number engines
		unsigned int bootstrap;

		bool divideByInit;// whether or not to divide simulation outputs by the value at t=0
		int logTransformSimData;// whether or not to log transform simulation data. this value acts as the base.
		bool standardizeSimData;// whether or not to standardize simulation data
		bool standardizeExpData;// whether or not to standardize experimental data

		bool deleteOldFiles; // whether or not to delete unneeded files during the fitting run

		// Genetic algorithm options
		unsigned int extraWeight;	// how much extra weight to add while breeding in genetic algorithm
		float swapRate;	// the rate at which to swap parent parameters during breeding
		bool forceDifferentParents;// whether or not to force difference parents when breeding
		unsigned int maxRetryDifferentParents;// how many times to attempt selection of different parents if forceDifferentParents is true
		unsigned int smoothing;
		unsigned int keepParents;

		string maxFitTime; //Raquel fix SLURM support
		//unsigned long maxFitTime;	// Maximum amount of time to let the fit run
		unsigned long maxNumSimulations; // Maximum number of simulations to run

		// PSO options
		float inertia; // 0.72
		float cognitive; // 1.49
		float social; // 1.49
		unsigned int nmax; // 20
		unsigned int nmin; // 80
		float inertiaInit; // 1
		float inertiaFinal; // 0.1
		float absTolerance; // 10E-4
		float relTolerance; // 10E-4
		bool mutateQPSO;
		float betaMax;
		float betaMin;

		std::string topology; // fullyconnected
		std::string psoType; // bbpso

		bool enhancedStop; // true
		bool enhancedInertia; // true

		// DE Options
		unsigned int numIslands;
		unsigned int mutateType;
		float cr;
		unsigned int migrationFrequency;
		unsigned int numToMigrate;

		// SA options
		double minTemp;
		double minRadius;
		float localSearchProbability;
		float randParamsProbability;

		unsigned int outputEvery; // In an asynchronous fit, output a fit summary every n simulations

		// Cluster options
		std::string clusterSoftware;// which cluster software to use
		std::string clusterAccount;	// user account to specify in cluster submission commands
		bool saveClusterOutput;		// whether or not to save output during a cluster fit
		std::string clusterQueue;	// The cluster queue to submit to
		bool emailWhenFinished;
		std::string emailAddress;
		std::string hostfile;

		template<class Archive>
		void serialize(Archive &ar, const unsigned int version)
		{
			//std::cout << " serializing options" << std::endl;

			ar & jobName;
			ar & fitType;
			ar & outputDir;
			ar & jobOutputDir;
			ar & bngCommand;
			ar & expFiles;
#ifdef VER2
			ar & freeParams_;
			ar & models;
			ar & defaultModel;
			ar & constraints_;
			ar & constraintWeight;
#else
			ar & model;
#endif
			ar & verbosity;
			ar & synchronicity;
			ar & maxGenerations;
			ar & swarmSize;
			ar & minFit;
			ar & bootstrap;
			ar & parallelCount;
			ar & objFunc;
			ar & usePipes;
			ar & useCluster;
			ar & seed;

			ar & divideByInit;
			ar & logTransformSimData;
			ar & standardizeSimData;
			ar & standardizeExpData;

			ar & deleteOldFiles;

			ar & extraWeight;
			ar & swapRate;
			ar & forceDifferentParents;
			ar & maxRetryDifferentParents;
			ar & smoothing;
			ar & keepParents;

			ar & maxFitTime;
			ar & maxNumSimulations;

			ar & inertia;
			ar & cognitive;
			ar & social;
			ar & nmax;
			ar & nmin;
			ar & inertiaInit;
			ar & inertiaFinal;
			ar & absTolerance;
			ar & relTolerance;
			ar & mutateQPSO;
			ar & betaMin;
			ar & betaMax;

			ar & topology;
			ar & psoType;

			ar & enhancedStop;
			ar & enhancedInertia;

			ar & numIslands;
			ar & mutateType;
			ar & cr;
			ar & migrationFrequency;
			ar & numToMigrate;

			ar & minTemp;
			ar & minRadius;
			ar & localSearchProbability;
			ar & randParamsProbability;

			ar & outputEvery;

			ar & clusterSoftware;
			ar & clusterAccount;
			ar & saveClusterOutput;
			ar & clusterQueue;
			ar & hostfile;
		}
	};
	SwarmOpts options;



private:
	friend class boost::serialization::access;

	void initFit();
	std::vector<std::vector<unsigned int> > generateTopology(unsigned int populationSize);
	void launchParticle(unsigned int pID, bool nextGen = false);
#ifdef VER2
	void launchSubParticle(unsigned int pID, unsigned int mid, bool nextGen);
	std::vector<unsigned int> update_finished_running_particles();
#endif
	void runGeneration();
	void runAsyncGeneration();

	void breedGenerationGA(std::vector<unsigned int> children = std::vector<unsigned int>());
	void runNelderMead(unsigned int receiver, unsigned int cpu);
	std::map<double, unsigned int> getNearestNeighbors(unsigned int it, unsigned int N);

	void cleanupFiles(const char * path);
	void finishFit();
	void getAllParticleParams();
	void outputRunSummary(std::string outputDir);
	void outputRunSummary();
	void outputBootstrapSummary();
	void killAllParticles(int tag);
	std::unordered_map<unsigned int, std::vector<double>> checkMasterMessagesDE();
	//std::vector<unsigned int> checkMasterMessages();
	void checkExternalMessages();
	void resetVariables();
	void generateBestFitModel(std::string outputDirectory);

	void initPSOswarm(bool resumeFit = false);
	void processParticlesPSO(std::vector<unsigned int> particles, bool nextFlight = false);
	void updateEnhancedStop();
	double getEuclidianNorm(double y, unsigned int n);
	void updateParticleWeights();
	double calcParticleWeight(unsigned int particle);
	double calcWeightedAveragePosition();
	void processParamsPSO(std::vector<double> &params, unsigned int pID, double fit);
	void processParamsDE(std::vector<double> &params, unsigned int pID, double fit);

	bool checkStopCriteria();
	void updateInertia();
	void updateContractionExpansionCoefficient();
	double calcMeanFit();

	void saveSwarmState();

	std::vector<double> calcParticlePosPSO(unsigned int particle);
	std::vector<double> calcParticlePosBBPSO(unsigned int particle, bool exp = false);
	std::vector<double> calcParticlePosQPSO(unsigned int particle, std::vector<double> mBests);
	std::vector<double> getNeighborhoodBestPositions(unsigned int particle);
	std::vector<double> calcQPSOmBests();

	unsigned int pickWeighted(double weightSum, std::multimap<double, unsigned int> &weights, unsigned int extraWeight);
	std::string mutateParamGA(FreeParam* fp, double paramValue);

	std::vector<double> mutateParticleDE(unsigned int particle, float mutateFactor = 0);
	std::vector<double> mutateParticleSA(unsigned int particle, float mutateFactor = 0);
	std::vector<double> crossoverParticleDE(unsigned int particle, std::vector<double> mutationSet, float cr = 0, bool normalize = false);
	void sendMigrationSetDE(unsigned int island, std::vector<std::vector<unsigned int>> islandTopology, std::map<unsigned int, std::vector<std::vector<double>>> &migrationSets);
	void recvMigrationSetDE(unsigned int island, std::map<unsigned int, std::vector<std::vector<double>>> &migrationSets);

	std::vector<double> generateParticleTemps();
	std::vector<double> generateParticleRadii();
	std::vector<float> generateParticleFs();
	std::vector<float> generateParticleCRs();
	unsigned int pickWeightedSA();
	bool metropolisSelection(unsigned int particle, double fit, float particleTemp);
	void swapTR(std::vector<double> particleRadii, std::vector<double> particleTemps);
	std::vector<double> generateTrialPointSA(unsigned int controller, unsigned int receiver, std::vector<double> particleRadii, std::vector<float>particleCRs, std::vector<float>particleFs, std::vector<std::vector<float>> &trialParams);
	double normalizeParam(double oldParam, double min, double max, bool log);
	double deNormalizeParam(double oldParam, double min, double max, bool log);

	void insertKeyByValue(std::multimap<double, unsigned int> &theMap, double key, unsigned int value);
	std::string exePath_;
	std::string configPath_;

#ifdef VER2	//razi added
	std::set<unsigned int> runningSubParticles_;
	std::set<unsigned int>::iterator runningSubParticlesIterator_;

	std::set<unsigned int> failedSubParticles_;
	std::set<unsigned int>::iterator failedSubParticlesIterator_;

	std::set<unsigned int> finishedSubParticles_;
	//std::set<unsigned int>::iterator finishedSubParticlesIterator_;

	std::set<unsigned int> finishedParticles_;
	//std::set<unsigned int>::iterator finishedParticlesIterator_;


	vector<pair<int,float>> subParRankFinal; //Raquel added this variable to store the model checging results as a global variable

	std::map<unsigned int, std::vector<std::string>> expPaths_;
	std::vector<std::string> sConf_;

	// razi: Counts how many iterations each particle/subparticle has performed
	std::map<unsigned int, std::map<unsigned int, unsigned int>> particleIterationCounter_; //std::map<PID, MID, Iteration>> particleIterationCounter_


	std::vector<std::map<unsigned int, unsigned int>> freeParamMapping; //razi this keeps the mapping from models free parameters to the full list of free parameters
	//Ex: for two models we may have <<1,1>,<2,2>,<3,3>, <4,5>, <1,2>,<2,3>, <3,4>, <4,6>> which means there are 4 free paramas in model1 and 4 free paramas in model 2 with the given ordering
	// model1: k1,k2,k3,k5, model2: k2,k3,k4,k6,   Full list :K1,K2,K3,K4,k5,k6


#else
	std::string sConf_;

	// Counts how many iterations each particle has performed
	std::map<unsigned int, unsigned int> particleIterationCounter_;


#endif
	std::set<unsigned int> runningParticles_;
	std::set<unsigned int>::iterator runningParticlesIterator_;

	std::set<unsigned int> failedParticles_;
	std::set<unsigned int>::iterator failedParticlesIterator_;

	std::vector<std::vector<unsigned int> > populationTopology_;

	// Maybe we can change them to vectors, too
	std::map<unsigned int, double> particleBestFits_;
	std::map<unsigned int, double> subparticleBestFits_; //Raquel added support to subparticles
	std::map<unsigned int, double> subparticleBestFits_copy; //Raquel added support to subparticles

	std::multimap<double, unsigned int> particleBestFitsByFit_;
	std::multimap<double, unsigned int> particleBestFitsByFit_copy;
	std::multimap<double, unsigned int> subparticleBestFitsByFit_; //Raquel added support to subparticles
	std::multimap<double, unsigned int> subparticleBestFitsByFit_copy; //Raquel added support to subparticles

	std::multimap<double, unsigned int> currentsubparticleBestFitsByFit_; //Raquel added support to subparticles
	std::multimap<double, unsigned int> currentsubparticleBestFitsByFit_copy; //Raquel added support to subparticles

	std::multimap<double, unsigned int> swarmBestFits_;
	std::multimap<double, unsigned int> subswarmBestFits_; //Raquel added support to subparticles
	std::multimap<double, unsigned int> currentsubswarmBestFits_; //Raquel added support to subparticles
	std::multimap<double, unsigned int> currentsubswarmBestFits_copy; //Raquel added support to subparticles

	std::map<unsigned int, std::vector<double>> particleParamVelocities_;

	// Holds the running best parameter set for each particle
	std::map<unsigned int, std::vector<double>> particleBestParamSets_;


	// Holds particle weights for use in enhancedStop stop criteria
	std::map<unsigned int, double> particleWeights_;


	// In DE, maps islands and particles together
	std::vector<unsigned int> particleToIsland_; // particle -> island
	std::vector<std::vector<unsigned int>> islandToParticle_; // island -> particle

	// Holds the current parameter set being used by each particle
	std::map<unsigned int, std::vector<double>> particleTrialParamSets_;

	unsigned int permanenceCounter_; // 0
	unsigned int flightCounter_; // 0
	double weightedAvgPos_; // 0
	double optimum_; // 0
	unsigned int inertiaUpdateCounter_; // 0;
	double beta_;
	double cauchyMutator_;


	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		//std::cout << " serializing swarm" << std::endl;

		ar & options;

		ar & allGenFits;
		ar & bootstrapMaps;

		ar & isMaster;
		ar & currentGeneration;
		ar & resumingSavedSwarm;
		ar & hasMutate;

#ifdef VER2	//razi added
		ar & runningSubParticles_;
		ar & failedSubParticles_;
		ar & finishedSubParticles_;
		ar & finishedParticles_;
		ar & freeParamMapping; //razi this keeps the mapping from models free parameters to the full list of free parameters
#endif
		ar & runningParticles_;
		ar & failedParticles_;

		ar & exePath_;
		ar & configPath_;
		ar & sConf_;

		ar & populationTopology_;
		ar & particleBestFits_;
		ar & particleBestFitsByFit_;
		ar & swarmBestFits_;
		ar & particleParamVelocities_;
		ar & particleBestParamSets_;
		ar & particleCurrParamSets_;
		ar & particleWeights_;
		ar & particleIterationCounter_;
		ar & particleToIsland_;
		ar & islandToParticle_;

		ar & permanenceCounter_;
		ar & flightCounter_;
		ar & weightedAvgPos_;
		ar & optimum_;
		ar & inertiaUpdateCounter_;
		ar & beta_;
		ar & cauchyMutator_;
	}

	Timer tmr_;
	Timer generationTime_; //Raquel adding walltime functionality

};



#endif /* SWARM_HH_ */
