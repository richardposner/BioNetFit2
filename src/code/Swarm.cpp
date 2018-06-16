
/*============================================================================
// Authors     : Brandon Thomas, Abolfazl Razi, Raquel Dias
// Version     : 2.0
// Lat Update: : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/


//razi: check particleIterationCounter_ on initfit and also check messages ....

#include "Swarm.hh"

using namespace std;
using namespace std::chrono;
namespace fs = boost::filesystem;

Swarm::Swarm() {

	// Whether or not we are master
	isMaster = false;
	resumingSavedSwarm = false;
	swarmComm = 0;
	fitCompareTolerance = 1e-6;
	bootstrapCounter = 0;
	commInit = false;
	verbose = false;  //razi added for debugging

	options.jobName = "";
	options.fitType = "";
	options.outputDir = "";
	options.bngCommand = "";
	options.outputEvery = 100;

	options.models.clear();
	options.defaultModel = -1;
	sConf_.clear();

	options.synchronicity = 0;
	options.maxGenerations = 0;
	options.swarmSize = 0;
	options.minFit = 0;
	options.bootstrap = 0;
	options.parallelCount = 0;
	options.seed = 0;

	options.useCluster = false;
	options.saveClusterOutput = false;
	options.emailWhenFinished = false;
	options.emailAddress = "";

	options.usePipes = false;
	options.divideByInit = false;
	options.logTransformSimData = false;
	options.standardizeSimData = false;
	options.standardizeExpData = false;

	options.deleteOldFiles = true;
	options.objFunc = 1;
	options.extraWeight = 0;
	options.swapRate = 0.5;
	options.forceDifferentParents = true;
	options.maxRetryDifferentParents = 100;
	options.smoothing = 1;
	options.keepParents = 0;

	//options.maxFitTime = MAX_LONG;
	options.maxFitTime = std::string("700:00:00");
	options.maxNumSimulations = MAX_LONG;

	// PSO options
	options.inertia = 0.72; // 0.72
	options.cognitive = 1.49; // 1.49
	options.social = 1.49; // 1.49
	options.nmax = 0; // 20
	options.nmin = 80; // 80
	options.inertiaInit = 1; // 1
	options.inertiaFinal = 0.1; // 0.1
	options.absTolerance = 10e-4; // 10E-4
	options.relTolerance = 10e-4; // 10E-4
	options.mutateQPSO = false;
	options.betaMin = 0.5;
	options.betaMax = 1.0;

	options.topology = "fullyconnected"; // fullyconnected
	options.psoType = "pso"; // pso

	options.enhancedStop = false; // true
	options.enhancedInertia = false; // true

	options.numIslands = 0;

	options.minTemp = pow(10, -10);
	options.minRadius = pow(10, -6);
	options.localSearchProbability = 0.01;
	options.randParamsProbability = 0.1;

	options.verbosity = 1;
	hasMutate = false;

	currentGeneration = 0;

	/*
	std::map<int, double> particleBestFits_;
	std::map<int, std::vector<double> > particleBestParamSets_;
	std::map<int, std::vector<double> > particleCurrParamSets_;
	std::map<int, double> particleWeights_;
	std::map<double, int> particleBestFitsByFit_;
	 */

	permanenceCounter_ = 0; // 0
	flightCounter_ = 0; // 0
	weightedAvgPos_ = 0; // 0
	optimum_ = 0; // 0
	inertiaUpdateCounter_ = 0; // 0;
	beta_ = 0.7;
	cauchyMutator_ = 0.2;

	auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	generalRand.seed(seed);
	generalRand.discard(700000);

	/*
	if (options.seed) {
		generalRand.seed(options.seed);
		generalRand.discard(700000);

		srand(options.seed);
	}
	else {
		// TODO: Make sure everything is being seeded properly and in the proper place. Also let's do away with rand()
		// Seed our random number engine
		auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
		generalRand.seed(seed);
		generalRand.discard(700000);

		srand (std::tr1::random_device{}());
	}
	 */

}

void Swarm::initRNGS(int seed) {
	if (seed) {
		srand(seed);
		if(options.verbosity >= 3){
			cout << "seeding rng with: " << seed << endl;
		}
	}
	else {
		// TODO: Make sure everything is being seeded properly and in the proper place. Also let's do away with rand()
		// Seed our random number engine
//razi changed was 		srand (std::tr1::random_device{}());
		srand (std::random_device{}());
	}
}

void Swarm::initComm() {
	// Create the communication class
	Pheromones *ph = new Pheromones();

	// Initialize the communication class
	ph->init(this);

	// Store our communication class
	swarmComm = ph;
}


double Swarm::normalizeParam(double oldParam, double min, double max, bool log) {
	if (log) {
		return (log10(oldParam) - log10(min)) / (log10(max) - log10(min));
		//return (oldParam - min) / (max - min);
	}
	else {
		return (oldParam - min) / (max - min);
	}
}


double Swarm::deNormalizeParam(double oldParam, double min, double max, bool log) {
	if (log) {
		return pow(10, log10(min) + oldParam * (log10(max) - log10(min)));
		//return min + oldParam * (max - min);
	}
	else {
		return min + oldParam * (max - min);
	}
}

































void Swarm::addExp(string path, int mid) {
//	Timer tmr;

//cout << "Swarm::addExp   1-path:"<< path<<endl; mypause();
	path = convertToAbsPath(path);
//cout << "Swarm::addExp   2-fullpath:"<< path<<endl; mypause();
	string basename = getFilename(path);
//cout << "Swarm::addExp   3-basename:"<< basename<<endl; mypause();

//cout << "Swarm::addExp inserting EXP?" << endl; mypause();
	this->options.expFiles.insert(make_pair(basename, new Data(path, this, true, mid)));
//cout << "Swarm::addExp exp inserted" << endl; mypause();

	//double t = tmr.elapsed();
	//cout << "Adding .exp took " << t << " seconds" << endl;
}


/*
//razi: this version is based on having different set of free parameters for each model
// the implementation is subject to change, lets keep the same set of free parameters for all models and use a mask vector
// to account for heterogenous set of free parameters, note: check what happens if there is a conflict between the free parameter and fixed parameters
*/

void Swarm::consolidate_model_params(){

	if(options.verbosity>=1){
		cout<<"consolidating " << this->getNumModels()<<" models."<<endl;
	}

	//Raquel: commented, unnecessary
	/*if (options.models.size() <1){
		if(options.verbosity >= 3) cout<<"No need to consolidate models. Less than two models exist"<<endl;
		return;
	}
 *
 * Bug Note by Raquel
 * There was a bug in the consolidation function bellow
 * When you have more then 2 models that share equal paramaters, the command:
 * delete options.models.at(j)->freeParams_[fj->first];
 * will try to delete the free parameter from the list more then once
 * this will cause an error of memory access violation
 * I solved this problem by creating a new map containing only the unique parameters
 */



	map<string,FreeParam*> uniqueList;
	int uniqueIndex = 0;
	//remove dulicate free parameters
	unsigned int i,j, k, cnt=0;


	if(options.models.size() > 1 ){

	for (i=0; i <options.models.size(); i++){  //k-1 first models
		//cout << "in the first loop" << endl;
		for (j=i+1; j <options.models.size(); j++){  //the subsequent model to the end
			//cout << "in the second loop" << endl;
			for (auto fi=options.models.at(i)->freeParams_.begin(); fi!=options.models.at(i)->freeParams_.end(); ++fi){
				//cout << "in the thrid loop" << endl;
				for (auto fj=options.models.at(j)->freeParams_.begin(); fj!=options.models.at(j)->freeParams_.end(); ++fj){
					//cout << "in the fourth loop" << endl;
					if(fi->first.compare(fj->first)==0){  //common free parameters
						//cout << "in the if if the fourth loop" << endl;
						cnt++;
						if(options.verbosity >= 4)
							cout<<"Common free parameters found.  Param:"<< fi->first <<"  model1:" <<options.models.at(i)->getName() <<"   model2:" << options.models.at(j)->getName() <<endl;
						//cout << "before deleting" << endl;

						if ( uniqueList.find(fj->first) == uniqueList.end() ) {
							cout << fj->first << endl;
							//delete options.models.at(j)->freeParams_[fj->first]; //delete free param objecyt in j
							uniqueList[fj->first] = fj->second;
							uniqueIndex++;
							//cout << "Found parameter " << fj->first << endl;
						  // found
						} else {
							  // not found
								//cout << "parameter already found: " << fj->first << endl;

						}
						//cout << "after delete" << endl;
						//options.models.at(j)->freeParams_[fj->first] = options.models.at(i)->freeParams_[fi->first];  //map to i
						fj->second = fi->second;
					}
				}
			}
		}
	}
	//cout << "Found " << uniqueIndex << " unique parameters" << endl;
	//Raquel changed the way that the function find unique parameters

	if(options.verbosity >= 3) cout<<"Removing duplicate free parameters is finished. "<<cnt<<" duplicate free parameters are deleted."<<endl;

	}else{

		uniqueList = options.models.at(0)->freeParams_;
	}



		for (i=0; i <options.models.size(); i++){
			options.models.at(i)->freeParams_ = uniqueList;
			if(options.verbosity >= 1){
				cout << "model " << i << " free parameters " << options.models.at(i)->freeParams_.size() << endl;
			}
		}


		options.freeParams_.clear(); cnt=0;
		options.freeParams_ = uniqueList;
		cnt = options.freeParams_.size();

	 //razi make full list of union of free parameters
	options.freeParams_.clear(); cnt=0;
	for (i=0; i <options.models.size(); i++){  //razi: all models, was k-1 first models, later test
		for (auto fi=options.models.at(i)->freeParams_.begin(); fi!=options.models.at(i)->freeParams_.end(); ++fi){
			if (options.freeParams_.count(fi->first) < 1){  //a new free parameters
				options.freeParams_.insert(make_pair(fi->first, fi->second)); cnt++;
				if(options.verbosity >= 3) cout<<"Free parameter:"<<fi->first <<" is added to the list"<<endl;
			}
		}
	}


	if(options.verbosity >= 3) cout<<"Free parameters found:"<< cnt << endl;


	//razi make mapping from models free parameters to the full list of free parameters


	std::map<unsigned int, unsigned int> dummay_map;
	for (i=0; i <options.models.size(); i++){  //razi: all models, was k-1 first models, later test
		j = 0;

		freeParamMapping.push_back(dummay_map);
		for (auto fj=options.models.at(i)->freeParams_.begin(); fj!=options.models.at(i)->freeParams_.end(); ++fj){
			k=0;
			for (auto fk=options.freeParams_.begin(); fk!=options.freeParams_.end(); ++fk){
//cout<<"Full List   FreeParam["<<k<<"] :" << fk->first<<endl; //mypause();
				if (fj->first == fk->first){
					freeParamMapping.at(i).insert(make_pair(j,k));
				}
				k++;
			}
			j++;
		}
	}


	if(options.verbosity >= 3)cout<<"Model consolidation is completed. "<<cnt<<" free parameters are found."<<endl;
}




void Swarm::setsConf(std::string sConf, int mid){  //razi added
	if (options.verbosity>=3) cout<<"Setting sConf id:" << mid << " File:" << sConf <<endl;
	int sz = sConf_.size();
	if (mid < 0){
		outputError("Error in setting sConf. Invalid id. Quitting ....");
	}else if ((mid >= 0) && (mid == sz)){
		sConf_.push_back(sConf);
	}else if ((mid >= 0) && (mid < sz)){
		sConf_.at(mid) = sConf;
	}else if ((mid >= 0) && (mid > sz)){

		for(int i = 0; i < mid - sz; i++){
			sConf_.push_back("");
		}
		sConf_.at(mid) = sConf;
	}

//	cout<<"set conf done, new size is:"<< sConf_.size() <<". Eneter a number to continue..."; cin>>i;
//	if( sConf_.size()>0)
//		for(i=0; i<sConf_.size(); i++)
//			cout<<" SConf[" << i <<"] is:" << sConf_.at(i)<<", ";
//
//	}
//	catch(string err){
//		cout<<"Error occurred in set conf Err:"<< err;
//	}
}

void Swarm::setExpPath(std::string prefixPath, std::string path, int mid){

/*   old ver, just one exp file for each model
	std::vector<string> paths;
	string expfile;

	expPaths_.clear();
	paths = split_string(path, ",");
	if (paths.size()>1)
		outputError("Can not set more than 1 exp file for each model at this version");

	for(int i=0; i< paths.size(); ++i){
		expfile = convertToAbsPath(paths.at(i));
//cout<<"Swarm::setExpPath AAA1  Enter a number..."; mypause();
		expPaths_.push_back(expfile);
//cout<<"Swarm::setExpPath AAA2, File:"<< expfile <<"  Enter a number..."; mypause();
		addExp(expfile, mid);
//cout<<"Swarm::setExpPath AAA3  Enter a number..."; mypause();

		if (options.verbosity >=3) cout<<"Exp["<<i<<"]  :"<<expfile <<" is added to the list of exp files.\n";
	}
*/

	std::vector<string> Paths;
	std::vector<string> absPaths;
	fs::path prefix(prefixPath);
	string expfile;
	Paths = split_string(path, ",");

	if (mid==-1){ //razi: means that each file is for one model [in the sqame order]

		for(unsigned int i=0; i< Paths.size(); ++i){
			fs::path path(Paths.at(i));
			fs::path abspath = path;
			
			if (!fs::exists(path.string())){
				abspath = prefix / path;
			}
			absPaths.clear();
			expfile = abspath.string();
			if(options.verbosity>=3){
				cout<<"*** exp file found"<<expfile<< " set for model:"<<i<<endl;
			}
			absPaths.push_back(expfile);
			expPaths_.insert(make_pair(i, absPaths));
			addExp(expfile, i);
			if (options.verbosity >=3) cout<<"Exp["<<i<<"]  :"<<expfile <<" is added to the list of exp files.\n";
		}


	}else{
		if((mid>=0)&&((unsigned) mid < options.models.size())){
			absPaths.clear();
			for(unsigned int i=0; i< Paths.size(); ++i){
				fs::path path(Paths.at(i));
				fs::path abspath = path;
				
				if (!fs::exists(path.string())){
					abspath = prefix / path;
				}
				expfile = abspath.string();
				absPaths.push_back(expfile);
				if(options.verbosity>=3){
					cout<<"*** exp file found"<<expfile<<endl;
				}
				if (options.verbosity >= 3){ cout<<"Swarm::setExpPath try to add  File:"<< expfile <<endl;}
				addExp(expfile, mid);
				if (options.verbosity >=3) cout<<"Exp["<<i<<"]  :"<<expfile <<" is added to the list of exp files for model:"<<mid<<".\n";
			}
			//cout<<"Number of exp files for this model: "<< absPaths.size()<<endl; cout<< "Number of expPaths_ before adding :"<< expPaths_.size()<<endl;
			expPaths_[mid]= absPaths; //razi: there is eeror for duplicate key, expPaths_.insert(make_pair(mid, absPaths));
		}
		else
			outputError("Invalid request to set exp files. I am quitting...");
	}
}



void Swarm::setModels(std::string pathPrefix, std::string path, bool overwrite){
	std::vector<string> paths;
	fs::path prefix(pathPrefix);

	if(overwrite)
		this->options.models.clear();

	if  (options.verbosity >= 3) cout<<"Swarm::setModels:  processing model files:"<<path<<endl;
	paths = split_string(path, ",");
	for(unsigned int i=0; i< paths.size(); ++i){
		fs::path path(paths.at(i));
		fs::path abspath = path;
		
		if (!fs::exists(path.string())){
		                abspath = prefix / path;
		}

		setModel(abspath.string(), i, overwrite);
		if (options.verbosity >=3) cout<<"Model["<<i<<"]  :"<<paths.at(i)<<" is added to the list of models.\n" << endl;
	}
	check_model_consistency(options.models);
	//this->options.model = this->options.models.at(0); //set the first model as the basic model for backward comaptibility
}

void Swarm::setModel(std::string path, int mid, bool overwrite){
	std::string modelfile = convertToAbsPath(path);
	unsigned int  nmodel=getNumModels();

	if ((unsigned) mid == nmodel){
		this->options.models.push_back(new Model(this, modelfile));  //add to the end
		if(options.verbosity >=4) cout<<" Model["<< mid<<"]: "<<path<<" is appended to the list of models."<<endl;
	}else if ((unsigned) mid < nmodel){ //model may exist
		if ((overwrite) || (! options.models.at(mid))){ //if no model is set or overwrite is allowed
			this->options.models[mid] = new Model(this, modelfile); //you may need to kill the overwritten model
			if(options.verbosity >=4) cout<<" Model["<< mid<<"]: "<<path<<" is added to the list of models."<<endl;
		}
	}else {    //include dummy
		for (unsigned int i=0; i < mid - nmodel; i++) //add dummy file
			this->options.models.push_back(0);
		this->options.models.push_back(new Model(this, modelfile));  //add to the end
		if(options.verbosity >=4) cout<<" Model["<< mid<<"]: "<<path<<" is added to the list of models after filling dummy models."<<endl;
	}
}



void Swarm::setDefaultModelIndex(int mid){
	if((mid>=0) &&((unsigned) mid < this->options.models.size())){
		//this->options.model = this->options.models.at(mid);
		options.defaultModel = mid;
		if (options.verbosity >= 3) cout<<"Default model not set, Index is set to " << mid << " # of models:"<< this->options.models.size()<<endl;
	}else outputError("Default model not set, Index is out of range:"+toString(mid)+" # of models:"+ toString((int)this->options.models.size()));
}




std::string Swarm::getExp(unsigned int modelId) {
	return getExp(modelId, 0);
}
std::string Swarm::getExp(unsigned int modelId, unsigned int expid) {
	std::string expPath = expPaths_[modelId][expid];
	return getFilename(expPath);
}



std::string Swarm::getModelName(unsigned int modelId, bool FullPath){
	std::string path = options.models.at(modelId)->getName();
	if (FullPath)
		return path;
	else
		return getFilename(path);
}

bool Swarm::check_model_exp_consistency(){
	unsigned int i,j;

	return true;

	//razi: this check is based on the assumption that there is exactly one exp file per model
	// this test can be avoided later
	if(this->expPaths_.size()!=this->options.models.size()){
		cout<<"Number of exp paths is : "<<this->options.models.size()<< endl;
		cout<<"Number of model paths is : "<<this->options.models.size()<< endl;

		for(i=0; i<this->options.models.size(); i++)
			cout<<"model ["<< i <<"] is: "<< this->options.models.at(i)->getName() <<endl;
		for(i=0; i<this->options.models.size(); i++)
			for(j=0; i<this->expPaths_[i].size(); j++)
			cout<<"exp["<< i <<","<<j<<"] is: "<< this->expPaths_[i].at(j)<<endl;

		return false;
		//this->outputError("The number of exp files is not equal to the number of exp files. quitting...");
	}else return true;
	//Todo, check other consistency checks
}

void Swarm::initFit () {  //razi: modified version
	int i, nModels=options.models.size();
	if (resumingSavedSwarm) {
		//razi: TODO modify this part
		cout << "Resuming fitting run. May need some modifcation to support multiple files. Enter a number..." << endl; cin>>i;

		/*
		for (auto particle = particleIterationCounter_.begin(); particle != particleIterationCounter_.end(); ++particle) {
			swarmComm->univMessageSender.push_back(toString(currentGeneration - 1));

			cout<<"Init fit, modify this line ...\n";
			for (unsigned int mid = 0; mid < options.models.size(); ++mid)
				swarmComm->sendToSwarm(0, particle->first, SEND_NUMFLIGHTS_TO_PARTICLE, false, swarmComm->univMessageSender);
			swarmComm->univMessageSender.clear();

			for (auto param = particleCurrParamSets_.at(particle->first).begin(); param != particleCurrParamSets_.at(particle->first).end(); ++param) {
				//cout << "Sending " << *param << " to " << particle->first << endl;
				swarmComm->univMessageSender.push_back(toString(*param));
			}
			swarmComm->sendToSwarm(0, particle->first, SEND_FINAL_PARAMS_TO_PARTICLE, false, swarmComm->univMessageSender);
			swarmComm->univMessageSender.clear();
		}

		// We need to set current generation to the iteration counter because the saved swarm
		// currentGeneration may not be accurate due to it changed between runGeneration() and
		// breedGenerationGA()
		if (options.fitType == "ga") {
			currentGeneration = currentGeneration - 1;
		}

		*/
	}
	else {
		if (options.verbosity >= 3) {
			cout << "Initializing fitting run" << endl;
		}

		// If we are running ODE, we want to generate a network and network reader.
		// The network file has parameters that are replaced in every generation,
		// and the network reader is used to read those network files

		//razi, generate subparticles, set particle id to 1;
		unsigned int ParID = 1; unsigned int subParId;

		Particle p = Particle(this, 1); bool genNetworkFlag = false;
		vector<subParticle *> sp(nModels);
		for (int mid=0; mid <nModels; mid++){
			if (options.models.at(mid)->getHasGenerateNetwork() && bootstrapCounter == 0) {

				genNetworkFlag = true;
				if (options.verbosity >= 3) {
					cout << "Creating dummy .bngl, running to get a .net file" << endl;
				}

				// First create a particle which will generate the network
				// and fill it with dummy params
				p.setModel(options.models[mid], mid);

				subParId= fcalcsubParID(ParID, mid, nModels); //razi: calc subparticle id from model index
				sp[mid] = new subParticle(&p, mid, subParId); //razi: generate subparticles
				//sp[mid]->setModel(mid);

			}
		}
		if (genNetworkFlag){
			p.generateParams();

			// Output model file with dummy parameters. This file will be used to generate
			// our initial network
			//razi, was options.model->outputModelWithParams(p.getParams(), options.jobOutputDir, "base.bngl", "", false, false, false, false, false);
			for (unsigned int mid = 0; mid < options.models.size(); ++mid){
				options.models[mid]->outputModelWithParams(p.getParams(mid), options.jobOutputDir, p.swarm_->getModelName(mid, false)+"_base.bngl", "", false, false, false, false, false);

				// Construct our simulation command and run the network generator

				string modelPath;
				if (options.jobOutputDir.substr(options.jobOutputDir.size()-1,1)=="/")   //razi added to avoid "//" in windows
					modelPath = options.jobOutputDir + p.swarm_->getModelName(mid, false)+"_base.bngl";
				else
					modelPath = options.jobOutputDir + "/" + p.swarm_->getModelName(mid, false)+"_base.bngl";

				string command = options.bngCommand + " --outdir " + options.jobOutputDir + " " + modelPath + " > " + options.jobOutputDir + "netgen_output 2>&1";
				if(options.verbosity>=1){
					cout << "Generating initial .net file with command: " << command << endl;
				}
				if (options.useCluster) {
					if(options.clusterSoftware == "slurm"){
						command = generateSlurmCommand(command, false, 1);
					}

				}

				int ret; ret= runCommand(command);
				if (ret != 0){
#ifdef PC_VER  //This is only for debugging in windows, the return value is 0 due to some error in BNG2.pl, later to be checked
					cout<<"Warning: Net file may not be generated properly: " <<command<<endl;
#else
					outputError("Error(32): Couldn't generate initial .net file with command: " + command + ". Quitting.");
#endif
				}else{
						// Now that we have a .net file, replace the full .bngl with a .bngl
					// containing ONLY action commands and a .net file loader. This is
					// used later to run our .net files.
					options.models[mid]->outputModelWithParams(p.getParams(mid), options.jobOutputDir, p.swarm_->getModelName(mid, false)+"_base.bngl", "", false, true, false, false, false);

					// Now store our .net file for later use
					//razi, was string netPath = options.jobOutputDir + "/base.net";
					string netPath = options.jobOutputDir + "/"  + p.swarm_->getModelName(mid, false)+"_base.net";
					//netPath=removeLastSlash(netPath);  //razi added for windows, already taken care of
					options.models[mid]->parseNet(netPath);
				}
			}
		}


/*
#if(defined(PC_VER))  //This is only for debugging in windows, the return value is 0 due to some error in BNG2.pl, later to be checked
		// random simulation is requested, just copy the exp file plus some noise into gdat files
		if (options.bngCommand.find("random")!=std::string::npos){
			ret=0;

		}
#endif
*/


		vector<double> params;
		int subPar = 0; //Raquel added to support subparticles
		for (unsigned int param = 0; param < getNumFreeParams(); ++param) {  //Razi: All free params tat appear in some models
			params.push_back(0);
		}

		for (int particle = 1; particle <= options.swarmSize; ++particle) {
			particleCurrParamSets_.insert(pair<int, vector<double> >(particle, params));

			//razi: init for each subparticle
			map<unsigned int, vector<double> > zmap;
			for (unsigned int mid = 0; mid < options.models.size(); ++mid){
				vector<double> zvec(options.models[mid]->getNumFreeParams());	//razi: dummy vector t contain free parameter values
				zmap.insert(make_pair(mid, zvec));
			}
			subparticleCurrParamSets_.insert(make_pair(particle, zmap));

		}


		if (options.fitType == "pso") {
			for (int particle = 1; particle <= options.swarmSize; ++particle) {
				particleParamVelocities_[particle] = params;
				particleBestParamSets_[particle] = params;
				particleWeights_[particle] = 0;
				for (int mid=0; mid<nModels; mid++)
					particleIterationCounter_[particle][mid] = 0;
			}
		}
		if (options.fitType == "pso" || "de") {
			int subPar = 0;

			for (int particle = 1; particle <= options.swarmSize; ++particle) {
				particleBestFits_[particle] = 0;

				//Raquel: added support for subparticles
				for(int mid = 0; (unsigned) mid < options.models.size(); mid++){
					subPar = fcalcsubParID(particle, mid, options.models.size());
					subparticleBestFits_[subPar] = 0;
				}//Raquel: added support for subparticles

			}
		}
		if(options.fitType == "de") {
			for (unsigned int island = 0; island <= options.numIslands; ++island) {
				islandToParticle_.push_back(vector<unsigned int>(options.swarmSize / options.numIslands,0));
			}
			for (int particle = 1; particle <= options.swarmSize; ++particle) {
				particleToIsland_.push_back(0);
				particleBestFitsByFit_.insert(pair<double, unsigned int>(0, particle));

				//Raquel added this loop to support subparticles
				for(unsigned int mid = 0; mid < options.models.size(); mid++){
					subPar = fcalcsubParID(particle, mid, options.models.size());
					subparticleBestFitsByFit_.insert(pair<double, unsigned int>(0, subPar));

				}//Raquel added this loop to support subparticles

			}
		}
		if (options.fitType == "sa") {
			for (int particle = 1; particle <= options.swarmSize; ++particle) {
				particleBestFitsByFit_.insert(pair<double, unsigned int>(0, particle));

				//Raquel added this loop to support subparticles
				for(unsigned int mid = 0; mid < options.models.size(); mid++){
					subPar = fcalcsubParID(particle, mid, options.models.size());
					subparticleBestFitsByFit_.insert(pair<double, unsigned int>(0, subPar));
					currentsubparticleBestFitsByFit_.insert(pair<double, unsigned int>(0, subPar));

				}//Raquel added this loop to support subparticles


			}


		}



	}
}

void Swarm::addMutate(std::string mutateString) {
	vector<string> mutComponents;
	split(mutateString, mutComponents);

	// Make sure we have three components to work with
	if (mutComponents.size() == 3) {
		// Make sure first parameter name exists as a free parameter
		if (options.freeParams_.count(mutComponents[0]) > 0 || mutComponents[0] == "default") {
			// Make sure 2nd and 3rd components are numeric
			if (isFloat(mutComponents[1]) && isFloat(mutComponents[2])) {
				if (mutComponents[0] != "default") {
					//cout << "not default" << endl;
					options.freeParams_.at(mutComponents[0])->setMutationRate(stof(mutComponents[1]));
					options.freeParams_.at(mutComponents[0])->setMutationFactor(stof(mutComponents[2]));
					options.freeParams_.at(mutComponents[0])->setHasMutation(true);

					//cout << "setting " << mutComponents[0] << " to " << mutComponents[1] << ":" << mutComponents[2] << endl;
				}
				else {
					for (map<string,FreeParam*>::iterator fp = options.freeParams_.begin(); fp != options.freeParams_.end(); ++fp) {
						//cout << "setting " << mutComponents[0] << " to " << mutComponents[1] << ":" << mutComponents[2] << endl;

						fp->second->setMutationRate(stof(mutComponents[1]));
						fp->second->setMutationFactor(stof(mutComponents[2]));
						fp->second->setHasMutation(true);
					}
				}
			}
			else {
				outputError("Error: Problem parsing the mutation option in your .conf file. The mutation rate and/or factor are non-numeric.");
			}
		}
		else {
			cout << "Warning: We found a mutation option '" << mutComponents[0] << "' in your .conf file, but don't see a matching free parameter specification in your model file. We will ignore this mutation factor." << endl;
		}
	}
	else {
		outputError("Error: Problem parsing a mutation option in your .conf file. Each mutation option requires three components: the parameter name, mutation rate, and mutation factor.");
	}
}




vector<double> Swarm::calcQPSOmBests() {
	// Eq 2b Liu et al
	//cout<<"Modification needed\n...."; mypause();

	vector<double> mBests;

	//for (unsigned int param = 0; param < options.model->getFreeParams_().size(); ++param) {
	//for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) { //for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
	unsigned int i = 0;



		double sum = 0;
		for (int particle = 1; particle <= options.swarmSize; ++particle) {
			i = 0;
			for (auto param = particleCurrParamSets_.at(particle).begin(); param != particleCurrParamSets_.at(particle).end(); ++param) {

//				for(int mid=0; mid < options.models.size(); mid++){

					sum += abs(particleBestParamSets_.at(particle)[i]);

	//			}
				++i;
			}



			double mean = sum / (double)options.swarmSize;

			// TODO: Adaptive mutation ala Liu 2005
			if (options.mutateQPSO) {
				boost::random::cauchy_distribution<double> dist(0, cauchyMutator_);
				double mutator = dist(generalRand);
				mean += mutator;
				//cout << mean << " mutated by: " << mutator << endl;
			}
			mBests.push_back(sum / (double)options.swarmSize);

		}


	//}

	return mBests;

}




void Swarm::setfitType(string type) {
	this->options.fitType = type;

	if (options.fitType == "swarm") {
		options.parallelCount = options.swarmSize;
	}
}

void Swarm::setJobOutputDir(string path) {
	options.jobOutputDir = path;

	if (!checkIfFileExists(options.outputDir)) {
		string cmd = "mkdir " + options.outputDir;
		//cout << "running: " << cmd << endl;
		if(runCommand(cmd) != 0) {
			outputError("Error: Couldn't create base output directory with command: " + cmd + ". Quitting.");
		}
	}

	if (checkIfFileExists(options.jobOutputDir)) {
		string input;
		string answer;
		while (1) {

			//Raquel fixing MPI support issues
			if(options.clusterSoftware != "mpi" && options.clusterSoftware != "BNF2mpi"){
				cout << "Warning: Your output directory " << options.jobOutputDir << " already exists. Overwrite? (Y or N) ";
				getline(cin, input);
				stringstream myInp(input);
				myInp >> answer;
			}

			if(options.clusterSoftware == "mpi" || options.clusterSoftware == "BNF2mpi"){
				answer = "Y";
			}

			if (answer == "Y" || answer == "y") {
				string cmd = "rm -r " + options.jobOutputDir;
				if (options.verbosity >= 1) {
					cout << "Deleting old output directory, this may take a few minutes..." << endl;
				}

				if(runCommand(cmd) != 0) {
					outputError("Error: Couldn't delete existing output directory with the command: " + cmd + ". Quitting.");
				}
				break;
			}
			else if (answer == "N" || answer == "n") {
				outputError("Error: Output directory already exists. Quitting.");
			}
		}
	}
	string cmd = "mkdir " + options.jobOutputDir;
	//cout << "else running: " << cmd << endl;
	//int ret = system(cmd.c_str());
	if(runCommand(cmd) != 0) {
		outputError("Error: Couldn't create output directory with command: " + cmd + ". Quitting.");
	}
}



void Swarm::doSwarm() {
	if(options.verbosity >=3) cout<<"doSwarm invoked: bootsrtarps:"<< options.bootstrap<<" \n";
	for (bootstrapCounter = 0; bootstrapCounter <= options.bootstrap; ++bootstrapCounter) {
		initFit();
		// Main fit loops

		if (options.synchronicity) {
			if (options.fitType == "ga") {
				runSGA();
				if(options.verbosity>=3){
					cout << "Finished GA" << endl;
				}
			}
			else if (options.fitType == "pso") {
				//cout << "RAQUEL BEFORE PSO" << endl;
				runSPSO();
				//cout << "RAQUEL After PSO" << endl;

			}
			else if (options.fitType == "de") {
				runSDE();
			}
			else if (options.fitType == "sa") {
				runSSA();
			}
		}

		// Asynchronous fit loops
		else {
			// Genetic fit
			if (options.fitType == "ga") {
				runAGA();
			}
			// PSO fit
			else if (options.fitType == "pso") {
				runAPSO();
			}
			else if (options.fitType == "de") {
				runADE();
			}
			else if (options.fitType == "sa") {
				runASA();
			}
		}
		if(options.verbosity>=4){

			cout << "RAQUEL: starting finishFit function" << endl;

		}
		finishFit();
		if(options.verbosity>=4){

			cout << "RAQUEL: done finishFit function" << endl;
		}

	}

}

bool Swarm::checkStopCriteria() {
	if (options.verbosity >= 3) {
		//cout << "Checking stop criteria. Flight count is " << flightCounter_ << " and max is " << options.maxNumSimulations << endl;
		//cout << "Fittype is " << options.fitType << endl;
	}

	if (options.fitType == "pso") {
		if (options.enhancedStop) {
			// Eq 4 from Moraes et al
			// Algorithm 1 from Moraes et al
			if (abs(particleBestFitsByFit_.begin()->first - optimum_) < options.absTolerance + (options.relTolerance * particleBestFitsByFit_.begin()->first) ) {
				if (options.verbosity >= 3) {
					cout << "Incrementing permanence counter" << endl;
				}

				++permanenceCounter_;
				++inertiaUpdateCounter_;

				double quotient = ((double)permanenceCounter_ / options.nmin);
				// Check to see if quotient is a real number
				if (floor(quotient) == quotient) {
					double oldAvgPos = weightedAvgPos_;
					updateEnhancedStop();
					if ( getEuclidianNorm(weightedAvgPos_ - oldAvgPos, this->getNumFreeParams()) < options.absTolerance ) {
						// Fit finished!

						if (options.verbosity >= 3) {
							cout << "Stopped according to enhanced stop criteria" << endl;
						}
						return true;
					}
				}
			}
			else {
				permanenceCounter_ = 0;
				if (particleBestFitsByFit_.begin()->first < optimum_) {
					optimum_ = particleBestFitsByFit_.begin()->first;
				}

			}

		}

		//cout << "RAQUEL: MAX GENERATION: " << options.maxGenerations << endl;
		if (options.maxGenerations && currentGeneration > options.maxGenerations) {

			cout << "Stopped because currentGeneration > maxGenerations" << endl;
			return true;
		}
		/*
		else if (options.nmax) {
			if (abs(particleBestFitsByFit_.begin()->first - optimum_) < options.absTolerance + (options.relTolerance * particleBestFitsByFit_.begin()->first) ) {
				if (options.verbosity >= 3) {
					cout << "Incrementing permanence counter. Optimum is " << optimum_ << " and f(g) is " << particleBestFitsByFit_.begin()->first << ". Difference is " << abs(particleBestFitsByFit_.begin()->first - optimum_) << " and conditional is " << options.absTolerance + (options.relTolerance * particleBestFitsByFit_.begin()->first) << endl;
				}

				++permanenceCounter_;
				++inertiaUpdateCounter_;

				if (permanenceCounter_ >= options.nmax) {

					if (options.verbosity >= 3) {
						cout << "Stopped according to permanence > n_max" << endl;
					}
					return true;
				}
			}
			else {
				permanenceCounter_ = 0;
				if (particleBestFitsByFit_.begin()->first < optimum_) {
					optimum_ = particleBestFitsByFit_.begin()->first;
				}
			}
		}
		 */
	}//	else if ((options.fitType == "ga" || options.fitType == "de") && options.synchronicity) {
	else if ( (options.fitType == "ga" || options.fitType == "de") ) {
		if (options.verbosity >= 3) {
			cout << "Checking if we've reached max generation. Current is " << (currentGeneration - 1) << " and max is " << options.maxGenerations << endl;
		}
		if (options.maxGenerations && currentGeneration > options.maxGenerations) {
			cout << "Stopped because currentGeneration > maxGenerations" << endl;
			return true;
		}
	}

	// If we've run more than the maximum number of simulations
	if (flightCounter_ >= options.maxNumSimulations) {
		if (options.verbosity >= 3) {
			cout << "Stopped according to flightCounter (" << flightCounter_ << ") >= maxNumSimulations (" << options.maxNumSimulations << ")" << endl;
		}
		return true;
	}

	if ( (swarmBestFits_.size() && swarmBestFits_.begin()->first <= options.minFit) || (particleBestFitsByFit_.size() && particleBestFitsByFit_.begin()->first <= options.minFit) ) {
		cout << "Stopped according to swarmBestFit (" << particleBestFitsByFit_.begin()->first << ") <= options.minFit (" << options.minFit << ")" << endl;
		return true;
	}

	// TODO: Check if all particles have converged to same solution

	/*
			// TIME
			if () {
			}
	 */
	//cout << "RAQUEL: Returning FALSE to Stop Criteria." << endl;
	return false;
}

void Swarm::updateEnhancedStop() {
	// First calculate particle weights

	if (options.verbosity >= 3) {
		cout << "Updating enhanced stop" << endl;
	}

	updateParticleWeights();
	weightedAvgPos_ = calcWeightedAveragePosition();
}

double Swarm::getEuclidianNorm(double y, unsigned int n) {

	// Eq 10 in Moraes at al
	double sum = 0;
	for (unsigned int i = 1; i <= n; ++i) {
		sum += pow(y, 2);
	}

	double norm = sqrt( (1/n) * sum);

	return norm;
}

void Swarm::updateParticleWeights() {

	if (options.verbosity >= 3) {
		cout << "Updating particle weights" << endl;
	}

	for (int p = 1; p <= options.swarmSize; ++p) {
		particleWeights_[p] = calcParticleWeight(p);
	}
	if (options.verbosity >= 3) {
		cout << "done here" << endl;
	}
}

double Swarm::calcWeightedAveragePosition() {

	if (options.verbosity >= 3) {
		cout << "Calculating weighted average position" << endl;
	}

	double sum = 0;
	for (int p = 1; p <= options.swarmSize; ++p) {

		if ((unsigned) p == particleBestFitsByFit_.begin()->second) {
			continue;
		}

		sum += particleWeights_.at(p) * particleBestFits_[p];
	}

	return sum;
}

double Swarm::calcParticleWeight(unsigned int particle) {
	// Get reciprocal of euclidian norm of the difference between swarm best fit and particle best fit
	double numerator = 1 / getEuclidianNorm( (particleBestFits_[particle] - particleBestFitsByFit_.begin()->first), this->getNumFreeParams());

	// Eq 8 in Moraes et al
	double sum = 0;
	for (int i = 1; i <= options.swarmSize; ++i) {
		// Make sure we're not using the particle with the best fit -- it will
		// result in a div_by_0
		if ((unsigned)i == particleBestFitsByFit_.begin()->second) {
			continue;
		}

		// Add 1/euclidian
		sum += 1 / numerator;
	}

	double weight = numerator / sum;

	return weight;
}

Particle * Swarm::createParticle(unsigned int pID) {
	Particle * p = new Particle(this, pID);

	return p;
}

vector<vector<unsigned int> > Swarm::generateTopology(unsigned int populationSize) {
	Timer tmr;

	vector<vector<unsigned int> > allParticles (populationSize + 1);

	if (options.fitType == "pso" || options.fitType == "de") {

		if (options.verbosity >= 3) {
			cout << "Generating initial particles with a " << options.topology << " topology" << endl;
		}

		if (options.topology == "fullyconnected") {
			for (unsigned int p = 1; p <= populationSize; ++p) {
				vector<unsigned int> connections;
				for (unsigned int c = 1; c <= populationSize; ++c) {
					if (c != p) {
						connections.push_back(c);
					}
				}
				allParticles[p] = connections;
			}
		}
		else if (options.topology == "ring") {

			// Connect the first particle manually
			unsigned int firstConnection[] = {2, populationSize};
			allParticles[1] = vector<unsigned int> (firstConnection, firstConnection + sizeof(firstConnection) / sizeof(firstConnection[0]));

			vector<unsigned int> connections;
			for (unsigned int p = 2; p <= populationSize - 1; ++p) {

				connections.clear();
				// Connect to particle before, and particle after
				connections.push_back(p-1);
				connections.push_back(p+1);
				allParticles[p] = connections;

			}

			// Connect the last particle manually
			unsigned int lastConnection[] = {populationSize - 1, 1};
			allParticles[populationSize] = vector<unsigned int> (lastConnection, lastConnection + sizeof(lastConnection) / sizeof(lastConnection[0]));
		}
		else if (options.topology == "star") {
			vector<unsigned int> connections;

			// First connect the central particle to all others
			for (unsigned int c = 2; c <= populationSize; ++c) {
				connections.push_back(c);
			}
			allParticles[1] = connections;

			// Then connect all particles to the central particles
			for (unsigned int p = 2; p <= populationSize; ++p) {
				connections.clear();
				connections.push_back(1);
				allParticles[p] = connections;
			}
		}
		else if (options.topology == "mesh") {

			unsigned int desiredArea = populationSize;
			unsigned int divisor = ceil(sqrt(desiredArea));
			if (divisor != 0) {
			  while(desiredArea % divisor != 0) {
				  ++divisor;
			  }
			}
			unsigned int length = divisor;
			unsigned int width = desiredArea / divisor;

			// Construct a matrix of dimension length x width
			// and fill it with particles
			int p = 0;
			vector<vector<unsigned int> > matrix(length, vector<unsigned int>(width));
			for (unsigned int x = 0; x < length; ++x) {
				for (unsigned int y = 0; y < width; ++y) {
					matrix[x][y] = ++p;
				}
			}

			// Make our connections
			vector<unsigned int> connections;
			for (unsigned int x = 0; x < length; ++x) {
				for (unsigned int y = 0; y < width; ++y) {
					connections.clear();
					p = matrix[x][y];

					// If we're on the left bound
					if (x == 0) {
						// Add the particle to the right of us
						connections.push_back(matrix[x+1][y]);
					}
					// If we're on the right bound
					else if (x == length - 1) {
						// Add the particle to the left of us
						connections.push_back(matrix[x-1][y]);
					}
					// If we're in a center
					else {
						// Add the particles on either side of us
						connections.push_back(matrix[x-1][y]);
						connections.push_back(matrix[x+1][y]);
					}

					if (y == 0) {
						connections.push_back(matrix[x][y+1]);
					}
					else if (y == width - 1) {
						connections.push_back(matrix[x][y-1]);
					}
					else {
						connections.push_back(matrix[x][y-1]);
						connections.push_back(matrix[x][y+1]);
					}
					allParticles[p] = connections;
				}
			}
		}
		else if (options.topology == "toroidal") {
			unsigned int desiredArea = populationSize;
			unsigned int divisor = ceil(sqrt(desiredArea));
			if (divisor != 0) {
			  while(desiredArea % divisor != 0) {
				  ++divisor;
			  }
			}
			unsigned int length = divisor;
			unsigned int width = desiredArea / divisor;

			// Construct a matrix of dimensions length x width
			// and fill it with particles
			unsigned int p = 0;
			vector<vector<unsigned int> > matrix(length, vector<unsigned int>(width));
			for (unsigned int x = 0; x < length; ++x) {
				for (unsigned int y = 0; y < width; ++y) {
					matrix[x][y] = ++p;
				}
			}

			// Make our connections
			vector<unsigned int> connections;
			for (unsigned int x = 0; x < length; ++x) {
				for (unsigned int y = 0; y < width; ++y) {
					connections.clear();
					p = matrix[x][y];

					if (x == 0) {
						connections.push_back(matrix[x+1][y]);
						connections.push_back(matrix[length-1][y]);
					}
					else if (x == (length - 1)) {
						connections.push_back(matrix[x-1][y]);
						connections.push_back(matrix[0][y]);
					}
					else {
						connections.push_back(matrix[x-1][y]);
						connections.push_back(matrix[x+1][y]);
					}

					if (y == 0) {
						connections.push_back(matrix[x][y+1]);
						connections.push_back(matrix[x][width-1]);
					}
					else if (y == (width - 1)) {
						connections.push_back(matrix[x][y-1]);
						connections.push_back(matrix[x][0]);
					}
					else {
						connections.push_back(matrix[x][y-1]);
						connections.push_back(matrix[x][y+1]);
					}
					allParticles[p] = connections;
				}
			}
		}
		else if (options.topology == "tree") {
			unsigned int usedParticles = 1;
			unsigned int numLevels = 1;
			unsigned int previousLevel = 1;
			unsigned int currentLevel;

			// Determine number of levels in tree
			while (usedParticles < populationSize) {
				currentLevel = previousLevel*2;
				usedParticles += currentLevel;
				previousLevel = currentLevel;
				++numLevels;
			}

			vector<vector<unsigned int> > tree(numLevels, vector<unsigned int>());
			usedParticles = 1;
			currentLevel = 1;
			unsigned int currentParticle = 2;
			bool doneFilling = false;

			// Construct tree and fill it with particles
			tree[0].push_back(1);

			// For each level in tree
			for (unsigned int level = 1; level <= numLevels; ++level) {
				// Current level's particle count is double that of previous level
				currentLevel = currentLevel * 2;

				// For each slot in current level
				for (unsigned int i = 0; i < currentLevel; ++i) {
					// Fill slot with particle
					tree[level].push_back(currentParticle);
					++currentParticle;

					// Make sure we're not filling past our swarm size
					if (currentParticle == populationSize + 1) {
						doneFilling = true;
						break;
					}
				}
				if (doneFilling) {
					break;
				}
			}

			// For each level in tree
			for (unsigned int level = 1; level < numLevels; ++level) {
				int prevGroupCounter = 0;
				int particleCounter = 0;

				// For each particle in level
				for (unsigned int p = 0; p < tree[level].size(); ++p) {

					// Connect current particle (p) with proper particle in last level
					allParticles[tree[level][p]].push_back(tree[level-1][prevGroupCounter]);

					// Connect particle in last level to current particle (p)
					allParticles[tree[level-1][prevGroupCounter]].push_back(tree[level][p]);

					// If our particle counter reaches two, we're in a new pair in (or new particle
					// in previous level)
					if (++particleCounter == 2) {
						particleCounter = 0;
						++prevGroupCounter;
					}
				}
			}
		}
		else {
			outputError("Error: BioNetFit did not recognize your specified topology of: " + options.topology);
		}
	}

	/*
	int p = 0;
	for (auto o = allParticles.begin(); o != allParticles.end(); ++o) {
		cout << p << " is connected to:" << endl;
		for (auto i = (*o).begin(); i != (*o).end(); ++i) {
			cout << *i << endl;
		}
		++p;
	}
	 */

	double t = tmr.elapsed();
	if (options.verbosity >= 3) {
		cout << "Particle creation took " << t << " seconds" << endl;
	}

	return allParticles;
}

void Swarm::processParticlesPSO(vector<unsigned int> particles, bool newFlight) {

	map<unsigned int, std::vector<double> > particleNewParamSets;

	// Create the output directory for the next generation
	if (!checkIfFileExists(options.jobOutputDir + toString(currentGeneration))) {
		string createDirCmd = "mkdir " + options.jobOutputDir + toString(currentGeneration);
		int retryCounter = 0;
		while (runCommand(createDirCmd) != 0 && !checkIfFileExists(options.jobOutputDir + toString(currentGeneration))) {
			if(++retryCounter >= 100) {
				outputError("Error: Couldn't create " + options.jobOutputDir + toString(currentGeneration) + " to hold next generation's output.");
			}
			//sleep(1);
			cout << "Trying again to create dir" << endl;
		}
	}
	else if (options.verbosity >= 4) cout<<"Swarm::breedGenerationGA-AAA4, dir already exist\n.";


	if (options.verbosity >= 3) {
		cout << "Processing " << particles.size() << " particles" << endl;
	}

	vector<double> mBests;
	if (options.psoType == "qpso") {
		mBests = calcQPSOmBests();
	}

	// For each particle in our particle set
	for (auto particle = particles.begin(); particle != particles.end(); ++particle) {
		if (options.verbosity >= 3) {
			cout << "Processing particle " << *particle << endl;
		}

		// Calculate the next iteration positions according
		// to user preference
		if (options.psoType == "pso") {
			if (options.enhancedInertia) {
				updateInertia();
			}
			//cout << "RAQUEL before calcParticlePosPSO" << endl;

			particleCurrParamSets_[*particle] = calcParticlePosPSO(*particle);

			//for(int mid = 0; mid < options.models.size(); mid++){
			//	update_cur_particle_params(*particle, mid, true);
			//}

			//cout << "RAQUEL after calcParticlePosPSO size: " << particleCurrParamSets_[*particle].size() << endl;

		}
		else if (options.psoType == "bbpso") {
			if (options.verbosity >= 3) {
				cout << "RAQUEL before calcParticlePosBBPSO" << endl;
			}
			particleCurrParamSets_[*particle] = calcParticlePosBBPSO(*particle);
			if (options.verbosity >= 3) {
				cout << "RAQUEL after calcParticlePosBBPSO" << endl;
			}
		}
		else if (options.psoType == "bbpsoexp") {
			particleCurrParamSets_[*particle] = calcParticlePosBBPSO(*particle, true);
		}
		else if (options.psoType == "qpso") {
			updateContractionExpansionCoefficient();
			particleCurrParamSets_[*particle] = calcParticlePosQPSO(*particle, mBests);
		}

		if (options.verbosity >= 3) {
			cout << "Got next positions for particle " << *particle << endl;
		}
		// Convert positions to string so they can be sent to the particle
		vector<string> nextPositionsStr;
		vector<double> nextPositionsDouble;
		for (auto param = particleCurrParamSets_.at(*particle).begin(); param != particleCurrParamSets_.at(*particle).end(); ++param) {
			//nextPositionsStr.push_back(toString(static_cast<long double>(*param)));
			nextPositionsStr.push_back(toString(*param));
			nextPositionsDouble.push_back(*param);
			//cout << "Param " << *param << " size: " << particleCurrParamSets_.at(*particle).size() << endl;
			//cout << "Total size: " << particleCurrParamSets_.size() << endl;
		}

		//particleCurrParamSets_[*particle] = nextPositionsDouble;

		//subparticleCurrParamSets_.insert(make_pair(*particle, zmap));

		if (options.verbosity >= 3) {
			cout << "Sending new positions to " << *particle << endl;
		}

		// Replace any changed param sets in the master set
	//	for (auto child = particleCurrParamSets_.begin(); child != particleNewParamSets.end(); ++child) {
	//		cout << "param " << child<-first << endl;
	//	}

		// Finally, send the parameters

		//swarmComm->sendToSwarm(0, *particle, SEND_FINAL_PARAMS_TO_PARTICLE, false, nextPositionsStr);
		//for(int mid = 0; mid < options.models.size(); mid++){
		//	update_cur_particle_params(*particle, mid, true);
		//}

		int sp = 0;
		//Raquel: had to change the line above to send the message to subparticles
		if (options.verbosity >= 3) {
			cout << "RAQUEL total models " << options.models.size() << endl;
		}
		for(unsigned int i = 0; i < options.models.size(); i++){

			sp = fcalcsubParID(*particle, i, options.models.size());
			if (options.verbosity >= 3) {
				cout << "RAQUEL: sending message to subPar: " << sp << " model " << i << " particle " << *particle << endl;
			}
			swarmComm->sendToSwarm(0, sp, SEND_FINAL_PARAMS_TO_PARTICLE, false, nextPositionsStr);
			if (options.verbosity >= 3) {
				cout << "RAQUEL: sent message to subPar: " << sp << endl;
			}

			subparticleCurrParamSets_[*particle][i] = nextPositionsDouble;

			//for(auto param = nextPositionsStr.begin(); param != nextPositionsStr.end(); ++param) {
			//	particleNewParamSets[sp].push_back(stod(*param));
			//}
		}



		// Replace any changed param sets in the master set
		//for (auto child = particleNewParamSets.begin(); child != particleNewParamSets.end(); ++child) {
		//	particleCurrParamSets_[child->first] = child->second;
		//}



	}
}

vector<double> Swarm::calcParticlePosPSO(unsigned int particle) {

	if (options.verbosity >= 3) {
		cout << "Calculating velocity and position of particle " << particle << endl;
	}

	// This vector holds the new positions to be sent to the particle
	vector<double> nextPositions(particleCurrParamSets_.at(particle).size());
	vector<double> nextVelocities(particleCurrParamSets_.at(particle).size());
	if (options.verbosity >= 4) {
		cout << "RAQUEL getting best neighbor inside calcParticlePosPSO" << endl;
	}
	// Get the best positions for particle's neighborhood
	vector<double> neighborhoodBestPositions = getNeighborhoodBestPositions(particle);
	if (options.verbosity >= 4) {
		cout << "RAQUEL done best neighbor inside calcParticlePosPSO" << endl;
	}

	int i = 0;

	//cout << "inertia: " << options.inertia << endl;
	//cout << "cognitive: " << options.cognitive << endl;
	//cout << "social: " << options.social << endl << endl;
	if (options.verbosity >= 4) {
		cout << "RAQUEL entering for loop particleCurrParamSets_ " << endl;
	}
	// For each parameter in the current parameter set
	for (auto param = particleCurrParamSets_.at(particle).begin(); param != particleCurrParamSets_.at(particle).end(); ++param) {
		if (options.verbosity >= 3) {
			cout << "before " << *param << endl;
		}
		// Set up formula variables
		double currVelocity = options.inertia * particleParamVelocities_.at(particle)[i];
		//cout << "cv: " << currVelocity << endl;
		double r1 = ((double) rand() / (RAND_MAX));
		//cout << "r1: " << r1 << endl;
		double r2 = ((double) rand() / (RAND_MAX));
		//cout << "r2: " << r2 << endl;
		double personalBestPos = particleBestParamSets_.at(particle)[i];
		//cout << "pb: " << personalBestPos << endl;
		double currPos = particleCurrParamSets_.at(particle)[i];
		if (options.verbosity >= 3) {
			cout << "cp: " << currPos << endl;
		}
		//cout << "nbp: " << neighborhoodBestPositions[i] << endl;
		// Set velocity

		// TODO: Look into constriction factor - Clerc
		// Also, according to Eberhart and Shi, constriction factor + velocity clamping
		// results in fastest convergence
		double nextVelocity = (options.inertia * currVelocity) + options.cognitive * r1 * (personalBestPos - currPos) + options.social * r2 * (neighborhoodBestPositions[i] - currPos);
		if (options.verbosity >= 3) {
			cout << "nv: " << nextVelocity << endl;
		}
		// Set velocity
		particleParamVelocities_.at(particle)[i] = nextVelocity;
		// Set position
		nextPositions[i] = currPos + nextVelocity;

		/*
		if (nextPosition <= 0) {
			cout << "NP less than 0. Setting to Very Small Number" << endl;
			nextPosition = 0.00000001;
			particleParamVelocities_.at(particle)[i] = 0;
		}
		 */
		if (options.verbosity >= 3) {
			cout << "after " << nextPositions[i] << endl << endl;
		}
		++i;
	}

	if (options.verbosity >= 3) {
		cout << "RAQUEL done for loop particleCurrParamSets_ " << endl;
	}

	return nextPositions;
}

vector<double> Swarm::calcParticlePosBBPSO(unsigned int particle, bool exp) {

	if (options.verbosity >= 3) {
		cout << "Calculating BBPSO for particle " << particle << endl;
	}

	// Get the best positions for particle's neighborhood
	vector<double> neighborhoodBestPositions = getNeighborhoodBestPositions(particle);

	/*
	cout << "nbps:" << endl;
	for (auto p = neighborhoodBestPositions.begin(); p != neighborhoodBestPositions.end(); ++p) {
		cout << "nbp: " << *p << endl;
	}
	 */

	if (options.verbosity >= 4) {
		cout << "RAQUEL before nextPositions" << endl;
	}
	// This vector holds the new positions to be sent to the particle
	vector<double> nextPositions(particleCurrParamSets_.size());
	if (options.verbosity >= 4) {
		cout << "RAQUEL after nextPositions" << endl;
	}
	bool usePersonalBest = false;
	if (options.verbosity >= 3) {
		cout << "RAQUEL: entering for" << endl;
	}
	// For each parameter in the current parameter set
	int i = 0;
	for (auto param = particleBestParamSets_.at(particle).begin(); param != particleBestParamSets_.at(particle).end(); ++param) {
		if (exp) {
			if ( ((float) rand() / (RAND_MAX)) < 0.5 ) {
				usePersonalBest = true;
			}
		}

		double personalBestPos = particleBestParamSets_.at(particle)[i];

		if (usePersonalBest) {
			nextPositions[i] = personalBestPos;
			usePersonalBest = false;
		}
		else {
			// Calculate our mean and std
			if (options.verbosity >= 3) {
				cout << "personalbest: " << personalBestPos << " neighborbest: " << neighborhoodBestPositions[i] << endl;
			}
			double mean = (abs(personalBestPos) + abs(neighborhoodBestPositions[i])) / 2;
			if (options.verbosity >= 3) {
				cout << "mean: " << mean << endl;
			}
			double std = abs(abs(personalBestPos) - abs(neighborhoodBestPositions[i]));
			if (options.verbosity >= 3) {
				cout << "std: " << std << endl;
			}
			// Create the gaussian distribution
			boost::random::normal_distribution<double> dist(mean, std);

			// Pick our next position randomly from distribution
			nextPositions[i] = dist(generalRand);
			if (options.verbosity >= 3) {
				cout << "picked: " << nextPositions[i] << endl;
			}
		}
		++i;
	}
	if (options.verbosity >= 4) {
		cout << "RAQUEL: exiting for" << endl;
	}

	return nextPositions;
}

vector<double> Swarm::calcParticlePosQPSO(unsigned int particle, vector<double> mBests) {

	vector<double> nextPositions;
	vector<double> neighborhoodBests = getNeighborhoodBestPositions(particle);
	for (unsigned int d = 0; d < mBests.size(); ++d) {
		//cout << particle << " before: " << particleCurrParamSets_[particle][d] << endl;
		//cout << particle << " mbest: " << mBests[d] << endl;
		//cout << particle << " best: " << particleBestParamSets_[particle][d] << endl;
		//cout << particle << " swarmbest: " << getNeighborhoodBestPositions(particle)[d] << endl;
		double fi1 = ((double) rand() / (RAND_MAX));
		//cout << "f1: " << fi1 << endl;
		//double fi2 = ((double) rand() / (RAND_MAX));
		//double p = ((fi1 * particleBestFits_.at(particle)) + (fi2 * swarmBestFits_.begin()->first)) / (fi1 + fi2);
		// TODO: Work out whether or not we should be using abs() for this and below
		double p = fi1 * abs(particleBestParamSets_[particle][d]) + (1 - fi1) * neighborhoodBests[d];
		double u = ((double) rand() / (RAND_MAX));

		// TODO: Linearly decrease beta from 1.0 to 0.5? See Liu et al.

		// Liu et all, 2005, eq 2a
		if (u > 0.5) {
			nextPositions.push_back(p - beta_ * abs(mBests[d] - abs(particleCurrParamSets_.at(particle)[d])) * log(1/u));
			//cout << particle << " p: " << p << endl;
			//cout << "mbest - curr: " << mBests[d] - particleCurrParamSets_.at(particle)[d] << endl;
			//cout << "log: " << log(1/u) << endl;
			//cout << particle << " after: " << p - beta_ * abs(mBests[d] - particleCurrParamSets_.at(particle)[d]) * log(1/u) << endl;
		}
		else {
			nextPositions.push_back(p + beta_ * abs(mBests[d] - abs(particleCurrParamSets_.at(particle)[d])) * log(1/u));
			//cout << particle << " after: " << p + beta_ * abs(mBests[d] - particleCurrParamSets_.at(particle)[d]) * log(1/u) << endl;
		}
	}

	return nextPositions;
}

void Swarm::updateInertia() {
	// Eq 3 from Moraes et al
	options.inertia = options.inertiaInit + (options.inertiaFinal - options.inertiaInit) * (inertiaUpdateCounter_ / (float)(options.nmax + inertiaUpdateCounter_));
	if (options.verbosity >= 3) {
		cout << "Setting inertia to " << options.inertia << endl;
	}
}

void Swarm::updateContractionExpansionCoefficient() {
	beta_ = options.betaMax - (((float)flightCounter_ / options.maxNumSimulations) * (options.betaMax - options.betaMin));
	if (options.verbosity >= 3) {
		cout << "updating beta_ to " << beta_ << endl;
	}
}

vector<double> Swarm::getNeighborhoodBestPositions(unsigned int particle) {

	if (options.verbosity >= 3) {
		cout << "Getting neighborhood best for particle " << particle << endl;
	}

	double currBestFit;

	int subPar;

	//Raquel, if we are using constraints, select fit and constraints combined
	if(options.constraints_.size()>=1){

		subPar = fcalcsubParID(particle, 0, options.models.size());
		currBestFit = subParRankFinal[subPar].second;

		//Raquel: also look at other models inside the same particle
		if(options.models.size()>1){

			for(unsigned int mid=1; mid < options.models.size(); mid++){

				subPar = fcalcsubParID(particle, mid, options.models.size());
				if(subParRankFinal[subPar].second < currBestFit){
					currBestFit = subParRankFinal[subPar].second;
				}

			}

		}

		if (options.verbosity >= 3) {
			cout << "############### currBestFit = " << subParRankFinal[subPar].second << " subpar" << subParRankFinal[subPar].first << endl;
		}
	}else{

		// Set the current best fit to our own best fit
		currBestFit = particleBestFits_.at(particle);
		//cout << "cbf: " << currBestFit << endl;
		//cout << "RAQUEL: defining currentBestNeighbor" << endl;


	}


	int currentBestNeighbor = particle;
	//cout << "set initial nb to " << currentBestNeighbor << endl;

	//cout << "allParticles size: " << populationTopology_.size() << endl;
	//cout << "allParticles neighbor size: " << populationTopology_[particle].size() << endl;

	//cout << "RAQUEL: entering for loop populationTopology_" << endl;

	// For every neighbor in this particle's neighborhood
	for (auto neighbor = populationTopology_[particle].begin(); neighbor != populationTopology_[particle].end(); ++neighbor) {


		auto it = particleBestFits_.find(*neighbor);




		// Skip this neighbor if it doesn't contain a fit value
		if (it->second == 0) {
			continue;
		}


		//Raquel adding support to constraints with the if below
		//cout << "checking if neighbor " << *neighbor << " has a better fit of " << it->second << " than " << currBestFit << endl;
		if(options.constraints_.size()>=1){
			//Raquel: support to constraints
			subPar = fcalcsubParID(*neighbor, 0, options.models.size());
			double it2 = subParRankFinal[subPar].second;

			if (it2 < currBestFit) {
				//cout << "it does! setting current best fit of particle " << *neighbor << " of " << it->second << endl;
				currentBestNeighbor = *neighbor;
				currBestFit = it2;
			}

			//Raquel: also look into other models from the same particle
			if(options.models.size()>1){

				for(unsigned int mid=1; mid < options.models.size(); mid++){

					subPar = fcalcsubParID(*neighbor, mid, options.models.size());
					it2 = subParRankFinal[subPar].second;

					if (it2 < currBestFit) {
						//cout << "it does! setting current best fit of particle " << *neighbor << " of " << it->second << endl;
						currentBestNeighbor = *neighbor;
						currBestFit = it2;
					}

				}

			}
		}else{
			// If Neighbor's best fit is better than ours, update the best
			// neighbor. Also, update the current best fit value
			if (it->second < currBestFit) {
				//cout << "it does! setting current best fit of particle " << *neighbor << " of " << it->second << endl;
				currentBestNeighbor = *neighbor;
				currBestFit = it->second;
			}

		}


	}

	//cout << "RAQUEL: exiting for loop populationTopology_" << endl;

	//cout << "cbn: " << currentBestNeighbor << endl;
	return particleBestParamSets_.at(currentBestNeighbor);

}

void Swarm::processParamsPSO(vector<double> &params, unsigned int subParID, double fit) {
	//unsigned int mid = fcalcMID(subParID,options.models.size());
	unsigned int pID = fcalcParID(subParID,options.models.size());

	if (options.verbosity >= 3) {
		cout << "Processing finished params for particle " << pID << " subParticle:" << subParID<<"  with fit of " << fit << endl;
	}

	unsigned int i = 0;

	//razi: TODO: modify, the params order is based on free parameters in the corresonding model, reorder based on the global list

	//outputError("Needs modifications: for models that have different numbers of free parameters, we need to check if the values are in the correct order.");

	// If this is the particle's first iteration, we need to store its params
	if (particleIterationCounter_[pID][0] == 1) {   //razi: later check consistency, models parameters
		//cout << "Storing " << params.size() << " params for particle " << pID << endl;
		for (auto param = params.begin(); param != params.end(); ++param) {
			if (options.verbosity >= 3) {
				cout << *param << endl;
			}
			particleCurrParamSets_[pID][i] = *param;

			//Raquel added: now parameters for subparticles are processed as well
			for(unsigned int mid=0; mid < options.models.size(); mid++){
				subparticleCurrParamSets_[pID][mid][i] = *param;
			}

			++i;
		}
	}
	// The the fit value of this param set is less than the particles best fit
	// we should update the particle's best fit, then store the best fit params
	if (particleBestFits_.at(pID) == 0 || fit < particleBestFits_.at(pID)) {
		if (options.verbosity >= 3) {
			cout << "Updating best fit and params for particle " << pID << endl;
		}
		// Insert into best fit lists, erasing in the case of the map with fits as keys
		particleBestFits_[pID] = fit;
		subparticleBestFits_[subParID] = fit; //Raquel added support to subparticles
		subparticleBestFits_copy[subParID] = fit;

		map<double, unsigned int>::iterator toDelIt = particleBestFitsByFit_.end();

		for (auto it = particleBestFitsByFit_.begin(); it != particleBestFitsByFit_.end(); ++it) {

			//cout << "loop: " << it->second << endl;
			if (it->second == pID) {
				//cout << "Erasing old best fit for particle " << pID << endl;
				toDelIt = it;
			}
		}

		if (toDelIt != particleBestFitsByFit_.end()) {
			particleBestFitsByFit_.erase(toDelIt);
		}

		particleBestFitsByFit_.insert(pair<double, unsigned int>(fit, pID));
		particleBestFitsByFit_copy.insert(pair<double, unsigned int>(fit, pID));

		unsigned int i = 0;
		for (auto param = params.begin(); param != params.end(); ++param) {
			//cout << "Updating best param for particle " << pID << ": " << *param << endl;
			particleBestParamSets_[pID][i] = *param;
			++i;
		}

		for(unsigned int mid = 0; mid < options.models.size(); mid++){
			update_cur_particle_params(pID, mid, true);
		}



	}


	//Raquel added to support subparticles
	if (subparticleBestFits_.at(subParID) == 0 || fit < subparticleBestFits_.at(subParID)) {
			if (options.verbosity >= 3) {
				cout << "Updating best fit and params for particle " << pID << endl;
			}
			// Insert into best fit lists, erasing in the case of the map with fits as keys
			subparticleBestFits_[subParID] = fit; //Raquel added support to subparticles
			subparticleBestFits_copy[subParID] = fit;

			map<double, unsigned int>::iterator toDelIt2 = subparticleBestFitsByFit_.end();

			for (auto it = subparticleBestFitsByFit_.begin(); it != subparticleBestFitsByFit_.end(); ++it) {
				//cout << "loop: " << it->second << endl;
				if (it->second == subParID) {
					//cout << "Erasing old best fit for particle " << pID << endl;
					toDelIt2 = it;
				}
			}

			if (toDelIt2 != subparticleBestFitsByFit_.end()) {
				subparticleBestFitsByFit_.erase(toDelIt2); //Raquel added support to subparticles
			}


			//cout << "Adding " << fit << " to the fit values of subpar " << subParID << endl;
			subparticleBestFitsByFit_.insert(pair<double, unsigned int>(fit, subParID)); //Raquel added support to subparticles
			currentsubparticleBestFitsByFit_.insert(pair<double, unsigned int>(fit, subParID)); //Raquel added support to subparticles
			currentsubparticleBestFitsByFit_copy.insert(pair<double, unsigned int>(fit, subParID));
			for(unsigned int mid = 0; mid < options.models.size(); mid++){
				update_cur_particle_params(pID, mid, true);
			}
	}



}







//Raquel added this function
void Swarm::processParamsDE(vector<double> &params, unsigned int subParID, double fit) {
	//unsigned int mid = fcalcMID(subParID,options.models.size());
	unsigned int pID = fcalcParID(subParID,options.models.size());

	if (options.verbosity >= 3) {
		cout << "Processing finished params for particle " << pID << " subParticle:" << subParID<<"  with fit of " << fit << endl;
	}

	unsigned int i = 0;

	//razi: TODO: modify, the params order is based on free parameters in the corresonding model, reorder based on the global list

	//outputError("Needs modifications: for models that have different numbers of free parameters, we need to check if the values are in the correct order.");

	// If this is the particle's first iteration, we need to store its params
	if (particleIterationCounter_[pID][0] == 1) {   //razi: later check consistency, models parameters
		//cout << "Storing " << params.size() << " params for particle " << pID << endl;
		for (auto param = params.begin(); param != params.end(); ++param) {
			if (options.verbosity >= 3) {
				cout << *param << endl;
			}
			particleCurrParamSets_[pID][i] = *param;

			//Raquel added: now parameters for subparticles are processed as well
			for(unsigned int mid=0; mid < options.models.size(); mid++){
				subparticleCurrParamSets_[pID][mid][i] = *param;
			}

			++i;
		}
	}

	// The the fit value of this param set is less than the particles best fit
	// we should update the particle's best fit, then store the best fit params
	if (particleBestFits_.at(pID) == 0 || fit < particleBestFits_.at(pID)) {
		if (options.verbosity >= 3) {
			cout << "Updating best fit and params for particle " << pID << endl;
		}

		// Insert into best fit lists, erasing in the case of the map with fits as keys
		particleBestFits_[pID] = fit;

		map<double, unsigned int>::iterator toDelIt = particleBestFitsByFit_.end();
		if (options.verbosity >= 4) {
			cout << "looping fits" << endl;
		}
		for (auto it = particleBestFitsByFit_.begin(); it != particleBestFitsByFit_.end(); ++it) {
			//cout << "loop: " << it->second << endl;
			if (it->second == pID) {
				//cout << "Erasing old best fit for particle " << pID << endl;
				toDelIt = it;
				//it = particleBestFitsByFit_.erase(it);
			}
		}
		if (options.verbosity >= 3) {
			cout << "setting new fit value" << endl;
		}
		if (toDelIt != particleBestFitsByFit_.end()) {
			particleBestFitsByFit_.erase(toDelIt);

		}
		if (options.verbosity >= 3) {
			cout << "done" << endl;
		}
		particleBestFitsByFit_.insert(pair<double, unsigned int>(fit, pID));

		//unsigned int i = 0;
		//for (auto param = params.begin(); param != params.end(); ++param) {
		//	cout << "Updating best param for particle " << pID << ": " << *param << endl;
		//	particleBestParamSets_[pID][i] = *param;
		//	++i;
		//}

		if (options.verbosity >= 3) {
			cout << "updating params" << endl;
		}
		for(unsigned int mid = 0; mid < options.models.size(); mid++){
			update_cur_particle_params(pID, mid, true);
		}
		if (options.verbosity >= 3) {
			cout << "done" << endl;
		}

	}



	if (subparticleBestFits_.at(subParID) == 0 || fit < subparticleBestFits_.at(subParID)) {
		if (options.verbosity >= 3) {
			cout << "Updating best fit and params for particle " << pID << endl;
		}

		// Insert into best fit lists, erasing in the case of the map with fits as keys
		subparticleBestFits_[subParID] = fit; //Raquel added support to subparticles
		subparticleBestFits_copy[subParID] = fit; //Raquel added support to subparticles

		map<double, unsigned int>::iterator toDelIt2 = subparticleBestFitsByFit_.end();
		if (options.verbosity >= 3) {
			cout << "looping fits" << endl;
		}
		for (auto it = subparticleBestFitsByFit_.begin(); it != subparticleBestFitsByFit_.end(); ++it) {
			//cout << "loop: " << it->second << endl;
			if (it->second == subParID) {
				//cout << "Erasing old best fit for particle " << pID << endl;
				toDelIt2 = it;
				//it = particleBestFitsByFit_.erase(it);
			}
		}
		if (options.verbosity >= 3) {
			cout << "seting new fit value" << endl;
		}
		if (toDelIt2 != subparticleBestFitsByFit_.end()) {
			subparticleBestFitsByFit_.erase(toDelIt2);

		}

		if (options.verbosity >= 3) {
			cout << "done" << endl;
		}
		subparticleBestFitsByFit_.insert(pair<double, unsigned int>(fit, subParID));//Raquel added support to subparticles
		subparticleBestFitsByFit_copy.insert(pair<double, unsigned int>(fit, subParID));//Raquel added support to subparticles
		currentsubparticleBestFitsByFit_.insert(pair<double, unsigned int>(fit, subParID));//Raquel added support to subparticles
		currentsubparticleBestFitsByFit_copy.insert(pair<double, unsigned int>(fit, subParID));//Raquel added support to subparticles

		//unsigned int i = 0;
		//for (auto param = params.begin(); param != params.end(); ++param) {
		//	cout << "Updating best param for particle " << pID << ": " << *param << endl;
		//	particleBestParamSets_[pID][i] = *param;
		//	++i;
		//}

		if (options.verbosity >= 3) {
			cout << "updating params" << endl;
		}
		for(unsigned int mid = 0; mid < options.models.size(); mid++){
			update_cur_particle_params(pID, mid, true);
		}
		if (options.verbosity >= 3) {
			cout << "done" << endl;
		}

	}
}

void Swarm::launchParticle(unsigned int pID, bool nextGen) {
	//run all subparticles associated with different models
	for(unsigned int mid=0; mid<options.models.size(); ++mid)
		launchSubParticle(pID, mid, nextGen);
}

void Swarm::launchSubParticle(unsigned int pID, unsigned int mid, bool nextGen) {

	//calculate new particle number
	unsigned int subParID = (pID - 1)* options.models.size()+ mid + 1;
	if (options.verbosity>= 3) {
		cout << "RAQUEL launchSUB suParID: " << subParID << endl;
	}
#ifdef TEST_SIMULATOR
//	if (subParID!=3){
//		cout<<"For test all subParticles are avoided except 3. Uncomment later.\n";
//		return;
//	}
#endif
	if (options.verbosity>= 1) {
		cout << "currentGeneration: " << currentGeneration << " options.useCluster: " << options.useCluster << " nextGen: " << nextGen << " bootstrapCounter: " << bootstrapCounter << endl;
	}
	if (currentGeneration == 1 && !options.useCluster && !nextGen && bootstrapCounter == 0) {

		// Construct command needed to run the particle
		//unsigned int mid = getParticleModelId(pID);

		string command = exePath_ + " -v -t particle -p " + toString(subParID) + " -a run -g " + toString(currentGeneration) + " -c " + sConf_[mid] + " -e " + expPaths_[mid][0];  //razi: consider more than 1 exp
		command = command + " -n " + toString(this->getNumModels());
command = command + " >> pOUT 2>&1";
		command = command + " &";

		//just for test, later uncomment
		if (options.verbosity >= 3) {
			cout << "Running SubParticle:"<<subParID<<"[" << pID <<"-"<< mid<<"] using command: "<< command << endl; //razi changed to show the command, was cout << "Running Particle " << pID <<endl;
		}

		if (runCommand(command) != 0) {
			cout << "Warning: Couldn't launch SubParticle:"<<subParID<<"[" << pID <<"-"<< mid<<"]  with command: " << command <<  endl;
			failedParticles_.insert(pID);
			failedSubParticles_.insert(subParID); //razi: add, later check consistency between the two lists
			return;
		}
	}
	if (options.verbosity >= 3) {
		cout << "Adding subparticle " << subParID << " to the list of running particles." << endl;
	}
	runningSubParticles_.insert(subParID); //razi: was //runningParticles_.insert(pID);
	//cout << "RAQUEL sending message to supar" << subParID << "message=NEXT_GEN" << endl;
	swarmComm->sendToSwarm(0, subParID, NEXT_GENERATION, false, swarmComm->univMessageSender);
	if (options.verbosity >= 3) cout << "Launching SubParticle:"<<subParID<<"[" << pID <<"-"<< mid<<"] finished, hence it is active now...\n\n";
}

void Swarm::runGeneration () {   //razi: modified to include subparticles
	// TODO: Implement walltime

	//Raquel added, master will handle the directory creation
	string path = options.jobOutputDir + toString(currentGeneration) + "/";

	if (!checkIfFileExists(path)) {
		runCommand("mkdir " + path);
	}


	unsigned int nModels = options.models.size();
	if(options.verbosity >= 1) {
		cout << "Running generation " << currentGeneration << " with " << options.swarmSize << " particles..." << endl;
	}

	unsigned int numFinishedParticles = 0;

//	vector<unsigned int, vector<bool> > runningParticles;
//vector<unsigned int, vector<bool> > finishedParticles;
//	vector<unsigned int> runningSubParticles;
//	vector<unsigned int> finishedSubParticles;

	unsigned int p, mid, sp;
	vector<unsigned int> newFinishedParticles;
	// razi: handle running particles that include subparticles, dont exceed parallel count
	//we need to track number of completed particles which is not simply equal to (nModels * number of subParticles)
	//a particle is considered completed of all its subparticles are finished. XXXX
	if(options.verbosity>=3){

		cout << "@@@@@@@@@@@@@@@@@ size before cleaning " << currentsubswarmBestFits_.size() << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

	}

	sp=0;
	currentsubparticleBestFitsByFit_.clear();
	currentsubswarmBestFits_.clear();
	finishedParticles_.clear();
	finishedSubParticles_.clear(); //Raquel added to solve problem of less and less result files as generations go
	unsigned int maxSubPar = fcalcsubParID(options.swarmSize, options.models.size()-1, options.models.size());
	int tries = 0;

	if(options.verbosity>=3){
		cout << "Maxsubpar = "  << maxSubPar << endl;

	}

	generationTime_.reset(); //Raquel implementing walltime

	while (finishedSubParticles_.size() < maxSubPar ) { //&& runningSubParticles_.size() < maxSubPar) { //razi: loop over particles, each particle includes nModels subPArticles

		while (runningSubParticles_.size()< options.parallelCount && sp < maxSubPar) {//razi: make sure the number of subparticles don't exceed parallel count limit
			sp++;
			p = fcalcParID(sp, nModels);
			mid = fcalcMID(sp, nModels);

			launchSubParticle(p, mid, false);
			//launchParticle(p, false);
			//cout << "RAQUEL runningSubParticles_.size() " << runningSubParticles_.size() << " options.parallelCount " << options.parallelCount << endl;

			tries = 0;
		}
		// Check for any messages from particles
		//usleep(250000);
		//cout << "RAQUEL Entering checkMasterMessages" << endl;
		newFinishedParticles = checkMasterMessages();
		//cout << "RAQUEL Done checkMasterMessages finishedSubParticles_.size() " << finishedSubParticles_.size() << endl;

		//cout << "finishedSubParticles_.size() " << finishedSubParticles_.size() << " runningSubParticles_.size() " << runningSubParticles_.size() << " maxSubPar " << maxSubPar << " sp " << sp << endl;
		numFinishedParticles = finishedParticles_.size();
		//cout << "RAQUEL rungeneration numFinishedParticles " << numFinishedParticles << endl;
		//tries++; // Enable this fix if the program stays stuck waiting for a message from the slave, even though the master already sent the message

		if(options.verbosity >= 5) {
			cout << numFinishedParticles << " finished particles" << endl;

		}
		if(tries > 50){

			cout << "RAQUEL " << runningSubParticles_.size() << " subparticles failed in generation " << currentGeneration << endl;
//			cout << "RAQUEL resending message NEXT_GENERATION to failed particles" << endl;
			for (auto stuckParticle = runningSubParticles_.begin(); stuckParticle != runningSubParticles_.end(); ++stuckParticle ){

				swarmComm->sendToSwarm(0, *stuckParticle, NEXT_GENERATION, false, swarmComm->univMessageSender);

			}

			tries = 0;
		}

	}

	if(options.verbosity >= 2) {

		cout << "Total time to run generation: " << generationTime_.elapsed() << " seconds." << endl;
		cout << "Max fit time: " << options.maxFitTime << endl;
		int convertedMaxFitTime = convertTimeToSeconds(options.maxFitTime);


		cout << "Max fit time converted: " << convertedMaxFitTime << endl;
	}
	finishedParticles_.clear();
	finishedSubParticles_.clear(); //Raquel added to solve problem of less and less result files as generations go



	/*
	while (numFinishedParticles < options.swarmSize) { //razi: loop over particles, each particle includes nModels subPArticles

		if ((runningSubParticles_.size()*nModels+nModels )<= options.parallelCount && numLaunchedParticles < options.swarmSize) {//razi: make sure subparticles don't exceed parallel count limit
			launchParticle(p);
			numLaunchedSubParticles += 1;
			++p;
		}

		// Check for any messages from particles
		usleep(10000);
		finishedSubParticles = checkMasterMessages();
		numFinishedSubParticles += finishedSubParticles.size();
		finishedSubParticles.clear();
	}
*/

	if ( options.verbosity >=3){ cout<< "running a swarm generation finished ....\n";}

	if (failedParticles_.size() > (unsigned) (options.swarmSize - 3) ) {
		finishFit();
		outputError("Error: You had too many failed runs. Check simulation output (.BNG_OUT files) or adjust walltime.");
	}
	currentGeneration += 1;



}



void Swarm::runAsyncGeneration() {   //Raquel: added new function to run assynchronous generations
	// TODO: Implement walltime
	unsigned int nModels = options.models.size();
	if(options.verbosity >= 1) {
		cout << "Running generation " << currentGeneration << " with " << options.swarmSize << " particles..." << endl;
	}

	//unsigned int numFinishedParticles = 0;

//	vector<unsigned int, vector<bool> > runningParticles;
//vector<unsigned int, vector<bool> > finishedParticles;
//	vector<unsigned int> runningSubParticles;
//	vector<unsigned int> finishedSubParticles;

	unsigned int p, mid, sp;
	vector<unsigned int> newFinishedParticles;
	// razi: handle running particles that include subparticles, dont exceed parallel count
	//we need to track number of completed particles which is not simply equal to (nModels * number of subParticles)
	//a particle is considered completed of all its subparticles are finished. XXXX


	sp=0;
	currentsubparticleBestFitsByFit_.clear();
	currentsubswarmBestFits_.clear();
	finishedParticles_.clear();
	finishedSubParticles_.clear(); //Raquel added to solve problem of less and less result files as generations go
	unsigned int maxSubPar = fcalcsubParID(options.swarmSize, options.models.size()-1, options.models.size());
	int tries = 0;

	if(options.swarmSize<4){

		finishFit();
		outputError("Error: If running an asynchronous run please set a swarm size of at least 4, otherwise run the synchronous method.");


	}

	while ( (finishedParticles_.size() <= (unsigned) (options.swarmSize/2) || finishedParticles_.size() <= 2 ) ) { //razi: loop over particles, each particle includes nModels subPArticles

		while (runningSubParticles_.size()< options.parallelCount && sp < maxSubPar) {//razi: make sure the number of subparticles don't exceed parallel count limit
			sp++;
			p = fcalcParID(sp, nModels);
			mid = fcalcMID(sp, nModels);

			launchSubParticle(p, mid, false);
			//launchParticle(p, false);
			//cout << "RAQUEL runningSubParticles_.size() " << runningSubParticles_.size() << " options.parallelCount " << options.parallelCount << endl;

			tries = 0;
		}
		// Check for any messages from particles
		//usleep(250000);
		//cout << "RAQUEL Entering checkMasterMessages" << endl;
		newFinishedParticles = checkMasterMessages();
		if(options.verbosity>=3){

			cout << "@@@@@@@@@@@@@@ SIZE AFTER CHECKING MESSAGES " <<  currentsubswarmBestFits_.size() << " @@@@@@@@@@@@@@@@@@@@@@@" <<endl;

		}
		flightCounter_ += newFinishedParticles.size();
		//cout << "RAQUEL Done checkMasterMessages newFinishedParticles.size() " << newFinishedParticles.size() << endl;
		//cout << "Flight counter = " << flightCounter_ << endl;
		//cout << "finishedSubParticles_.size() " << finishedSubParticles_.size() << " runningSubParticles_.size() " << runningSubParticles_.size() << " maxSubPar " << maxSubPar << " sp " << sp << endl;
		//numFinishedParticles = finishedParticles_.size();
		//cout << "RAQUEL rungeneration numFinishedParticles " << numFinishedParticles << endl;
		//tries++; // Enable this fix if the program stays stuck waiting for a message from the slave, even though the master already sent the message
		if(tries > 50){

			cout << "RAQUEL " << runningSubParticles_.size() << " subparticles failed in generation " << currentGeneration << endl;
//			cout << "RAQUEL resending message NEXT_GENERATION to failed particles" << endl;
			for (auto stuckParticle = runningSubParticles_.begin(); stuckParticle != runningSubParticles_.end(); ++stuckParticle ){

				swarmComm->sendToSwarm(0, *stuckParticle, NEXT_GENERATION, false, swarmComm->univMessageSender);

			}

			tries = 0;
		}
		if(options.verbosity >= 3){
			cout << "FINISHED SUBPARTICLES: " << finishedSubParticles_.size() << endl;
			cout << "FINISHED PARTICLES: " << finishedParticles_.size() << endl;
		}

	}



	/*
	while (numFinishedParticles < options.swarmSize) { //razi: loop over particles, each particle includes nModels subPArticles

		if ((runningSubParticles_.size()*nModels+nModels )<= options.parallelCount && numLaunchedParticles < options.swarmSize) {//razi: make sure subparticles don't exceed parallel count limit
			launchParticle(p);
			numLaunchedSubParticles += 1;
			++p;
		}

		// Check for any messages from particles
		usleep(10000);
		finishedSubParticles = checkMasterMessages();
		numFinishedSubParticles += finishedSubParticles.size();
		finishedSubParticles.clear();
	}
*/

	if ( options.verbosity >=3){ cout<< "running a swarm generation finished .... Particles finished: " << finishedParticles_.size() << endl;}

	if (failedParticles_.size() > (unsigned) (options.swarmSize - 3) ) {
		finishFit();
		outputError("Error: You had too many failed runs. Check simulation output (.BNG_OUT files) or adjust walltime.");
	}
	currentGeneration += 1;



		//finishedParticles_.clear();
		finishedSubParticles_.clear(); //Raquel added to solve problem of less and less result files as generations go




}



void Swarm::processLateParticles(map<string, double> simParams_, int subParID, bool breed, int currGen) {

		//int found = 0;
			//for(auto i = 1; i<=options.swarmSize; i++){
			//	found = 0;
				//for(auto k = 0; k<options.models.size(); k++){
				//	for(auto j = finishedSubParticles_.begin(); j != finishedSubParticles_.end(); ++j){
						//subParID = fcalcsubParID(i, k, options.models.size());
					//	if(subParID == *j){
							//found = 1;

//						}

	//				}

		if(currGen==0){
			currGen = currentGeneration;

		}
	//
					//if(found==0){
		int pID = fcalcParID(subParID, options.models.size());
		int mid = fcalcMID(subParID, options.models.size());


						map<string, double> paramSet;

							vector<string> paramVals;
							vector<string> paramVecStr;

							//auto best = allGenFits.begin();

							//allGenFits.insert(pair<double, string>(fitCalc, order_params(paramsString, mid)));
							//swarmBestFits_.insert(pair<double, unsigned int>(fitCalc, pID));
							//subswarmBestFits_.insert(pair<double, unsigned int>(fitCalc, subParID));
							//currentsubswarmBestFits_.insert(pair<double, unsigned int>(fitCalc, subParID));

							int bestSP = 0;
							int bestMD = -1;
							int bestP = 0;

							//int mysize = subparticleBestFits_.size();
							//cout << "SIZE1 " << mysize << endl;
							//mysize = swarmBestFits_.size();
							//cout << "SIZE2 " << mysize << endl;

							if(simParams_.size()==0){

								//Raquel: we gotta fix this, make this loop work and discover why the variables are getting erased.
								for(auto pi = subparticleBestFits_.begin(); pi != subparticleBestFits_.end(); pi++){

									bestMD = fcalcMID(pi->second, options.models.size());
									if(bestMD == mid){
										bestSP = pi->second;
										bestP = fcalcParID(bestSP, options.models.size());
										break;
									}

								}


								//split(best->second, paramVals);

								auto fp = options.models.at(mid)->freeParams_.begin();
								for (unsigned int i = 1; i < subparticleCurrParamSets_[pID][mid].size(); i++) {
									//paramSet.insert(pair<string, double> (fp->first, subparticleCurrParamSets_[bestP][mid][i])); //make this work later, to clone the best particle
									paramSet.insert(pair<string, double> (fp->first, subparticleCurrParamSets_[pID][mid][i]));

									++fp;

									paramVecStr.push_back(toString(subparticleCurrParamSets_[bestP][mid][i]));
								}




							}else{

								for (auto i = simParams_.begin(); i != simParams_.end(); ++i) {
									//paramSet.insert(pair<string, double> (fp->first, subparticleCurrParamSets_[bestP][mid][i])); //make this work later, to clone the best particle
									paramSet.insert(pair<string, double> (i->first, i->second));

									paramVecStr.push_back(toString(i->second));
								}
							}

							if(breed){
								swarmComm->sendToSwarm(0, subParID, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramVecStr);
							}
							if(options.verbosity>=3){
								cout << "@@@@@@@@@@@@@@@@@@@@@@@@ LATE SUBPARTICLE " << subParID << " @@@@@@@@@@@@@@@@@@@@@@@@@@@" << " paramvals size: " << paramVals.size() << " freeparam names size: " << options.models.at(mid)->freeParams_.size() << endl;
							}
							for (unsigned int r = 1; r <= options.smoothing; ++r) {
								string outputDir = options.jobOutputDir + toString(currGen) + "/";
								string bnglFilename = getModelName(mid, false)+"_" + toString(pID) + "_" + toString(r) + ".bngl";
								string suffix = toString(pID) + "_" + toString(r) ;
								//options.models[mid]->outputModelWithParams(paramSet, outputDir, bnglFilename, suffix, false, false, true, false, false);

								if (options.models[mid]->getHasGenerateNetwork()){
									// If we're using ODE solver, output .net file
									bnglFilename = boost::regex_replace(bnglFilename, boost::regex("bngl$"), string("net"));
									options.models[mid]->outputModelWithParams(simParams_, outputDir, bnglFilename, suffix,  false, false, true, false, false);
								}
								else {
									// If we're using network free simulation, output .bngl
									options.models[mid]->outputModelWithParams(simParams_, outputDir, bnglFilename, suffix, false, false, false, false, false);
								}


							}
						/*if(breed){
							int numFinishedBreeding = 0;
							while (numFinishedBreeding == 0) {
									if(options.verbosity>=3){

										cout << "checking for DONEBREED FOR LATE PARTICLE " << subParID << endl;
									}
									unsigned int numMessages = swarmComm->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm->univMessageReceiver, true);

							//		int Pheromones::recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage, int messageID) {

									numFinishedBreeding += numMessages;
									//cout << numFinishedBreeding << endl;
									if(options.verbosity>=3){

										cout << "RAQUEL numFinishedBreeding " << numFinishedBreeding << endl;
									}
							}
						}*/
					//}
				//}

			//}

}


void Swarm::cleanupFiles(const char * path) {
	// TODO: Should we fork to do this? ...yes.
	DIR* dp;
	dirent* de;
	errno = 0;

	vector<string> filesToDel;

	dp = opendir(path);

	if (dp)
	{
		while (true)
		{
			//cout << "about to readdir" << endl;
			errno = 0;
			de = readdir( dp );
			if (de == NULL) break;
			//cout << "read " << de->d_name << endl;
			if (boost::regex_match(string( de->d_name ),boost::regex(".+xml$|.+species$|.+cdat$|.+gdat$|.+net$|.+BNG_OUT"))) {
				//cout << "found a file";
				filesToDel.push_back( string(de->d_name) );
			}
		}
		closedir( dp );
	}
	else {
		cout << "Warning: Couldn't open " << path << " to delete un-needed simulation files." << endl;
	}

	//for (auto i: filesToDel) {
	for (auto i = filesToDel.begin(); i != filesToDel.end(); ++i) {
		string fullPath = string(path) + "/" + *i;
		remove(fullPath.c_str());
	}
}


void Swarm::finishFit() {

	if (options.verbosity >=3) cout<<"Fitting finished. Summarizing the results.\n";
	string outputDir;
	if (options.bootstrap) {
		outputDir = options.outputDir + "/" + options.jobName + "_bootstrap/";

		if (!checkIfFileExists(outputDir)) {
			string command = "mkdir " + outputDir + " && cp " + configPath_ + " " + outputDir;

			int errorCounter = 0;
			while (runCommand(command) != 0 && errorCounter < 10) {
				//sleep(1);
				++errorCounter;
			}

			if (errorCounter > 10) {
				cout << "Error: Couldn't create directory: " + outputDir + " to contain final bootstrap results.";
				//outputBootstrapSummary();
			}
		}

		// Output to bootstrap file
		outputBootstrapSummary();

		// Output fit summary
		string outputFilePath = outputDir + toString(bootstrapCounter + 1) + "_all_fits.txt";
		outputRunSummary(outputFilePath);

		// Reset variables
		resetVariables();

		if ((bootstrapCounter + 1) < options.bootstrap) {
			cout << "Deleting last fitting run and beginning next fit.." << endl;

			string command = "cd " + options.jobOutputDir + " && find -mindepth 1 -maxdepth 1 -type d -exec rm -r {} \\;";
			runCommand(command);

			killAllParticles(NEW_BOOTSTRAP);
		}
		else {
			killAllParticles(FIT_FINISHED);
			exit(EXIT_SUCCESS);
		}
	}
	else {
		outputDir = options.jobOutputDir + "Results/";   //was outputDir = options.jobOutputDir + "Results/";
		string command = "mkdir " + outputDir + " && cp " + configPath_ + " " + outputDir;

		int errorCounter = 0;
		while (runCommand(command) != 0 && errorCounter < 10) {
			//sleep(1);
			++errorCounter;
		}

		if (errorCounter > 10) {
			cout << "Error: Couldn't create directory: " + outputDir + " to contain final fitting results. Outputting results to screen.";
			outputRunSummary();
		}

		string outputFilePath = outputDir + "all_fits.txt";
		outputRunSummary(outputFilePath);

		generateBestFitModel(outputDir);

		if(options.verbosity>=3){
			cout << "Killing all particles, sending FIT_FINISHED message" << endl;
		}
		killAllParticles(FIT_FINISHED);
		if(options.verbosity>=3){
			cout << "Done killing all particles." << endl;
		}
		cout << "Finished fitting in " << tmr_.elapsed() << " seconds. Results can be found in " << options.jobOutputDir << "Results" << endl;
	}
	//copyBestFitToResults(outputDir);
}

void Swarm::resetVariables() {
	allGenFits.clear();
	currentGeneration = 1;
	flightCounter_ = 0;
}





void Swarm::getClusterInformation() {

	// If user didn't specify cluster platform, let's figure it out ourself
	if (options.clusterSoftware.size() == 0) {
		// Test for slurm
		string output;
		string output2;
		string output3;

		runCommand("which srun", output);
		if (output.length() > 0) {
			options.clusterSoftware = "slurm";
		}
		// Test for PBS-type
		else{
			runCommand("which qsub", output);
			if(output.length() > 0) {
				// Test for Torque/PBS
				runCommand("which maui", output);
				runCommand("which moab", output2);
				if(output.length() > 0 || output2.length() > 0) {
					options.clusterSoftware = "torque";
				}
				// Test for SGE
				else {
					string output3;
					runCommand("which sge_execd", output);
					runCommand("which qconf", output2);
					runCommand("which qmon", output3);
					if (output.length() > 0 || output2.length() > 0 || output3.length() > 0) {
						outputError("Error: BioNetFit doesn't support GridEngine clusters. If you are not running on a GridEngine cluster, specify the cluster platform in the .conf file using the 'cluster_software' option.");
					}
				}
			}
			// Test for mpi
			else {
				output.clear();
				runCommand("which mpirun", output);

				if (output.length() > 0) {
					options.clusterSoftware = "mpi";

					if (options.hostfile.size() <= 0) {
						cout << "Warning: You are using mpi, but failed to specify a hostfile in your .conf file.." << endl;
					}
				}
			}
		}
	}
	else {
		if (options.clusterSoftware != "slurm" && options.clusterSoftware != "torque" && options.clusterSoftware != "mpi" && options.clusterSoftware != "BNF2mpi") {
			outputError("You specified an unrecognized cluster type in your .conf file. BioNetFit only supports 'torque' or 'slurm' cluster types.");
		}
	}

	// If we still don't know the cluster type, let's ask the user.
	if (options.clusterSoftware.size() == 0) {
		string input;
		string clusterType;

		while (1) {
			cout << "BioNetFit couldn't determine which type of cluster software you are using. Specify (T) for Torque, or (S) for Slurm" << endl;
			getline(cin, input);
			stringstream myInp(input);
			myInp >> clusterType;
			
			if (clusterType == "S" || clusterType == "s") {
				options.clusterSoftware = "slurm";
				break;
			}
			else if (clusterType == "T" || clusterType == "t") {
				options.clusterSoftware = "torque";
				break;
			}
			
			
			if (options.fitType == "de") {
			//	processParamsDE(paramsVec, subParID, fitCalc);
			}
		}
	}
}



vector<unsigned int> Swarm::checkMasterMessages() {  //razi:  modified version, particles replaced with subparticles
	unsigned int pID, mid, NumNewlyFinishedParticles=0;
	std::vector<unsigned int> NewlyFinishedParticles;

	if (options.verbosity >= 3) {
		cout << "Checking messages" << endl;
	}

	// First let's check interswarm communication
	int numMessages = swarmComm->recvMessage(-1, 0, -1, true, swarmComm->univMessageReceiver);
	//if we set false in the boolean above, the program still works for multiple models, but doesn't work for single models, who knows why...
	if (numMessages >= 1) {
		if (options.verbosity >= 3) {
			cout << "Found " << numMessages << " messages" << endl;
		}

		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		smhRange = swarmComm->univMessageReceiver.equal_range(SIMULATION_END);

		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {

			unsigned int subParID;
			subParID = sm->second.sender;

			pID = fcalcParID(subParID, options.models.size());
			mid= fcalcMID(subParID, options.models.size());

			fixRunningParticle(subParID);


			if (options.verbosity >= 3) {
				cout << "subParticle [" << subParID << ":"<< pID <<"-"<< mid <<"] finished simulation" << endl;
			}
			// Get an iterator to the particle in our list of running particles
			runningSubParticlesIterator_ = runningSubParticles_.find(subParID);
			// Then remove it
			if (runningSubParticlesIterator_ == runningSubParticles_.end()) {
				string errMsg = "Error: Couldn't remove subParticle " + toString(subParID) + " from runningParticle list.";

				if(options.verbosity >= 3) {

					cout << "Size of runningSubParticles_: " << runningSubParticles_.size() << endl;
					cout << "Subparticles in running list are: " << endl;

					for(auto myIt = runningSubParticles_.begin(); myIt != runningSubParticles_.end(); myIt++){
						cout << "first: " << *myIt << endl;

					}
					cout << "Done listing running particles. " << endl;
				}
				outputError(errMsg);
			}

			if(options.verbosity >= 3) {

				cout << "Master messages Added supar " << subParID << " to finnished list and removed it from running list" << endl;
			}
			runningSubParticles_.erase(runningSubParticlesIterator_);
			finishedSubParticles_.insert(subParID); //razi: add to the list
			++particleIterationCounter_[pID][mid];
			if(options.verbosity >= 5) {

				cout << "Master messages running sub" << runningSubParticles_.size() << " finished size " << finishedSubParticles_.size() << endl;
			}
			update_cur_particle_params(pID, mid, true); //razi: TODO  the same procedure should be done in other places, the slave sends theliust of free parameters for subparticle that should be reordered for particles

			//finishedsubParticles.push_back(subParID);


			//razi: was ++flightCounter_;    //razi increment the flight counter if this subparticle completes a particle
			NewlyFinishedParticles = update_finished_running_particles(); //update particle lists
			NumNewlyFinishedParticles = NewlyFinishedParticles.size(); //update particle lists
			flightCounter_=flightCounter_+ NumNewlyFinishedParticles;
			if (options.verbosity == 3) {
				cout << "Master NumNewlyFinishedParticles " <<  NumNewlyFinishedParticles << " Flightcounter " << flightCounter_ << endl;
			}
			unsigned int gen = currentGeneration;
			string paramsString;
			if (options.synchronicity == 1 || options.synchronicity == 0) { //Raquel testing
				if (options.fitType == "ga") {
					paramsString = "gen" + toString(gen) + "perm" + toString(subParID) + " ";
				}
				else if (options.fitType == "pso") {
					paramsString = "flock" + toString(gen) + "SubParticle" + toString(subParID) + " ";
				}

			}
			/*else if (options.synchronicity == 0) {
				// Increment our flight counter

				// TODO: This run summary contains 1 less particle than it should.
				// Output a run summary every outputEvery flights
				if (flightCounter_ % options.outputEvery == 0) {
					string outputPath = options.jobOutputDir + toString(flightCounter_) + "_summary.txt";
					outputRunSummary(outputPath);
				}

				paramsString = toString(flightCounter_) + " ";
			}
	*/
			// Store the parameters given to us by the particle


//cout<< "check messages: this part may need modifications ...."<<endl;

//razi: just for debugging, later del
/*int i1, i2;
cout<<"size 1:"<< subparticleCurrParamSets_.size()<<endl;
i1=0;
for (auto ii=subparticleCurrParamSets_.begin(); ii!=subparticleCurrParamSets_.end(); ii++){
	i2=0;
	for (auto jj=ii->second.begin(); jj!=ii->second.end(); jj++){
		cout<<"size [p-sp]: ["<< i1<<","<<i2 <<"]: "<< jj->second.size()<<endl;
		i2++;
	}i1++;
}
*/

			vector<double> paramsVec;
			int i = 0;
			for (vector<string>::iterator m = sm->second.message.begin() + 1; m != sm->second.message.end(); ++m) {
				paramsString += *m + " ";
//cout<<"paramsString["<<i<<"]= "<<paramsString<<endl; //mypause();
				paramsVec.push_back(stod(*m));
				if (options.fitType == "ga") {
					subparticleCurrParamSets_[pID][mid][i] = stod(*m);   //razi: the list and order of free params can be different for  different models
					if (options.verbosity>= 3) {
						cout << "pID " << pID << " mid " << mid << " i " << i << " *m " << *m << endl;
					}
					particleCurrParamSets_[pID][i] = stod(*m);
				}
				++i;
			}
			if(options.verbosity >= 5) {

cout<<i <<" Params collected for one subparticle: "<<paramsString<<endl; //mypause();
subparticleCurrParamSets_copy = subparticleCurrParamSets_;
particleCurrParamSets_copy = particleCurrParamSets_;

			}
			update_cur_particle_params(pID, mid, true); //razi: TODO  the same procedure should be done in other places, the slave sends theliust of free parameters for subparticle that should be reordered for particles

			if (options.verbosity>= 3) {cout << "params are stored" << endl;} //mypause();}

			double fitCalc = stod(sm->second.message[0]);
			if (options.verbosity>= 3) {cout << "fitCalc: " << fitCalc <<endl;} //mypause();}


			update_fitCalcs(); //razi: TODO: later develop and test this function
			if(options.verbosity >= 3) {

cout << "add to allGenFits fitCalc: " << fitCalc <<" params:" << paramsString<<endl; //mypause();
			}
			allGenFits.insert(pair<double, string>(fitCalc, order_params(paramsString, mid)));
			swarmBestFits_.insert(pair<double, unsigned int>(fitCalc, pID));
			subswarmBestFits_.insert(pair<double, unsigned int>(fitCalc, subParID));
			currentsubswarmBestFits_.insert(pair<double, unsigned int>(fitCalc, subParID));
			currentsubswarmBestFits_copy.insert(pair<double, unsigned int>(fitCalc, subParID));

			if(options.verbosity>=3){

				cout << "@@@@@@@@@@@@@@@@@ currentsubswarmBestFits_ size A " << currentsubswarmBestFits_.size() << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

			}
//cout << "add to swarmBestFits_ fitCalc: " << fitCalc <<endl;mypause();

			if (options.fitType == "pso") {
				processParamsPSO(paramsVec, subParID, fitCalc);
			}

			if (options.fitType == "de") {
				processParamsDE(paramsVec, subParID, fitCalc);
			}
		}

		// TODO: When sending NEXT_GENERATION, make sure failed particles have actually run again. If not, they need re-launched.
		if(options.verbosity >= 5) {

			cout << "RAQUEL: before univMessageReceiver.equal_range(SIMULATION_FAIL)" << endl;
			cout << "Number of runnning subparticles before" << runningSubParticles_.size() << endl;
		}
		smhRange = swarmComm->univMessageReceiver.equal_range(SIMULATION_FAIL);
		if(options.verbosity >= 5) {

			cout << "RAQUEL: after univMessageReceiver.equal_range(SIMULATION_FAIL);" << endl;
		}
		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {

			unsigned int subParID = sm->second.sender;
			pID = fcalcParID(subParID, options.models.size());
			mid= fcalcMID(subParID, options.models.size());

			fixRunningParticle(subParID);

			// Get an iterator to the particle in our list of running particles
			runningSubParticlesIterator_ = runningSubParticles_.find(subParID);
			if(options.verbosity >= 5) {

				cout << "Number of runnning subparticles before erasing" << runningSubParticles_.size() << endl;
			}


			if (runningSubParticlesIterator_ == runningSubParticles_.end()) {
				string errMsg = "Error: Couldn't remove particle " + toString(pID) +  "  SubParticle " + toString(subParID) + " from runningParticle list.";
				outputError(errMsg);
			}else{
				runningSubParticles_.erase(runningSubParticlesIterator_);
				if(options.verbosity >= 5) {

					cout << "RAQUEL Erased one subparticle fro running list " << endl;
					cout << "New size: " << runningSubParticles_.size() << endl;
				}
			}
			runningSubParticlesIterator_ = failedSubParticles_.find(subParID);
			if (runningSubParticlesIterator_== failedSubParticles_.end())  //add to the list of failed subparticles
				failedParticles_.insert(subParID);
			else{
				string errMsg = "Error: subparticle failed, but failure was already reported for particle " + toString(pID) +  "  SubParticle " + toString(subParID);
				outputError(errMsg);
			}
			update_finished_running_particles();


			cout << "Particle " << pID << " SubParticle:" << subParID << "  failed in gen " << currentGeneration << endl;

			// Store particle ID in our list of failed particles

		}
		if(options.verbosity >= 5) {

			cout << "RAQUEL: After for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {" << endl;
		}
		int messageCount = 0;
		if(options.verbosity >= 5) {

			cout << "RAQUEL: Before for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) { {" << endl;
		}
		smhRange = swarmComm->univMessageReceiver.equal_range(GET_RUNNING_PARTICLES);
		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {

			vector<string> intParts;
			//for (auto i: runningParticles_) {
			for (auto i = runningSubParticles_.begin(); i != runningSubParticles_.end(); ++i) { //razi: get list of subparticles
				intParts.push_back(toString(*i));
			}
			messageCount++;
			//cout << "RAQUEL: Sending " << messageCount << " messages." << endl;
			//cout << "RAQUEL: sending message to subpar " << sm->second.sender << " message=Running_particles" << endl;
			swarmComm->sendToSwarm(0, sm->second.sender, SEND_RUNNING_PARTICLES, true, intParts);
			//cout << "RAQUEL: Sent " << messageCount << " messages." << endl;
		}

		if(options.verbosity >= 5) {

			cout << "RAQUEL: After for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) { {" << endl;
		}

		swarmComm->univMessageReceiver.clear();

	}

	// Now let's check for any external messages
	//checkExternalMessages();
	if (options.verbosity>= 4) {

		cout<<"checking message is finished.\n";
		cout << "RAQUEL: NewlyFinishedParticles = " << toString(NewlyFinishedParticles.size()) << endl;
	}

	return NewlyFinishedParticles;

}


void Swarm::addRunningParticle(int subParID) {
	//int pID = fcalcParID(subParID, options.models.size());
	runningSubParticles_.insert(subParID);

	vector<unsigned int> newPar = update_finished_running_particles();
	//runningParticles_.insert(pID);

}

void Swarm::fixRunningParticle(int subParID){


	runningSubParticlesIterator_ = runningSubParticles_.find(subParID);
	if (runningSubParticlesIterator_ == runningSubParticles_.end()) {
		addRunningParticle(subParID);
	}

}


unordered_map<unsigned int, vector<double> > Swarm::checkMasterMessagesDE() {

	unsigned int pID, subParID, mid=0; //Removed newlyFinishedParticles
	if (options.verbosity >= 3) { //razi: TODO test for multiple models
		//cout << "Checking messages DE, may need some modifications/test to support sbParticle" << endl;
	}

	unordered_map<unsigned int, vector<double> > subParticleParams;

	// First let's check interswarm communication
	int numMessages = swarmComm->recvMessage(-1, 0, -1, false, swarmComm->univMessageReceiver);

	if (numMessages >= 1) {
		if (options.verbosity >= 3) {
			//cout << "Found " << numMessages << " messages" << endl;
		}

		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		smhRange = swarmComm->univMessageReceiver.equal_range(SIMULATION_END);
		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
			subParID = sm->second.sender;
			pID = fcalcParID(subParID, options.models.size());
			mid= fcalcMID(subParID, options.models.size());

			if (options.verbosity >= 3) {
				//cout << "Particle " << pID << " finished simulation" << endl;
			}

			// Get an iterator to the particle in our list of running particles
			runningSubParticlesIterator_ = runningSubParticles_.find(subParID);
			// Then remove it
			if (runningSubParticlesIterator_ == runningSubParticles_.end()) {
				string errMsg = "Error: Couldn't remove subParticle " + toString(subParID) + " from runningParticle list.";
				cout << "Size of runningSubParticles_(2): " << runningSubParticles_.size() << endl;
				cout << "Subparticles in running list are: " << endl;
				for(auto myIt = runningSubParticles_.begin(); myIt != runningSubParticles_.end(); myIt++){
					cout << "first: " << *myIt << endl;

				}
				cout << "Done listing running particles. " << endl;
				outputError(errMsg);
			}
			runningSubParticles_.erase(runningSubParticlesIterator_);

			++particleIterationCounter_[pID][mid];
			update_finished_running_particles();

			/*
			// Only increment flight counter if we're finishing a trial set
			if (trial == true) {
				++flightCounter_;
			}
			 */

			/*
			unsigned int gen = currentGeneration;
			string paramsString;
			// Only output a run summary after a trial set has finished
			if (options.synchronicity == 1 && trial == true) {
				paramsString = "gen" + toString(static_cast<long long int>(gen)) + "perm" + toString(static_cast<long long int>(pID)) + " ";
			}
			else if (options.synchronicity == 0 && trial == true) {
				// TODO: This run summary contains 1 less particle than it should.
				// Output a run summary every outputEvery flights
				if (flightCounter_ % options.outputEvery == 0) {
					string outputPath = options.jobOutputDir + toString(static_cast<long long int>(flightCounter_)) + "_summary.txt";
					outputRunSummary(outputPath);
				}

				paramsString = toString(static_cast<long long int>(flightCounter_)) + " ";
			}
			 */

			double fitCalc = stod(sm->second.message[0]);

			/*
			bool replaceParams = false;
			if (trial == false) {
				particleBestFits_[pID] = fitCalc;
				insertKeyByValue(particleBestFitsByFit_, fitCalc, pID);
				replaceParams = true;
			}
			else if (trial == true && fitCalc < particleBestFits_.at(pID)) {
				cout << "replacing" << endl;
				particleBestFits_[pID] = fitCalc;
				insertKeyByValue(particleBestFitsByFit_, fitCalc, pID);
				replaceParams = true;
			}
			 */

			// Store the parameters given to us by the particle
			vector<double> paramsVec;
			paramsVec.push_back(fitCalc);
			//int i = 0;
			for (vector<string>::iterator m = sm->second.message.begin() + 1; m != sm->second.message.end(); ++m) {

				paramsVec.push_back(stod(*m));

				/*
				if (replaceParams) {
					//cout << "stored " << stod(*m) << " for particle " << pID << " as " << particleCurrParamSets_[pID][i] << endl;
					particleCurrParamSets_[pID][i] = stod(*m);
					paramsString += *m + " ";
				}
				else if (trial) {
					paramsString += toString(static_cast<long double>(particleCurrParamSets_[pID][i])) + " ";
				}
				 */

				//++i;
			}

			/*
			if (trial == true) {
				allGenFits.insert(pair<double, string>(particleBestFits_[pID], paramsString));
				swarmBestFits_.insert(pair<double, unsigned int>(fitCalc, pID));
			}
			 */

			subParticleParams.insert(pair<unsigned int, vector<double> >(subParID, paramsVec)); //razi: params for subparticle, then translate to particle params
			//cout << "done processing particle " << pID << endl;
		}

		// TODO: When sending NEXT_GENERATION, make sure failed particles have actually run again. If not, they need re-launched.
		smhRange = swarmComm->univMessageReceiver.equal_range(SIMULATION_FAIL);
		for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {

			subParID = sm->second.sender;
			pID = fcalcParID(subParID, options.models.size());
			mid= fcalcMID(subParID, options.models.size());

			// Get an iterator to the particle in our list of running particles
			runningSubParticlesIterator_ = runningSubParticles_.find(subParID);



			// Get an iterator to the particle in our list of running particles
			runningSubParticlesIterator_ = runningSubParticles_.find(subParID);
			if (runningSubParticlesIterator_ == runningSubParticles_.end()) {
				string errMsg = "Error: Couldn't remove particle " + toString(pID) +  "  SubParticle " + toString(subParID) + " from runningParticle list.";
				outputError(errMsg);
			}else
				runningSubParticles_.erase(runningSubParticlesIterator_);

			runningSubParticlesIterator_ = failedSubParticles_.find(subParID);
			if (runningSubParticlesIterator_== failedSubParticles_.end())  //add to the list of failed subparticles
				failedParticles_.insert(subParID);
			else{
				string errMsg = "Error: subParticle failed, but failure was already reported for particle " + toString(pID) +  "  SubParticle " + toString(subParID);
				outputError(errMsg);
			}
			update_finished_running_particles();



			// Increment our finished counter
			//finishedParticles.push_back(pID);
			subParticleParams.insert(pair<unsigned int, vector<double> >(subParID, vector<double>()));  //razi: params for subparticle, then translate to particle params
		}

		swarmComm->univMessageReceiver.clear();

	}

	return subParticleParams;
}

void Swarm::checkExternalMessages() {
	string messagePath = options.jobOutputDir + ".req";
	if (checkIfFileExists(messagePath)) {
		ifstream messageFile(messagePath);

		std::string line;

		if (messageFile.is_open()) {
			while (getline(messageFile, line)) {
				if (line == "output results") {
					string outputDir = options.jobOutputDir + "Results";
					string command = "mkdir " + outputDir + " && cp " + configPath_ + " " + outputDir;

					int errorCounter = 0;
					while (!runCommand(command) && errorCounter < 10) {
						//sleep(1);
						++errorCounter;
					}

					if (errorCounter > 10) {
						cout << "Error: Couldn't create directory: " + outputDir + " to contain final fitting results. Outputting results to screen.";
						outputRunSummary();
					}

					outputRunSummary(outputDir);
				}
			}
		}
		else {
			string errMsg = "Error: Couldn't open data file " + messagePath + " for parsing.";
			outputError(errMsg);
		}
		messageFile.close();
		string delCmd = "rm " + messagePath;
		runCommand(delCmd);
	}
}

string Swarm::getClusterCommand(string runCmd) {
	if (options.clusterSoftware == "slurm") {
		return generateSlurmMultiProgCmd(runCmd);
	}
	else if (options.clusterSoftware == "torque") {
		return generateTorqueBatchScript(runCmd);
	}
	else if (options.clusterSoftware == "mpi") {
		cout << "runCmd is: " << runCmd << endl;
		return generateMPICommand(runCmd);
	}
	else if (options.clusterSoftware == "BNF2mpi") {
			if(options.verbosity >= 3){
				cout << "runCmd is: " << runCmd << endl;
			}
			return generateBNF2MPICommand(runCmd);
	}
	else {
		return 0;
	}
}

string Swarm::generateMPICommand(string cmd) {

	//sbatch << "mpirun -prepend-rank -np 1 " << runCmd << " load " << sConf_ << " : -np " << options.swarmSize << " " << runCmd << " particle 0 run " << sConf_ << endl;
	string newCmd = "mpirun -nooversubscribe ";
	if (options.hostfile.size()) {
		newCmd += "-hostfile " + options.hostfile + " ";
	}
	cout<<"Swarm::generateMPICommand may need some modifications 654\n";
	newCmd += "-tag-output -np 1 " + cmd + " -a load -c " + sConf_[0] + " : -nooversubscribe ";
	if (options.hostfile.size()) {
		newCmd += "-hostfile " + options.hostfile + " ";
	}
	cout<<"Swarm::generateMPICommand may need some modifications 654\n";
	newCmd += "-tag-output -np " + toString(options.swarmSize) + " " + cmd + " -t particle -p 0 -a run -c" + sConf_[0];

	if (options.saveClusterOutput) {
		newCmd += " > " + options.outputDir + "/" + options.jobName + "_cluster_output/" + options.jobName + " 2>&1";
	}
	return newCmd;
}


string Swarm::generateBNF2MPICommand(string cmd) {

	if(options.verbosity >=4){
		cout << "Entering generateBNF2MPICommand()" << endl;
	}
	int maxSubPar = 0;
	maxSubPar = fcalcsubParID(options.swarmSize, options.models.size()-1, options.models.size());


	//sbatch << "mpirun -prepend-rank -np 1 " << runCmd << " load " << sConf_ << " : -np " << options.swarmSize << " " << runCmd << " particle 0 run " << sConf_ << endl;
	string newCmd = "mpirun -oversubscribe ";
	if (options.hostfile.size()) {
		newCmd += "-hostfile " + options.hostfile + " ";
	}
	newCmd += "-tag-output -np 1 " + cmd + " -a load -c " + sConf_[0] + " : -oversubscribe ";
	if (options.hostfile.size()) {
		newCmd += "-hostfile " + options.hostfile + " ";
	}
	newCmd += "-tag-output -np " + toString(maxSubPar) + " " + cmd + " -t particle -p 0 -a run -c " + sConf_[0];

	if (options.saveClusterOutput) {
		newCmd += " > " + options.outputDir + "/" + options.jobName + "_cluster_output/" + options.jobName + " 2>&1";
	}
	if(options.verbosity >=3){
	cout << "MPI command built from generateBNF2MPICommand: " << endl << newCmd << endl;
	}
	return newCmd;
}

string Swarm::generateSlurmMultiProgCmd(string runCmd) {
	string multiProgConfPath = options.jobOutputDir + "multiprog.conf";
	ofstream multiprog(multiProgConfPath, ios::out);

	if (multiprog.is_open()) {

		vector<string> paths;
		paths.push_back(sConf_[0]);
		paths.push_back(expPaths_[0][0]);
		paths.push_back(exePath_);

		string mycommonPath = commonPath(paths, '/');
		string expsubPath;
		string exesubPath;
		string sconfsubPath;

		multiprog << "#SBATCH --workdir=" << mycommonPath << endl;

		string command = exePath_ + " -a load -c " + sConf_[0];
		//multiprog << "0 " << runCmd << " load " << sConf_[0] << endl;
		multiprog << "0 " << command << endl;
		int subPID;
		for (int id = 1; id <= options.swarmSize; ++id) {
			for (unsigned int mid = 0; mid < options.models.size(); ++mid){

				subPID = fcalcsubParID(id, mid, options.models.size());

				exesubPath = exePath_.substr(mycommonPath.length(), mycommonPath.length()+exePath_.length());
				expsubPath = expPaths_[mid][0].substr(mycommonPath.length(), mycommonPath.length()+expPaths_[mid][0].length());
				sconfsubPath = sConf_[mid].substr(mycommonPath.length(), mycommonPath.length()+sConf_[mid].length());

				//multiprog << toString(subPID) + " " << runCmd << " particle " << toString(id) << " run " << sConf_[mid] << endl;
				//				multiprog << toString(id) + " " << runCmd << " particle " << toString(id) << " run " << sConf_[mid] << endl;
				command = "." + exesubPath + " -t particle -p " + toString(subPID) + " -a cluster -g " + toString(currentGeneration) + " -c ." + sconfsubPath + " -e ." + expsubPath;  //razi: consider more than 1 exp
				command = command + " -n " + toString(this->getNumModels());
				//command = command + " >> pOUT 2>&1";
				multiprog << toString(subPID) + " " << command << endl;



			}
		}
		multiprog.close();
	}
	return generateSlurmCommand(multiProgConfPath);
}

string Swarm::generateTorqueBatchScript(string cmd) {
	string batchScriptPath = options.jobOutputDir + "batch_script.sh";
	ofstream batchScript(batchScriptPath, ios::out);

	if (batchScript.is_open()) {
		batchScript << "#!/bin/bash" << endl << endl;

		batchScript << "#PBS -N " << options.jobName << endl;

		if (options.saveClusterOutput) {
			batchScript << "#PBS -o " + options.outputDir + "/" + options.jobName + "_cluster_output/" + options.jobName << endl;
			batchScript << "#PBS -j oe" << endl;
		}
		else {
			batchScript << "#PBS -o /dev/null" << endl;
		}

		if (!options.clusterQueue.empty()) {
			batchScript << "#PBS -p	" + options.clusterQueue << endl;
		}

		// Specify the cluster account if needed
		if (!options.clusterAccount.empty()) {
			batchScript << "#PBS -A " + options.clusterAccount << endl;
		}

		batchScript << "#PBS -l procs=" << options.swarmSize + 1 << ",";

		//if (options.maxFitTime != MAX_LONG) {
		//	batchScript << " walltime=00:00:" << toString(options.maxFitTime) << ",";
		//}
		if (!options.maxFitTime.empty()) { //Raquel fixed torque support
			batchScript << " walltime=" << options.maxFitTime << ","; //Raquel fixed torque support
		}

		batchScript << endl << endl;

		cout<<"Swarm::generateTorqueBatchScript Might need some modifications ....642\n";
		batchScript << "mpirun -np1" << cmd << " load " << sConf_[0] << " : cmd ";
		batchScript.close();
	}

	return cmd;
}

string Swarm::generateSlurmCommand(string cmd, bool multiProg, unsigned int nCPU) {
	string command;

	//TODO: Need to display terminal output from cluster jobs

	// srun submits the job to the cluster
	command += "srun";

	// Add the job name
	command += " -J " + options.jobName;

	// Specify the cluster account if needed
	if (!options.clusterAccount.empty()) {
		command += " -A " + options.clusterAccount;
	}

	if (!options.clusterQueue.empty()) {
		command += " -p	" + options.clusterQueue;
	}

	if (options.emailWhenFinished) {
		command += " --mail-type=END,FAIL --mail-user=" + options.emailAddress;
	}

	//if (options.maxFitTime != MAX_LONG) {
	if (!options.maxFitTime.empty()) { //Raquel fixed slurm support
		//command += " -t 00:00:" + toString(options.maxFitTime);
		command += " -t " + options.maxFitTime; //Raquel fixed slurm support
	}

	// Specify output directory if needed
	if (options.saveClusterOutput) {
		command += " -o " + options.outputDir + "/" + options.jobName + "_cluster_output/" + options.jobName;
	}
	else {
		command += " -o /dev/null";
	}

	if (nCPU == 0) {
		int subParID;
		subParID = fcalcsubParID(options.parallelCount, options.models.size(), options.models.size());
		//command += " -n" + toString(options.parallelCount + 1);
		command += " -n" + toString(subParID);
	}
	else {
		command += " -n" + toString(nCPU);
	}

	command += " -l";


	if (multiProg) {
		command += " --multi-prog " + cmd;
	}
	else {
		command+= " " + cmd;
	}

	return command;
}

string Swarm::generateSlurmBatchFile(string runCmd) {
	string sbatchPath = options.jobOutputDir + "slurm_script.sh";
	ofstream sbatch(sbatchPath, ios::out);

	if (sbatch.is_open()) {
		//TODO: Need to display terminal output from cluster jobs

		sbatch << "#!/bin/sh" << endl << endl;

		// Add the job name
		sbatch << "#SBATCH -J " + options.jobName << endl;

		// Specify the cluster account if needed
		if (!options.clusterAccount.empty()) {
			sbatch << "#SBATCH -A " + options.clusterAccount << endl;
		}

		if (!options.clusterQueue.empty()) {
			sbatch << "#SBATCH -p " + options.clusterQueue << endl;
		}

		// Specify output directory if needed
		if (options.saveClusterOutput) {
			sbatch << "#SBATCH -o " + options.outputDir + "/" + options.jobName + "_cluster_output/" + options.jobName << endl;
		}
		else {
			sbatch << "#SBATCH -o /dev/null" << endl;
		}

		sbatch << "#SBATCH -n" + toString(options.parallelCount + 1) << endl;

		sbatch << endl;

		sbatch << "module load openmpi" << endl;
		//sbatch << "module load intel << endl;"
		sbatch << "echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH" << endl;
		sbatch << "echo PATH: $PATH" << endl;
		sbatch << "echo which mpirun:" << endl;
		sbatch << "which mpirun" << endl;
		sbatch << "echo pwd:" << endl;
		sbatch << "pwd"	<< endl;
		sbatch << "echo ldd GenFit2:" << endl;
		sbatch << "ldd /home/bt285/BioNetFit/bin/GenFit2" << endl;;

		cout<<"Swarm::generateSlurmBatchFile Might need some modifications ....582\n";
		sbatch << "mpirun -mca mca_component_show_load_errors 10 -v --tag-output -np 1 " << runCmd << " load " << sConf_[0] << " : -np " << options.swarmSize << " " << runCmd << " particle 0 run " << sConf_[0] << endl;

		sbatch.close();
	}
	else {
		outputError("Error: Couldn't generate slurm batch file: " + sbatchPath + ". Quitting.");
	}

	runCmd = "sbatch " + sbatchPath;

	return runCmd;
}

unsigned int Swarm::pickWeighted(double weightSum, multimap<double, unsigned int> &weights, unsigned int extraWeight) {
	double lowerBound = 0;
	double upperBound = weightSum;
	if(options.verbosity>=4){
		cout << "lower limit = " << lowerBound << "; upper limit = " << upperBound << endl;
	}
	// TODO: Better error handling here?
	if (upperBound <= 0) {
		return 0;
	}

	if(options.verbosity>=4){

		cout << "lower limit = " << lowerBound << "; upper limit = " << upperBound << endl;
	}
	boost::random::uniform_real_distribution<double> unif(lowerBound, upperBound);

	double random = unif(generalRand);
	double chosen = random * ( 1 - (extraWeight / 10 ));

	//cout << "chosen: " << chosen << endl;

	double currentSum = 0;
	unsigned int indexCounter = 0;
	for (multimap<double, unsigned int>::reverse_iterator w = weights.rbegin(); w != weights.rend(); ++w) {
		currentSum += w->first;
		//cout << "adding " << w->first << endl;

		if (currentSum >= chosen) {
			return indexCounter;
		}
		++indexCounter;
	}
	//cout << "fell off the end" << endl;
	return 0;
}

void Swarm::insertKeyByValue(multimap<double, unsigned int> &theMap, double key, unsigned int value) {
	//cout << "looking for " << value << endl;
	for (auto thePair = theMap.begin(); thePair != theMap.end(); ++thePair) {
		if (thePair->second == value) {
			//cout << "found " << value << " inserting " << key << endl;
			theMap.erase(thePair);
			theMap.insert(pair<double, unsigned int>(key, value));
			return;
		}
	}
}

string Swarm::mutateParamGA(FreeParam* fp, double paramValue) {
	//Timer tmr;
	//uniform_real_distribution<double> unif(0,1);
	if(options.verbosity>=3){
		cout << "generating random dist" << endl;
	}
	boost::random::uniform_real_distribution<double> unif(0,1);
	if(options.verbosity>=3){
		cout << "done" << endl;
	}
	// Generate a random number and see if it's less than our mutation rate.  If is, we mutate.
	if (unif(generalRand) < fp->getMutationRate()) {
		// Store our mutation factor
		if(options.verbosity>=3){
			cout << "setting maxrange" << endl;
		}
		float maxChange = abs(paramValue) * fp->getMutationFactor();
		if(options.verbosity>=3){
			cout << "done" << endl;
		}
		// Generate a new distribution between 0 and double maxChange
		//using param_t = uniform_real_distribution<>::param_type;
		//param_t p{0.0, maxChange * 2};
		//unif.param(p);

		boost::random::uniform_real_distribution<double> unif(0.0, maxChange * 2);

		// Ger our new random number between 0 and maxChange, subtract maxChange.
		double change = unif(generalRand) - maxChange;

		// Add/subtract the value from the parameter
		paramValue+= change;
		//cout << "factor: " << fp->getMutationFactor() << " max change: " << maxChange << " rand: " << rand << endl << endl;
	}
	//double t = tmr.elapsed();
	//cout << "Mutate took " << t << " seconds" << endl;
	return toString(paramValue);
}

void Swarm::saveSwarmState() {
	string serializedSwarmPath = options.jobOutputDir + "swarmState.sconf";

	std::ofstream ofs(serializedSwarmPath);
	if (ofs.is_open()) {

		boost::archive::binary_oarchive ar(ofs);
		ar & *this;
		ofs.close();
	}
	else {
		cout << "Warning: Couldn't save swarm state to " << serializedSwarmPath << endl;
	}
}

void Swarm::initPSOswarm(bool resumeFit) {
	vector<unsigned int> finishedParticles;
	//unsigned int numFinishedParticles = 0;

	if (options.verbosity >= 3) {
		cout << "Launching swarm" << endl;
	}

	// Launch all particles
	//for (unsigned int p = 1; p <= options.swarmSize; ++p) {
	//	launchParticle(p);
	//}

	//if (options.verbosity >= 3) {
	//	cout << "Waiting for particles to finish" << endl;
	//}

	// Wait until they are all finished so we can initialize
	// stuff for subsequent iterations
	/*while(numFinishedParticles < options.swarmSize) {
		// Check for our finished particles

		cout << "RAQUEL swarm size is " << options.swarmSize << " finished: " << numFinishedParticles << endl;
		finishedParticles = checkMasterMessages();
		numFinishedParticles += finishedParticles.size();

		// Sleep for a second
		usleep(1000000);
	}

	*/
	//Raquel replaced the loop above by a synchronous generation
	runGeneration();



	// Fill a vector with all pID's to send for processing
	vector<unsigned int> allParticles;
	for (int p = 1; p <= options.swarmSize; ++p) {
		allParticles.push_back(p);
	}

	// Send all pID's in for processing (update velocities and positions)
	processParticlesPSO(allParticles, false);

	// Set our optimum
	optimum_ = particleBestFitsByFit_.begin()->first;

	// If we're using enhanced stop criteria, update the enhanced stop variables
	if (options.enhancedStop) {
		// Initialize our particle weights and weighted average position
		updateEnhancedStop();
	}


}

void Swarm::outputError(string errorMessage) {
	cout << errorMessage << endl;

	if (commInit) {
		// Tell all particles to die
		killAllParticles(FIT_FINISHED);

		// Deconstruct our communicator
		swarmComm->~Pheromones();
	}

	exit (EXIT_FAILURE);
}





void Swarm::sendMigrationSetDE(unsigned int island, vector<vector<unsigned int> > islandTopology, map<unsigned int, vector<vector<double> > > &migrationSets) {
	if (options.verbosity >= 3) {
		cout << "Sending migration set from island " << island << endl;
	}

	// Generate a random number between 1 and the number of neighbors
	// in the island topology
	boost::random::uniform_int_distribution<int> unif(0, islandTopology[island].size() - 1);
	int parID = 0;
	// Fill a vector with particles from this island, starting with best fits and ending with worst
	vector<unsigned int> particlesToSend;
	//Raquel: if the user provided constraints, then use the ranks calculated from the fit values and constraints combined in order to select the best particle
	if(options.constraints_.size()>0){

		//the ranks are sorted from the best particle to the worst, so the for loop will stop in the first time it finds a particle that belongs to the current iland
		for (auto fitIt = subParRankFinal.begin(); fitIt != subParRankFinal.end(); ++fitIt) {
			parID = fcalcParID(fitIt->first, options.models.size());

			if (particleToIsland_[parID] == island) {
				particlesToSend.push_back(parID);
				if (particlesToSend.size() == options.numToMigrate) {
					break;
				}
			}
		}

	}else{

		for (auto fitIt = particleBestFitsByFit_.begin(); fitIt != particleBestFitsByFit_.end(); ++fitIt) {
			if (particleToIsland_[fitIt->second] == island) {
				particlesToSend.push_back(fitIt->second);
				if (particlesToSend.size() == options.numToMigrate) {
					break;
				}
			}
		}

	}

	// Choose the index of our receiver
	unsigned int receivingIslandIndex = unif(generalRand);

	// Loop through particles to send
	vector<double> migrationSet;
	//for (auto particle = particlesToSend.begin(); particle != particlesToSend.end(); ++particle) {
	for (unsigned int i = 0; i < options.numToMigrate; ++i) {
		if (options.verbosity >= 3) {
			cout << "sending set " << i + 1 << " from particle " << particlesToSend[i] << endl;
		}
		// Fill migrationSet with this particles current parameters
		for (auto param = particleCurrParamSets_[particlesToSend[i]].begin(); param != particleCurrParamSets_[particlesToSend[i]].end(); ++param) {
			migrationSet.push_back(*param);
		}
		// Add that migration set to the list of the receiving island's
		// migration sets
		migrationSets[islandTopology[island][receivingIslandIndex]].push_back(migrationSet);
		migrationSet.clear();
		if (options.verbosity >= 3) {
			cout << island << " sending migration set to island " << islandTopology[island][receivingIslandIndex] << endl;
		}
	}
}

void Swarm::recvMigrationSetDE(unsigned int island, map<unsigned int, vector<vector<double> > > &migrationSets) {
	// If we have any sets to receive..
	if (migrationSets.find(island) != migrationSets.end() && migrationSets.at(island).size()) {

		if (options.verbosity >= 3) {
			cout << "Receiving " << migrationSets.at(island).size() << " migration sets for island " << island << endl;
		}

		// Create a list of particles from the receiving island that will
		// receive the migration set
		vector<unsigned int> particlesToRecv;
		int parID = 0;

		//Raquel: if the user provided constraints, then use the ranks calculated from the fit values and constraints combined in order to select the best particle
		if(options.constraints_.size()>0){
			//Raquel: this will receive the particles with the worst fit to migrate
			for (vector<pair<int,float> >::reverse_iterator fitIt = subParRankFinal.rbegin(); fitIt != subParRankFinal.rend(); ++fitIt) {
				if (options.verbosity >= 3) {
					cout << fitIt->first << endl;
				}
				parID = fcalcParID(fitIt->first, options.models.size());

				if (particleToIsland_[parID] == island) {
					particlesToRecv.push_back(parID);
				}
			}

		}else{

			for (map<double, unsigned int>::reverse_iterator fitIt = particleBestFitsByFit_.rbegin(); fitIt != particleBestFitsByFit_.rend(); ++fitIt) {
				if (options.verbosity >= 3) {
					cout << fitIt->second << endl;
				}
				if (particleToIsland_[fitIt->second] == island) {
					particlesToRecv.push_back(fitIt->second);
				}
			}


		}




		// Iterates through the list of receiving particles
		auto recvIt = particlesToRecv.begin();
		//cout << 999 << endl;
		// For each migration set destined for this island
		unsigned int replacementCounter = 0;
		int subParID = 0;

		for (auto migrationSet = migrationSets[island].begin(); migrationSet != migrationSets[island].end();) {
			if (options.verbosity >= 3) {
				cout << "Replacing " << migrationSet->size() << " params for particle " << *recvIt << " with fit value of " << particleBestFits_.at(*recvIt) << endl;
			}
			vector<string> paramStr;
			unsigned int i = 0;
			// For each parameter in the migration set
			for (auto param = migrationSet->begin(); param != migrationSet->end(); ++param) {
				// Insert parameter to the currentParamSets tracker
				if (options.verbosity >= 3) {
					cout << "param: " << *param << endl;
				}
				particleCurrParamSets_.at(*recvIt)[i] = *param;
				paramStr.push_back(toString(*param));
				++i;
			}
			// Send the parameters to the particle
			//Raquel updating to support subparticles
			for(unsigned int mid = 0; mid < options.models.size(); mid++){

				subParID = fcalcsubParID(*recvIt, mid, options.models.size());
				swarmComm->sendToSwarm(0, subParID, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramStr);

			}
			//swarmComm->sendToSwarm(0, *recvIt, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramStr);
			// Next migration set, choose a new receiver from the island
			++recvIt;
			++replacementCounter;
			migrationSet = migrationSets[island].erase(migrationSet);

			if (replacementCounter >= options.swarmSize / options.numIslands) {
				break;
			}
		}
	}
}

void Swarm::runSGA() {
	if (options.verbosity >= 3) {
		cout << "Running synchronous genetic fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + "1")) {
		string createDirCmd = "mkdir " + options.jobOutputDir + "1";
		if (runCommand(createDirCmd) != 0) {
			outputError("Error: Couldn't create first generation output directory with command: " + createDirCmd + ". Quitting.");
		}
	}

	bool stopCriteria = false;
	while (!stopCriteria){
		//cout << "RAQUEL: starting runGeneration" << endl;
		runGeneration();
		//cout << "RAQUEL: finished runGeneration" << endl;
		//cout << "RAQUEL: starting checkStopCriteria" << endl;
		//cout << "RAQUEL: STOP CRITERIA value" << stopCriteria << endl;
		//cout << "RAQUEL: starting saveSwarmState" << endl;

		saveSwarmState();
		//cout << "RAQUEL: finished saveSwarmState" << endl;

		string currentDirectory = options.jobOutputDir + toString(currentGeneration);
		if (options.deleteOldFiles) {
			//cout << "RAQUEL: starting cleanupFiles" << endl;

			cleanupFiles(currentDirectory.c_str());
			//cout << "RAQUEL: finished cleanupFiles" << endl;

		}

		string outputPath = options.jobOutputDir + toString(currentGeneration - 1 ) + "_summary.txt";
		cout << "Generation: " << toString(currentGeneration - 1 ) << endl;
		//cout << "RAQUEL: starting outputRunSummary from runSGA" << endl;

		outputRunSummary(outputPath);
		if(options.verbosity>=3){
			cout << "RAQUEL: finished outputRunSummary from runSGA" << endl;
			//cout << "RAQUEL: starting breedGenerationGA from runSGA" << endl;
		}
		stopCriteria = checkStopCriteria();

		if (!stopCriteria) {

			breedGenerationGA();
			if(options.verbosity>=3){

				cout << "Finished breedGenerationGA from runSGA" << endl;
			}
		}
	}
	//cout << "RAQUEL: got out from WHILE STOP CRITERIA" << endl;
}

void Swarm::runSPSO() {
	if (options.verbosity >= 3) {
		cout << "Running a synchronous PSO fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + toString(currentGeneration))) {//Raquel: changed, subdir name was hardcoded to 1
		string createDirCmd = "mkdir " + options.jobOutputDir + toString(currentGeneration);//Raquel: changed, subdir name was hardcoded to 1
		runCommand(createDirCmd);
	}

	if (!resumingSavedSwarm) {
		// Generate all particles that will be present in the swarm
		populationTopology_ = generateTopology(options.swarmSize);
	}

	//saveSwarmState();
	string path; //Raquel added
	vector<unsigned int> finishedParticles;
	vector<unsigned int> particlestoProcess;

	bool stopCriteria = false;

	//int numFinishedParticles = 0;

	//unsigned int p = 1;
	//int mid; //Raquel added
	//int sp = 0; //Raquel added
	//int nModels = options.models.size(); //Raquel added


	while (!stopCriteria) {

		if (options.verbosity >= 0) {
			cout << "Generation (flock): " << currentGeneration << endl;
		}


/*
		while (finishedParticles.size() < options.swarmSize && p <= options.swarmSize) {
			if (runningParticles_.size() < options.parallelCount) {

				//Raquel Changed this to run the correct number of subparticles
				sp++; //Raquel added
				p = fcalcParID(sp, nModels);//Raquel added
				mid = fcalcMID(sp, nModels);//Raquel added
				launchSubParticle(p, mid, false);//Raquel added
				//launchParticle(p++);
			}
			cout << "RAQUEL: parallel count " << options.parallelCount << endl;
			cout << "RAQUEL: runningParticles_.size() " << runningParticles_.size() << endl;
			cout << "RAQUEL: options.swarmSize " << options.swarmSize << endl;
			cout << "RAQUEL : p " << p << endl;
			cout << "RAQUEL: finished " << finishedParticles.size() << endl;

			usleep(250000);

			vector<unsigned int> currFinishedParticles;
			cout << "Raquel before checkmastermessages" << endl;
			currFinishedParticles = checkMasterMessages();
			cout << "Raquel after checkmastermessages" << endl;


			if (currFinishedParticles.size()) {
				finishedParticles.insert(finishedParticles.end(), currFinishedParticles.begin(), currFinishedParticles.end());
			}




		}

	*/


		//sp=0;

	/*	finishedParticles_.clear();
		finishedSubParticles_.clear(); //Raquel added to solve problem of less and less result files as generations go
		particlestoProcess.clear();
		numFinishedParticles = 0;

		while (numFinishedParticles < options.swarmSize) { //razi: loop over particles, each particle includes nModels subPArticles

			if (runningSubParticles_.size()< options.parallelCount) {//razi: make sure the number of subparticles don't exceed parallel count limit
				sp++;
				p = fcalcParID(sp, nModels);
				mid = fcalcMID(sp, nModels);

				launchSubParticle(p, mid, false);
				//launchParticle(p, false);

			}
			// Check for any messages from particles
			usleep(10000);
			checkMasterMessages();

			numFinishedParticles = finishedParticles_.size();
		}

*/
			//numFinishedParticles = options.swarmSize;

				//cout << "RAQUEL: starting runGeneration" << endl;
				runGeneration();
				//usleep(250000);

				//cout << "RAQUEL: finished runGeneration" << endl;
				//cout << "RAQUEL: starting checkStopCriteria" << endl;
				stopCriteria = checkStopCriteria();
				//cout << "RAQUEL: STOP CRITERIA value" << stopCriteria << endl;
				//cout << "RAQUEL: starting saveSwarmState" << endl;

				saveSwarmState();
				//cout << "RAQUEL: finished saveSwarmState" << endl;

				string currentDirectory = options.jobOutputDir + toString(currentGeneration);
				if (options.deleteOldFiles) {
					//cout << "RAQUEL: starting cleanupFiles" << endl;

					cleanupFiles(currentDirectory.c_str());
					//cout << "RAQUEL: finished cleanupFiles" << endl;

				}



				// Check for stop criteria
				stopCriteria = checkStopCriteria();
				if (options.verbosity >= 3) {
					cout << "Stop Criteria value: " << stopCriteria << endl;
				}
				string outputPath = options.jobOutputDir + toString(currentGeneration - 1 ) + "_summary.txt";

				outputRunSummary(outputPath);


		if (!stopCriteria) {


			//for(int i=1; i <= options.swarmSize; i++){//
			//	for(int j=0; j<options.models.size(); j++){
			//		subparticleCurrParamSets_[i][j] = particleCurrParamSets_[i];
			//	}
			//}



			for(int i=1; i <= options.swarmSize; i++){//
							particlestoProcess.push_back(i);
			}

			// Process particles
			//cout << "Raquel before processParticlesPSO " << particlestoProcess.size() << " finished particles" << endl;
			processParticlesPSO(particlestoProcess, true);


			//cout << "Raquel After processParticlesPSO in the main loop" << endl;

			//cout << "RAQUEL FALSE STOP CRITERIA IN RUNSPSO" << endl;

		}




	}
}

void Swarm::runSDE() {
	if (options.verbosity >= 3) {
		cout << "Running a synchronous DE fit" << endl;
	}

	// Create an output directory for each island
	if (!checkIfFileExists(options.jobOutputDir + toString(1))) {
		string createDirCmd = "mkdir " + options.jobOutputDir + toString(1);
		runCommand(createDirCmd);
	}

	if (options.verbosity >= 3) {
		cout << "Generating island topology" << endl;
	}
	vector<vector<unsigned int> > islandTopology;
	if (!resumingSavedSwarm) {
		// Generate all particles that will be present in the swarm
		islandTopology = generateTopology(options.numIslands);
		// Map all particles and islands
		unsigned int particlesPerIsland = options.swarmSize / options.numIslands;
		unsigned int currParticle = 1;
		for (unsigned int i = 1; i <= options.numIslands; ++i) {
			for (unsigned int p = 0; p < particlesPerIsland; ++p) {
				particleToIsland_[currParticle] = i;
				islandToParticle_[i][p] = currParticle;
				if (options.verbosity >= 4) {
					cout << i << " " << p << " " << islandToParticle_[i][p] << endl;
				}
				if(options.verbosity >= 3){
					cout << "Particle to island " << particleToIsland_[currParticle] << endl;
					cout << "currPar " << currParticle << endl;
				}
				++currParticle;
			}
		}
	}

	saveSwarmState();

	//if (options.verbosity >= 3) {
	//	cout << "Launching first generation" << endl;
	//}

	unordered_map<unsigned int, vector<double> > finishedParticles;
	bool stopCriteria = false;
	vector<unsigned int> islandFinishedParticles(options.numIslands + 1, 0);

	map<unsigned int, vector<vector<double> > > migrationSets;

	bool trialLoop = false;




	//saveSwarmState();
	string path; //Raquel added
	//vector<unsigned int> finishedParticles;
	//vector<unsigned int> particlestoProcess;


	//int numFinishedParticles = 0;

	//unsigned int p = 1;
	//int mid; //Raquel added
	int sp = 0; //Raquel added
	//int nModels = options.models.size(); //Raquel added




	while(!stopCriteria) {


		//numFinishedParticles = options.swarmSize;
		string currentDirectory = options.jobOutputDir + toString(currentGeneration);
		if (options.deleteOldFiles) {
			//cout << "RAQUEL: starting cleanupFiles" << endl;

			cleanupFiles(currentDirectory.c_str());
			//cout << "RAQUEL: finished cleanupFiles" << endl;

		}


		//cout << "RAQUEL: starting runGeneration" << endl;
		runGeneration();
		//usleep(250000);

		// Create the output directory for the next generation
		if (!checkIfFileExists(options.jobOutputDir + toString(currentGeneration))) {
			string createDirCmd = "mkdir " + options.jobOutputDir + toString(currentGeneration);
			int retryCounter = 0;
			while (runCommand(createDirCmd) != 0 && !checkIfFileExists(options.jobOutputDir + toString(currentGeneration))) {
				if(++retryCounter >= 100) {
					outputError("Error: Couldn't create " + options.jobOutputDir + toString(currentGeneration) + " to hold next generation's output.");
				}
				//sleep(1);
				if(options.verbosity>=3){

					cout << "Trying again to create dir" << endl;
				}
			}
		}

		saveSwarmState();
		//cout << "RAQUEL: finished saveSwarmState" << endl;

		string outputPath = options.jobOutputDir + toString(currentGeneration - 1 ) + "_summary.txt";

		cout << "Generation: " << toString(currentGeneration - 1 ) << endl;

		outputRunSummary(outputPath);

		unordered_map<unsigned int, vector<double> > finished;


		for (int i = 1; i <= options.swarmSize; i++){
			//for (int mid=0; mid < options.models.size(); mid++){
				finished.insert(make_pair(i, subparticleCurrParamSets_[i][0]));


			//}

		}

		if(options.verbosity >= 4){
			cout << "RAQUEL FINISHED SIZE " << finished.size() << endl;
			cout << "vector size " << islandFinishedParticles.size() << endl;
		}
		if (finished.size()) {
			finishedParticles.insert(finished.begin(), finished.end());
			//numFinishedParticles = finished.size();
		}
		int counter = 0;


//		for(auto index = particleToIsland_.begin(); index != particleToIsland_.end(); )

		for (auto particle = finishedParticles.begin(); particle != finishedParticles.end(); ++particle) {
			// Increment the counter that tracks the number of particles finished
			// for a given island
			counter++;
			if(options.verbosity >= 4){
				cout << "counter" << counter << endl;
				cout << "Particle first " << particle->first  << endl;
				cout << "Particle to island " << particleToIsland_[particle->first] << endl;
				cout << "Particle to isalnd size " << particleToIsland_.size() << endl;
			}
			islandFinishedParticles[particleToIsland_[particle->first]] += 1;
		}

		finishedParticles.clear();
		//numFinishedParticles = 0;

		// Done processing params. Now we need to...
		// Check each island to see if it has completed its generation
		for (unsigned int island = 1; island <= options.numIslands; ++island) {
			// If the number finished in this island is equal to the total size of the island
			if (islandFinishedParticles.at(island) == (options.swarmSize / options.numIslands)) {
				if(options.verbosity >= 3){
					cout << "Island " << island << " finished" << endl;
				}
				// Loop through the particles in the island
				for (auto particle = islandToParticle_.at(island).begin(); particle != islandToParticle_.at(island).end(); ++particle) {

					// Will hold mutated and crossed over parameter sets
					vector<double> newParamSet;
					// If we're not in a trial loop, we should mutate and crossover
					if (trialLoop == false) {
						// Create a mutation set for the particle
						newParamSet = mutateParticleDE(*particle);

						// Run crossover for the particle
						newParamSet = crossoverParticleDE(*particle, newParamSet);
					}
					// If we're in a main loop, we just send the current parameter
					// set back to the particle
					else {
						newParamSet = particleCurrParamSets_.at(*particle);
					}

					// Convert our param set to string for sending
					vector<string> paramVecStr;
					for (auto param = newParamSet.begin(); param != newParamSet.end(); ++param) {
						paramVecStr.push_back(toString(*param));
					}

					// Send new param sets to particles for next generation
//					swarmComm->sendToSwarm(0, *particle, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramVecStr);

					//Raquel updated the interprocess communication, added the loop bellow
					for(unsigned int mid = 0; mid < options.models.size(); mid++){
						sp = fcalcsubParID(*particle, mid, options.models.size());
						swarmComm->sendToSwarm(0, sp, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramVecStr);
					}


					// Empty our finished particle counter for the next loop
					islandFinishedParticles[island] = 0;
				}
			}
		}


			// Check for stop criteria
			stopCriteria = checkStopCriteria();
			if(options.verbosity >= 3){
				cout << "Stop Criteria value: " << stopCriteria << endl;
			}
			if (!stopCriteria) {
				// Send/receive migration sets
				if ((options.migrationFrequency != 0) && (currentGeneration % options.migrationFrequency == 0)) {

					for (unsigned int island = 1; island <= options.numIslands; ++island) {
						sendMigrationSetDE(island, islandTopology, migrationSets);
					}

					for (unsigned int island = 1; island <= options.numIslands; ++island) {
						recvMigrationSetDE(island, migrationSets);
					}
				}

//				string outputPath = options.jobOutputDir + toString(currentGeneration) + "_summary.txt";
	//			outputRunSummary(outputPath);
				//++currentGeneration;
				//cout << "Switching to main loop" << endl;
				//trialLoop = false;

			}
		//}
	}
}

void Swarm::runAGA() {
	if (options.verbosity >= 3) {
		cout << "Running an asynchronous genetic fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + "1")) {
		string createDirCmd = "mkdir " + options.jobOutputDir + "1";
		if (runCommand(createDirCmd) != 0) {
			outputError("Error: Couldn't create first generation output directory with command: " + createDirCmd + ". Quitting.");
		}
	}

	//int sp = 0;

	// Holds finished particles
	vector<unsigned int> finishedParticles;

	// Run the first generation
	runGeneration();

	// Save swarm state
	saveSwarmState();


	string outputPath = options.jobOutputDir + toString(currentGeneration - 1 ) + "_summary.txt";

	outputRunSummary(outputPath);

	cout << "Generation: " << toString(currentGeneration - 1 ) << endl;


	// Breed generation
	breedGenerationGA();


	// Re-launch the particles in to the Swarm proper
	/*for (int p = 1; p <= options.swarmSize; ++p) {


		//for(int mid = 0; mid < options.models.size(); mid++){

		//	launchSubParticle(p, mid, false);

		//}

		launchParticle(p, false);
	}

*/

	//int newFinishedPar = 0;

	//bool stopCriteria = false;
	//int iteration = 0;

	/*
	cout << "==================================ENTERING LOOP==================================" << endl;
	while (!stopCriteria){
		usleep(250000);
		cout << "RAQUEL Iteration " << iteration;
		iteration++;
		// Check for any messages from particles and store finished particles

		//finishedParticles = checkMasterMessages();
		finishedParticles = checkMasterMessages();
		newFinishedPar += finishedParticles.size();

		cout << "new finished particles " << newFinishedPar << endl;
		cout << "Running particles " << runningParticles_.size() << endl;
		cout << "Global finished " << finishedParticles_.size() << endl;
		cout << "Global finished sub " << finishedSubParticles_.size() << endl;

		// Check stop criteria
		cout << "Checking stop criteria" << endl;
		stopCriteria = checkStopCriteria();
		cout << "Done, stop criteria value = " << stopCriteria << endl;

		cout << "flightcounter = " << flightCounter_ << endl;
		cout << "max flightcounter" << options.maxNumSimulations << endl;

		// Cleanup old files if needed
		string currentDirectory = options.jobOutputDir + toString(currentGeneration);
		if (options.deleteOldFiles) {
			cleanupFiles(currentDirectory.c_str());
		}

		// If we haven't reached stop criteria, breed and re-launch finished particles
//		if (!stopCriteria && finishedParticles.size() >= options.swarmSize/2) {
		if (!stopCriteria && finishedParticles_.size() >= options.swarmSize/2) {
			cout << "RAQUEL preparing for asynchronous breed, 1/2 of particles are done" << endl;
			// Get new params for any finished particles
			//breedGenerationGA(finishedParticles);
			copy(finishedParticles_.begin(), finishedParticles_.end(), inserter(finishedParticles, finishedParticles.begin()));
			cout << "Async Particles done: " << finishedParticles.size() << endl;

			breedGenerationGA(finishedParticles);
			newFinishedPar = 0;

		//for(auto i = finishedParticles.begin(); i != finishedParticles.end(); ++i){
		//		for(int j=0; j < options.models.size(); j++){
		//			sp = fcalcsubParID(*i, j, options.models.size());
		//			cout << "sending nex gen to " << sp << endl;
		//			swarmComm->sendToSwarm(0, sp, NEXT_GENERATION, false, swarmComm->univMessageSender);
		//		}
		//	}
			//finishedParticles.clear();
			//finishedParticles_.clear();
			//finishedSubParticles_.clear(); //Raquel added to solve problem of less and less result files as generations go

			// Re-launch any finished particles
			for (auto pID = finishedParticles_.begin(); pID != finishedParticles_.end(); ++pID) {

				launchParticle(*pID, false);

			}

			currentGeneration++;
			flightCounter_++;



		}
	}
	*/
	bool stopCriteria = false;
	while (!stopCriteria){
		if (options.verbosity >= 4) {
			cout << "RAQUEL: Started runGeneration();" << endl;
		}
		runAsyncGeneration();
		if (options.verbosity >= 4) {
			cout << "RAQUEL: Finished runGeneration();" << endl;

			cout << "RAQUEL: Started checkStopCriteria();" << endl;
		}

		saveSwarmState();
		if (options.verbosity >= 4) {
			cout << "RAQUEL: Finished saveSwarmState();" << endl;
		}


		string currentDirectory = options.jobOutputDir + toString(currentGeneration);
		if (options.deleteOldFiles) {
			if (options.verbosity >= 4) {
				cout << "RAQUEL: Started cleanupFiles();" << endl;
			}
			cleanupFiles(currentDirectory.c_str());
			if (options.verbosity >= 4) {
				cout << "RAQUEL: Finished cleanupFiles();" << endl;
			}

		}

		string outputPath = options.jobOutputDir + toString(currentGeneration - 1 ) + "_summary.txt";

		outputRunSummary(outputPath);

		cout << "Generation: " << toString(currentGeneration - 1 ) << endl;

		stopCriteria = checkStopCriteria();
		if (options.verbosity >= 4) {
			cout << "RAQUEL: Stop criteria value is " << stopCriteria << endl;

			cout << "RAQUEL: Finished checkStopCriteria();" << endl;
			cout << "RAQUEL: Started saveSwarmState();" << endl;
		}

		if (!stopCriteria) {

			if(options.verbosity>=3){

				cout << "RAQUEL: started breedGenerationGA();" << endl;
			}
			finishedParticles.clear();
			for(auto i = finishedParticles_.begin(); i!=finishedParticles_.end(); ++i){
				//cout << "RAQUEL: processing particle number: " << *i << endl;
				finishedParticles.push_back(*i);
			}
			breedGenerationGA(finishedParticles);
			if(options.verbosity>=3){

				cout << "RAQUEL: finished breedGenerationGA();" << endl;
			}
		}

		currentsubswarmBestFits_.clear();

	}

}

void Swarm::runAPSO() {
	if (options.verbosity >= 3) {
		cout << "Running an asynchronous PSO fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + "1")) {
		string createDirCmd = "mkdir " + options.jobOutputDir + "1";
		runCommand(createDirCmd);
	}

	if (!resumingSavedSwarm) {
		// Generate all particles that will be present in the swarm
		populationTopology_ = generateTopology(options.swarmSize);

		// Run first flight
		initPSOswarm();
	}


	if (options.verbosity >= 0) {
		cout << "Generation (flock): " << toString(currentGeneration - 1 ) << endl;
	}
		vector<unsigned int> particlestoProcess;


	bool stopCriteria = false;
	while (!stopCriteria){
		if (options.verbosity >= 3) {
			cout << "RAQUEL: Started runasyncGeneration();" << endl;
		}
		runAsyncGeneration();
		if (options.verbosity >= 3) {
			cout << "RAQUEL: Finished runasyncGeneration();" << endl;

			cout << "RAQUEL: Started checkStopCriteria();" << endl;
		}

		if (options.verbosity >= 0) {
			cout << "Generation (flock): " << toString(currentGeneration - 1 ) << endl;
		}

		stopCriteria = checkStopCriteria();
		if (options.verbosity >= 3) {
			cout << "RAQUEL: Stop criteria value is " << stopCriteria << endl;

			cout << "RAQUEL: Finished checkStopCriteria();" << endl;
			cout << "RAQUEL: Started saveSwarmState();" << endl;
		}
		saveSwarmState();
		if (options.verbosity >= 3) {

			cout << "RAQUEL: Finished saveSwarmState();" << endl;
		}

		string currentDirectory = options.jobOutputDir + toString(currentGeneration);
		if (options.deleteOldFiles) {
			if (options.verbosity >= 3) {
				cout << "RAQUEL: Started cleanupFiles();" << endl;
			}
			cleanupFiles(currentDirectory.c_str());
			if (options.verbosity >= 3) {
				cout << "RAQUEL: Finished cleanupFiles();" << endl;
			}
		}

		string outputPath = options.jobOutputDir + toString(currentGeneration - 1 ) + "_summary.txt";
		if (options.verbosity >= 3) {
			cout << "RAQUEL: Started outputRunSummary();" << endl;
		}
		outputRunSummary(outputPath);
		if (options.verbosity >= 3) {
			cout << "RAQUEL: Finished outputRunSummary();" << endl;
		}

		if (!stopCriteria) {
			if (options.verbosity >= 3) {
				cout << "RAQUEL: started processParticlesPSO();" << endl;
			}
			particlestoProcess.clear();

			for(auto i = finishedParticles_.begin(); i!=finishedParticles_.end(); ++i){
				if (options.verbosity >= 3) {
					cout << "RAQUEL: processing particle number: " << *i << endl;
				}
				particlestoProcess.push_back(*i);
			}
			processParticlesPSO(particlestoProcess, true);

			if (options.verbosity >= 3) {
				cout << "RAQUEL: finished processParticlesPSO();" << endl;
			}
		}
	}
}

void Swarm::runADE() {

	if (options.verbosity >= 3) {
		cout << "Running asynchronous DE fit" << endl;
	}

	// Create an output directory for each island
	if (!checkIfFileExists(options.jobOutputDir + toString(1))) {
		string createDirCmd = "mkdir " + options.jobOutputDir + toString(1);
		runCommand(createDirCmd);
	}

	if (options.verbosity >= 3) {
		cout << "Generating island topology" << endl;
	}
	vector<vector<unsigned int> > islandTopology;
	if (!resumingSavedSwarm) {
		// Generate all particles that will be present in the swarm
		islandTopology = generateTopology(options.numIslands);
		// Map all particles and islands
		unsigned int particlesPerIsland = options.swarmSize / options.numIslands;
		unsigned int currParticle = 1;
		for (unsigned int i = 1; i <= options.numIslands; ++i) {
			for (unsigned int p = 0; p < particlesPerIsland; ++p) {
				particleToIsland_[currParticle] = i;
				islandToParticle_[i][p] = currParticle;
				if (options.verbosity >= 3) {
					cout << i << " " << p << " " << islandToParticle_[i][p] << endl;

					cout << "Particle to island " << particleToIsland_[currParticle] << endl;
					cout << "currPar " << currParticle << endl;
				}
				++currParticle;
			}
		}
	}

	saveSwarmState();

	//if (options.verbosity >= 3) {
	//	cout << "Launching first generation" << endl;
	//}

	unordered_map<unsigned int, vector<double> > finishedParticles;
	bool stopCriteria = false;
	vector<unsigned int> islandFinishedParticles(options.numIslands + 1, 0);

	map<unsigned int, vector<vector<double> > > migrationSets;

	bool trialLoop = false;




	//saveSwarmState();
	string path; //Raquel added
	//vector<unsigned int> finishedParticles;
	//vector<unsigned int> particlestoProcess;


	//int numFinishedParticles = 0;

	//unsigned int p = 1;
	//int mid; //Raquel added
	int sp = 0; //Raquel added
	//int nModels = options.models.size(); //Raquel added


	bool first = true;

	while(!stopCriteria) {


		//numFinishedParticles = options.swarmSize;

		//cout << "RAQUEL: starting runGeneration" << endl;
		if(first == true){
			//Raquel added: if this is the first geration, run one generation synchronously, otherwaise run it asynchronously
			runGeneration();
			first = false;
		}else{
			runAsyncGeneration();

		}
		// Create the output directory for the next generation
		if (!checkIfFileExists(options.jobOutputDir + toString(currentGeneration))) {
			string createDirCmd = "mkdir " + options.jobOutputDir + toString(currentGeneration);
			int retryCounter = 0;
			while (runCommand(createDirCmd) != 0 && !checkIfFileExists(options.jobOutputDir + toString(currentGeneration))) {
				if(++retryCounter >= 100) {
					outputError("Error: Couldn't create " + options.jobOutputDir + toString(currentGeneration) + " to hold next generation's output.");
				}
				//sleep(1);
				if(options.verbosity>=3){

					cout << "Trying again to create dir" << endl;
				}
			}
		}
		//usleep(250000);

		saveSwarmState();
		//cout << "RAQUEL: finished saveSwarmState" << endl;

		string currentDirectory = options.jobOutputDir + toString(currentGeneration);
		if (options.deleteOldFiles) {
			//cout << "RAQUEL: starting cleanupFiles" << endl;

			cleanupFiles(currentDirectory.c_str());
			//cout << "RAQUEL: finished cleanupFiles" << endl;

		}


		string outputPath = options.jobOutputDir + toString(currentGeneration - 1 ) + "_summary.txt";


		if (options.verbosity >= 0) {
			cout << "Generation: " << toString(currentGeneration - 1 ) << endl;
		}

		outputRunSummary(outputPath);

		unordered_map<unsigned int, vector<double> > finished;


		for (int i = 1; i <= options.swarmSize; i++){
			//for (int mid=0; mid < options.models.size(); mid++){
				finished.insert(make_pair(i, subparticleCurrParamSets_[i][0]));


			//}

		}

		if (options.verbosity >= 4) {
			cout << "RAQUEL FINISHED SIZE " << finished.size() << endl;
			cout << "vector size " << islandFinishedParticles.size() << endl;
		}

		if (finished.size()) {
			finishedParticles.insert(finished.begin(), finished.end());
			//numFinishedParticles = finished.size();
		}
		int counter = 0;


//		for(auto index = particleToIsland_.begin(); index != particleToIsland_.end(); )

		for (auto particle = finishedParticles.begin(); particle != finishedParticles.end(); ++particle) {
			// Increment the counter that tracks the number of particles finished
			// for a given island
			counter++;
			if (options.verbosity >= 4) {
				cout << "counter" << counter << endl;
				cout << "Particle first " << particle->first  << endl;
				cout << "Particle to island " << particleToIsland_[particle->first] << endl;
				cout << "Particle to isalnd size " << particleToIsland_.size() << endl;
			}
			islandFinishedParticles[particleToIsland_[particle->first]] += 1;
		}

		finishedParticles.clear();
		//numFinishedParticles = 0;


		// Done processing params. Now we need to...
		// Check each island to see if it has completed its generation
		for (unsigned int island = 1; island <= options.numIslands; ++island) {
			// If the number finished in this island is equal to the total size of the island
			if (islandFinishedParticles.at(island) == (options.swarmSize / options.numIslands)) {
				if (options.verbosity >= 3) {
					cout << "Island " << island << " finished" << endl;
				}
				// Loop through the particles in the island
				for (auto particle = islandToParticle_.at(island).begin(); particle != islandToParticle_.at(island).end(); ++particle) {

					// Will hold mutated and crossed over parameter sets
					vector<double> newParamSet;
					// If we're not in a trial loop, we should mutate and crossover
					if (trialLoop == false) {
						// Create a mutation set for the particle
						newParamSet = mutateParticleDE(*particle);

						// Run crossover for the particle
						newParamSet = crossoverParticleDE(*particle, newParamSet);
					}
					// If we're in a main loop, we just send the current parameter
					// set back to the particle
					else {
						newParamSet = particleCurrParamSets_.at(*particle);
					}

					// Convert our param set to string for sending
					vector<string> paramVecStr;
					for (auto param = newParamSet.begin(); param != newParamSet.end(); ++param) {
						paramVecStr.push_back(toString(*param));
					}

					// Send new param sets to particles for next generation
//					swarmComm->sendToSwarm(0, *particle, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramVecStr);

					//Raquel updated the interprocess communication, added the loop bellow
					for(unsigned int mid = 0; mid < options.models.size(); mid++){
						sp = fcalcsubParID(*particle, mid, options.models.size());
						swarmComm->sendToSwarm(0, sp, SEND_FINAL_PARAMS_TO_PARTICLE, false, paramVecStr);
					}


					// Empty our finished particle counter for the next loop
					islandFinishedParticles[island] = 0;
				}
			}
		}


			// Check for stop criteria
			stopCriteria = checkStopCriteria();
			if (options.verbosity >= 3) {
				cout << "Stop Criteria value: " << stopCriteria << endl;
			}
			//Raquel printing the output summary even if we didn't reach the stop criteria


			if (!stopCriteria) {
				// Send/receive migration sets
				if ((options.migrationFrequency != 0) && (currentGeneration % options.migrationFrequency == 0)) {

					for (unsigned int island = 1; island <= options.numIslands; ++island) {
						sendMigrationSetDE(island, islandTopology, migrationSets);
					}

					for (unsigned int island = 1; island <= options.numIslands; ++island) {
						recvMigrationSetDE(island, migrationSets);
					}
				}


//				string outputPath = options.jobOutputDir + toString(currentGeneration) + "_summary.txt";
	//			outputRunSummary(outputPath);
				//++currentGeneration;
				//cout << "Switching to main loop" << endl;
				//trialLoop = false;
			}
		//}
	}



}

void Swarm::runASA() {
	if (options.verbosity >= 3) {
		cout << "Running an asynchronous PSADE fit" << endl;
	}

	if (!checkIfFileExists(options.jobOutputDir + "1")) {
		string createDirCmd = "mkdir " + options.jobOutputDir + "1";
		if (runCommand(createDirCmd) != 0) {
			outputError("Error: Couldn't create first generation output directory with command: " + createDirCmd + ". Quitting.");
		}
	}

	// Initialize and/or fill our parameter vectors
	vector<unsigned int> isLocal = vector<unsigned int>(options.swarmSize + 1);
	vector<double> particleTemps = vector<double>(options.swarmSize + 1, 0);
	vector<double> particleRadii = vector<double>(options.swarmSize + 1, 0);
	vector<float> particleFs = generateParticleFs();
	vector<float> particleCRs = generateParticleCRs();
	vector<float> cpuToParticle = vector<float>(options.swarmSize + 1, 0);
	vector<vector<float> > trialParams = vector<vector<float> >(options.swarmSize + 1, vector<float>(2, 0));
	map<unsigned int, unsigned int> particleToController;

	// Launch the initialization population and map
	// CPUs to particles
	for (int p = 1; p <= options.swarmSize; ++p) {
		launchParticle(p);
		cpuToParticle[p] = p;
		particleToController[p] = p;
	}

	unordered_map<unsigned int, vector<double> > finishedParticles;
	unordered_map<unsigned int, vector<double> > initFinishedParticles;
	int numFinishedParticles = 0;

	// Wait for initialization population to finish simulations
	while (numFinishedParticles < options.swarmSize) {
		initFinishedParticles = checkMasterMessagesDE();
		// Save finished particles in another map, we still need to process them in the main loop later on
		if (initFinishedParticles.size()) {
			finishedParticles.insert(initFinishedParticles.begin(), initFinishedParticles.end());
			numFinishedParticles += initFinishedParticles.size();
		}

		// Process each finished particle
		for (auto particle = initFinishedParticles.begin(); particle != initFinishedParticles.end(); ++particle) {
			//++flightCounter_;
			cout << "particle " << cpuToParticle[particle->first] << " on cpu " << particle->first << " finished" << endl;
			string paramString = toString(flightCounter_) + " ";

			// Update parameter values
			unsigned int i = 0;
			for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
				particleCurrParamSets_[cpuToParticle[particle->first]][i++] = *param;
				cout << "updating particle " << cpuToParticle[particle->first] << " to " << *param << " at " << i << endl;
				paramString += toString(*param) + " ";
			}

			// Update fit calcs
			particleBestFits_[cpuToParticle[particle->first]] = particle->second[0];
			insertKeyByValue(particleBestFitsByFit_, particle->second[0], cpuToParticle[particle->first]);
			//allGenFits.insert(pair<double, string>(particle->second[0], paramString));
			cout << "updating particle " << cpuToParticle[particle->first] << " calc to " << particle->second[0] << endl;
			//cout << "count: " << particleBestFitsByFit_.count(particle->second[0]) << endl;

			// Update F and CR
			trialParams[cpuToParticle[particle->first]][0] = particleFs[cpuToParticle[particle->first]];
			trialParams[cpuToParticle[particle->first]][1] = particleCRs[cpuToParticle[particle->first]];
		}
	}

	// Generate our initial temps and radii
	particleTemps = generateParticleTemps();
	particleRadii = generateParticleRadii();

	// Main loop
	cout << "entering main loop" << endl;
	while (!checkStopCriteria()) {

		// Check to make sure we aren't processing init swarm before checking for new
		if (!finishedParticles.size()) {
			finishedParticles = checkMasterMessagesDE();
		}

		// Process any finished particles
		for (auto particle = finishedParticles.begin(); particle != finishedParticles.end(); ++particle) {
			cout << "particle " << cpuToParticle[particle->first] << " on cpu " << particle->first << " finished" << endl;

			string paramString = toString(flightCounter_ + 1) + " ";

			// If the particle finished their local search
			if (isLocal[particle->first]) {
				cout << "particle just finished local search" << endl;

				cout << "greedy acceptance of fit of " << particle->second[0] << endl;
				// Greedy acceptance
				particleBestFits_[cpuToParticle[particle->first]] = particle->second[0];
				insertKeyByValue(particleBestFitsByFit_, particle->second[0], cpuToParticle[particle->first]);

				unsigned int i = 0;
				for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
					particleCurrParamSets_[cpuToParticle[particle->first]][i++] = *param;
					cout << "updating particle " << cpuToParticle[particle->first] << " to " << *param << " at " << i << endl;
					paramString += toString(*param) + " ";
				}

				particleFs[cpuToParticle[particle->first]] = trialParams[cpuToParticle[particle->first]][0];
				particleCRs[cpuToParticle[particle->first]] = trialParams[cpuToParticle[particle->first]][1];

				allGenFits.insert(pair<double, string>(particle->second[0], paramString));

				isLocal[particle->first] = 0;
			}
			else {
				// If we pass metropolis selection, store the params for assignment to receiver
				if (particle->second[0] < particleBestFits_.at(particleToController[cpuToParticle[particle->first]]) || metropolisSelection(particleToController[cpuToParticle[particle->first]], particle->second[0], particleTemps[particleToController[cpuToParticle[particle->first]]])) {

					cout << "particle was accepted through metropolis selection" << endl;

					unsigned int i = 0;
					for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
						particleCurrParamSets_[cpuToParticle[particle->first]][i++] = *param;
						cout << "updating particle " << cpuToParticle[particle->first] << " to " << *param << " at " << i << endl;
						paramString += toString(*param) + " ";
					}

					// Update fitness
					particleBestFits_[cpuToParticle[particle->first]] = particle->second[0];
					insertKeyByValue(particleBestFitsByFit_, particle->second[0], cpuToParticle[particle->first]);
					cout << "accepting params and fit of " << particle->second[0] << endl;

					// Update CR and F
					cout << "updating F to " << trialParams[cpuToParticle[particle->first]][0] << " and CR to " << trialParams[cpuToParticle[particle->first]][1] << endl;
					particleFs[cpuToParticle[particle->first]] = trialParams[cpuToParticle[particle->first]][0];
					particleCRs[cpuToParticle[particle->first]] = trialParams[cpuToParticle[particle->first]][1];
				}
				else {
					for (auto param = particle->second.begin() + 1; param != particle->second.end(); ++param) {
						paramString += toString(*param) + " ";
					}
				}

				// Save params for output no matter what
				allGenFits.insert(pair<double, string>(particle->second[0], paramString));
				cout << "inserting " << paramString << endl;
				// Do a local search at some probability
				// Or if particle is best in swarm
				if (particle->second[0] < particleBestFits_.at(particleToController[cpuToParticle[particle->first]]) || options.localSearchProbability >= ((float)rand() / (float)RAND_MAX)) {
					cout << "going to do a local search" << endl;
					// This receiver will go on to do a local search
					isLocal[particle->first] = cpuToParticle[particle->first];
				}
			}

			if ((options.outputEvery != 0) && (++flightCounter_ && flightCounter_ % options.outputEvery == 0)) {
				string outputPath = options.jobOutputDir + toString(flightCounter_) + "_summary.txt";
				outputRunSummary(outputPath);
				//cout << "fc is " << flightCounter_ << ", outputting" << endl;
			}

			unsigned int it = rand() % options.swarmSize + 1;
			cout << "chose receiver of " << it << endl;
			vector<double> newParams;

			// Swap temps and radii
			swapTR(particleRadii, particleTemps);

			// If we aren't going to do a local search, generate a new trial point
			if (!isLocal[particle->first]) {
				// Pick the controller
				unsigned int controller = pickWeightedSA();
				cout << "not doing local search. chose controller of " << controller << endl;
				// Generate a new trial vector
				particleToController[cpuToParticle[particle->first]] = controller;
				newParams = generateTrialPointSA(controller, it, particleRadii, particleCRs, particleFs, trialParams);
				cout << "generated new trial points" << endl;

				// Make sure we know this CPU will be running this Receiver
				cpuToParticle[particle->first] = it;
				cout << "setting cpu " << particle->first << " as receiver " << it << endl;

				// Save new params for sending to particle
				vector<string> newParamsStr;
				for (auto param = newParams.begin(); param != newParams.end(); ++param) {
					newParamsStr.push_back(toString(*param));
				}

				cout << "running " << cpuToParticle[particle->first] << " with new params on cpu " << particle->first << endl;

				// Send the params to the CPU
				swarmComm->sendToSwarm(0, particle->first, SEND_FINAL_PARAMS_TO_PARTICLE, false, newParamsStr);
				launchParticle(particle->first);
			}
			else {
				cout << "running local search on cpu " << particle->first << endl;
				// Do local search using PREVIOUS receiver
				runNelderMead(isLocal[particle->first], particle->first);
			}
		}
		finishedParticles.clear();
	}
}

void Swarm::runSSA() {
	if (options.verbosity >= 3) {
		cout << "Running an asynchronous PSADE fit" << endl;
	}

	// Initialize and/or fill our parameter vectors
	vector<unsigned int> isLocal = vector<unsigned int>(options.swarmSize, 0);
	vector<float> particleTemps = vector<float>(options.swarmSize + 1, 0);
	vector<float> particleRadii = vector<float>(options.swarmSize + 1, 0);
	vector<float> particleFs = generateParticleFs();
	vector<float> particleCRs = generateParticleCRs();
	vector<float> cpuToParticle = vector<float>(options.swarmSize + 1, 0);
	vector<vector<float> > trialParams = vector<vector<float> >(options.swarmSize + 1, vector<float>(2, 0));

	// Launch the initialization population and map
	// CPUs to particles
	for (int p = 1; p <= options.swarmSize; ++p) {
		launchParticle(p);
		cpuToParticle[p] = p;
	}

	unordered_map<unsigned int, vector<double> > finishedParticles;
	int numFinishedParticles = 0;

	// Wait for initialization population to finish simulations
	while (numFinishedParticles < options.swarmSize) {


	}
}

vector<double> Swarm::generateParticleTemps() {
	// Eq 4 in Olensek et al

	// maxTemp is the difference between highest and lowest fits
	double maxTemp = particleBestFitsByFit_.rbegin()->first - particleBestFitsByFit_.begin()->first;
	//cout << "maxtemp is " << maxTemp << endl;
	//cout << "minTemp is " << options.minTemp << endl;
	//cout << "ss is " << options.swarmSize - 1 << endl;
	double ct = (1 / ((double)options.swarmSize - 1)) * log(maxTemp / options.minTemp);


	//cout << fixed << setprecision(10) << "ct is " << ct << endl;

	vector<double> temps = vector<double>(options.swarmSize + 1, 0);
	for (int p = 1; p <= options.swarmSize; ++p) {
		temps[p] = exp(0 - (ct * (p - 1)));
		//cout << fixed << setprecision(10) << "temp of " << p << " is " << temps[p] << endl;
	}

	return temps;
}

vector<double> Swarm::generateParticleRadii() {
	// Eq 4 in Olensek et al

	float rMax = 1;
	float cr = (1 / ((double)options.swarmSize - 1)) * log(rMax / options.minRadius);

	vector<double> radii = vector<double>(options.swarmSize + 1, 0);
	for (int p = 1; p <= options.swarmSize; ++p) {
		radii[p] = exp(0 - (cr * (p - 1)));
		//cout << "radius of " << p << " is " << radii[p] << endl;
	}

	return radii;
}

vector<float> Swarm::generateParticleFs() {
	float min = 0.5;
	float max = 1.5;

	vector<float> Fs = vector<float>(options.swarmSize + 1, 0);

	for (int p = 1; p <= options.swarmSize; ++p) {
		float r = (float)rand() / (float)RAND_MAX;
		float diff = max - min;

		Fs[p] = (r * diff) + min;
	}

	return Fs;
}

vector<float> Swarm::generateParticleCRs() {
	float min = 0.1;
	float max = 0.9;

	vector<float> CRs = vector<float>(options.swarmSize + 1, 0);

	for (int p = 1; p <= options.swarmSize; ++p) {
		float r = (float)rand() / (float)RAND_MAX;
		float diff = max - min;

		CRs[p] = (r * diff) + min;
	}

	return CRs;
}

unsigned int Swarm::pickWeightedSA() {

	//Eq 6 Olensek et al

	// Calculate the reimann sum
	float sum = 0;
	for (int r = 1; r <= options.swarmSize; ++r) {
		sum += exp(0 - (signed)r);
		//cout << "0 - r: " << 0 - (signed)r << endl;
	}

	//cout << "sum: " << sum << endl;
	// Fill vector with probabilities of particles being selected
	map<double, unsigned int> probabilities;
	unsigned int r = 0;
	double S = 0;
	for (auto particle = particleBestFitsByFit_.begin(); particle != particleBestFitsByFit_.end(); ++particle) {
		S += exp(0 - (signed)++r) / sum;
		probabilities[S] = particle->second;
	}

	float rnd = ((double) rand() / (RAND_MAX));

	int count = 0;
	// Go through and choose the particle
	for (auto probability = probabilities.begin(); probability != probabilities.end(); ++probability) {
		cout << "rnd: " << rnd << " count: " << ++count << " prob: " << probability->first << endl;
		if (rnd < probability->first) {
			return probability->second;
		}

		//rnd -= probability->first;
	}

	// Fell off the end..return the worst particle.
	return probabilities.begin()->second;
}

bool Swarm::metropolisSelection(unsigned int particle, double fit, float particleTemp) {

	float p = min <float>(1, exp( (0 - (fit - particleBestFits_.at(particle))) / particleTemp ) );
	float r = (float)rand() / (float)RAND_MAX;

	if (r < p) {
		return true;
	}
	else {
		return false;
	}
}

void Swarm::swapTR(vector<double> particleRadii, vector<double> particleTemps) {
	// Pick two random individuals
	unsigned int r1 = rand() % options.swarmSize + 1;
	unsigned int r2 = rand() % options.swarmSize + 1;

	// Make sure they're different
	while (r1 == r2) {
		r2 = rand() % options.swarmSize + 1;
	}

	// Choose probability and a random number [0,1]
	float p = min <float>(1, exp( (particleBestFits_.at(r1) - particleBestFits_.at(r2) ) * ((1/particleRadii[r1]) - (1/particleRadii[r2]))));
	float r = (float)rand() / (float)RAND_MAX;

	// Swap
	if (r < p) {
		float r1Temp = particleTemps[r1];
		float r1Radius = particleRadii[r1];

		particleTemps[r1] = particleTemps[r2];
		particleRadii[r1] = particleRadii[r2];

		particleTemps[r2] = r1Temp;
		particleRadii[r2] = r1Radius;
	}
}
vector<double> Swarm::generateTrialPointSA(unsigned int controller, unsigned int receiver, vector<double> particleRadii, vector<float>particleCRs, vector<float>particleFs, vector<vector<float> > &trialParams) {
	vector<double> currParams = normalizeParams(particleCurrParamSets_.at(controller));
	vector<double> testParams = deNormalizeParams(currParams);

	float cr;
	float f;

	// Choose either random F and CR or take from receiver
	if ( ((float)rand() / (float)RAND_MAX) < options.randParamsProbability) {
		cout << "choosing random F and CR" << endl;
		float min = 0.5;
		float max = 1.5;

		float r = (float)rand() / (float)RAND_MAX;
		float diff = max - min;
		f = (r * diff) + min;

		min = 0.1;
		max = 0.9;
		r = (float)rand() / (float)RAND_MAX;
		diff = max - min;
		cr = (r * diff) + min;

		//cout << "f is " << f << " and cr is " << cr << endl;
	}
	else {
		f = particleFs[receiver];
		cr = particleCRs[receiver];
		//cout << "using receiver f and cr of " << f << " and " << cr << endl;
	}

	trialParams[receiver][0] = f;
	trialParams[receiver][1] = cr;

	cout << "before mutate: " << currParams[0] << endl;
	// Mutation functions pulls parameters from the global list
	currParams = mutateParticleSA(controller, f);
	cout << "before crossover: " << currParams[0] << endl;
	// Crossover function uses params returned from mutation function
	currParams = crossoverParticleDE(controller, currParams, cr, true);

	//cout << "mutated and crossed over params" << endl;

	vector<double> newParams;
	for (auto param = currParams.begin(); param != currParams.end(); ++param) {
		float r = (float)rand() / (float)RAND_MAX;
		cout << "before jiggle: " << *param << endl;
		float newParam = *param + particleRadii[controller] * tan(3.141592654 * (r - 0.5));
		cout << "after jiggle: " << newParam << endl;

		if (newParam <= 0 || newParam > 1) {
			newParam = (double)rand() / (double)RAND_MAX;
			//cout << "new: " << newParam << endl;
		}

		newParams.push_back(newParam);
	}

	return deNormalizeParams(newParams);
}

void Swarm::generateBootstrapMaps(vector<map<string, map<string, map<double,unsigned int> > > > &bootStrapMaps) {
	for (unsigned int i = 0; i < options.bootstrap; ++i) {
		//cout << "i: " << i << endl;
		// For each .exp file
		map<string, map<string, map<double, unsigned int> > > maps;
		for (auto dataSet = options.expFiles.begin(); dataSet != options.expFiles.end(); ++dataSet) {
			//cout << "set: " << dataSet->first << endl;
			// For each column
			map<string, map<double, unsigned int> > bsMap;
			for (auto col = dataSet->second->dataCurrent->begin(); col != dataSet->second->dataCurrent->end(); ++col) {
				//cout << "col: " << col->first << endl;
				// Fill the map with 0's
				map<double, unsigned int> colVals;
				for (auto tp = col->second.begin(); tp != col->second.end(); ++tp) {
					//cout << "tp: " << tp->first << endl;
					colVals[tp->first] = 0;
				}

				// Select datapoints at random. If a datapoint is selected,
				// increment it's integer value in the vals map
				for (unsigned int o = 0; o < col->second.size(); ++o) {
					int i;
					if (col->second.size() > 1) {
					  i = rand() % (col->second.size() - 1);
					}
					else
					{
					  i = 0;
					}
					auto it = col->second.begin();
					advance(it, i);
					colVals[it->first] += 1;
					//cout << "chose " << i << " (" << it->first << ")" << endl;
				}

				// Insert this column into the map
				bsMap.insert(pair<string, map<double, unsigned int> >(col->first, colVals));
			}
			// Insert this dataset into the maps set
			maps.insert(pair<string, map<string, map<double, unsigned int> > >(dataSet->first, bsMap));
		}
		// Insert this set of maps into the master set
		bootStrapMaps.push_back(maps);
	}

	/*
	unsigned int setCounter = 1;
	for (auto bsSet = bootStrapMaps.begin(); bsSet != bootStrapMaps.end(); ++bsSet) {
		cout << "Set " << setCounter << endl;
		for (auto dataSet = bsSet->begin(); dataSet != bsSet->end(); ++dataSet) {
			cout << "Dataset " << dataSet->first << endl;
			for (auto col = dataSet->second.begin(); col != dataSet->second.end(); ++col) {
				cout << "Col " << col->first << endl;
				for (auto tp = col->second.begin(); tp != col->second.end(); ++tp) {
					cout << tp->first << " " << tp -> second << endl;
				}
			}
		}
	}
	 */
}






map<double, unsigned int> Swarm::getNearestNeighbors(unsigned int it, unsigned int N) {
	map<double, unsigned int> neighbors;

	// Normalize it params
	vector<double> itParams = normalizeParams(particleCurrParamSets_.at(it));

	for (auto paramSet = particleCurrParamSets_.begin(); paramSet != particleCurrParamSets_.end(); ++paramSet) {
		// Skip this particle if it is the 'it', or if it has the same fit as 'it'
		if (paramSet->first == it || particleCurrParamSets_.at(it) == particleCurrParamSets_.at(paramSet->first)) {
			continue;
		}

		// Calculate sum of squares
		double sum = 0;
		vector<double> pParams = normalizeParams(particleCurrParamSets_.at(paramSet->first));
		auto itIt = itParams.begin();
		for (auto param = pParams.begin(); param != pParams.end(); ++param) {
			sum += pow((*itIt - *param), 2);
			++itIt;
		}

		// Insert sum and id to map. It is automatically ordered.
		neighbors.insert(pair<double, unsigned int>(sum, paramSet->first));
	}

	// We only need the top N, so delete the rest
	auto nIt = neighbors.begin();
	advance(nIt, N);
	neighbors.erase(nIt, neighbors.end());

	return neighbors;
}






int Swarm::printDetails(){ //razi added for debugging

 try{
 cout <<endl<< "================================================================================================"<<endl;
 cout << "=                              DETAILS OF SWARM OBJECT                                         ="<<endl;
 cout << "================================================================================================"<<endl;
 cout<<"general details for swarm object:..."<<endl;
 cout<<"       isMaster:"<<isMaster<<endl;
 cout<<"       currentGeneration:"<<currentGeneration<<endl;
 cout<<"       bootstrapCounter:"<<bootstrapCounter<<endl;
 cout<<"       resumingSavedSwarm:"<<resumingSavedSwarm<<endl;
 cout<<"       hasMutate:"<<hasMutate<<endl;
 cout<<"       fitCompareTolerance:"<<fitCompareTolerance<<endl;
 cout<<"       commInit:"<<commInit<<endl<<endl;


 cout<<"Swarm options ...."<<endl;
 cout<<"       jobName:"<<options.jobName <<endl;
 cout<<"       fitType:"<<options.fitType <<endl;
 cout<<"       outputDir:"<<options.outputDir <<endl;
 cout<<"       jobOutputDir:"<<options.jobOutputDir <<endl;
 cout<<"       bngCommand:"<<options.bngCommand <<endl<<endl;

 cout<<"Swarm options ...."<<endl;
 cout<<"       verbosity:"<<options.verbosity <<endl;
 cout<<"       synchronicity:"<<options.synchronicity <<endl;
 cout<<"       maxGenerations:"<<options.maxGenerations <<endl;
 cout<<"       swarmSize:"<<options.swarmSize <<endl;
 cout<<"       minFit:"<<options.minFit <<endl;
 cout<<"       parallelCount:"<<options.parallelCount <<endl;
 cout<<"       objFunc:"<<options.objFunc <<endl;
 cout<<"       usePipes:"<<options.usePipes <<endl;
 cout<<"       useCluster:"<<options.useCluster <<endl;
 cout<<"       seed:"<<options.seed <<endl;
 cout<<"       bootstrap:"<<options.bootstrap <<endl;
 cout<<"       divideByInit:"<<options.divideByInit <<endl;
 cout<<"       logTransformSimData:"<<options.logTransformSimData <<endl;
 cout<<"       standardizeSimData:"<<options.standardizeSimData <<endl;
 cout<<"       standardizeExpData:"<<options.standardizeExpData <<endl;
 cout<<"       deleteOldFiles:"<<options.deleteOldFiles <<endl<<endl;

 cout<<"Genetic Algorithm parameters..."<<endl;
 cout<<"       extraWeight:"<<options.extraWeight <<endl;
 cout<<"       swapRate:"<<options.swapRate <<endl;
 cout<<"       forceDifferentParents:"<<options.forceDifferentParents <<endl;
 cout<<"       maxRetryDifferentParents:"<<options.maxRetryDifferentParents <<endl;
 cout<<"       smoothing:"<<options.smoothing <<endl;
 cout<<"       keepParents:"<<options.keepParents <<endl;
 cout<<"       maxFitTime:"<<options.maxFitTime <<endl;
 cout<<"       maxNumSimulations:"<<options.maxNumSimulations <<endl<<endl;


 cout<<"PSO parameters..."<<endl;
 cout<<"       inertia:"<<options.inertia <<endl;
 cout<<"       cognitive:"<<options.cognitive <<endl;
 cout<<"       social:"<<options.social <<endl;
 cout<<"       nmax:"<<options.nmax <<endl;
 cout<<"       nmin:"<<options.nmin <<endl;
 cout<<"       inertiaInit:"<<options.inertiaInit <<endl;
 cout<<"       inertiaFinal:"<<options.inertiaFinal <<endl;
 cout<<"       absTolerance:"<<options.absTolerance <<endl;
 cout<<"       relTolerance:"<<options.relTolerance <<endl;
 cout<<"       mutateQPSO:"<<options.mutateQPSO <<endl;
 cout<<"       betaMax:"<<options.betaMax <<endl;
 cout<<"       betaMin:"<<options.betaMin <<endl;
 cout<<"       topology:"<<options.topology <<endl;
 cout<<"       psoType:"<<options.psoType <<endl;
 cout<<"       enhancedStop:"<<options.enhancedStop <<endl;
 cout<<"       enhancedInertia:"<<options.enhancedInertia <<endl;


 cout<<"DE parameters..."<<endl;
 cout<<"       numIslands:"<<options.numIslands <<endl;
 cout<<"       mutateType:"<<options.mutateType <<endl;
 cout<<"       cr:"<<options.cr <<endl;
 cout<<"       migrationFrequency:"<<options.migrationFrequency <<endl;
 cout<<"       numToMigrate:"<<options.numToMigrate <<endl;


 cout<<"DE parameters..."<<endl;
 cout<<"       minTemp:"<<options.minTemp <<endl;
 cout<<"       minRadius:"<<options.minRadius <<endl;
 cout<<"       localSearchProbability:"<<options.localSearchProbability <<endl;
 cout<<"       randParamsProbability:"<<options.randParamsProbability <<endl;
 cout<<"       outputEvery:"<<options.outputEvery <<endl;

 cout<<"Cluster options..."<<endl;
 cout<<"       clusterSoftware:"<<options.clusterSoftware <<endl;
 cout<<"       clusterAccount:"<<options.clusterAccount <<endl;
 cout<<"       saveClusterOutput:"<<options.saveClusterOutput <<endl;
 cout<<"       clusterQueue:"<<options.clusterQueue <<endl;
 cout<<"       emailWhenFinished:"<<options.emailWhenFinished <<endl;
 cout<<"       emailAddress:"<<options.emailAddress <<endl;
 cout<<"       hostfile:"<<options.hostfile <<endl;

 cout<<"List of Free Parameters ..."<<endl;
	 for (auto frp = options.freeParams_.begin(); frp != options.freeParams_.end(); ++frp)
		cout<<" Free Parameter: " << frp->first <<"  GenMethod:" << frp->second->getGenerationMethod() <<" Name:" << frp->second->getParameterName()<< " Range:["<<frp->second->getGenMin()<<"..."<< frp->second->getGenMax() <<"]"<<  std::endl;
 cout<<endl;


 cout<<"Modles and Actions..."<<endl;
 for (unsigned int i=0; i< options.models.size(); i++ ){
	 cout<<" Model ["<<i<<"]: " << options.models.at(i)->getName()<<endl<<"              ";
	 for (auto act = options.models.at(i)->actions.begin(); act != options.models.at(i)->actions.end(); ++act)
			cout<<"Action: " << act->first  <<", ";
 	 cout<<endl;

 	 for (auto frp = options.models.at(i)->freeParams_.begin(); frp != options.models.at(i)->freeParams_.end(); ++frp)
			cout<<" Free Parameter: " << frp->first <<"  GenMethod:" << frp->second->getGenerationMethod() <<" Name:" << frp->second->getParameterName()<< " Range:["<<frp->second->getGenMin()<<"..."<< frp->second->getGenMax() <<"]"<<  std::endl;
	 cout<<endl;
 }


 cout<<"EXP files ..."<<endl;
	 for (unsigned int i=0; i< options.models.size(); ++i){
		 cout<<" EXP Files for model:"<< i<<"   ";
		 for (unsigned int j=0; j< expPaths_[i].size(); ++j)
	 		 cout<<expPaths_[i][j]<<",  ";
		 cout<<endl;
	 }


 cout << "================================================================================================"<<endl<<endl;

 }
 catch (exception& e)
 {
    cout << "Error occurred when printing details of a swarm object: Standard exception: " << e.what() << endl;
 }
 return 1;
}


void Swarm::generateBestFitModel(string outputDir) {
	map<string, double> paramSet;

	vector<string> paramVals;
	auto best = allGenFits.begin();
	split(best->second, paramVals);

	auto fp = this->getFreeParams_().begin();
	for (unsigned int i = 1; i < paramVals.size(); i++) {
		paramSet.insert(pair<string, double> (fp->first, stod(paramVals[i])));
		++fp;
	}
	if (options.verbosity>=3) cout<<"Swarm: Generating Best Fit Model in path:"<< outputDir<< " File:"<<options.jobName<<".bngl .\n";

	//razi TODO:   adjust paramSet based on the model
	for (unsigned int mid =0; mid < options.models.size(); mid++){
		options.models[mid]->outputModelWithParams(paramSet, outputDir, (options.jobName + ".bngl"), "", false, false, false, false, false);
	}

	if(options.verbosity>=3){
		cout << "Finished saving " << options.jobName << ".bngl" << endl;
	}
}

void Swarm::outputRunSummary(string outputPath) {
	int counter; //Raquel: this will count the index of the parameter value for a given parameter
	int particleCounter; //Raquel: This will count the total number of particles
	vector<string> paramList; //Raquel: This is the global list of parameters
	int paramCounter; //Raquel: this is the index of the parameter name

	//int parID = 0;

	string paramString;

	map<unsigned int, map<unsigned int, map<string, double> > > alignedResults;   //Raquel: Particle -> Model -> Parameter -> value

	//cout<<"Swarm::outputRunSummary(string outputPath) may need some modifications ...."<<endl; //mypause(); //razi TODO: later test

	if (options.verbosity >=3) cout<<"Fitting finished: Writing Summary into File: " << outputPath<<endl;

	ofstream outputFile;
	if (options.verbosity >=5){
		cout << "RAQUEL opening file" << endl;
	}
	//outputFile.open(outputPath, ofstream::out);
	//outputFile.open(outputPath);
	//outputFile.open(outputPath, ios::out | ios::app);
	//outputFile.open(outputPath.c_str());
	outputFile.open(outputPath.c_str(), ios::out | ios::app);
	if (options.verbosity >=5){
		cout << "entering if" << endl;
	}
	//int paramIndex = 0;




	if (outputFile.is_open()){

		outputFile.precision(8);
		outputFile.setf(ios::scientific);


	if(options.fitType == "ga" || options.fitType == "pso" || options.fitType == "de") {


		paramCounter = 0;
		 for (unsigned int i=0; i< options.models.size(); i++ ){//loops through models
	 			if (options.verbosity >=1){

	 				cout<<" Model ["<<i<<"]: " << options.models.at(i)->getName()<<endl<<"              ";
	 			}
			 counter = 0;
		 	 for (auto it1 = options.models.at(i)->freeParams_.begin(); it1 != options.models.at(i)->freeParams_.end(); ++it1){//loop through parameters
		 		if (options.verbosity >=3){
		 		    cout<<" Free Parameter: " << it1->first << endl;
		 		}
		 		    if((find(paramList.begin(), paramList.end(), it1->first)) == paramList.end()){
		 		    	paramList.insert(paramList.end(), it1->first);
		 		    	paramCounter++;

		 		    }
		 		    particleCounter = 0;

		 		    for(auto it2 = subparticleCurrParamSets_.begin(); it2!=subparticleCurrParamSets_.end(); ++it2){//loops through subparticle IDs
		 		    	particleCounter++;
		 		    	//cout << "subparticleCurrParamSets_ first " << it2->first << endl;
		 		      	for(auto it3 = it2->second.begin(); it3 != it2->second.end(); ++it3){//loops through model IDs
		 		      		//cout << "model i " << i << endl;
		 		    		if(it3->first==i){
		 		    			if (options.verbosity >=3){
		 		    				cout << "parameter for index " <<  counter  << " : "  << it3->second[counter] << endl;
		 		    			}
		 		    			alignedResults[it2->first][i][it1->first] = it3->second[counter];
		 		    			// Particle -> Model -> Parameter -> value
		 		    		}
		 		    	}

		 		    }
					counter++; //free param counter
		 	 }
			 cout << endl;
		 }


	}

		 //outputFile << left << setw(16) << "Model" << left << setw(16) << "Particle" << left << setw(16) << "Fit";
 		 //cout << left << setw(16) << "Model" << left << setw(16) << "Particle" << left << setw(16) << "Fit";
		 outputFile << left << setw(16) << "Particle";
		 if (options.verbosity >=1){
			 cout << left << setw(16) << "Particle";
		 }
 		 for(unsigned int i = 0; i < options.models.size(); i++){
 			if (options.verbosity >=1){
 			 cout << "Fit_M";
 			 cout << left << setw(11) << i;
 			}
 			 outputFile << "Fit_M";
 			 outputFile << left << setw(11) <<  i;
 		 }


		 for(unsigned int i = 0; i < paramList.size(); i++){
			 if (options.verbosity >=1){
				 cout << left << setw(16) << paramList[i];
			 }
			 outputFile << left << setw(16) << paramList[i];


		 }
		 if (options.verbosity >=1){
			 cout << endl;
		 }
		 outputFile << endl;

		// int foundModel;

int foundParam;
int pcounter = 0; //index of the particles
int fcounter; //index of the fit results
//Raquel: this will print the output in the correct order, aligned
//Raquel TODO: save this printed output to the output file, and calculate overall fit, merging the fits for each parameter
//int fmid;
int fpid;
//int pmid;
//int psub;



multimap<double, unsigned int> subparticleFitIDMap;

if(options.fitType == "ga" || options.fitType == "pso" || options.fitType == "de") {

	 subparticleFitIDMap = currentsubswarmBestFits_;
	 currentsubswarmBestFits_copy = currentsubswarmBestFits_;
	 if(options.verbosity >= 3){
		 cout << "@@@@@@@@@@@@@@@@@@@@@ FIT SIZE @@@@@@@@@@@@@@@@@" << currentsubswarmBestFits_.size() << endl;
	 }

}

		 for(auto ri1 = alignedResults.begin(); ri1 != alignedResults.end(); ++ri1){//loop through particles


			 for (auto ri2 = ri1->second.begin(); ri2 != ri1->second.end(); ++ri2){//loop through models
				 pcounter++;// Raquel: subparticle id counter

				 //str =
				 //str = ;
				 //cout << "Model_";
				 //cout << left << setw(10) << ri2->first;
		 			if (options.verbosity >=1){

		 				cout << "Particle_";

		 				cout << left << setw(7) << ri1->first;
		 			}

				 outputFile << "Particle_";
				 outputFile << left << setw(7) << ri1->first;

				 //psub = fcalcsubParID(pcounter, pmid, options.models.size());
				 //pmid = fcalcMID(psub,options.models.size());

				 //Raquel: this fcounter will be equal to the subparticle ID
				 fcounter = 0;
				 for (auto fits = subparticleFitIDMap.begin(); fits != subparticleFitIDMap.end(); ++fits){//loops through the fit resuts
					 fcounter++;

					 //fmid = fcalcMID(fcounter,options.models.size());
					 fpid = fcalcParID(fits->second, options.models.size());

					 if(fpid==pcounter){
				 			if (options.verbosity >=1){

				 				cout << left << setw(16) << fits->first;
				 			}
						 outputFile << left << setw(16) << fits->first;

					 }

					 //cout << "fcounter " << fcounter << " pcounter " << pcounter << endl;

				 }


				 for(unsigned int i = 0; i < paramList.size(); i++){
					 foundParam = 0;
					 for(auto ri3 = ri2->second.begin(); ri3 != ri2->second.end(); ++ri3){//loop through parameters

						 if(paramList[i] == ri3->first){
					 			if (options.verbosity >=1){

					 				cout << left << setw(16) << ri3->second;
					 			}
							 outputFile << left << setw(16) << ri3->second;

							 foundParam = 1;
						 }

					 }
					 if(foundParam == 0){
				 			if (options.verbosity >=1){

				 				cout << left << setw(16) << "NA";
				 			}
						 outputFile << left << setw(16) << "NA";

					 }
				 }
		 			if (options.verbosity >=1){

		 				cout << endl;
		 			}
				 outputFile << endl;

				 break;
			 }

		 }


		 if(options.constraintWeight!=0){
	 			if (options.verbosity >=4){

	 				cout << "CONSTRAINT WEIGHT DIFFERENT OF ZERO: " << options.constraintWeight << endl;
	 			}
		 //this variable will store the number of constraints fulfilled per subparticle
		 vector<pair<int,float> > constraintsCount;
		 int mid = 0;
		 constraintsCount = resultChecking();

		 //sort results from the larger number of constraints to the smaller
		 sort(constraintsCount.begin(), constraintsCount.end(), [](const pair<int,float> &left, const pair<int,float> &right){return left.second > right.second;});


		 cout << "SubParID\tConstraints" << endl;
		 int constraintRank=0;
		 int previous=-1;

		 vector<pair<int,float> > subParRankCons;

		 //rank the results based on the number of constraints
		 for (unsigned int i = 0 ; i < constraintsCount.size(); i++){
		     cout << constraintsCount[i].first << " " << constraintsCount[i].second << "\n";
		     if(previous!=constraintsCount[i].second){
		    	 constraintRank++;
		     }
		     subParRankCons.push_back(make_pair(constraintsCount[i].first, constraintRank));
		     previous=constraintsCount[i].second;
		 }

		 cout << "SubParID\tConstraintRank" << endl;

		 for (unsigned int i = 0 ; i < subParRankCons.size(); i++){
		     cout << subParRankCons[i].first << " " << subParRankCons[i].second << "\n";
		 }

		 vector<pair<int,float> > fitCount;

		 fcounter = 0;





		 for (auto itr = subparticleFitIDMap.begin(); itr != subparticleFitIDMap.end(); ++itr){
			 fcounter++;
			 fitCount.push_back(make_pair(itr->second,itr->first));

		 }

		 //sort results from the smaller fit to the larger fit
		 sort(fitCount.begin(), fitCount.end(), [](const pair<int,float> &left, const pair<int,float> &right){return left.second < right.second;});


		 cout << "SubParID\tFit" << endl;
		 int fitRank=0;
		 float previousFit=-1;

		 vector<pair<int,float> > subParRankFit;

		 //rank the results based on the Fit value
		 for (unsigned int i = 0 ; i < fitCount.size(); i++){
		     cout << fitCount[i].first << " " << fitCount[i].second << "\n";
		     if(previousFit!=fitCount[i].second){
		    	 fitRank++;
		     }
		     subParRankFit.push_back(make_pair(fitCount[i].first, fitRank));
		     previousFit=fitCount[i].second;
		 }


		 cout << "SubParID\tFitRank" << endl;

		 for (unsigned int i = 0 ; i < subParRankFit.size(); i++){
		     cout << subParRankFit[i].first << " " << subParRankFit[i].second << "\n";
		 }


		 vector<pair<int,float> > finalScore;
		 float consWeight = options.constraintWeight; //Here it defines how important the consraints will be
		 //put a large consWeight number = the constraint ranks will prevail
		 //put a small consWeight number = the fit ranks will prevail
		 float fitWeight = 1-options.constraintWeight; //Here it defines how important the fits will be

		 float scoreTmp = 0;

		 //this loop will combine both scores based on fit ranks and fulfilled constraints ranks
/*		 for (int i = 0 ; i < subParRankCons.size(); i++){


			 for (int j = 0 ; j < subParRankFit.size(); j++){



				 if(subParRankCons[i].first == subParRankFit[j].first){
					 //Here it defines how important the cnsraints will be
					 scoreTmp = ((subParRankCons[i].second*consWeight) + (subParRankFit[j].second * fitWeight));
					 finalScore.push_back(make_pair(subParRankCons[i].first, scoreTmp));

				 }else{
					 //if the model that we are looking at belongs to the wild type
					 mid = fcalcMID(subParRankFit[j].first,options.models.size());

					 //if(mid == 0 && (subParRankCons[i].first-1) == subParRankFit[j].first){
					 if( find(subParRankCons.begin(), subParRankCons.end(), subParRankFit[j].first) == subParRankCons.end() ){
						 //cout << "mid=0 subParRankFit[j].first= " << subParRankFit[j].first << endl;
						 cout << " subpar subParRankFit[j].first " << subParRankFit[j].first << " mid " << mid << endl;
						 scoreTmp = subParRankFit[j].second*fitWeight;
						 finalScore.push_back(make_pair(subParRankFit[j].first, scoreTmp));
					 }

				 }

			 }


		 }
*/
		 int found;

		 for (unsigned int i = 0 ; i < subParRankFit.size(); i++){

			 found = 0;
			 for (unsigned int j = 0 ; j < subParRankCons.size(); j++){

				 if(subParRankCons[j].first == subParRankFit[i].first){
					 //Here it defines how important the cnsraints will be
					 scoreTmp = ((subParRankCons[j].second*consWeight) + (subParRankFit[i].second * fitWeight));
					 finalScore.push_back(make_pair(subParRankCons[j].first, scoreTmp));

					 found = 1;

				 }

			 }
			 //if the subparticle is not in the list of constraints, keep the original rank
			 if(found==0){

				 scoreTmp = subParRankFit[i].second;
				 //scoreTmp = subParRankFit[j].second*fitWeight;
				 finalScore.push_back(make_pair(subParRankFit[i].first, scoreTmp));


			 }



		 }

		 //sort results from the smaller fit to the larger fit
		 sort(finalScore.begin(), finalScore.end(), [](const pair<int,float> &left, const pair<int,float> &right){return left.second < right.second;});


		 cout << "SubParID\tFinal_Rank_Score" << endl;
		 fitRank=0;
		 previousFit=-1;

		 //vector<pair<int,float> > subParRankFinal;
		 subParRankFinal.clear();
		 //rank the results based on the Fit value
		 for (unsigned int i = 0 ; i < finalScore.size(); i++){
		     cout << finalScore[i].first << " " << finalScore[i].second << "\n";
		     if(previousFit!=finalScore[i].second){
		    	 fitRank++;
		     }
		     subParRankFinal.push_back(make_pair(finalScore[i].first, fitRank));
		     previousFit=finalScore[i].second;
		 }

		 cout << "SubParID\tFinal_Rank" << endl;

		 for (unsigned int i = 0 ; i < subParRankFinal.size(); i++){
		     cout << subParRankFinal[i].first << " " << subParRankFinal[i].second << "\n";
		 }



		 //Print the reranked results

		 string outdir = options.outputDir + "/" + options.jobName + "/" + toString(currentGeneration-1) + "/";

		 string outname = outdir + "Ranked_results.txt";


		 ofstream outFile;

		 outFile.open(outname);

		 cout << "Particle\tModel\tFit_Rank\tConstraint_Rank\tFinal_Rank\tConstraints" << endl;
		 outFile << "Particle\tModel\tFit_Rank\tConstraint_Rank\tFinal_Rank\tConstraints" << endl;

		 //int maxPar = fcalcParID(subparticleFitIDMap.size(), options.models.size());

		 //int mid;
		 int pid;

		 for (unsigned int i = 0 ; i < subParRankFinal.size(); i++){

			 mid = fcalcMID(subParRankFinal[i].first,options.models.size());
			 pid = fcalcParID(subParRankFinal[i].first, options.models.size());

			 for (unsigned int j = 0 ; j < subParRankFit.size(); j++){


				 found = 0;
				 for (unsigned int k = 0 ; k < subParRankCons.size(); k++){

					 if(subParRankFinal[i].first==subParRankFit[j].first && subParRankFinal[i].first==subParRankCons[k].first){

						 for(unsigned int r=0; r < constraintsCount.size(); r++){

							 if(subParRankFinal[i].first ==  constraintsCount[r].first ){


								 cout << pid << "\t" << mid << "\t" << subParRankFit[j].second << "\t" << subParRankCons[k].second << "\t" << subParRankFinal[i].second << "\t" <<constraintsCount[r].second <<  endl;
								 outFile << pid << "\t" << mid << "\t" << subParRankFit[j].second << "\t" << subParRankCons[k].second << "\t" << subParRankFinal[i].second << "\t" << constraintsCount[r].second << endl;

								 found = 1;


							 }


						 }

					 }

				 }


				 if(found==0 && subParRankFinal[i].first==subParRankFit[j].first){

					 cout << pid << "\t" << mid << "\t" << subParRankFit[j].second << "\t" << "-" << "\t" << subParRankFinal[i].second << "\t-" << endl;
					 outFile << pid << "\t" << mid << "\t" << subParRankFit[j].second << "\t" << "-" << "\t" << subParRankFinal[i].second << "\t-" << endl;


				 }


			 }



		 }

		 }
		 /*
		 //rank results
		 int maxPar = fcalcParID(allGenFits.size(), options.models.size());


		 vector<pair<int,float> > sortedFits;
		 for (auto itr = allGenFits.begin(); itr != allGenFits.end(); ++itr){
			 sortedFits.push_back(*itr);
		 }


		 sort (sortedFits.begin(),sortedFits.end());
		 cout << "sorted result" << endl;
		 for (int i = 0 ; i < sortedFits.size(); i++){
		     cout << sortedFits[i].first << " " << sortedFits[i].second << "\n";
		 }

		 for(int i=1; i <= maxPar; i++ ){

			 fcounter = 0;
			 for (auto fits = allGenFits.begin(); fits != allGenFits.end(); ++fits){//loops through the fit resuts
				 fcounter++;

				 fmid = fcalcMID(fcounter,options.models.size());
				 fpid = fcalcParID(fcounter, options.models.size());

				 if(i==fmid){
					 cout << left << setw(16) << fits->first;
					 outputFile << left << setw(16) << fits->first;

				 }

				 //cout << "fcounter " << fcounter << " pcounter " << pcounter << endl;

			 }


		 }
		*/

//Raquel: here ends the block of code to fix the alignment problem with the parameter values

		outputFile.close();
	}
	else {
		cout << "Warning: couldn't open: " << outputPath << " to write fit summary." << endl;
	}



}


//Raquel: result checking function, will call Evaluate.cpp
vector<pair<int,float> > Swarm::resultChecking(){
float result = 0;
	string outdir = options.outputDir + "/" + options.jobName + "/" + toString(currentGeneration-1) + "/";

	string inputFile1;
	string inputFile2;

	string path;
	string basename1;
	string basename2;

	string outname;

	vector<pair<int,float> > constraintsCount;
	map<int,string> constraint;

	int pnumber;
	int mnumber;
	int subnumber;
	float iteration1;
	float iteration2;

	vector<string> constraintParams;

	/* for(auto it = subparticleCurrParamSets_.begin(); it!=subparticleCurrParamSets_.end(); ++it){
		 //particle, model, parameters
		 for(int i = 0; i < options.models.size(); i++){
			 subnumber = fcalcsubParID(it->first, i, options.models.size());

			 constraintsCount.push_back(make_pair(subnumber,0));
		 }
	 }
*/
	 for(auto it = subparticleCurrParamSets_.begin(); it!=subparticleCurrParamSets_.end(); ++it){
		 cout << "particle :" << it->first << endl;

		 for(unsigned int i = 0; i < options.models.size(); i++){

			 //cout <<  "model: " << i << endl;
			 for(unsigned int j= 0; j < options.models.size(); j++){

				 //if(i!=j && i<j){

					 path = options.models.at(i)->getName();
					 basename1 = getFilename(path);
					 inputFile1 = outdir + basename1 + "_" + toString(it->first) + "_" + "1" + ".gdat";

					 path = options.models.at(j)->getName();
					 basename2 = getFilename(path);
					 inputFile2 = outdir + basename2 + "_" + toString(it->first) + "_" + "1" + ".gdat";

					 cout << "infile 1: " << inputFile1 << endl;
					 cout << "infile 2: " << inputFile2 << endl;
					 //		std::map<int,string> constraints_; //Raquel: added constraint options support
					 outname = outdir + basename1 + "_vs_" + basename2 + "_" + toString(it->first) + ".txt";
					 //add outname as the fourth input
					 //Raquel: adding support to receive different model numbers and time points as input
					 //cout << "Starting model checking, testing for " << options.constraints_.size() << " contraints." << endl;

					 for(auto constraintIt = options.constraints_.begin(); constraintIt != options.constraints_.end(); ++constraintIt){

						 constraintParams = split_string(constraintIt->second, " ");

						 //cout << "consParam 0 " << constraintParams[0] << "consParam 1 " <<  constraintParams[1] << "consParam 2 " <<  constraintParams[2] << "consParam 3 " << constraintParams[3] << "consParam 4 " << endl;
						 constraint.insert(make_pair(constraintIt->first,constraintParams[0]));
						 //cout << "COUNSTRAINT size from swarm: " << constraint.size() << endl;

						 if((i == (unsigned) atoi(constraintParams[1].c_str()) && j== (unsigned) atoi(constraintParams[2].c_str()))){
						 //if((i==atoi(constraintParams[1].c_str()) && j==atoi(constraintParams[2].c_str())) || (i==atoi(constraintParams[1].c_str()) && i==atoi(constraintParams[2].c_str())) || (j==atoi(constraintParams[1].c_str()) && j==atoi(constraintParams[2].c_str()))){

							 iteration1 = atof(constraintParams[3].c_str());
							 iteration2 = atof(constraintParams[4].c_str());


							 result += evaluateResults(inputFile1,inputFile2,constraint, outname, iteration1, iteration2);


							 pnumber = it->first;
							 mnumber = i;
							 subnumber = fcalcsubParID(pnumber, mnumber, options.models.size());

							 //for(auto conscountIt = constraintsCount.begin(); conscountIt != constraintsCount.end(); ++conscountIt){
								 //if(conscountIt->first == subnumber){
									 //constraintsCount.erase(conscountIt);
									 //cout << "adding result " << toString(conscountIt->second+result) << " to constraint count to subpar" << subnumber << endl;
									 //constraintsCount.push_back(make_pair(conscountIt->first, conscountIt->second+result));
								// }


							 //}
							 //constraintsCount[subnumber] = constraintsCount[subnumber]+(float)result;
							 //constraintsCount.push_back(make_pair(subnumber,(float)result));
						 }
						 constraint.clear();



					 }





				 //}

			 }

			 pnumber = it->first;
			 mnumber = i;
			 subnumber = fcalcsubParID(pnumber, mnumber, options.models.size());

			 constraintsCount.push_back(make_pair(subnumber, result));

			 result = 0;
			 //Raquel: added this break so only the first model(wild type), is compared versus the others for now
			//break;

		 }


	 }

	 //for(auto it2 = particleIterationCounter_.begin(); it2!=particleIterationCounter_.end(); it2++){
	//	 cout << "Iteration: " << it2->first << endl;

	 //}


	// for(auto i = constraintsCount.begin(); i != constraintsCount.end(); ++i){
	//	 cout << "constraintCount values: ";
	//	 cout << i->first << " ";
	//	 cout << i->second << endl;
//
//	 }

	 return constraintsCount;

}


void Swarm::outputRunSummary() {
	// Output first two fields of header
	cout<<"Swarm::outputRunSummary() may need some modifications ...."<<endl; //mypause(); //razi TODO: later test
	//Raquel

	cout << left << setw(16) << "Fit" << left << setw(16) << "Iteration" << endl;

	cout << left << setw(16) << "Fit" << left << setw(16) << "Iteration";

	// Output parameter names
	//for (autoDE i: options.model->freeParams_) {
	for (auto i = options.freeParams_.begin(); i != options.freeParams_.end(); ++i) {

		cout  << left << setw(16) << i->first << endl;

		cout << left << setw(16) << i->first;
	}

	cout << endl;

	vector<string> paramVals;

	//for (auto i: allGenFits) {
	for (auto i = allGenFits.begin(); i != allGenFits.end(); ++i) {   //razi: TODO later check if allGenFits filled correctly when different model have different list of free parameters
		split(i->second, paramVals);

		cout << left << setw(16) << i->first << left << setw(16) << paramVals[0];

		for (unsigned int i = 1; i < paramVals.size(); i++) {
			cout << left << setw(16) << stod(paramVals[i]);
		}

		paramVals.clear();

		cout << endl;
	}
}





vector<double> Swarm::normalizeParams(vector<double> oldParams) {
	vector<double> newParams;

	outputError("needs modifications to support multiple files");
/*	unsigned int d = 0;
	for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
		//cout << "denormalized: " << oldParams[d] << endl;
		//cout << "max: " << param->second->getGenMax() << " min: " << param->second->getGenMin() << " diff: " << (param->second->getGenMax() - param->second->getGenMin()) << endl;
		//newParams.push_back((oldParams[d++] - param->second->getGenMin()) / (param->second->getGenMax() - param->second->getGenMin()));
		newParams.push_back(normalizeParam(oldParams[d], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog()));
		//cout << "normalized: " << *(newParams.end() - 1) << endl;
		++d;
	}
	*/
	return newParams;
}

vector<double> Swarm::deNormalizeParams(vector<double> oldParams) {
	vector<double> newParams;

/*	outputError("needs modifications to support multiple files");
	unsigned int d = 0;
	for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
		//cout << "normalized: " << oldParams[d] << endl;
		//newParams.push_back(pow(10, log10(param->second->getGenMin()) + oldParams[d++] * (log10(param->second->getGenMax()) - log10(param->second->getGenMin() ))));
		//newParams.push_back( param->second->getGenMin() + oldParams[d++] * (param->second->getGenMax() - param->second->getGenMin()));
		newParams.push_back(deNormalizeParam(oldParams[d], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog()));
		//cout << "denormalized: " << *(newParams.end() - 1) << endl;
		++d;
	}
*/
	return newParams;
}



void Swarm::outputBootstrapSummary() {
	outputError("Swarm::outputBootstrapSummary()");
/*
	string outputPath = options.outputDir + "/" + options.jobName + "_bootstrap/bootstrap_results.txt";
	ofstream outputFile;
	outputFile.open(outputPath, ofstream::out | std::ofstream::app);

	if (outputFile.is_open()) {
		outputFile.precision(8);
		outputFile.setf(ios::scientific);

		if (bootstrapCounter == 0) {
			// Output first two fields of header
			outputFile << left << setw(16) << "Fit";

			// Output parameter names
			for (auto i = options.model->freeParams_.begin(); i != options.model->freeParams_.end(); ++i) {
				outputFile << left << setw(16) << i->first;
			}
			outputFile << endl;
		}

		vector<string> paramVals;
		split(allGenFits.begin()->second, paramVals);
		outputFile << left << setw(16) << allGenFits.begin()->first;
		for (unsigned int i = 1; i < paramVals.size(); i++) {
			outputFile << left << setw(16) << stod(paramVals[i]);
		}
		outputFile << endl;
		outputFile.close();
	}
	else {
		cout << "Warning: couldn't open: " << outputPath << " to write fit summary." << endl;
	}
	*/
}

void Swarm::killAllParticles(int tag) {
	int p, sp;
	unsigned int mid;
	for (p = 1; p <= options.swarmSize; ++p) {
		//cout << "killing " << p << " with tag: " << tag << endl;
		for(mid=0; mid < options.models.size(); mid++){//Raquel modified this
			//sp = (p-1)*getNumModels()+ mid + 1; //Raquel replaced by code bellow
			sp = fcalcsubParID(p, mid, options.models.size());

			if(options.verbosity>=3){
				cout << "RAQUEL sending message to subPar " << sp << "message=" << tag << endl;
			}
			swarmComm->sendToSwarm(0, sp, tag, false, swarmComm->univMessageSender);
		}
	}
	swarmComm->univMessageSender.clear();
	//cout << "killed all particles" << endl;
}


vector<double> Swarm::mutateParticleDE(unsigned int particle, float mutateFactor) {

	if (options.verbosity >= 3) {
		cout << "Mutating particle " << particle << endl;
	}

	int parID = 0;
	vector<double> mutatedParams;
	unsigned int currIsland = 0;
	unsigned int particlesPerIsland = 0;

	if (options.fitType == "de") {
		currIsland = particleToIsland_[particle]; //Raquel changed, was .at(particle)
		particlesPerIsland = options.swarmSize / options.numIslands;
	}

	// f between 0.1 and 0.9
	//float f = 0.1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.9-0.1)));
	float f;
	if (mutateFactor) {
		f = mutateFactor;
	}
	else {
		f = 0.9;
	}

	boost::random::uniform_int_distribution<int> unif(0, particlesPerIsland - 1);
	if(options.verbosity>=3){
		cout << "done setting uniform_int_distribution" << endl;
	}
	//Raquel: adding support to use constraints when selecting the best particles in DE
	if (options.mutateType == 1) {
		int bestParticle = 0;

		if(options.constraints_.size()>0){//If the user provided any constraint
			//use the combined rank between fit values and constraint values to select the best particle
			//The ranks are sorted from the best fit to the worst fit values, so the loop will stop at the first particle that it finds, with the best fit as possible per iland
			for (auto fitVal = subParRankFinal.begin(); fitVal != subParRankFinal.end(); ++fitVal) {
				//cout << "Raquel setting particle to island value" << endl;

				parID = fcalcParID(fitVal->first, options.models.size());

				if (fitVal->second > 0 && particleToIsland_[parID] == currIsland) {
					//cout << "Raquel setting particle to island value done" << endl;

					bestParticle = fitVal->second;
					//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
					break;
				}

			}


		}else{

			for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
				if(options.verbosity>=4){
					cout << "Raquel setting particle to island value" << endl;
				}
				if (fitVal->first > 0 && particleToIsland_[fitVal->second] == currIsland) {
					if(options.verbosity>=4){
						cout << "Raquel setting particle to island value done" << endl;
					}
					bestParticle = fitVal->second;
					if(options.verbosity>=4){
						cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
					}
					break;
				}

			}


		}



		vector<double> bestParamSet = particleCurrParamSets_.at(bestParticle);
		unsigned int pi = 0;
		for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) { //for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;

			while (p1 == p2) {
				p1 = unif(generalRand);
				p2 = unif(generalRand);
			}

			double p1Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p1])[pi];
			double p2Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p2])[pi];

			mutatedParams.push_back(bestParamSet[pi] + f * (p1Param - p2Param));
			++pi;
		}
	}
	else if (options.mutateType == 2) {
		unsigned int pi = 0;
		int bestParticle = 0;
		//for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
		//	if (fitVal->first > 0 && particleToIsland_.at(fitVal->second) == currIsland) {
		//		bestParticle = fitVal->second;
		//	}
		//}
		if(options.constraints_.size()>0){//If the user provided any constraint
					//use the combined rank between fit values and constraint values to select the best particle
					//mutate type 2 will take the worst particle for the current island and mutate it
					for (auto fitVal = subParRankFinal.begin(); fitVal != subParRankFinal.end(); ++fitVal) {
						if(options.verbosity>=4){
							cout << "Raquel setting particle to island value" << endl;
						}
						parID = fcalcParID(fitVal->first, options.models.size());

						if (fitVal->second > 0 && particleToIsland_.at(parID) == currIsland) {
							if(options.verbosity>=4){
								cout << "Raquel setting particle to island value done" << endl;
							}
							bestParticle = fitVal->second;
							//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
							//break;
						}

					}


				}else{


					for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
						//cout << "Raquel setting particle to island value" << endl;
						if (fitVal->first > 0 && particleToIsland_[fitVal->second] == currIsland) {
							//cout << "Raquel setting particle to island value done" << endl;

							bestParticle = fitVal->second;
							//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
							//break;
						}

					}


				}


		vector<double> bestParamSet = particleCurrParamSets_[bestParticle];
		for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) {//for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			int p4 = 0;
			while (p1 == p2 || p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4 || p3 == p4) {
				p1 = unif(generalRand);
				p2 = unif(generalRand);
				p3 = unif(generalRand);
				p4 = unif(generalRand);
			}

			double p1Param = particleCurrParamSets_[islandToParticle_[currIsland][p1]][pi];
			double p2Param = particleCurrParamSets_[islandToParticle_[currIsland][p2]][pi];
			double p3Param = particleCurrParamSets_[islandToParticle_[currIsland][p3]][pi];
			double p4Param = particleCurrParamSets_[islandToParticle_[currIsland][p4]][pi];

			mutatedParams.push_back(bestParamSet[pi] + f * (p1Param - p2Param) + f * (p3Param - p4Param));
			++pi;
		}
	}
	else if (options.mutateType == 3) {
		unsigned int pi = 0;
		int bestParticle = 0;
		//for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
		//	if (fitVal->first > 0 && particleToIsland_.at(fitVal->second) == currIsland) {
		//		bestParticle = fitVal->second;
		//	}
		//}

		if(options.constraints_.size()>0){//If the user provided any constraint
					//use the combined rank between fit values and constraint values to select the best particle
					//mutate type 3 will take the worst particle for the current island and mutate it
					for (auto fitVal = subParRankFinal.begin(); fitVal != subParRankFinal.end(); ++fitVal) {
						if(options.verbosity>=4){
							cout << "Raquel setting particle to island value" << endl;
						}
						parID = fcalcParID(fitVal->first, options.models.size());

						if (fitVal->second > 0 && particleToIsland_[parID] == currIsland) {
							if(options.verbosity>=4){
								cout << "Raquel setting particle to island value done" << endl;
							}
							bestParticle = fitVal->second;
							//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
							//break;
						}

					}


				}else{


					for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
						//cout << "Raquel setting particle to island value" << endl;
						if (fitVal->first > 0 && particleToIsland_[fitVal->second] == currIsland) {
							//cout << "Raquel setting particle to island value done" << endl;

							bestParticle = fitVal->second;
							//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
							//break;
						}

					}


				}


		vector<double> bestParamSet = particleCurrParamSets_[bestParticle];
		for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) { //for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			while (p1 == p2 || p1 == p3 || p2 == p3) {
				p1 = unif(generalRand);
				p2 = unif(generalRand);
				p3 = unif(generalRand);
			}

			double p1Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p1])[pi];
			double p2Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p2])[pi];
			double p3Param = particleCurrParamSets_.at(islandToParticle_.at(currIsland)[p3])[pi];

			mutatedParams.push_back(p1Param + f * (p2Param - p3Param));
			++pi;
		}
	}
	else if (options.mutateType == 4) {
		unsigned int pi = 0;
		int bestParticle = 0;
		//for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
		//	if (fitVal->first > 0 && particleToIsland_.at(fitVal->second) == currIsland) {
		//		bestParticle = fitVal->second;
		//	}
		//}

		if(options.constraints_.size()>0){//If the user provided any constraint
					//use the combined rank between fit values and constraint values to select the best particle
					//mutate type 4 will take the worst particle for the current island and mutate it
					for (auto fitVal = subParRankFinal.begin(); fitVal != subParRankFinal.end(); ++fitVal) {
						if(options.verbosity>=4){
							cout << "Raquel setting particle to island value" << endl;
						}
						parID = fcalcParID(fitVal->first, options.models.size());

						if (fitVal->second > 0 && particleToIsland_[parID] == currIsland) {
							if(options.verbosity>=4){
								cout << "Raquel setting particle to island value done" << endl;
							}
							bestParticle = fitVal->second;
							//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
							//break;
						}

					}


				}else{


					for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
						//cout << "Raquel setting particle to island value" << endl;
						if (fitVal->first > 0 && particleToIsland_[fitVal->second] == currIsland) {
							//cout << "Raquel setting particle to island value done" << endl;

							bestParticle = fitVal->second;
							//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
							//break;
						}

					}


				}



		vector<double> bestParamSet = particleCurrParamSets_[bestParticle];
		for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) {//for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			int p4 = 0;
			int p5 = 0;
			// This (and above instances) should be optimized
			while (p1 == p2 || p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4 || p3 == p4 || p1 == p5 || p2 == p5 || p3 == p5 || p4 == p5) {
				p1 = unif(generalRand);
				p2 = unif(generalRand);
				p3 = unif(generalRand);
				p4 = unif(generalRand);
				p5 = unif(generalRand);
			}

			double p1Param = particleCurrParamSets_[islandToParticle_[currIsland][p1]][pi];
			double p2Param = particleCurrParamSets_[islandToParticle_[currIsland][p2]][pi];
			double p3Param = particleCurrParamSets_[islandToParticle_[currIsland][p3]][pi];
			double p4Param = particleCurrParamSets_[islandToParticle_[currIsland][p4]][pi];
			double p5Param = particleCurrParamSets_[islandToParticle_[currIsland][p5]][pi];

			mutatedParams.push_back(p1Param + f * (p2Param - p3Param) + f * (p4Param - p5Param));
			++pi;
		}
	}

	return mutatedParams;
}

vector<double> Swarm::mutateParticleSA(unsigned int particle, float mutateFactor) {

	if (options.verbosity >= 3) {
		cout << "Mutating particle " << particle << endl;
	}

	vector<double> mutatedParams;

	// f between 0.1 and 0.9
	//float f = 0.1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.9-0.1)));
	float f;
	if (mutateFactor) {
		f = mutateFactor;
	}
	else {
		f = 0.9;
	}

	boost::random::uniform_int_distribution<int> unif(1, options.swarmSize);

	if (options.mutateType == 1) {
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0 ) {
				bestParticle = fitVal->second;
				//cout << "setting particle " << fitVal->second << " as best with fit of " << fitVal -> first << endl;
				break;
			}
		}
		vector<double> bestParamSet = normalizeParams(particleCurrParamSets_.at(bestParticle));
		unsigned int pi = 0;
		for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) {//for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;

			while (p1 == p2) {
				p1 = unif(generalRand);
				p2 = unif(generalRand);
			}

			//cout << "p1: " << p1 << " p2: " << p2 << " pi: " << pi << endl;
			double p1Param = normalizeParam(particleCurrParamSets_.at(p1)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p2Param = normalizeParam(particleCurrParamSets_.at(p2)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			mutatedParams.push_back(bestParamSet[pi] + f * (p1Param - p2Param));
			++pi;
		}
	}
	else if (options.mutateType == 2) {
		unsigned int pi = 0;
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0) {
				bestParticle = fitVal->second;
			}
		}
		vector<double> bestParamSet = normalizeParams(particleCurrParamSets_.at(bestParticle));
		for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) { //for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			int p4 = 0;
			while (p1 == p2 || p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4 || p3 == p4) {
				p1 = unif(generalRand);
				p2 = unif(generalRand);
				p3 = unif(generalRand);
				p4 = unif(generalRand);
			}

			double p1Param = normalizeParam(particleCurrParamSets_.at(p1)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p2Param = normalizeParam(particleCurrParamSets_.at(p2)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p3Param = normalizeParam(particleCurrParamSets_.at(p3)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p4Param = normalizeParam(particleCurrParamSets_.at(p4)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());

			mutatedParams.push_back(bestParamSet[pi] + f * (p1Param - p2Param) + f * (p3Param - p4Param));
			++pi;
		}
	}
	else if (options.mutateType == 3) {
		unsigned int pi = 0;
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0) {
				bestParticle = fitVal->second;
			}
		}
		vector<double> bestParamSet = normalizeParams(particleCurrParamSets_.at(bestParticle));
		for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) { //for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			while (p1 == p2 || p1 == p3 || p2 == p3) {
				p1 = unif(generalRand);
				p2 = unif(generalRand);
				p3 = unif(generalRand);
			}

			double p1Param = normalizeParam(particleCurrParamSets_.at(p1)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p2Param = normalizeParam(particleCurrParamSets_.at(p2)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p3Param = normalizeParam(particleCurrParamSets_.at(p3)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());

			mutatedParams.push_back(p1Param + f * (p2Param - p3Param));
			++pi;
		}
	}
	else if (options.mutateType == 4) {
		unsigned int pi = 0;
		int bestParticle = 0;
		for (auto fitVal = particleBestFitsByFit_.begin(); fitVal != particleBestFitsByFit_.end(); ++fitVal) {
			if (fitVal->first > 0) {
				bestParticle = fitVal->second;
			}
		}
		vector<double> bestParamSet = normalizeParams(particleCurrParamSets_.at(bestParticle));
		for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) {  //for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
			int p1 = 0;
			int p2 = 0;
			int p3 = 0;
			int p4 = 0;
			int p5 = 0;
			// This (and above instances) should be optimized
			while (p1 == p2 || p1 == p3 || p1 == p4 || p2 == p3 || p2 == p4 || p3 == p4 || p1 == p5 || p2 == p5 || p3 == p5 || p4 == p5) {
				p1 = unif(generalRand);
				p2 = unif(generalRand);
				p3 = unif(generalRand);
				p4 = unif(generalRand);
				p5 = unif(generalRand);
			}

			double p1Param = normalizeParam(particleCurrParamSets_.at(p1)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p2Param = normalizeParam(particleCurrParamSets_.at(p2)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p3Param = normalizeParam(particleCurrParamSets_.at(p3)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p4Param = normalizeParam(particleCurrParamSets_.at(p4)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());
			double p5Param = normalizeParam(particleCurrParamSets_.at(p5)[pi], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog());

			mutatedParams.push_back(p1Param + f * (p2Param - p3Param) + f * (p4Param - p5Param));
			++pi;
		}
	}

	return mutatedParams;
}

vector<double> Swarm::crossoverParticleDE(unsigned int particle, vector<double> mutationSet, float crossoverRate, bool normalize) {
	if (options.verbosity >= 3) {
		cout << "Crossing over particle " << particle << endl;
	}

	// Uses binomial crossover
	boost::random::uniform_int_distribution<int> intDist(0, this->getNumFreeParams() - 1);//boost::random::uniform_int_distribution<int> intDist(0, options.model->getNumFreeParams() - 1);
	unsigned int randParam = intDist(generalRand);
	boost::random::uniform_real_distribution<float> floatDist(0, 1);

	float cr;
	if (crossoverRate) {
		cr = crossoverRate;
	}
	else {
		cr = options.cr;
	}

	//cout << "crossing over particle " << particle << ". randParam is " << randParam << endl;
	vector<double> newParamSet;
	unsigned int p = 0;
	for (auto param = this->getFreeParams_().begin(); param != this->getFreeParams_().end(); ++param) { //for (auto param = options.model->getFreeParams_().begin(); param != options.model->getFreeParams_().end(); ++param) {
		//for (unsigned int p = 0; p < options.model->getNumFreeParams(); ++p) {
		float rand = floatDist(generalRand);
		//cout << "p: " << p << ". rand: " << rand << endl;

		if (rand < cr || p == randParam) {
			newParamSet.push_back(mutationSet[p]);
			//cout << "pushing back " << mutationSet[p] << " from mutation set" << endl;
		}
		else {
			if (normalize) {
				newParamSet.push_back(normalizeParam(particleCurrParamSets_.at(particle)[p], param->second->getGenMin(), param->second->getGenMax(), param->second->getIsLog()));
			}
			else {
				newParamSet.push_back(particleCurrParamSets_.at(particle)[p]);
			}
			//cout << "pushing back " << particleCurrParamSets_.at(particle)[p] << " from current set" << endl;
		}
		++p;
	}

	return newParamSet;
}


void Swarm::breedGenerationGA(vector<unsigned int> children) {
	//check all messages that are sent by the master
	//send double messages in case you have subparticles
	//one message per subparticle, one subparticle per model
	//use functions in the Util
	//function parse model, creates the file qith the contents of the model
	//children
	unsigned int numFinishedBreeding = 0;
	if(options.verbosity>=3){

		cout << "RAQUEL number of children received: " << children.size() << endl;
	}
	//later:
	//where the parameters are comming from, full list of consolidate

	// TODO: Need to take into consideration failed particles? Definitely need to in first generation.
	if (options.verbosity >= 3) {
		cout << "Breeding generation, generation:" << currentGeneration <<"  swarm size:"<<options.swarmSize<<"   keep size:"<< options.keepParents<<"  children size:"<< children.size()<< " " <<endl;

		for (std::map<std::string, FreeParam*>::iterator it= this->options.freeParams_.begin(); it!=this->options.freeParams_.end(); it++) //for (std::map<std::string, FreeParam*>::iterator it= options.model->freeParams_.begin(); it!=options.model->freeParams_.end(); it++)
			cout<<"Free parameter is: "<< it->first<<"  range:"<< it->second->getGenMin()<<" to "<< it->second->getGenMax() <<endl;

	}

	if (children.size() == 0) {
		for (int i = 1; i <= options.swarmSize; ++i) {
			children.push_back(i);
		}
	}

	if(options.verbosity>=3){

		cout << "RAQUEL swarmsize: " << options.swarmSize << endl;
		cout << "RAQUEL children.size: " << children.size() << endl;
	}
	std::map<unsigned int, std::vector<double> > particleNewParamSets;

	// Create the output directory for the next generation
	if (!checkIfFileExists(options.jobOutputDir + toString(currentGeneration))) {
		string createDirCmd = "mkdir " + options.jobOutputDir + toString(currentGeneration);
		int retryCounter = 0;
		while (runCommand(createDirCmd) != 0 && !checkIfFileExists(options.jobOutputDir + toString(currentGeneration))) {
			if(++retryCounter >= 100) {
				outputError("Error: Couldn't create " + options.jobOutputDir + toString(currentGeneration) + " to hold next generation's output.");
			}
			//sleep(1);
			if(options.verbosity>=3){

				cout << "Trying again to create dir" << endl;
			}
		}
	}
	else if (options.verbosity >= 4) cout<<"Swarm::breedGenerationGA-AAA4, dir already exist\n.";

	unsigned int parentPoolSize = options.swarmSize;
	parentPoolSize = parentPoolSize - options.keepParents;

	double maxWeight;
	multimap<double, unsigned int>::iterator w;

	//use model checking results if user provided constraints
	if(options.constraints_.size()>=1){
		int last = (options.swarmSize*options.models.size()) - 1;

		maxWeight = subParRankFinal[last].second;
		if(options.verbosity>=4){

			cout << "WWWWWWWWWWWWWWWWWWWW subParRankFinal[last]     " << subParRankFinal[last].second << "    WWWWWWWWWWWWWWWWWWW     " << subParRankFinal[last].first << endl;

		}

	}else{
		// Create an iterator to our fit list, use it to get our maximum fit value, then reset the iterator to the beginning of the list
		//Raquel updated to use subparticles

		w = subswarmBestFits_.begin();
		advance(w, (options.swarmSize*options.models.size()) - 1);
		maxWeight = w->first;
		w = subswarmBestFits_.begin();
		if(options.verbosity>=4){

			cout << "WWWWWWWWWWWWWWWWWWWW      " << w->first << "    WWWWWWWWWWWWWWWWWWW     " << w->second << endl;
		}
	}


	// Fill the second element of the weight map with difference between maxWeight and fit value
	multimap<double, unsigned int> weightDiffs;
	double diff;
	double weightSum = 0;

	if(options.constraints_.size()>=1){
		//Raquel updated to use subparticles
		for (unsigned int i = 0; i < options.swarmSize*options.models.size(); ++i) {
			diff = maxWeight - subParRankFinal[i].second;
			weightSum += diff;
			weightDiffs.insert(pair<double, unsigned int>(diff, subParRankFinal[i].first));
		}

	}else{
		//Raquel updated to use subparticles
		for(auto w = subswarmBestFits_.begin(); w!=subswarmBestFits_.end();++w){
		//for (unsigned int i = 0; i < options.swarmSize*options.models.size(); ++i) {
			diff = maxWeight - w->first;
			if(options.verbosity>=4){

				cout << "diff= " << diff << endl;
				cout << "maxWeight = " << maxWeight << endl;
				cout << "w->first= " << w->first << endl;
			}
			weightSum += diff;
			weightDiffs.insert(pair<double, unsigned int>(diff, w->second));
			//++w;
		}

	}

	if (weightSum == 0) {
		finishFit();
		outputError("Your population has converged. Quitting.");
	}

	if (options.verbosity >= 3){
		cout << "max: " << maxWeight << endl;
		cout << "weight sum: " << weightSum << endl;
		for (auto i = weightDiffs.begin(); i != weightDiffs.end(); ++i) {
			cout << "weight diff: " << i->first << " " << i->second << endl;
		}

		for (auto i = swarmBestFits_.begin(); i != swarmBestFits_.end(); ++i) {
			cout << "sbf: " << i->first << " " << i->second << endl;
		}

		for (auto i = allGenFits.begin(); i != allGenFits.end(); ++i) {
			cout << "agf: " << i->first << " " << i->second << endl;
		}
	}


	// If we want to keep any parents unchanged, send unchanged param sets to children.
	// We start with the global best fit params and iterate through the fit list from there.
	unsigned int childCounter = 1;
	int subParID = 0;
	int nModels = options.models.size();
	int pID = 0;
	// Only keep parents if we're doing an entire generation at once
	if (children.size() == (unsigned) options.swarmSize) {
		auto parent = allGenFits.begin();
		for (unsigned int p = 1; p <= options.keepParents; ++p) {

			vector<string> params;

			//Raquel, if we are using constraints, select parents with the best final rank (Fit and constraints combined)
			if(options.constraints_.size()>=1){
				//vector<pair<int,float> > subParRankFinal;
				int index = 0;
				int bestSubparID = subParRankFinal[index].first;
				index++;
				pID = fcalcParID(bestSubparID, nModels);
				vector<double> nparams = particleCurrParamSets_[pID];

				params.clear();
				for (auto param = nparams.begin(); param != nparams.end(); ++param) {
					particleNewParamSets[childCounter].push_back((*param));
					params.push_back(toString(*param));
				}

			}else{//Raquel, otherwise, use the original fit values
				split(parent->second, params);
				params.erase(params.begin());

				for (auto param = params.begin(); param != params.end(); ++param) {
					particleNewParamSets[childCounter].push_back(stod(*param));
				}

			}
			if(options.verbosity>=3){

				cout << "Sending unchanged params to " << children[childCounter - 1] << " counter: " << childCounter << endl;
			}
			//swarmComm->sendToSwarm(0, childCounter, SEND_FINAL_PARAMS_TO_PARTICLE, false, params);
			for (unsigned int mid = 0; mid < options.models.size(); ++mid){


				subParID = 1+mid+(children[childCounter - 1]-1)*nModels;
				if(options.verbosity>=3){

					cout << "Sending to SubPar: " << subParID << endl;
					cout << "parameters to send: " << endl;

					for(unsigned int i = 0; i<params.size(); i++){
						cout << "'" << params[i] << "'" << endl;

					}
				}
				swarmComm->sendToSwarm(0, subParID, SEND_FINAL_PARAMS_TO_PARTICLE, false, params);
				//sleep(1);


/*				while (numFinishedBreeding < (unsigned) subParID) {
					if(options.verbosity>=3){

						cout << "checking for DONEBREED" << endl;
					}
					unsigned int numMessages = swarmComm->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm->univMessageReceiver, true);

			//		int Pheromones::recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage, int messageID) {

					numFinishedBreeding += numMessages;

					//cout << numFinishedBreeding << endl;
					if(options.verbosity>=3){

						cout << "RAQUEL numFinishedBreeding " << numFinishedBreeding << endl;
					}
				}
*/
			}
			++parent;
			++childCounter;
			if(options.verbosity>=3){

				cout << "Raquel: (keep parents loop) Sent " << childCounter << " children to Swarm with SEND_FINAL_PARAMS_TO_PARTICLE" << endl;
			}
		}

	}


	//cout << "Raquel: (keep parents loop) DONE" << endl;

	float parentPairs;
	//cout << "children.size(): " << children.size() << endl;
	if (children.size() == (unsigned) options.swarmSize) {
		parentPairs = (float)parentPoolSize / 2;
	}
	else {
		parentPairs = (float)children.size() / 2;
	}

	//cout << "we have " << parentPairs << " parent pairs" << endl;

	//cout<<"Swarm::breedGenerationGA-AAA8: parentPairs:"<<parentPairs<<endl;

	boost::random::uniform_int_distribution<long int> unif(1, 100);
	for (unsigned int i = 0; i < parentPairs; ++i) {

		unsigned int p1;
		unsigned int p2;
		// Pick the fit values (particle parents) used in breeding
		//do {
		p1 = pickWeighted(weightSum, weightDiffs, options.extraWeight);
		p2 = pickWeighted(weightSum, weightDiffs, options.extraWeight);

		// Quit if we try too many times to select suitable parents
		//if (++maxFitCounter >= 10000) {
		//	outputError("Error: Tried too many times to select parents that didn't exceed the specified max_fit value of " + toString(static_cast<long double>(options.maxFit)) + ". Quitting.");
		//}
		// Make sure we don't exceed max_fit
		//} while (options.maxFit != 0 && particleBestFits_.at(p1) >= options.maxFit && particleBestFits_.at(p2) >= options.maxFit);

		// If we want different parents used in breeding, make sure that happens
		unsigned int retryCount = 0;
		while (p1 == p2 && options.forceDifferentParents) {
			retryCount++;
			if (retryCount > options.maxRetryDifferentParents) {
				if (options.verbosity >= 1) {
					cout << "Tried too many times to select different parents for breeding. Selecting the first two." << endl;
				}

				// Get reverse iterator to the weight map (best parents are at the end)
				multimap<double, unsigned int>::reverse_iterator w = weightDiffs.rbegin();

				// The weight map is sorted, so the first element will be the best fit
				p1 = w->second;
				if(options.verbosity>=3){

					cout << "1: " << w->second;
				}
				// Increment the map iterator until we find a fit value that isn't ours
				while (p1 == p2 && w != weightDiffs.rend()) {
					++w;
					p2 = w->second;

					if(options.verbosity>=3){

						cout << "tried to set p2 as " << w->second;
					}
				}

				break;
			}
			p2 = pickWeighted(weightSum, weightDiffs, options.extraWeight);
			if(options.verbosity>=3){

				cout << "retry2: " << p2 << endl;
			}
		}
		if(options.verbosity>=3){

			cout << "chose " << p1 << " and " << p2 << endl;
		}
		auto p1It = allGenFits.begin();
		advance(p1It, p1);
		vector<string> p1Vec;
		split(p1It->second, p1Vec);
		p1Vec.erase(p1Vec.begin());

		auto p2It = allGenFits.begin();
		advance(p2It, p2);
		vector<string> p2Vec;
		split(p2It->second, p2Vec);
		p2Vec.erase(p2Vec.begin());

		vector<string> c1Vec, c2Vec;
		unsigned int pi = 0;

		for (auto p = this->getFreeParams_().begin(); p != this->getFreeParams_().end(); ++p) {//for (auto p = options.model->getFreeParams_().begin(); p != options.model->getFreeParams_().end(); ++p) {
			string p1Param, p2Param;
			if(options.verbosity>=3){

				cout << "start working on free param["<<pi<<"]: " << p->first << endl;
			}
			p1Param = p1Vec[pi];
			p2Param = p2Vec[pi];
			if(options.verbosity>=3){

				cout << "p1: " << p1Param << endl;
				cout << "p2: " << p2Param << endl;
			}
			//Raquel: added
			if(p1Param.find_first_not_of("-01234567890.") != string::npos && p2Param.find_first_not_of("-01234567890.") != string::npos){
				if(options.verbosity>=3){
					cout << "no more parameters, done" << endl;
				}
				break;
			}

			if (unif(generalRand) < (options.swapRate * 100) ) {

				// TODO: Make sure individual mutation rates work
				if (hasMutate && p->second->isHasMutation()) {

					if(options.verbosity>=3){
						cout << "about to mutate" << endl;
					}
				//Raquel: fixed problem when one model has less parameters then de other, the stod function returned an error
				if(p1Param.find_first_not_of("-01234567890.") != string::npos){
					if(options.verbosity>=3){

						cout << "p1: " << p1Param << " is not a number, skiping..." << endl;
					}
				}else{
					//Raquel:  conversion from string to number worked
					p1Param = mutateParamGA(p->second, stod(p1Param));
					if(options.verbosity>=3){

						cout << p->second << " " << stod(p1Param) << endl;
					}
				}
				if(p2Param.find_first_not_of("-01234567890.") != string::npos){
					if(options.verbosity>=3){

						cout << "p2: " << p2Param << " is not a number, skiping..." << endl;
					}
				}else{
					//Raquel: conversion from string to number worked
					p2Param = mutateParamGA(p->second, stod(p2Param));

				}
				//Raquel: end fix

				//p1Param = mutateParamGA(p->second, stod(p1Param));
				//p2Param = mutateParamGA(p->second, stod(p2Param));

					if(options.verbosity>=3){

						cout << "mutated" << endl;
					}

				}
				if(options.verbosity>=3){

					cout << "swapping p1: " << p1Param << " with p2: " << p2Param << endl;
				}
					if(p1Param.find_first_not_of("-01234567890.") == string::npos && p2Param.find_first_not_of("01234567890.") == string::npos ){

						c1Vec.push_back(p2Param);
						particleNewParamSets[children[childCounter - 1]].push_back(stod(p2Param));
						if(options.verbosity>=3){

							cout << "c1 final: " << p2Param << endl;
						}
						c2Vec.push_back(p1Param);
						if (childCounter < children.size())
						{
							particleNewParamSets[children[childCounter]].push_back(stod(p1Param));
						}
					}else{cout << "Error p2param or p1Param are not numbers" << endl;}
			}
			else {
				if(p1Param.find_first_not_of("-01234567890.") == string::npos){

				c1Vec.push_back(p1Vec[pi]);
				particleNewParamSets[children[childCounter - 1]].push_back(stod(p1Vec[pi]));
				if(options.verbosity>=3){

					cout << "c1 final: " << p1Vec[pi] << endl;
				}
				}else{cout << "Error p2param or p1Param are not numbers" << endl;}
				if(p2Param.find_first_not_of("-01234567890.") == string::npos){

				c2Vec.push_back(p2Vec[pi]);
				if (childCounter < children.size())
				{
					particleNewParamSets[children[childCounter]].push_back(stod(p2Vec[pi]));
				}
				}else{cout << "Error p2param or p1Param are not numbers" << endl;}
				if(options.verbosity>=3){

					cout << "c1 final: " << p1Vec[pi] << "pushed back."<<endl;
				}
			}
			++pi;
			if(options.verbosity>=3){

				cout << "free param id:" << pi << "."<<endl;
			}
		}
		if(options.verbosity>=3){

			cout << "sending to " << children[childCounter - 1] << endl;
		}
	for (unsigned int mid = 0; mid < options.models.size(); ++mid){


		subParID = 1+mid+(children[childCounter - 1]-1)*nModels;
		if(options.verbosity>=3){

			cout << "Sending to SubPar: " << subParID << endl;

			for(unsigned int i = 0; i<c1Vec.size(); i++){
				cout << "'" << c1Vec[i] << "'" << endl;

			}
		}

		swarmComm->sendToSwarm(0, subParID, SEND_FINAL_PARAMS_TO_PARTICLE, false, c1Vec);
		//sleep(1);


/*		while (numFinishedBreeding < (unsigned) subParID) {
			if(options.verbosity>=3){

				cout << "checking for DONEBREED" << endl;
			}
			unsigned int numMessages = swarmComm->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm->univMessageReceiver, true);

	//		int Pheromones::recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage, int messageID) {

			numFinishedBreeding += numMessages;
			//cout << numFinishedBreeding << endl;
			if(options.verbosity>=3){

				cout << "RAQUEL numFinishedBreeding " << numFinishedBreeding << endl;
			}

		}
*/
	}

//		swarmComm->sendToSwarm(0, children[childCounter - 1], SEND_FINAL_PARAMS_TO_PARTICLE, false, c1Vec);
		++childCounter;
		if(options.verbosity>=3){

			cout << "Raquel: (mutate loop1)  Sent " << childCounter << " children to Swarm with SEND_FINAL_PARAMS_TO_PARTICLE" << endl;
		}
		// Make sure we don't send to too many parents (only relevant in last breeding with odd number of parents)
		if ( !( (fmod(parentPairs * 2, 2)) == 1 && parentPairs - i == 0.5 ) ) {
			if(options.verbosity>=3){

				cout << "sending to " << children[childCounter - 1] << endl;
			}
	for (unsigned int mid = 0; mid < options.models.size(); ++mid){


		subParID = 1+mid+(children[childCounter - 1]-1)*nModels;
		if(options.verbosity>=3){

			cout << "Sending to SubPar: " << subParID << endl;

			for(unsigned int i = 0; i<c2Vec.size(); i++){
				cout << "'" << c2Vec[i] << "'" << endl;

			}
		}
		swarmComm->sendToSwarm(0, subParID, SEND_FINAL_PARAMS_TO_PARTICLE, false, c2Vec);
		//sleep(1);

		if(options.verbosity>=3){

			cout << "waiting for " << subParID << endl;
		}
		//cout << "RAQUEL: children size " << children.size() << endl;
		//cout << "RAQUEL: swarm size" << options.swarmSize << endl;
		//cout << "RAQUEL: failed particles " << failedParticles_.size() << endl;
		//subParID = 1+mid+(children[childCounter - 1]-1)*nModels;

/*			while (numFinishedBreeding < (unsigned) subParID) {
				if(options.verbosity>=3){

					cout << "checking for DONEBREED" << endl;
				}
				unsigned int numMessages = swarmComm->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm->univMessageReceiver, true);

		//		int Pheromones::recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage, int messageID) {

				numFinishedBreeding += numMessages;
				//cout << numFinishedBreeding << endl;
				if(options.verbosity>=3){

					cout << "RAQUEL numFinishedBreeding " << numFinishedBreeding << endl;
				}
			}
*/
	}
	//		swarmComm->sendToSwarm(0, children[childCounter - 1], SEND_FINAL_PARAMS_TO_PARTICLE, false, c2Vec);
			++childCounter;
			if(options.verbosity>=3){

				cout << "Raquel: (mutate loop2) Sent " << childCounter-1 << " children to Swarm with SEND_FINAL_PARAMS_TO_PARTICLE" << endl;
			}


		}

		if(options.verbosity>=3){

			cout << "The loop completed for parentpairs:" << i+1<<"/"<<parentPairs<<"."<< endl;


			cout << "c1Vec:" << endl;
			for(unsigned int veci = 0; veci < c1Vec.size(); veci++){
				cout << c1Vec[veci] << endl;

			}

			cout << "c2Vec:" << endl;
			for(unsigned int veci = 0; veci < c2Vec.size(); veci++){
				cout << c2Vec[veci] << endl;

			}
		}
	}//	for (unsigned int i = 0; i < parentPairs; ++i) {


	//particleCurrParamSets_ = particleNewParamSets;
//cout<<"Swarm::breedGenerationGA-AAA9\n.";

//	for (auto it = particleNewParamSets.begin(); it != particleNewParamSets.end(); ++it) {
//		cout << "RAQUEL NEW PARAMS " << it->first << "sec " << it->second << endl;
//	}

	// Replace any changed param sets in the master set
	for (auto child = particleNewParamSets.begin(); child != particleNewParamSets.end(); ++child) {
		particleCurrParamSets_[child->first] = child->second;
	}

	while (numFinishedBreeding < (unsigned) (childCounter-1)) {
		if(options.verbosity>=3){

			cout << "checking for DONEBREED" << endl;
		}
		unsigned int numMessages = swarmComm->recvMessage(-1, 0, DONE_BREEDING, true, swarmComm->univMessageReceiver, true);

//		int Pheromones::recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage, int messageID) {

		numFinishedBreeding += numMessages;
		//cout << numFinishedBreeding << endl;
		if(options.verbosity>=3){

			cout << "RAQUEL numFinishedBreeding " << numFinishedBreeding << endl;
		}
	}

	if(options.verbosity>=3){

		cout << "done?" << endl;
	}

//cout<<"Swarm::breedGenerationGA-AAA10\n.";
}

void Swarm::runNelderMead(unsigned int it, unsigned int cpu) {
	if (options.verbosity >= 3) {
		cout << "Running Nelder-Meade local search for particle " << it << " on cpu " << cpu << endl;
	}
	// On the master side, we just need to construct the simplex, serialize it,
	// then send it to the CPU..

	// calc -> params
	map<double, vector<double> > simplex {pair<double, vector<double> >(particleBestFits_.at(it), normalizeParams(particleCurrParamSets_.at(it)))};

	// Get our nearest neighbors in parameter space
	map<double, unsigned int> neighbors = getNearestNeighbors(it, this->getNumFreeParams());//map<double, unsigned int> neighbors = getNearestNeighbors(it, options.model->getNumFreeParams());

	// Fill our simplex
	for (auto particle = neighbors.begin(); particle != neighbors.end(); ++particle) {
		cout << "chose neighbor: " << particle->second << " with fit of " <<  particleBestFits_.at(particle->second) << endl;
		simplex[particleBestFits_.at(particle->second)] = normalizeParams(particleCurrParamSets_.at(particle->second));
		for (auto param = simplex[particleBestFits_.at(particle->second)].begin(); param != simplex[particleBestFits_.at(particle->second)].end(); ++ param) {
			cout << "p: " << *param << endl;
		}
	}

	std::stringstream oss;
	boost::archive::text_oarchive oa(oss);
	oa << simplex;
	std::string serializedSimplex(oss.str());

	vector<string> message;
	message.push_back(serializedSimplex);

	runningParticles_.insert(cpu);
	swarmComm->sendToSwarm(0, cpu, BEGIN_NELDER_MEAD, false, message);
}

std::vector<unsigned int> Swarm::update_finished_running_particles(){

	//unsigned int nFinishedParticles = runningParticles_.size();
	unsigned int i , j, pID, subPID, maxSubPID, maxPID, nModels;
	std::set<unsigned int>::iterator it;
	std::vector<unsigned int> NewlyFinishedParticles;
	nModels=options.models.size(); NewlyFinishedParticles.clear(); //razi: init

	if (!runningSubParticles_.empty())
		for (it = runningSubParticles_.begin(); it != runningSubParticles_.end(); ++it){
			pID = fcalcParID(*it, nModels);
			if (runningParticles_.find(pID) == runningParticles_.end())  //razi: add if already does not exist
				runningParticles_.insert(pID);
		}

	if (!failedSubParticles_.empty())
		for (it = failedSubParticles_.begin(); it != failedSubParticles_.end(); ++it){
			pID = fcalcParID(*it, nModels);
			if (failedParticles_.find(pID) ==  failedParticles_.end())  //razi: add if already does not exist
				failedParticles_.insert(pID);
		}


	if (finishedSubParticles_.empty()){
		if(options.verbosity>=1){
			cout<<"finished subPartcile list is empty ...."<<endl;

		}
	}else{



		maxSubPID=0;
		for (std::set<unsigned int>::iterator maxIndex = finishedSubParticles_.begin(); maxIndex!=finishedSubParticles_.end(); ++maxIndex){
			if (options.verbosity >= 3) {
				cout<<*maxIndex<<", ";
			}
			if ((*maxIndex) > maxSubPID)
				maxSubPID = * maxIndex;
		}

		if (maxSubPID== (unsigned) -1){
			cout<<"finished subPartcile list does not include valid max id ...."<<endl;
			return NewlyFinishedParticles;
		}else{
			maxPID = fcalcParID(maxSubPID, nModels);
			if (options.verbosity >= 3) {
				cout<<"Updating finished particles. Max subParticle: " <<  maxSubPID<<"  Max Particle:" << maxPID <<endl;
			}

			bool allfinished;
			for (i = 1; i <= maxPID; i++){
				allfinished = true;
				for (j=0; j< nModels; j++){
					subPID = fcalcsubParID(i, j, nModels);
					if (finishedSubParticles_.find(subPID) == finishedSubParticles_.end())
						allfinished = false;
				}
				if (allfinished){
					if (options.verbosity >= 3) {
						cout<< "Particle: " << i<<" is completed. Adding to the list"<<endl;
					}
					if (finishedParticles_.find(i) == finishedParticles_.end()){
						finishedParticles_.insert(i);
						NewlyFinishedParticles.push_back(i);
					}
					else{
						if (options.verbosity >= 3) {
							cout<< "Particle: " << i<<" already included in the list of finished particles."<<endl;
						}
					}
				}

			}
		}
	}
	if (options.verbosity >= 3) {
		cout << "updated number of finished susb" << finishedSubParticles_.size();
		cout << "updated number of new finished subs" << NewlyFinishedParticles.size();
	}

	return NewlyFinishedParticles;
}


void Swarm::update_cur_particle_params(unsigned int pID, unsigned int mid, bool overwrite){ //razi: this need to be tested


	std::vector<unsigned int> mids;
	unsigned int i, mapfrom, mapto, nF;

	if(options.verbosity >=4){cout<<"update_cur_particle_params started."<<endl;} //mypause();

	if ((mid < (unsigned) -1) && (mid >= options.models.size()))
		outputError("updating current particle parameters, wrong model number" + toString(mid)+". Quitting....");

	if (mid == (unsigned) -1){
		for (i=0; i< options.models.size(); i++)
			mids.push_back(i);  //razi: do for all models if mid=-1
	}else
		mids.push_back(mid);

	for (i = 0; i <mids.size(); i++){
		mid=mids[i]; nF=options.models[mid]->getNumFreeParams();
		vector <double> InVec = subparticleCurrParamSets_[pID][mid];
		if (!InVec.empty()){
			vector <double> OutVec(nF);
			if ((!overwrite) && (particleCurrParamSets_[pID].size()>= nF)){
				for (i=0; i<nF; i++)
					OutVec[i] = particleCurrParamSets_[pID][i];
			}
			for (auto fp = freeParamMapping[mid].begin(); fp != freeParamMapping[mid].end(); fp++){
				mapfrom = fp->first;
				mapto = fp->second;
				if(options.verbosity >=4) cout<<"copying value for free parames from model:"<<mid<< " to the full list index1:"<< mapfrom <<"  index2:"<<mapto<<endl;
				if (overwrite || (OutVec[mapto]==0))
					OutVec[mapto] = InVec[mapfrom];
				particleCurrParamSets_[pID] = OutVec;
			}

		}
	}
	if(options.verbosity >=4){cout<<"update_cur_particle_params finished."<<endl;} //mypause();
}


//razi: later develop to calculate the particle fitting based on subparticle results
void Swarm::update_fitCalcs(){

	//Raquel
	if(options.verbosity >=5){
		cout<<"later develop: calculate the particle fitting based on subparticle results"<<endl; //mypause();
	}


}

string Swarm::order_params(string paramsString, unsigned int mid){
	if(options.verbosity >=5){
		cout<<"Not developed yet. Reorder the string containing the Free parameter values based on the order of free parameters in model:"<<mid<<endl<<endl;
	}
	return paramsString;
}

