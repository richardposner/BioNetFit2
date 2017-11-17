/*============================================================================
// Name        : Particle.cpp
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Lat Update: : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/



// TODO: Let's replace all the static_cast's with a custom function that converts ints/doubles to strings using stringstream

#include "Particle.hh"
#include "Setting.hh"

using namespace std;
using namespace std::chrono;

Particle::Particle(Swarm * swarm, int id) {
	id_ = id;
	swarm_ = swarm;
#ifdef VER2

	//razi: init vectors corresponding to multiple models
	unsigned int nModel = swarm_->getNumModels();
	if(nModel < 1) outputError("Can not create particle without a model.");
	models.clear();
	subParticles.clear();
	dataFiles_.clear();
	simParams_.clear();
	fitCalcs.clear();

	//razi: fill with dummay fileds
	std::map<std::string, std::map<int, Data*>> dummyFile;
	std::map<std::string, double> dummyParam;
	std::map<int, double> dummy_fit;
	for (unsigned int i=0; i< nModel; i++){
		subParticles.push_back(0);
		dataFiles_.push_back(dummyFile);
		models.push_back(0);
		simParams_.push_back(dummyParam);
		fitCalcs.push_back(dummy_fit);
	}


#else
	model_ = 0;
#endif
	objFuncPtr = 0;
	currentGeneration_ = 1;
	island_ = 1;
}

void Particle::setID(int id) {
	id_ = id;
}


// #1
double Particle::objFunc_sumOfSquares(double sim, double exp, double dummyvar) {
	return pow((abs(sim) - exp), 2);
}

// #2
double Particle::objFunc_chiSquare(double sim, double exp, double stdev) {
	return pow(((abs(sim) - exp) / stdev), 2);

	// TODO: Missing stdev?
}

// #3
double Particle::objFunc_divByMeasured(double sim, double exp, double dummyvar) {
	if(swarm_->options.verbosity >=4){
		cout << "sim: " << sim << endl;
		cout << "exp: " << exp << endl;
	}
	if(sim!=0){
		return pow(((abs(sim) - exp) / sim), 2);
	}else{
		return pow(((abs(sim) - exp)), 2); //Raquel solving division by zero problem

	}

}

// #4
double Particle::objFunc_divByMean(double sim, double exp, double mean) {
	return pow(((abs(sim) - exp) / mean), 2);
}

vector<double> Particle::getCentroid(vector<vector<double>> centroidVectors) {
	vector<double> centroid;

	for (unsigned int d = 0; d < centroidVectors[0].size(); ++d) {
		double sum = 0;
		for (unsigned int set = 0; set < centroidVectors.size(); ++set) {
			sum += centroidVectors[set][d];
		}
		centroid.push_back(sum / centroidVectors[0].size());
	}

	return centroid;
}

void Particle::generateParams() {
	// for a possibly better way to generate numbers
	// freeParams_ is a map with the first element being a parameter name, and the second element being a pointer to a FreeParam object
	//razi: in this new version we have vectors of these maps. Each vector element belongs to one of the parallel models indexed by mid

//cout<<"Particle::generateParams XXX-1"<<endl;
#ifdef VER2
//	for (unsigned int mid = 0; mid < models.size(); ++mid)
//	if (models.at(mid))
//	for (map<string,FreeParam*>::iterator i = models.at(mid)->freeParams_.begin(); i != models.at(mid)->freeParams_.end(); ++i) {
	if (models.at(0))
	for (map<string,FreeParam*>::iterator i = models.at(0)->freeParams_.begin(); i != models.at(0)->freeParams_.end(); ++i) {

#else
	for (map<string,FreeParam*>::iterator i = model_->freeParams_.begin(); i != model_->freeParams_.end(); ++i) {
#endif

		//vector<string> values;

		// Split our generation method/range into parts.
		//split(i->second, values);

		// Make sure our array contains the generation method
		/*if (values[0].empty()){
			string errMsg = "Problem parsing initial parameter generator for param: " + i->first + ".";
			outputError (errMsg);
		}*/

		// Store the free parameter name that we're currently generating
		string paramName;
		paramName = i->first;

		// Store the generation method
		//string genType = values[0];
		string genType = i->second->getGenerationMethod();

		// Generate a number on a linear scale
		if (genType == "random_var"){
			// Store our min and max values
			double min = i->second->getGenMin();
			double max = i->second->getGenMax();

			std::random_device rd;
			std::mt19937 mt(rd());
			std::uniform_real_distribution<double> dist(min, max);
			
			float myrand = dist(mt);

			pair<string,double> paramPair = make_pair(paramName, myrand);

//cout<<"Particle::generateParams XXX-2"<<endl;

#ifdef VER2
			//code for forcing the subparticles to have the same values
			for (unsigned int mid = 0; mid < models.size(); mid++){
				if (models.at(mid)){
					setParam(paramPair, mid);

				}
			}
#else
			setParam(paramPair);
#endif
//cout<<"Particle::generateParams XXX-3"<<endl;

		}
		else if (genType == "loguniform_var") {

			// Store our min and max values
			double min = i->second->getGenMin();
			double max = i->second->getGenMax();

			//cout << "generating with " << min << ":" << max << endl;

			// Generate a random double between 0 and rand_max
			double myrand = static_cast <double> (rand()) / static_cast<double> (RAND_MAX);
			//cout << id_ << " rand is " << myrand << endl;

			//
			double exp = log10(min) + myrand * (log10(max) - log10(min) );
			//cout << id_ << " base is " << exp << endl;
			//cout << id_ << " final " << pow(10, exp) << endl << endl;

			pair<string,double> paramPair = make_pair(paramName, pow(10, exp));

//cout<<"Particle::generateParams XXX-4"<<endl; mypause();
			//	unsigned int subParID = parParticle->calcsubParID(mid);

#ifdef VER2
			//code for forcing the subparticles to have the same values
			for (unsigned int mid = 0; mid < models.size(); mid++){
				if (models.at(mid)){
					setParam(paramPair, mid);
				}
			}
			//setParam(paramPair, mid);
#else
			setParam(paramPair);
#endif
//cout<<"Particle::generateParams XXX-5"<<endl; mypause();

		}
//cout<<"Particle::generateParams XXX-6"<<endl; mypause()
		//TODO: Add the rest of the init parm generation options
	}
//cout<<"Particle::generateParams XXX-7"<<endl; mypause();
}




#ifdef VER2
void Particle::runNelderMead(map<double, vector<double>> simplex, unsigned int mid) {

	std::map<std::string, double> simParams = simParams_.at(mid);
	auto curPar = getSubParticle(mid);

#else //VER2
void Particle::runNelderMead(map<double, vector<double>> simplex) {
	std::map<std::string, double> simParams = simParams_;
	auto curPar = this;
#endif //VER2

	// The transformation coefficients
	float reflection = 1.0;
	float expansion = 2.0;
	float contraction = 0.5;
	float shrink = 0.5;
	unsigned int simulationCount = 0;

	while (simulationCount < 20) {
		// Get our important vertices
		auto sIt = simplex.end();
		advance(sIt, - 1); // Last element
		cout << "worst is: " << sIt->first << endl;
		auto worst = sIt;
		advance(sIt, - 1); // Second to last element
		cout << "second worst is: " << sIt->first << endl;
		auto good = sIt;
		auto best = simplex.begin();
		cout << "best is: " << best->first << endl;
		bool invalid = false;

		vector<vector<double>> centroidVectors;
		for (auto params = simplex.begin(); params != worst; ++params) {
			centroidVectors.push_back(swarm_->normalizeParams(params->second));
		}

		// Calculate the centroid
		vector<double> centroid = getCentroid(centroidVectors);

		// Reflect
		cout << "reflecting" << endl;
		vector<double> R; // The transformation vector
		for (unsigned int d = 0; d < centroid.size(); ++d) {
			double r = centroid[d] + reflection * (centroid[d] - worst->second[d]);
			if (r <= 0 || r > 1) {
				invalid = true;
			}
			cout << "creating new pt with eq " << centroid[d] << " + " << reflection << " * (" << centroid[d] << " - " << worst->second[d] << "): " << r << endl;
			R.push_back(r);
		}

		vector<double> deNormalizedReflection = swarm_->deNormalizeParams(R);

		auto tIt = deNormalizedReflection.begin();
		for (auto p = simParams.begin(); p != simParams.end(); ++p) {
			p->second = *tIt;
			++tIt;
		}

		double rCalc;
		if (!invalid) {
			for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
				curPar->runModel(i, true);
				++simulationCount;
			}
			if (swarm_->options.smoothing > 1) {
				curPar->smoothRuns();

			}
			curPar->calculateFit(true);
#ifdef VER2
			rCalc = fitCalcs.at(mid)[-1];
#else
			rCalc = fitCalcs[-1];
#endif
		}
		else {
			invalid = false;
			rCalc = pow(10, 10);
		}
		cout << "reflection: " << rCalc << endl;

		// R Better than good, but worse than best
		if (rCalc > best->first && rCalc < good->first) {
			if (swarm_->options.verbosity >= 3) {
				cout << "reflection was worse than best and better than good. erasing worst of " << worst->first << endl;
			}
			simplex.erase(worst);
			if (swarm_->options.verbosity >= 3) {
				cout << "sim size: " << simplex.size() << endl;
			}
			simplex.insert(pair<double, vector<double>> (rCalc, R));
			if (swarm_->options.verbosity >= 3) {
				cout << "inserting rCalc of " << rCalc << endl;
			}
			continue;
		}
		// R Better than best
		else if (rCalc < best->first) {
			if (swarm_->options.verbosity >= 3) {
				cout << "reflection was better than best. expanding" << endl;
			}
			// Expand

			vector<double> E;
			for (unsigned int d = 0; d < centroid.size(); ++d) {
				double e = centroid[d] + expansion * (R[d] - centroid[d]);
				if (e <= 0 || e > 1) {
					invalid = true;
				}
				E.push_back(e);
				if (swarm_->options.verbosity >= 3) {
				cout << "creating new pt with eq " << centroid[d] << " + " << expansion << " * (" << R[d] << " - " << centroid[d] << "): " << e << endl;
				}
			}

			vector<double> deNormalizedExpansion = swarm_->deNormalizeParams(E);

			auto tIt = deNormalizedExpansion.begin();
			for (auto p = simParams.begin(); p != simParams.end(); ++p) {
				p->second = *tIt;
				++tIt;
			}

			double eCalc;
			if (!invalid) {
				for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
					curPar->runModel(i, true);
					++simulationCount;
				}
				if (swarm_->options.smoothing > 1) {
					curPar->smoothRuns();
				}
				curPar->calculateFit(true);

#ifdef VER2
				eCalc = fitCalcs.at(mid)[-1];
#else
				eCalc = fitCalcs[-1];
#endif
				cout << "expansion: " << eCalc << endl;
			}
			else {
				eCalc = pow(10,10);
				invalid = false;
			}

			if (eCalc < rCalc) {
				cout << "expansion was better than reflection. looping." << endl;
				simplex.erase(worst);
				simplex.insert(pair<double, vector<double>> (eCalc, E));
				continue;
			}
			else {
				cout << "expansion was worse than reflection. looping." << endl;
				simplex.erase(worst);
				simplex.insert(pair<double, vector<double>> (rCalc, R));
				continue;
			}
		}
		// R worse than good
		else if (rCalc > good->first) {
			cout << "reflection was worse than good. contracting." << endl;
			// Contraction
			vector<double> C;
			for (unsigned int d = 0; d < centroid.size(); ++d) {
				double c = centroid[d] + contraction * (worst->second[d] - centroid[d]);
				if (c <= 0 || c > 1) {
					invalid = true;
				}
				C.push_back(c);
				cout << "creating new pt with eq " << centroid[d] << " + " << contraction << " * (" << worst->second[d] << " - " << centroid[d] << "): " << c << endl;
			}

			vector<double> deNormalizedContraction = swarm_->deNormalizeParams(C);

			auto tIt = deNormalizedContraction.begin();
			for (auto p = simParams.begin(); p != simParams.end(); ++p) {
				p->second = *tIt;
				++tIt;
			}

			double cCalc;
			if (!invalid) {
				for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
/*
#ifdef VER2
					getSubParticle(mid)->runModel(i, true);
#else
					runModel(i, true);
#endif
*/
					curPar->runModel(i, true);
					++simulationCount;
				}
				if (swarm_->options.smoothing > 1) {
#ifdef VER2
					smoothRuns(mid); //razi: run for subparticle associated with model mid
#else
					smoothRuns();
#endif
					}
#ifdef VER2
				calculateFit(true, mid); //razi: run for subparticle associated with model mid
				cCalc = fitCalcs.at(mid)[-1];
#else
				calculateFit(true);
				cCalc = fitCalcs[-1];
#endif
				cout << "contraction: " << cCalc << endl;
			}
			else {
				invalid = false;
				cCalc = pow(10, 10);
			}

			if (cCalc < worst->first) {
				cout << "contraction was better than worst. looping." << endl;
				simplex.erase(worst);
				simplex.insert(pair<double, vector<double>> (cCalc, C));
				continue;
			}
			else {
				//cout << "contraction was worse than worst. shrinking" << endl;
				// Shrink
				for (auto sIt = ++simplex.begin(); sIt != simplex.end(); ++sIt) {
					for (unsigned int d = 0; d < centroid.size(); ++d) {
						double s = simplex.begin()->second[d] + (shrink * (sIt->second[d] - best->second[d]));
						sIt->second[d] = s;
					}
				}
			}
		}
	}

#ifdef VER2
	fitCalcs.at(mid)[-1] = simplex.begin()->first;
#else
	fitCalcs[-1] = simplex.begin()->first;
#endif
	unsigned int d = 0;
	for (auto p = simParams.begin(); p != simParams.end(); ++p) {
		p->second = simplex.begin()->second[d++];
	}
}







/*
bool Particle::checkNelderMeadTerminationCriteria(vector<vector<double>> simplex, unsigned int numEvaluations, unsigned int numIterations) {
	if (swarm_->options.maxLocalEvaluations && numEvaluations >= swarm_->options.maxLocalEvaluations) {
		return true;
	}
	else if (swarm_->options.maxLocalIterations && numEvaluations >= swarm_->options.maxLocalIterations) {
		return true;
	}
	else {

	}

	return false;
}
 */

























































#ifdef VER2
void Particle::doParticle(unsigned int mid) {
	int subParID = calcsubParID(mid);
	int id = subParID;

	if (swarm_->options.verbosity >= 3) {cout << "In doParticle(), waiting for message from master. subParticle:["<<subParID <<":"<< id_  <<"-"<<mid <<"]  model id:"<< mid<< endl; mypause();}

// razi: this is to run the program in slave mode without using master-slave configuration
// this should be commented later on.
#ifdef STANDALONE_STANDALONE_PARTICLE_TEST_ACTIVE

	if (swarm_->options.verbosity >= 3){ cout << "doParticle. Generate random parameters instead of waiting for master to send the parameters. Results are not valid !!!" << endl; mypause();}
	currentGeneration_ =1; int rnd; srand(subParID); //razi random seed

	//razi generate random values for free parameters only for test purpose
	for (auto frp = swarm_->options.models.at(mid)->freeParams_.begin(); frp != swarm_->options.models.at(mid)->freeParams_.end(); ++frp) {
		rnd=rand()%100; simParams_.at(mid).insert(pair<string, double>(frp->first, rnd));
		//cout<<"virtual free param:"<<frp->first<<" is set to "<<rnd<<"."<<endl;
	}

#else //STANDALONE_STANDALONE_PARTICLE_TEST_ACTIVE
	if (swarm_->resumingSavedSwarm) {
		if (swarm_->options.verbosity >= 3)  cout << "Waiting for flight count from master" << endl;
		swarm_->swarmComm->recvMessage(0, id, SEND_NUMFLIGHTS_TO_PARTICLE, true, swarm_->swarmComm->univMessageReceiver);
		currentGeneration_ = stoi(*(swarm_->swarmComm->univMessageReceiver.find(SEND_NUMFLIGHTS_TO_PARTICLE)->second.message.begin()));

		swarm_->swarmComm->univMessageReceiver.clear();
		if (swarm_->options.verbosity >= 3) cout << "Received flight count from master. Waiting for param set from master...." << endl;
		swarm_->swarmComm->recvMessage(0, id, SEND_FINAL_PARAMS_TO_PARTICLE, true, swarm_->swarmComm->univMessageReceiver);

		auto freeParam = swarm_->options.models.at(mid)->getFreeParams_().begin();
		for (auto paramVal = swarm_->swarmComm->univMessageReceiver.find(SEND_FINAL_PARAMS_TO_PARTICLE)->second.message.begin(); paramVal != swarm_->swarmComm->univMessageReceiver.find(SEND_FINAL_PARAMS_TO_PARTICLE)->second.message.end(); ++paramVal) {
			if (swarm_->options.verbosity >= 3) cout << "updating " << freeParam->first << " to " << *paramVal << endl;
			simParams_.at(mid).insert(pair<string, double>(freeParam->first, stod(*paramVal)));
			++freeParam;
		}
		swarm_->swarmComm->univMessageReceiver.clear();
	}
#endif //STANDALONE_STANDALONE_PARTICLE_TEST_ACTIVE


	if (swarm_->options.fitType == "de") {
		// Need to figure out which island I'm on
		for (unsigned int island = 1; island <= swarm_->options.numIslands; ++island) {
			for (unsigned int particle = 1; particle <= (swarm_->options.swarmSize / swarm_->options.numIslands) * island; ++particle) {
				if (particle == id_) {
					island_ = island;
					goto theEnd;
				}
			}
		}
		theEnd:;
	}


#ifndef STANDALONE_PARTICLE_TEST_ACTIVE

//cout << "Particle " << id_ << " waiting to begin" << endl;
//	swarm_->swarmComm->recvMessage(0, id_, NEXT_GENERATION, true, swarm_->swarmComm->univMessageReceiver);
	if (swarm_->options.verbosity >= 3) {
		cout << "subParticle " << id << " waiting to begin" << endl;
	}
	swarm_->swarmComm->recvMessage(0, id, NEXT_GENERATION, true, swarm_->swarmComm->univMessageReceiver); //Master should send message with subParID
	swarm_->swarmComm->univMessageReceiver.clear();
#endif //STANDALONE_PARTICLE_TEST_ACTIVE
	if (swarm_->options.verbosity >= 3) {
		cout << "Particle " << id_ << " starting" << endl;
	}
	if (swarm_->options.verbosity >= 3) {
		cout << "In doParticle(), entering main run loop" << endl;
	}

	bool doContinue = true;
	bool doRunModel = true;
	while(doContinue) {
		if (swarm_->bootstrapCounter > 0 && currentGeneration_ == 1) {
			swarm_->swarmComm->recvMessage(0, id, NEXT_GENERATION, true, swarm_->swarmComm->univMessageReceiver);
			swarm_->swarmComm->univMessageReceiver.clear();
		}
		if (swarm_->options.verbosity >= 3) {
			cout << "DoParticle: Generation is " << currentGeneration_ << endl;
		}
		if (doRunModel) {
			for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
				if (swarm_->options.verbosity >= 3) cout<<"doParticle-1x try to run model for i:"<<i<<endl;

				if (swarm_->isMaster){
						for (unsigned int j=0; j<=swarm_->getNumModels(); j++){   //razi:run all subparticles
						if (getSubParticle(j)){
							getSubParticle(j)->runModel(i, false);   //razi: check later due to false
							if (swarm_->options.verbosity >= 3) {
								cout<<"Particle::doParticle. Running subParticle:["<<calcsubParID(j) <<":"<< id_  <<"-"<<j <<"] for iteration:"<<i << " is completed...\n";
							}
						}
					}
				}else{
					if (getSubParticle(mid)){
						if (swarm_->options.verbosity >= 3) {
							cout<<"Particle::doParticle. Running subParticle:["<<subParID <<":"<< id_  <<"-"<<mid <<"] for iteration:"<<i << " is started...\n";
						}
						getSubParticle(mid)->runModel(i, false);   //razi: check later due to false
						if (swarm_->options.verbosity >= 3) {
							cout<<"Particle::doParticle. Running subParticle:"<< mid <<" for iteration:"<<i << " is completed...\n";
						}
					}else outputError("SubParticle number error subParticle:["+ toString(subParID)+":"+ toString(id_) + "-:"+ toString(mid)+"] particleID:"+ toString(id_)+ " model id:" + toString(mid) + "  subParID:" + toString(subParID));
				}
				if (swarm_->options.verbosity >= 3) cout<<"doParticle-1 ran model for i:"<<i<<endl;
			}

			if (swarm_->isMaster){
				for (unsigned int j=0; j<=swarm_->getNumModels(); j++){
					if (getSubParticle(j))
						if (swarm_->options.smoothing > 1)
							smoothRuns(j);
						finalizeSim(j);
				}
			}else{
				if (getSubParticle(mid))
					smoothRuns(mid);
				finalizeSim(mid);
			}
		}
		else{
			//cout<<"do Runmodel was not active\n";
			doRunModel = true;
		}

if (swarm_->options.verbosity >= 3) cout<<"doParticle-4 now lets check messages..."<<endl;
		if (swarm_->options.verbosity >= 5) {
			cout << "RAQUEL before if" << endl;
		}
		if (swarm_->options.fitType == "ga") {
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL Entering check messages GA" << endl;
			}
			//checkMessagesGenetic(); //Raquel: changed this function to fix the problem of not proceeding to the next generation
			checkMessagesGenetic(mid);
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL done check messages" << endl;
			}
		}
		else if (swarm_->options.fitType == "pso") {
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL Entering check messages PSO" << endl;
			}
			checkMessagesPSO(mid);
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL Done check messages PSO" << endl;
			}

		}
		else if (swarm_->options.fitType == "de") {

			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL Entering check messages DE" << endl;
			}

			checkMessagesDE(mid);

			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL Done check messages DE" << endl;
			}

		}
		else if (swarm_->options.fitType == "sa") {
			doRunModel = checkMessagesDE();
		}
if (swarm_->options.verbosity >= 3) cout<<"doParticle-4 check messages done..."<<endl;

		if (doRunModel == true) {
			// Next generation
			++currentGeneration_;
		}
	}
	if (swarm_->options.verbosity >= 3) {
		cout << "In doParticle(), finished ...." << endl;
	}
}





void Particle::setModel(Model * model, unsigned int mid) {
	//cout<<"press a key to set model. Model size:"<<models.size()<<" mid:"<<mid;	//mypause();

	/* razi: moved to constructor
	int gap = mid+1-models.size();
	if(gap>0)
		for (int i=0; i<gap; i++)
			models.push_back(0);
	*/
	if (mid>=models.size())  outputError("setModel failed. Invalid model id:"+toString(mid)+". I am quitting....");
	models.at(mid)=model;   //the list is not long enough
	//cout<<"Model is set. press a key  Model size:"<<models.size()<<endl; //mypause();
}

void Particle::setParam(pair<std::string, double> myParams, unsigned int mid){

	/* razi: moved to constructor
	unsigned int gap= mid+1-simParams_.size();
	if (gap>0){   //fill dummy
		std::map<std::string, double> dummyParam;		//dummyParam[""]=0;
		for (int i=0; i < gap; i++)
			simParams_.push_back(dummyParam);//simParams_.push_back(std::map<std::string, double> dummyParam);
	}*/

	if (mid>=simParams_.size())  outputError("setModel failed. Invalid model id:"+toString(mid)+". I am quitting....");
		//Raquel: this code check with the master whether the parameter exists in model
	for (auto p = swarm_->options.models.at(mid)->getFreeParams_().begin(); p != swarm_->options.models.at(mid)->getFreeParams_().end(); ++p) {
		if(p->first == myParams.first) {
				//cout << "Raquel myParams first: " << myParams.first << " second: " << myParams.second << endl;
				simParams_.at(mid).insert(myParams);//Raquel: if it does exist, set the param value
		}//Raquel: otherwise, do nothing
	}

}

void subParticle::runModel(unsigned int iteration, bool localSearch) {

	unsigned int pID = parParticle->id_;
	if (parParticle->swarm_->options.verbosity >= 3) {
		cout << "RRR parParticle->id_: " << pID << endl;
	}
	Swarm * swp = parParticle->swarm_;  if(!swp){cout<<"subParticle::runModel: Invalid subParticle with no parent..."; return;}
	Model * mdl = parParticle->models.at(this->mid_);  if(!swp){cout<<"subParticle::runModel: Invalid subParticle with no model..."; return;}
	unsigned int mid = this->mid_;
	unsigned int subParID = parParticle->calcsubParID(mid);

	if ((swp->getNumModels()<1) || (mid<0) || (mid >=swp->getNumModels())){
		cout<<"subParticle::runModel. Number of models:" << swp->getNumModels() <<"subParticle id:" << mid<< endl;
		outputError("Something must be wrong. Particle can not be run. Quitting...");
	}


	// First get our path and filename variables set up for use in model generation, sim command, etc
	string bnglFilename = swp->getModelName(this->mid_,false)+"_"+toString(pID) + "_" + toString(iteration) + ".bngl";

	if (parParticle->swarm_->options.verbosity>= 3) {cout<<"Started running model for subParticle:["<<subParID <<":"<< pID  <<"-"<<mid_ <<"] model file is: " << bnglFilename<<".\n"; }

	string path = swp->options.jobOutputDir + toString(parParticle->currentGeneration_) + "/";
	string bnglFullPath = path + bnglFilename;

	string suffix = toString(pID) + "_" + toString(iteration);

	string pipePath;
	// Only need to generate files if we're in the first generation. In subsequent generations
	// the model generation is handled by breeding parents
	if (parParticle->currentGeneration_ == 1 && !localSearch) {
		if (mdl->getHasGenerateNetwork()){
			if (parParticle->swarm_->options.verbosity>= 3) cout<<"generating network...\n";
			if (swp->bootstrapCounter == 0) {
				//string netFilename = "base_"+toString(mid)+".net"; //Raquel: commented to fix problem with .net file name
				//string netFullPath = swp->options.jobOutputDir + netFilename; //Raquel: commented to fix problem with .net file name
				string path2 = swp->options.models.at(mid)->getName(); //Raquel: fixed problem with .net file name
				string basename1 = getFilename(path2); //Raquel: fixed problem with .net file name
				string netFilename = basename1 + "_base.net"; //Raquel: fixed problem with .net file name
				string netFullPath = swp->options.jobOutputDir + netFilename; //Raquel: fixed problem with .net file name
				if ((swp->options.jobOutputDir.substr(swp->options.jobOutputDir.size()-1,1)!="/") && (swp->options.jobOutputDir.substr(swp->options.jobOutputDir.size()-1,1)!="\\"))
					netFullPath = swp->options.jobOutputDir + "/" + netFilename;
				if (iteration == 1) {
					mdl->parseNet(netFullPath);
				}
			}
			//razi: check later
			mdl->outputModelWithParams(parParticle->simParams_.at(mid), path, bnglFilename, suffix, false, false, true, false, false);
		}
		else {
			//razi: check later
			mdl->outputModelWithParams(parParticle->simParams_.at(mid), path, bnglFilename, suffix, false, false, false, false, false);
		}
	}
	else if (localSearch) {
		cout << "regenerating models.." << endl;
		// And generate our models
		if (mdl->getHasGenerateNetwork()){
			// If we're using ODE solver, output .net file
			bnglFilename = boost::regex_replace(bnglFilename, boost::regex("bngl$"), string("net"));
			mdl->outputModelWithParams(parParticle->simParams_.at(mid), path, bnglFilename, suffix, false, false, false, false, true);
		}
		else {
			// If we're using network free simulation, output .bngl
			mdl->outputModelWithParams(parParticle->simParams_.at(mid), path, bnglFilename, suffix, false, false, false, false, false);
		}
	}

	ifstream bngfile(bnglFullPath);

	if(!bngfile){
		swp->processLateParticles(subParID, false, parParticle->currentGeneration_);
		swp->fixRunningParticle(subParID);


	}
	if (swp->options.verbosity>= 3) {
		cout << "RAQUEL GENERATION: " << parParticle->currentGeneration_ << " path" << path << endl;
	}
	// Generate .gdat pipes if we're using pipes

#ifdef PC_VER //pipe are defined only in linux operating systems
	if (swp->options.usePipes)
		outputError("Using pipes are not allowed in windows, disable use_pipes option in the config [.conf] file!!!");
#endif

	string outputSuffix;
	if (swp->options.usePipes) {
		for (std::map<std::string,Model::action>::iterator i = mdl->actions.begin(); i != mdl->actions.end(); ++i) {
			if (i->second.scanParam.size() > 0) {
				outputSuffix = ".scan";
			}
			else {
				outputSuffix = ".gdat";
			}

			//raquel: added line below to fix problem with real simulator file name
			//pipePath = swp->getModelName(this->mid_,false)+"_"+toString(pID) + "_" + toString(iteration) + ".gdat";
			//raquel commented below
			pipePath = path + i->first+ toString(pID) + "_" + toString(iteration) + outputSuffix; //new format to be consistent with bngl file
			if (parParticle->swarm_->options.verbosity >= 3) {
				cout << "Raquel: pipepath " << pipePath << endl;
				cout << "subParticle::runModel pp is: " << pipePath << endl;
			}
			createParticlePipe(pipePath.c_str());
		}
	}




	int ret = -1;
	if (swp->options.verbosity>=3) cout<<"subParticle::runModel: bngl command is:" << swp->options.bngCommand<<endl;
	//razi: generate random gdat if the simulation is set to random, only for test
	std::size_t found = swp->options.bngCommand.find("random");
	if (found!=std::string::npos){   //run artificial simulation
#ifndef TEST_SIMULATOR
		outputError("You are using the random simulator, but the macro is not activated. Recompile the code with TEST_SIMULATOR on !!!");
#endif
		if (parParticle->swarm_->options.verbosity >= 3) {
			cout<<"\nsubParticle::runModel: Random Simulation is chosen for test, generating random result file !!!\n\n"; //mypause();
		}
//swp->printDetails();
//if (mid==1) {return;}
		string exp_file = swp->getExpPath(mid, 0);
//if (mid==1) return; //razi for test uncomment
		if (parParticle->swarm_->options.verbosity >= 3) {
			cout<<"\nsubParticle::runModel: Tries to generate file from file :" << exp_file<<"\n."; //mypause();
		}
		string gdat_path = path +  swp->getExp(mid_,0)+ "_" + toString(pID) + "_" + toString(iteration) + ".gdat";
		if (parParticle->swarm_->options.verbosity >= 3) {
			cout<<"\nsubParticle::runModel: Finished generating file: "<< gdat_path << "   from file :" << exp_file<<"\n.";// mypause();
		}
		ret = generate_gdat_file(exp_file, gdat_path, pID);
		if (parParticle->swarm_->options.verbosity >= 3) {
			cout << "RAQUEL ret: " << ret << endl;
		}
	}else
	{
//if (mid==1) return; //razi for test uncomment
		// Construct our simulation command
		string command = swp->options.bngCommand + " --outdir " + path + " " + bnglFullPath + " >> " + path + toString(pID) + ".BNG_OUT 2>&1";
		if (swp->options.usePipes) {
			command += " &";
		}

		if (!checkIfFileExists(path)) {
			runCommand("mkdir " + path); //Raquel added this block
		}

		if (swp->options.verbosity >= 3) {
			cout << "subParticle::runModel: Running model with command: " << command << endl;
		}

		// Run simulation command
		//int ret = system(command.c_str());
		ret = system(command.c_str()); //Raquel: uncomented to run real Simulators, and removed the redeclaration of ret
	}


//if (mid==1) return; //razi for test uncomment

	// Check for simulation command success
	//Raquel: changed this temporarily
	if (ret == 0 || ret==-1 ) { // TODO: Need to check for simulation status when using pipes. Going by return code doesn't work there because we're using the & operator

		// really not sure how we can do this easily
		if (parParticle->swarm_->options.verbosity >= 3) {
			cout<<"subParticle::runModel: simulation completed successfully. Lets collect data ...\n"; //mypause();
		}
		string outputSuffix;
		// Save our simulation outputs to data objects
		for (std::map<std::string,Model::action>::iterator action = mdl->actions.begin(); action != mdl->actions.end(); ++action) {
			if (action->second.scanParam.size() > 0) {
				outputSuffix = ".scan";
				//cout<<"Particle::runModel: A1\n";
			}
			else {
				//cout<<"Particle::runModel: A2\n";
				outputSuffix = ".gdat";
			}
//cout<<"Particle::runModel: A3\n"; //mypause();
			//razi: this is the gdatfile including the simulation results for [ParID,subParID,iteration]
			string dataPath = path + action->first + "_" + toString(pID) + "_" + toString(iteration) + outputSuffix;  //was string dataPath = path + action->first + "_" + toString(id_) + "_" + toString(iteration) + outputSuffix;
			//cout<<"Particle::runModel: A3 dataPath :"<<dataPath<<endl; //mypause();
#ifdef VER2
			if (swp->options.verbosity>=3) cout<<"subParticle::runModel: Adding the results as a data object.  action:"<< action->first<< " datafile:"<<dataPath<<"   model id:"<< mid << "  iteration:" << iteration << endl;
			parParticle->dataFiles_.at(mid)[action->first].insert(pair<int, Data*>(iteration, new Data(dataPath, swp, false, mid)));

			if (parParticle->swarm_->options.verbosity >= 3) {
				cout << "RAQUEL COUNT AT MID " << mid << " = " << parParticle->dataFiles_.at(mid).count(action->first) << endl;
			}

			if (swp->options.verbosity>=3) cout<<"subParticle::runModel: Lets check if done properly. File:"<< parParticle->dataFiles_.at(mid)[action->first][iteration]->getPath() << "exists.\n";
#else
			parParticle->dataFiles_[action->first].insert(pair<int, Data*>(iteration, new Data(dataPath, swp, false)));
#endif

			//cout<<"Particle::runModel: A4\n";
		}
	}
	else {
		// If our return code is not 0, tell the master that the simulation failed
		if (parParticle->swarm_->options.verbosity >= 3) {
			cout << "RAQUEL ret " << ret << endl;

			cout << "subParticle::runModel: I failed completing runmodel ..." << endl; //mypause();
		}
		swp->swarmComm->sendToSwarm(int(subParID), 0, SIMULATION_FAIL, true, swp->swarmComm->univMessageSender);
	}
	if (parParticle->swarm_->options.verbosity >= 3) {

		cout<<"subParticle::runModel: runmodel finished ...\n";
	}
}


void Particle::checkMessagesGenetic(unsigned int mid) {
unsigned int subParID = calcsubParID(mid);
if ((mid<0)|| (mid >= swarm_->getNumModels()) || (!models.at(mid)))
    outputError("checkMessagesGenetic Failed. Invalid SubParticle.  PID:"+ toString(id_) + "  model id:" + toString(mid)+ "  subParID:" + toString(subParID)+". Quitting !!!" );
//if ((mid<0)|| (mid >= swarm_->getNumModels()) || (!getSubParticle(mid)) || (!models.at(mid)))
//outputError("checkMessagesGenetic Failed. Invalid SubParticle.  PID:"+ toString(id_) + "  model id:" + toString(mid)+ "  subParID:" + toString(subParID)+". Quitting !!!");

Model * mdl= models.at(mid);
	// Holds iterator ranges when finding items in the message holder
	pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;
	if (swarm_->options.verbosity >= 3) {
		cout << "RAQUEL entering while from function that takes mid as input" << endl;
	}
	string path; //Raquel: declaring it outside the while loop

	while (1) {
		// Retrieve any messages
		int numCheckedMessages = 0;
		int numMessages = swarm_->swarmComm->recvMessage(-1, subParID, -1, true, swarm_->swarmComm->univMessageReceiver, true);
		if (swarm_->options.verbosity >= 4) {
			cout << "RAQUEL inside LOOP ITERATION in checkMessagesGenetic" << endl;

			cout << "RAQUEL entering univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE)" << endl;
		}
		smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
		if (swarm_->options.verbosity >= 3) {
			cout << "RAQUEL passed univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE)" << endl;
		}

		if (smhRange.first != smhRange.second) {
			//cout << id_ << " found final " << endl;
			for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
				//Timer tmr;

				int messageIndex = 0;
				for (auto p = simParams_.at(mid).begin(); p != simParams_.at(mid).end(); ++p) {
					if (swarm_->options.verbosity >= 3) {
						cout << id_ << " updating parameter " << p->first << " to " << sm->second.message[messageIndex] << endl;
					}
					p->second = stod(sm->second.message[messageIndex]);
					if (swarm_->options.verbosity >= 3) {
						cout << "updated param " << p->first << ": " << p->second << endl;
					}
					++messageIndex;
				}

				path = swarm_->options.jobOutputDir + toString(currentGeneration_ + 1) + "/";

				if (!checkIfFileExists(path)) {
					if (swarm_->options.verbosity >= 3) {
						cout << "trying to create path " << path << endl;
					}
					runCommand("mkdir " + path);
				}

				for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
					// Construct our filenames
					string bnglFilename = swarm_->getModelName(mid, false)+"_" + toString(id_) + "_" + toString(i) + ".bngl";
					string bnglFullPath = path + bnglFilename;
					string suffix = toString(id_) + "_" + toString(i);
					if (swarm_->options.verbosity >= 3) {
						cout << "RAQUEL CREATING BNGL: " << bnglFullPath << endl;
					}
					// And generate our models
					if (swarm_->options.models.at(mid)->getHasGenerateNetwork()){
						// If we're using ODE solver, output .net and .bngl
						mdl->outputModelWithParams(simParams_.at(mid), path, bnglFilename, suffix, false, false, true, false, false);
					}
					else {
						// If we're using network free simulation, output .bngl
						mdl->outputModelWithParams(simParams_.at(mid), path, bnglFilename, suffix, false, false, false, false, false);
					}
					if (swarm_->options.verbosity >= 3) {
						cout << "RAQUEL BNGL generated: " <<  bnglFullPath << " id_ " << toString(id_) << " i " << toString(i) << " subParID " << subParID << endl;
					}
				}
				if (swarm_->options.verbosity >= 3) {
				// Tell the master we have our new params and are ready for the next generation
					cout << id_ << " telling master we're finished " << endl;
				}
				swarm_->swarmComm->sendToSwarm(subParID, 0, DONE_BREEDING, false, swarm_->swarmComm->univMessageSender);

				//double t = tmr.elapsed();
				//cout << "SEND_FINAL_PARAMS took " << t << " seconds" << endl;
				++numCheckedMessages;
			}
		}

		if (swarm_->options.verbosity >= 3) {
			cout << "RAQUEL BEFORE IF numCheckedMessages >= numMessages" << endl;
		}

		if (numCheckedMessages >= numMessages) {
			swarm_->swarmComm->univMessageReceiver.clear();
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL INSIDE IF numCheckedMessages >= numMessages = " << numCheckedMessages << ":" << numMessages << endl;
			}
			continue;
			//return 0; //Raquel: added
		}
		if (swarm_->options.verbosity >= 3) {
			cout << "RAQUEL PASSED IF numCheckedMessages >= numMessages" << endl;
		}

		if (swarm_->swarmComm->univMessageReceiver.find(FIT_FINISHED) != swarm_->swarmComm->univMessageReceiver.end()) {
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL entering ~pheromones" << endl;
			}
			swarm_->swarmComm->~Pheromones();
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL exiting ~pheromones" << endl;
			}
			//return;
			exit(0);
		}

		if (swarm_->swarmComm->univMessageReceiver.find(NEW_BOOTSTRAP) != swarm_->swarmComm->univMessageReceiver.end()) {
			if (swarm_->options.verbosity >= 3) {
				cout << "SubParticle " << subParID <<": Starting a new bootstrapping run" << endl;
			}

			++swarm_->bootstrapCounter;
			currentGeneration_ = 0;
			simParams_.at(mid).clear();
			generateParams();
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL INSIDE IF swarm_->swarmComm->univMessageReceiver.find(NEW_BOOTSTRAP)" << endl;
			}
			swarm_->swarmComm->univMessageReceiver.clear();
			return;
		}

		if (swarm_->swarmComm->univMessageReceiver.find(NEXT_GENERATION) != swarm_->swarmComm->univMessageReceiver.end()) {
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL INSIDE IF swarm_->swarmComm->univMessageReceiver.find(NEXT_GENERATION)" << endl;
			}
			return;
		}

		swarm_->swarmComm->univMessageReceiver.clear();
	}
	if (swarm_->options.verbosity >= 3) {
		cout << "RAQUEL existing while that takes mid as input" << endl;
	}
}

void Particle::checkMessagesPSO(unsigned int mid) {

	unsigned int subParID = calcsubParID(mid);
	if ((mid<0)|| (mid >= swarm_->getNumModels()) || (!models.at(mid)))
	    outputError("checkMessagesGenetic Failed. Invalid SubParticle.  PID:"+ toString(id_) + "  model id:" + toString(mid)+ "  subParID:" + toString(subParID)+". Quitting !!!" );
	//if ((mid<0)|| (mid >= swarm_->getNumModels()) || (!getSubParticle(mid)) || (!models.at(mid)))
	//outputError("checkMessagesGenetic Failed. Invalid SubParticle.  PID:"+ toString(id_) + "  model id:" + toString(mid)+ "  subParID:" + toString(subParID)+". Quitting !!!");
	if (swarm_->options.verbosity >= 3) {
		cout << "RAQUEL: I am subparticle " << subParID << " my model is " << mid << " generation " << currentGeneration_ << endl;
	}
	Model * mdl= models.at(mid);
		// Holds iterator ranges when finding items in the message holder
		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;
		if (swarm_->options.verbosity >= 3) {
			cout << "RAQUEL entering while from function that takes mid as input" << endl;
		}
		string path; //Raquel: declaring it outside the while loop


		while (1) {
			// Retrieve any messages
			int numCheckedMessages = 0;
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL I am subPar " << subParID << " model " << mid << " waiting for message." << endl;
			}
			int numMessages = swarm_->swarmComm->recvMessage(-1, subParID, -1, true, swarm_->swarmComm->univMessageReceiver, true);
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL inside LOOP ITERATION in checkMessagesGenetic" << endl;

				cout << "RAQUEL entering univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE)" << endl;
			}
			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL passed univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE)" << endl;
			}

			if (smhRange.first != smhRange.second) {
				if (swarm_->options.verbosity >= 3) {
					cout << subParID << " found final params" << endl;
				}
				for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
					//Timer tmr;

					int messageIndex = 0;
					for (auto p = simParams_.at(mid).begin(); p != simParams_.at(mid).end(); ++p) {
						//cout << id_ << " updating parameter " << p->first << " to " << sm->second.message[messageIndex] << endl;
						p->second = stod(sm->second.message[messageIndex]);
						//cout << "updated param " << p->first << ": " << p->second << endl;
						++messageIndex;
					}

					path = swarm_->options.jobOutputDir + toString(currentGeneration_ + 1) + "/";

					if (!checkIfFileExists(path)) {
						runCommand("mkdir " + path);
					}

					for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
						// Construct our filenames
						string bnglFilename = swarm_->getModelName(mid, false)+"_" + toString(id_) + "_" + toString(i) + ".bngl";
						string bnglFullPath = path + bnglFilename;
						string suffix = toString(id_) + "_" + toString(i);
						if (swarm_->options.verbosity >= 3) {
							cout << "RAQUEL CREATING BNGL: " << bnglFullPath << endl;
						}
						// And generate our models
						if (swarm_->options.models.at(mid)->getHasGenerateNetwork()){
							// If we're using ODE solver, output .net and .bngl
							mdl->outputModelWithParams(simParams_.at(mid), path, bnglFilename, suffix, false, false, true, false, false);
						}
						else {
							// If we're using network free simulation, output .bngl
							mdl->outputModelWithParams(simParams_.at(mid), path, bnglFilename, suffix, false, false, false, false, false);
						}
						if (swarm_->options.verbosity >= 3) {
							cout << "RAQUEL BNGL generated: " <<  bnglFullPath << " id_ " << toString(id_) << " i " << toString(i) << " subParID " << subParID << endl;
						}
					}

					// Tell the master we have our new params and are ready for the next generation
					//cout << id_ << " telling master we're finished " << endl;
					//swarm_->swarmComm->sendToSwarm(subParID, 0, DONE_BREEDING, false, swarm_->swarmComm->univMessageSender);

					//double t = tmr.elapsed();
					//cout << "SEND_FINAL_PARAMS took " << t << " seconds" << endl;
					++numCheckedMessages;
				}
			}

			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL BEFORE IF numCheckedMessages >= numMessages" << endl;
			}
			if (numCheckedMessages >= numMessages) {
				swarm_->swarmComm->univMessageReceiver.clear();
				if (swarm_->options.verbosity >= 3) {
					cout << "RAQUEL INSIDE IF numCheckedMessages >= numMessages = " << numCheckedMessages << ":" << numMessages << endl;
				}
				continue;
				//return 0; //Raquel: added
			}
			if (swarm_->options.verbosity >= 3) {
				cout << "RAQUEL PASSED IF numCheckedMessages >= numMessages" << endl;

			}

			if (swarm_->swarmComm->univMessageReceiver.find(FIT_FINISHED) != swarm_->swarmComm->univMessageReceiver.end()) {
				if (swarm_->options.verbosity >= 3) {
					cout << "RAQUEL entering pheromones" << endl;
				}
				swarm_->swarmComm->~Pheromones();
				if (swarm_->options.verbosity >= 3) {
					cout << "RAQUEL exiting pheromones" << endl;
				}
				exit(0);
				//return;
			}

			if (swarm_->swarmComm->univMessageReceiver.find(NEW_BOOTSTRAP) != swarm_->swarmComm->univMessageReceiver.end()) {
				if (swarm_->options.verbosity >= 3) {
					cout << "SubParticle " << subParID <<": Starting a new bootstrapping run" << endl;
				}

				++swarm_->bootstrapCounter;
				currentGeneration_ = 0;
				simParams_.at(mid).clear();
				generateParams();
				if (swarm_->options.verbosity >= 3) {
					cout << "RAQUEL INSIDE IF swarm_->swarmComm->univMessageReceiver.find(NEW_BOOTSTRAP)" << endl;
				}
				swarm_->swarmComm->univMessageReceiver.clear();
				return;
			}

			if (swarm_->swarmComm->univMessageReceiver.find(NEXT_GENERATION) != swarm_->swarmComm->univMessageReceiver.end()) {
				if (swarm_->options.verbosity >= 3) {
					cout << "RAQUEL INSIDE IF swarm_->swarmComm->univMessageReceiver.find(NEXT_GENERATION)" << endl;
				}
				swarm_->swarmComm->univMessageReceiver.clear();
				
				return;
			}

			swarm_->swarmComm->univMessageReceiver.clear();
		}
		if(swarm_->options.verbosity >=4){
			cout << "RAQUEL exiting while that takes mid as input" << endl;
		}


}

bool Particle::checkMessagesDE(unsigned int mid) {
	unsigned int subParID = calcsubParID(mid);
	if ((mid<0)|| (mid >= swarm_->getNumModels()) || (!getSubParticle(mid)) || (!models.at(mid)))
		outputError("checkMessagesDE Failed. Invalid SubParticle.  PID:"+ toString(id_) + "  model id:" + toString(mid)+ "  subParID:" + toString(subParID)+". Quitting !!!");

	Model * mdl= models.at(mid);

	while(1) {

		if (swarm_->options.verbosity >= 3) {
			cout << "Checking for messages from master" << endl;
		}

		int numCheckedMessages = 0;
		int numMessages = swarm_->swarmComm->recvMessage(-1, subParID, -1, true, swarm_->swarmComm->univMessageReceiver, true);


		if (swarm_->options.verbosity >= 3) {
			//cout << "Found " << numMessages << " messages" << endl;

			/*
			for (auto sm = swarm_->swarmComm->univMessageReceiver.begin(); sm != swarm_->swarmComm->univMessageReceiver.end(); ++sm) {
				cout << "tag: " << sm->second.tag << endl;
			}
			 */
		}

		// Holds iterator ranges when finding items in the message holder
		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		while (numCheckedMessages < numMessages) {
			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
			if (smhRange.first != smhRange.second) {
				if (swarm_->options.verbosity >= 3) {
					cout << "Receiving parameter list from master mid: " << mid << " simParams size: " << simParams_.size() << " SEND_FINAL_PARAMS_TO_PARTICLE " << SEND_FINAL_PARAMS_TO_PARTICLE << endl;
				}
				for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
					if (swarm_->options.verbosity >= 4){
						cout << "sm->second.message size: " << sm->second.message.size() << endl;
					}
					//cout << "sm->second.message[0]: " << sm->second.message[0] << endl;
					int messageIndex = 0;
					if(sm->second.message.size()>0){


					for (auto p = simParams_[mid].begin(); p != simParams_[mid].end(); ++p) {
						if (swarm_->options.verbosity >= 3) {
							cout << id_ << " updating parameter " << p->first << " to " << sm->second.message[messageIndex] << endl;
						}
						p->second = stod(sm->second.message[messageIndex]);
						++messageIndex;
					}

					for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {

						// Construct our filenames
						string bnglFilename = swarm_->getModelName(mid, false)+"_" + toString(id_) + "_" + toString(i) + ".bngl";
						string path = swarm_->options.jobOutputDir + toString(currentGeneration_ + 1) + "/";
						string bnglFullPath = path + bnglFilename;
						string suffix = toString(id_) + "_" + toString(i);

						if (!checkIfFileExists(path)) {

							runCommand("mkdir " + path);
						}

						// And generate our models
						if (mdl->getHasGenerateNetwork()){
							// If we're using ODE solver, output .net and .bngl
							if (swarm_->options.verbosity >= 3) {
								cout << "Using ODE solver, output .net and .bngl" << endl;
							}
							mdl->outputModelWithParams(simParams_[mid], path, bnglFilename, suffix, false, false, true, false, false);
						}
						else {
							// If we're using network free simulation, output .bngl
							if (swarm_->options.verbosity >= 3) {
								cout << " using network free simulation, output .bngl" << endl;
							}
							mdl->outputModelWithParams(simParams_[mid], path, bnglFilename, suffix, false, false, false, false, false);
						}
					}

					}
					++numCheckedMessages;

				}
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(FIT_FINISHED);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << "SubParticle " << subParID << ": Master told me to die" << endl;
				}
				if (swarm_->options.verbosity >= 3) {
					cout << "RAQUEL entering ~pheromones" << endl;
				}
				swarm_->swarmComm->~Pheromones();
				if (swarm_->options.verbosity >= 3) {
					cout << "RAQUEL exiting ~pheromones" << endl;
				}
				exit(0);
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEW_BOOTSTRAP);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << "SubParticle " << subParID << ": Starting a new bootstrapping run" << endl;
				}

				++swarm_->bootstrapCounter;
				currentGeneration_ = 1;

				swarm_->swarmComm->univMessageReceiver.clear();
				return true;
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(BEGIN_NELDER_MEAD);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << "Beginning Nelder-Mead local search" << endl;

					cout << "deserializing simplex" << endl;
				}
				std::stringstream iss;
				iss.str(smhRange.first->second.message[0]); // Simplex is serialized in the first indice of the message vector
				boost::archive::text_iarchive ia(iss);
				map<double, vector<double>> simplex;
				ia >> simplex;
				if (swarm_->options.verbosity >= 3) {
					cout << "running search" << endl;
				}
				runNelderMead(simplex, mid);
				if (swarm_->options.verbosity >= 3) {
					cout << "nelder mead finished. old calc: " << simplex.begin()->first << " new calc: " << fitCalcs.at(mid)[-1] << endl;
				}
				vector<string> paramsStr;
				paramsStr.push_back(toString(fitCalcs.at(mid)[-1]));
				for (auto p = simParams_.at(mid).begin(); p != simParams_.at(mid).end(); ++p) {
					paramsStr.push_back(toString(p->second));
					if (swarm_->options.verbosity >= 3) {
						cout << "new param: " << p->second << endl;
					}
				}
				swarm_->swarmComm->sendToSwarm(subParID, 0, SIMULATION_END, false, paramsStr);
				swarm_->swarmComm->univMessageReceiver.clear();

				return false;
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEXT_GENERATION);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout  << "SubParticle " << subParID <<  ": Master told me to move to the next iteration" << endl;
				}

				swarm_->swarmComm->univMessageReceiver.clear();
				return true;
			}
		}
		swarm_->swarmComm->univMessageReceiver.clear();
	}
}

void Particle::calculateFit(bool local, unsigned int mid) {

	//razi TODO: check if the fit calculation is corrent.

	//swarm_->outputError("Particle::calculateFit needs modification");

	bool usingSD = false;
	bool usingMean = false;

	// Construct pointers to the relevant fitting function.
	// This lets us avoid a switch/case or excess conditionals
	// within the calculation loop
	if (swarm_->options.objFunc == 1) {
		objFuncPtr = &Particle::objFunc_sumOfSquares;
	}
	else if (swarm_->options.objFunc == 2) {
		objFuncPtr = &Particle::objFunc_chiSquare;
		usingSD = true;
	}
	else if (swarm_->options.objFunc == 3) {
		objFuncPtr = &Particle::objFunc_divByMeasured;
	}
	else if (swarm_->options.objFunc == 4) {
		objFuncPtr = &Particle::objFunc_divByMean;
		usingMean = true;
	}

	double colSum;
	double setSum;
	double totalSum;
	double divisor;
	std::vector <double> totalSumAll[subParticles.size()]; //razi: o store sums for each subparticle


	//razi: loop over all subparticles associated with models

	unsigned int subParID = calcsubParID(mid);
	if ((mid<0)|| (mid >= swarm_->getNumModels()) || (!getSubParticle(mid)) || (!models.at(mid)) || (dataFiles_.at(mid).empty())){
		outputError("calculateFit Failed. Invalid SubParticle.  PID:"+ toString(id_) + "  model id:" + toString(mid)+ "  subParID:" + toString(subParID)+". Quitting !!!");
	}else{

		totalSum = 0;
		divisor = 0;

		// Loop through .exp files. Iterator points to string/dataset pair
		for (auto e = swarm_->options.expFiles.begin(); e != swarm_->options.expFiles.end(); ++e){
			if(swarm_->options.verbosity >=4){
				cout << "DATA FILE FOR MID " << mid << " " << e->first << endl;
			}
			//razi: we first need to check if the exp file belongs to this subParticle
			if (dataFiles_.at(mid).count(e->first)!=1){
				cout<<"The exp file "<< e->first<< " does not belong to model is:"<<mid<<". try other models."<<endl;
				continue; //Raquel: was a break, making the fit for MID 1 always become zero
			}
			if(swarm_->options.verbosity >=4){
				cout << "RAQUEL PASSED MID " << mid << endl;
			}
			// Loop through .exp columns. Iterator points to column/map pair
			//cout << "exp loop " << e->first << endl;
			setSum = 0;
			//cout << "first: " << e->second->dataCurrent->begin()->first << endl;
			//cout << "second " <<  e->second->dataCurrent->begin()->second.begin()->first << endl;
			//cout << " third " << e->second->dataCurrent->begin()->second.begin()->second << endl;
			for (std::map<std::string, std::map<double,double> >::iterator exp_col = e->second->dataCurrent->begin(); exp_col != e->second->dataCurrent->end(); ++exp_col) {
				// Loop through timepoints of column.  Iterator points to a timepoint/value pair

				if (usingMean) {
					//cout << "trying to set mean" << endl;
					divisor = e->second->colAverages.at(exp_col->first);
					//cout << "mean set" << endl;
					if(swarm_->options.verbosity >=4){
						cout << "RAQUEL divisor(colmean): " << divisor << endl;
					}
				}

				//cout << "col loop " << exp_col->first << endl;
				colSum = 0;
				for (std::map<double,double>::iterator timepoint = exp_col->second.begin(); timepoint != exp_col->second.end(); ++timepoint) {
					if (std::isnan(timepoint->second)) {
						continue;
					}

					//double exp = timepoint->second;
					//cout << "exp: " << exp << endl;

					if (swarm_->options.smoothing == 1) {
						//cout << "sim: " << dataFiles_.at(e->first).at(swarm_->options.smoothing)->dataCurrent->at(exp_col->first).at(timepoint->first) << endl;
					}
					else {
						//cout << "sim: " << dataFiles_.at(e->first).at(swarm_->options.smoothing + 1)->dataCurrent->at(exp_col->first).at(timepoint->first) << endl;
					}

					if (usingSD) {
						divisor = e->second->standardDeviations.at(exp_col->first).at(timepoint->first);
						//cout << "divisor: " << divisor << endl;
					}

					double sim;
					// TODO: Introduce fudge tolerance to account for precision loss in simulation control column
					if (swarm_->options.smoothing == 1) {
						if(swarm_->options.verbosity >=4){
							cout << "trying smoothing" << endl;
						}
						sim = dataFiles_.at(mid).at(e->first).at(swarm_->options.smoothing)->dataCurrent->at(exp_col->first).at(timepoint->first);
						//cout << dataCurrent->at(exp_col->first).at(timepoint->first) << endl;
						if(swarm_->options.verbosity >=4){
							cout << " [mid] " << mid << " [e->first] " << e->first << " [swarm_->options.smoothing] " << swarm_->options.smoothing << " exp_col->first " << exp_col->first << " timepoint->first " << timepoint->first << endl;
							cout << "done " << endl;
						}
					}
					else {
						if(swarm_->options.verbosity >=4){
							cout << "trying smoothing2 " << endl;
						}
						sim = dataFiles_.at(mid).at(e->first).at(swarm_->options.smoothing+1)->dataCurrent->at(exp_col->first).at(timepoint->first);
						if(swarm_->options.verbosity >=4){
							cout << "done" << endl;
						}
					}
					double sum = 0;
					if(sim != timepoint->second){ //if simulation and experimental values are exactly equal, the difference between them will be zero, so we can skip the objective function
						sum = (this->*objFuncPtr)(sim, timepoint->second, divisor);
					}
					if(swarm_->options.verbosity >= 3){
						cout << "sum = " << sum << " sim =" << sim << " timepoint->second = " << timepoint->second << " divisor = "<< divisor << endl;
						cout << "SUM: " << sum << "MID: " << mid << endl;
					}
					if (swarm_->options.bootstrap) {
						//cout << "all sets: "<< swarm_->bootstrapMaps.size() << endl;
						//cout << "files: " << swarm_->bootstrapMaps[swarm_->bootstrapCounter].size() << endl;
						//cout << "cols: " << (swarm_->bootstrapMaps[swarm_->bootstrapCounter].begin())->first << endl;
						//cout << "tps: " << swarm_->bootstrapMaps[swarm_->bootstrapCounter][exp_col->first][timepoint->first] << endl;
						//cout << "multiplying " << colSum << " by " << swarm_->bootstrapMaps[swarm_->bootstrapCounter - 1][e->first][exp_col->first][timepoint->first] << endl;
						//cout << "m: " << colSum << " by " << swarm_->bootstrapMaps[swarm_->bootstrapCounter].at(e->first).at(exp_col->first).at(timepoint->first) << endl;
						sum *= swarm_->bootstrapMaps[swarm_->bootstrapCounter][e->first][exp_col->first][timepoint->first];
					} //if (swarm_->options.bootstrap)
					colSum += sum;
				} //for through values
				setSum += colSum;
			} //for col loop
			totalSum += setSum;
		}//for (auto e = swarm_->options.expFiles.begin(); e != swarm_->options.expFiles.end(); ++e)
		// Erase our data sets
		dataFiles_.at(mid).clear();
		// Store our fit calc
		if(swarm_->options.verbosity >=4){
			cout << " totalSum " << totalSum << endl;
		}
		if (!local) {
			fitCalcs[mid][currentGeneration_] = pow(totalSum, 0.5);
			if(swarm_->options.verbosity >=4){
				cout << "RAQUEL LOCAL FIT for mid " << mid << ": " << fitCalcs[mid][currentGeneration_] << endl;
			}
		}
		else {
			fitCalcs[mid][-1] = pow(totalSum, 0.5);
			if(swarm_->options.verbosity >=4){
				cout << "RAQUEL GLOBALFIT FIT for mid " << mid << ": " << fitCalcs[mid][-1]<< endl;
			}
		}

		if (swarm_->options.verbosity >= 3) {
			cout << "Fit calculation["<<mid<<"]: " << pow(totalSum, 0.5) << endl;
		}

	} //if(!getSubParticle(mid))
}



void Particle::finalizeSim(unsigned int mid) {
	unsigned int subParID = calcsubParID(mid);
	if ((mid<0)|| (mid >= swarm_->getNumModels()) || (!getSubParticle(mid)) || (!models.at(mid)) || (dataFiles_.at(mid).empty())){
		outputError("finalizeSim Failed. Invalid SubParticle.  PID:"+ toString(id_) + "  model id:" + toString(mid)+ "  subParID:" + toString(subParID)+". Quitting !!!");
	}
	if (swarm_->options.verbosity >= 3) {cout << "Finalizing simulation. Started" << endl;}

	// Calculate our fit
	calculateFit(false, mid);  //razi: was calculateFit(), check if local should be false


	// Put our simulation params into the message vector
	if ((!this->getSubParticle(mid)) || (fitCalcs.at(mid).empty())){
		cout<<"can  not finalize simulation for subParticle with model id:"<<mid<<endl;
		return;
	}

		// Put our fit calc into the message vector
	swarm_->swarmComm->univMessageSender.push_back(toString(fitCalcs.at(mid).at(currentGeneration_)));
		//cout << "stored fit calc of " << fitCalcs.at(currentGeneration_) << " as " << swarm_->swarmComm->univMessageSender[0] << endl;


	//razi: here the slave fills in the values for free parameters [the list and order may be different for each model, should be taken care of by Master]
	for (map<string,double>::iterator i = simParams_.at(mid).begin(); i != simParams_.at(mid).end(); ++i){
		//cout << "stored param of " << i->second << endl;
		swarm_->swarmComm->univMessageSender.push_back(toString(i->second));
	}

	// Tell the swarm master that we're finished
	if (swarm_->options.verbosity >= 3) {
		cout << "Telling master that my simulation is finished" << endl;
	}

	swarm_->swarmComm->sendToSwarm(subParID, 0, SIMULATION_END, false, swarm_->swarmComm->univMessageSender);
	// Reset the message vector
	swarm_->swarmComm->univMessageSender.clear();
	if (swarm_->options.verbosity >= 3) {cout << "Finalizing simulation. Completed" << endl;}
}


void Particle::smoothRuns(unsigned int mid) {
	unsigned int subParID = calcsubParID(mid);
	//if ((mid<0)|| (mid >= swarm_->getNumModels()) || (!getSubParticle(mid)) || (!models.at(mid)) || (dataFiles_.at(mid).empty())){
	if ((mid<0)|| (mid > swarm_->getNumModels()) || (!models.at(mid)) || (dataFiles_.at(mid).empty())){ //Raquel: the test !getSubParticle(mid) is not working correctly


	if (mid<0)  cout<<"E1\n";
	if (mid >= swarm_->getNumModels()) cout<<"E2\n";
	if (!getSubParticle(mid)) cout<<"E3\n";
	if (!models.at(mid))  cout<<"E4\n";
	if (dataFiles_.at(mid).empty()) cout<<"E5\n";


		outputError("smoothRun Failed. Invalid SubParticle.  PID:"+ toString(id_) + "  model id:" + toString(mid)+ "  subParID:" + toString(subParID)+". Quitting !!!");
	}

	if (swarm_->options.verbosity >= 3) {
		cout << "Smoothing simulation outputs for particle:"<<id_<<" subParticle: "<<subParID<<endl;
	}

	map<string, map<double, double>> dataSet;
	// For each action/prefix/exp file
	for (auto action = dataFiles_.at(mid).begin(); action != dataFiles_.at(mid).end(); ++action) {
		// For each column
		// Insert a new iteration into the prefix set
		//cout << "action loop " << action->first << endl;
		for (auto col = action->second.at(1)->dataCurrent->begin(); col != action->second.at(1)->dataCurrent->end(); ++col) {
			// For each timepoint
			// This map holds time/param value pairs
			//cout << "col loop " << col->first << endl;
			map<double, double> timePairs;
			for (auto time = col->second.begin(); time != col->second.end(); ++time) {
				double sum = 0;
				int i = 0;
				// For each iteration
				for (unsigned int iteration = 1; iteration <= swarm_->options.smoothing; ++iteration) {
					//cout << "it loop " << iteration << endl;
					sum += dataFiles_.at(mid).at(action->first).at(iteration)->dataCurrent->at(col->first).at(time->first);
					++i;
				}
				double average = sum / (double)i;
				//cout << "average: " << average << endl;
				pair<double, double> timePair;
				timePair = make_pair(time->first, average);
				timePairs.insert(timePair);
			}
			dataSet.insert(pair<string, map<double, double>>(col->first, timePairs));
		}
		action->second.insert(pair<int, Data*>(swarm_->options.smoothing + 1, new Data(dataSet)));
		dataSet.clear();
	}


	/*
	std::cout.precision(18);

	for (auto action = dataFiles_.begin(); action != dataFiles_.end(); ++action) {
		cout << action->first << endl;
		Data * data = action->second.at(swarm_->options.smoothing + 1);
		for (auto col = data->dataCurrent->begin(); col != data->dataCurrent->end(); ++col) {
			cout << col->first << endl;
			for (auto tp = col->second.begin(); tp != col->second.end(); ++tp) {
				cout << tp->second << endl;
			}
		}
	}
	 */
	//cout << id_ << " done smoothing" << endl;
}


/*
subParticle::subParticle(Swarm * swarm, int Pid, Particle * p, unsigned int modelId): Particle::Particle(swarm, Pid){
	parParticle = p;
	mid_ = modelId;   //model id
}
*/
subParticle::subParticle(Particle * p, unsigned int modelId, unsigned int subParID) {
	parParticle = p;
	mid_ = modelId;   //razi: model id startuing from 0
	subParID_=subParID;   //razi: global subPar Id

	p->addSubParticle(this, modelId, true);
}


void subParticle::setModel(unsigned int mid){
	if (!parParticle)
		cout<<"Particle::setModel Error: Parent Particle is not set yet. I'm quitting ....\n"; exit(0);

	if (!parParticle->swarm_)
		cout<<"Particle::setModel Error: Swarm object is not set yet. I'm quitting ....\n"; exit(0);

	if (mid_!=mid)
		parParticle->swarm_->outputError("Particle::setModel Error: Trying to set invalid model id. I'm quitting ....");

	if ((mid<0) || (mid >= parParticle->swarm_->options.models.size()))
		parParticle->swarm_->outputError("Particle::setModel Error: Particle model out of range. I'm quitting ....");
	model = parParticle->swarm_->options.models.at(mid);
	//parParticle->swarm_->setParticleModelId(id_, mid);
}


void Particle::addSubParticle(subParticle* subparticle, unsigned int mid, bool overwrite){

	/* razi: moved to constructor
	unsigned int gap;
	gap = mid+1 - subParticles.size();
	if (gap>0){
		for (int i=0; i < gap; i++)
			subParticles.push_back(0);  //fill dummy
	}
	if(overwrite || !(subParticles.at(mid)))
		subParticles.at(mid)=subparticle;//	subParticles.at(mid)=.insert(subparticle);

	std::map<std::string, std::map<int, Data*>> dummyFile;
	gap= mid+1 - dataFiles_.size();
	if (gap>0)
		for (int i=0; i < gap; i++)
			dataFiles_.push_back(dummyFile);

	gap = mid+1-models.size();
	if(gap>0)
		for (int i=0; i<gap; i++)
			models.push_back(0);

	std::map<std::string, double> dummyParam;
	gap= mid+1 - simParams_.size();
	if (gap>0)
		for (int i=0; i < gap; i++)
			simParams_.push_back(dummyParam);


	std::map<int, double> dummy_fit;
	gap= mid+1 - fitCalcs.size();
	if (gap>0)
		for (int i=0; i < gap; i++)
			fitCalcs.push_back(dummy_fit);
	*/
	if((mid<0) || ( (mid>=models.size()) || (mid>=subParticles.size()) || (mid>=dataFiles_.size()) || (mid>=simParams_.size()) || (mid>=fitCalcs.size())) )
		 outputError("Add subParticle failed. Invalid model id:"+toString(mid)+". I am quitting....");
	subParIDs.push_back(subparticle->subParID_);
	if(overwrite || !(subParticles.at(mid)))
		subParticles.at(mid)=subparticle;//	subParticles.at(mid)=.insert(subparticle);
	else
		cout<<"SubParticle is not added. No empty room found!!!";

}


//void Particle::addSubParticle(int){}




unsigned int Particle::getNumSubParticle(){
	/* razi: this is also ok
	unsigned int cnt;
	for (unsigned int i =0; i < subParticles.size(); i++)
		if (subParticles.at(i)) cnt++
	return cnt;
	*/
	return subParIDs.size(); //razi: make sure to reomve the particle index when removing the particle
}



bool Particle::subParExist(unsigned int subParID){
	for(unsigned int i=0; i< subParIDs.size(); i++)
		if(subParIDs.at(i)==subParID)
			return true;
	return false;
}


unsigned int Particle::calcsubParID(unsigned int mid){
	return 1+mid+(id_-1)* swarm_->getNumModels();
}
unsigned int Particle::calcParID(unsigned int subParID, unsigned int mid){
	unsigned int pID = (int)((subParID-1-mid) / swarm_->getNumModels())+1;
	if (pID != id_)
		outputError("Error when calculating particle id: "+toString(pID)+" <>" + toString(id_) + "!!!\n");
	return pID;
}



void subParticle::runNelderMead(std::map<double, std::vector<double>> simplex){
	parParticle->runNelderMead(simplex, mid_);
}
void subParticle::doParticle() {
	parParticle->doParticle(mid_);
}
void subParticle::setModel(Model * model) {
	parParticle->setModel(model,mid_);
}
void subParticle::setParam(std::pair<std::string, double> myParams){
	parParticle->setParam(myParams, mid_);
}
void subParticle::calculateFit(bool local){
	parParticle->calculateFit(local, mid_);
}
void subParticle::finalizeSim(){
	parParticle->finalizeSim(mid_);
}
void subParticle::smoothRuns(){
	parParticle->smoothRuns(mid_);
}
























































































#else

void Particle::doParticle() {
	int id = id_;

	if (swarm_->options.verbosity >= 3) {
		cout << "In doParticle(), waiting for message from master" << endl;
	}



	if (swarm_->resumingSavedSwarm) {
		cout << "Waiting for flight count from master" << endl;
		swarm_->swarmComm->recvMessage(0, id, SEND_NUMFLIGHTS_TO_PARTICLE, true, swarm_->swarmComm->univMessageReceiver);
		currentGeneration_ = stoi(*(swarm_->swarmComm->univMessageReceiver.find(SEND_NUMFLIGHTS_TO_PARTICLE)->second.message.begin()));

		swarm_->swarmComm->univMessageReceiver.clear();
		cout << "Received flight count from master" << endl;

		cout << "Waiting for param set from master" << endl;
		swarm_->swarmComm->recvMessage(0, id, SEND_FINAL_PARAMS_TO_PARTICLE, true, swarm_->swarmComm->univMessageReceiver);

		auto freeParam = swarm_->options.model->getFreeParams_().begin();
		for (auto paramVal = swarm_->swarmComm->univMessageReceiver.find(SEND_FINAL_PARAMS_TO_PARTICLE)->second.message.begin(); paramVal != swarm_->swarmComm->univMessageReceiver.find(SEND_FINAL_PARAMS_TO_PARTICLE)->second.message.end(); ++paramVal) {
			cout << "updating " << freeParam->first << " to " << *paramVal << endl;
			simParams_.insert(pair<string, double>(freeParam->first, stod(*paramVal)));
			++freeParam;
		}
		swarm_->swarmComm->univMessageReceiver.clear();
	}



	if (swarm_->options.fitType == "de") {
		// Need to figure out which island I'm on
		for (unsigned int island = 1; island <= swarm_->options.numIslands; ++island) {
			for (unsigned int particle = 1; particle <= (swarm_->options.swarmSize / swarm_->options.numIslands) * island; ++particle) {
				if (particle == id_) {
					island_ = island;
					goto theEnd;
				}
			}
		}
		theEnd:;
	}


//cout << "Particle " << id_ << " waiting to begin" << endl;
//	swarm_->swarmComm->recvMessage(0, id_, NEXT_GENERATION, true, swarm_->swarmComm->univMessageReceiver);
cout << "subParticle " << id << " waiting to begin" << endl;
	swarm_->swarmComm->recvMessage(0, id, NEXT_GENERATION, true, swarm_->swarmComm->univMessageReceiver); //Master should send message with subParID
	swarm_->swarmComm->univMessageReceiver.clear();
cout << "Particle " << id_ << " starting" << endl;




	if (swarm_->options.verbosity >= 3) {
		cout << "In doParticle(), entering main run loop" << endl;
	}

	bool doContinue = true;
	bool doRunModel = true;
	while(doContinue) {
		if (swarm_->bootstrapCounter > 0 && currentGeneration_ == 1) {
			swarm_->swarmComm->recvMessage(0, id, NEXT_GENERATION, true, swarm_->swarmComm->univMessageReceiver);
			//int Pheromones::recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage, int messageID) {

			swarm_->swarmComm->univMessageReceiver.clear();
		}

		cout << "DoParticle: Generation is " << currentGeneration_ << endl;

		if (doRunModel) {
			for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
if (swarm_->options.verbosity >= 3) cout<<"doParticle-1 try to run model for i:"<<i<<endl;
				runModel(i);
			}

			if (swarm_->options.smoothing > 1) {
				if (swarm_->options.verbosity >= 3) cout<<"doParticle-2 try to smooth run "<<endl;
					smoothRuns();
			}else cout<<"Smoothing is not required\n";

			finalizeSim();
		}
		else{
			//cout<<"do Runmodel was not active\n";
			doRunModel = true;
		}

		if (swarm_->options.fitType == "ga") {
			checkMessagesGenetic();
		}
		else if (swarm_->options.fitType == "pso") {
			//cout << "RAQUEL before checkMessagesPSO" << endl;
			checkMessagesPSO(mid);
			//cout << "RAQUEL After checkMessagesPSO" << endl;

		}
		else if (swarm_->options.fitType == "de") {

			cout << "RAQUEL Entering check messages DE" << endl;

			checkMessagesDE(mid);
			cout << "RAQUEL Done check messages DE" << endl;

		}
		else if (swarm_->options.fitType == "sa") {
			doRunModel = checkMessagesDE();
		}

		if (doRunModel == true) {
			// Next generation
			++currentGeneration_;
		}
	}
	if (swarm_->options.verbosity >= 3) {
		cout << "In doParticle(), finished ...." << endl;
	}
}


void Particle::setModel(Model * model) {
	model_ = model;
}

void Particle::setParam(pair<string,double> myParams) {
	simParams_.insert(myParams);
}

void Particle::runModel(unsigned int iteration, bool localSearch) {

	// First get our path and filename variables set up for use in model generation, sim command, etc
	string bnglFilename = toString(id_) + "_" + toString(iteration) + ".bngl";
	string path = swarm_->options.jobOutputDir + toString(currentGeneration_) + "/";
	string bnglFullPath = path + bnglFilename;
	cout << "Particle::runModel: path is: " << path << endl;

	string suffix = toString(id_) + "_" + toString(iteration);

	string pipePath;

	// Only need to generate files if we're in the first generation. In subsequent generations
	// the model generation is handled by breeding parents
	if (currentGeneration_ == 1 && !localSearch) {
		if (swarm_->options.model->getHasGenerateNetwork()){
			if (swarm_->bootstrapCounter == 0) {
				string netFilename = toString(id_) + "_" + toString(iteration) + ".net";
				string netFullPath = path + netFilename;
				if (iteration == 1) {
					swarm_->options.model->parseNet(netFullPath);
				}
			}
			model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, true, false, false);
		}
		else {
			model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
		}
	}
	else if (localSearch) {
		cout << "regenerating models.." << endl;
		// And generate our models
		if (swarm_->options.model->getHasGenerateNetwork()){
			// If we're using ODE solver, output .net file
			bnglFilename = boost::regex_replace(bnglFilename, boost::regex("bngl$"), string("net"));
			model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, true);
		}
		else {
			// If we're using network free simulation, output .bngl
			model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
		}
	}

	// Generate .gdat pipes if we're using pipes

	string outputSuffix;
	if (swarm_->options.usePipes) {
		for (std::map<std::string,Model::action>::iterator i = model_->actions.begin(); i != model_->actions.end(); ++i) {
			if (i->second.scanParam.size() > 0) {
				outputSuffix = ".scan";
			}
			else {
				outputSuffix = ".gdat";
			}

			pipePath = path + i->first + "_" + toString(id_) + "_" + toString(iteration) + outputSuffix;
			cout << "Particle::runModel pp is: " << pipePath << endl;
			createParticlePipe(pipePath.c_str());
		}
	}




	int ret = -1;
	{
		// Construct our simulation command
		string command = swarm_->options.bngCommand + " --outdir " + path + " " + bnglFullPath + " >> " + path + toString(id_) + ".BNG_OUT 2>&1";
		if (swarm_->options.usePipes) {
			command += " &";
		}

		if (swarm_->options.verbosity >= 3) {
			cout << "Particle::runModel: Running model with command: " << command << endl;
		}

		// Run simulation command
		//int ret = system(command.c_str());
		ret = runCommand(command);
	}
	cout<<"Particle::runModel: simulation completed. collect data ...\n";

	// Check for simulation command success
	//Raquel: changed this temporarily
	if (ret == 0 || ret==-1 || ret==256) { // TODO: Need to check for simulation status when using pipes. Going by return code doesn't work there because we're using the & operator
		// really not sure how we can do this easily
		cout<<"Particle::runModel: completing runmodel: success ...\n";
		string outputSuffix;
		// Save our simulation outputs to data objects
		for (std::map<std::string,Model::action>::iterator action = model_->actions.begin(); action != model_->actions.end(); ++action) {
			if (action->second.scanParam.size() > 0) {
				outputSuffix = ".scan";
				//cout<<"Particle::runModel: A1\n";
			}
			else {
				//cout<<"Particle::runModel: A2\n";
				outputSuffix = ".gdat";
			}
			//cout<<"Particle::runModel: A3\n";
			string dataPath = path + action->first + "_" + toString(id_) + "_" + toString(iteration) + outputSuffix;
			dataFiles_[action->first].insert(pair<int, Data*>(iteration, new Data(dataPath, swarm_, false)));
			//cout<<"Particle::runModel: A4\n";
		}
	}
	else {
		// If our return code is not 0, tell the master that the simulation failed
		cout << "Particle::runModel: I failed completing runmodel ..." << endl;
		cout << " dataPath " << dataPath;
		swarm_->swarmComm->sendToSwarm(int(id_), 0, SIMULATION_FAIL, true, swarm_->swarmComm->univMessageSender);
	}
cout<<"Particle::runModel: runmodel finished ...\n";
}


void Particle::checkMessagesGenetic() {

	// Holds iterator ranges when finding items in the message holder
	pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

	while (1) {
		// Retrieve any messages
		int numCheckedMessages = 0;
		int numMessages = swarm_->swarmComm->recvMessage(-1, id_, -1, true, swarm_->swarmComm->univMessageReceiver, true);

		smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
		if (smhRange.first != smhRange.second) {
			cout << "RAQUEL id_ " << id_ << "found final " << endl;
			for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {
				//Timer tmr;

				int messageIndex = 0;
				for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
					cout << id_ << " updating parameter " << p->first << " to " << sm->second.message[messageIndex] << endl;
					p->second = stod(sm->second.message[messageIndex]);
					//cout << "updated param " << p->first << ": " << p->second << endl;
					++messageIndex;
				}

				string path = swarm_->options.jobOutputDir + toString(currentGeneration_ + 1) + "/";

				if (!checkIfFileExists(path)) {
					runCommand("mkdir " + path);
				}

				for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
					// Construct our filenames
					string bnglFilename = toString(id_) + "_" + toString(i) + ".bngl";
					string bnglFullPath = path + bnglFilename;
					string suffix = toString(id_) + "_" + toString(i);

					// And generate our models
					if (swarm_->options.model->getHasGenerateNetwork()){
						// If we're using ODE solver, output .net and .bngl
						model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, true, false, false);
					}
					else {
						// If we're using network free simulation, output .bngl
						model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
					}

				}

				// Tell the master we have our new params and are ready for the next generation
				cout << id_ << "RAQUEL: telling master we're finished id_: " << id_ << endl;
				swarm_->swarmComm->sendToSwarm(id_, 0, DONE_BREEDING, false, swarm_->swarmComm->univMessageSender);

				//double t = tmr.elapsed();
				//cout << "SEND_FINAL_PARAMS took " << t << " seconds" << endl;
				++numCheckedMessages;
			}
		}

		if (numCheckedMessages >= numMessages) {
			swarm_->swarmComm->univMessageReceiver.clear();
			continue;
		}


		if (swarm_->swarmComm->univMessageReceiver.find(FIT_FINISHED) != swarm_->swarmComm->univMessageReceiver.end()) {
			swarm_->swarmComm->~Pheromones();
			exit(0);
		}

		if (swarm_->swarmComm->univMessageReceiver.find(NEW_BOOTSTRAP) != swarm_->swarmComm->univMessageReceiver.end()) {
			if (swarm_->options.verbosity >= 3) {
				cout << id_ << ": Starting a new bootstrapping run" << endl;
			}

			++swarm_->bootstrapCounter;
			currentGeneration_ = 0;
			simParams_.clear();
			generateParams();

			swarm_->swarmComm->univMessageReceiver.clear();
			return;
		}

		if (swarm_->swarmComm->univMessageReceiver.find(NEXT_GENERATION) != swarm_->swarmComm->univMessageReceiver.end()) {
			return;
		}

		swarm_->swarmComm->univMessageReceiver.clear();
	}
}

void Particle::checkMessagesPSO() {

	while(1) {

		if (swarm_->options.verbosity >= 3) {
			cout << "Checking for messages from master" << endl;
		}

		int numCheckedMessages = 0;
		int numMessages = swarm_->swarmComm->recvMessage(-1, id_, -1, true, swarm_->swarmComm->univMessageReceiver, true);


		if (swarm_->options.verbosity >= 3) {
			cout << "Found " << numMessages << " messages" << endl;

			/*
			for (auto sm = swarm_->swarmComm->univMessageReceiver.begin(); sm != swarm_->swarmComm->univMessageReceiver.end(); ++sm) {
				cout << "tag: " << sm->second.tag << endl;
			}
			 */
		}

		// Holds iterator ranges when finding items in the message holder
		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		while (numCheckedMessages < numMessages) {
			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
			if (smhRange.first != smhRange.second) {
				//cout << "Receiving parameter list from master" << endl;
				for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {

					int messageIndex = 0;
					for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
						cout << id_ << " updating parameter " << p->first << " to " << sm->second.message[messageIndex] << endl;
						p->second = stod(sm->second.message[messageIndex]);
						++messageIndex;
					}

					string path = swarm_->options.jobOutputDir + toString(currentGeneration_ + 1) + "/";

					if (!checkIfFileExists(path)) {
						runCommand("mkdir " + path);
					}


					for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {
						// Construct our filenames
						string bnglFilename = toString(id_) + "_" + toString(i) + ".bngl";
						string bnglFullPath = path + bnglFilename;
						string suffix = toString(id_) + "_" + toString(i);

						// And generate our models
						if (swarm_->options.model->getHasGenerateNetwork()){
							// If we're using ODE solver, output .net and .bngl
							model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, true, false, false);
						}
						else {
							// If we're using network free simulation, output .bngl
							model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
						}

					}

					++numCheckedMessages;
				}
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(FIT_FINISHED);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Master told me to die" << endl;
				}

				exit(0);
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEW_BOOTSTRAP);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Starting a new bootstrapping run" << endl;
				}

				++swarm_->bootstrapCounter;
				currentGeneration_ = 1;

				swarm_->swarmComm->univMessageReceiver.clear();
				return;
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEXT_GENERATION);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Master told me to move to the next iteration" << endl;
				}

				swarm_->swarmComm->univMessageReceiver.clear();
				return;
			}
		}
		swarm_->swarmComm->univMessageReceiver.clear();
	}
}

bool Particle::checkMessagesDE() {

	while(1) {

		if (swarm_->options.verbosity >= 3) {
			cout << "Checking for messages from master" << endl;
		}

		int numCheckedMessages = 0;
		int numMessages = swarm_->swarmComm->recvMessage(-1, id_, -1, true, swarm_->swarmComm->univMessageReceiver, true);


		if (swarm_->options.verbosity >= 3) {
			//cout << "Found " << numMessages << " messages" << endl;

			/*
			for (auto sm = swarm_->swarmComm->univMessageReceiver.begin(); sm != swarm_->swarmComm->univMessageReceiver.end(); ++sm) {
				cout << "tag: " << sm->second.tag << endl;
			}
			 */
		}

		// Holds iterator ranges when finding items in the message holder
		pair <Pheromones::swarmMsgHolderIt, Pheromones::swarmMsgHolderIt> smhRange;

		while (numCheckedMessages < numMessages) {
			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(SEND_FINAL_PARAMS_TO_PARTICLE);
			if (smhRange.first != smhRange.second) {
				cout << "Receiving parameter list from master" << endl;
				for (Pheromones::swarmMsgHolderIt sm = smhRange.first; sm != smhRange.second; ++sm) {

					int messageIndex = 0;
					for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
						//cout << id_ << " updating parameter " << p->first << " to " << sm->second.message[messageIndex] << endl;
						p->second = stod(sm->second.message[messageIndex]);
						++messageIndex;
					}

					for (unsigned int i = 1; i <= swarm_->options.smoothing; ++i) {

						// Construct our filenames
						string bnglFilename = toString(id_) + "_" + toString(i) + ".bngl";
						string path = swarm_->options.jobOutputDir + toString(currentGeneration_ + 1) + "/";
						string bnglFullPath = path + bnglFilename;
						string suffix = toString(id_) + "_" + toString(i);

						if (!checkIfFileExists(path)) {
							runCommand("mkdir " + path);
						}

						// And generate our models
						if (swarm_->options.model->getHasGenerateNetwork()){
							// If we're using ODE solver, output .net and .bngl
							model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, true, false, false);
						}
						else {
							// If we're using network free simulation, output .bngl
							model_->outputModelWithParams(simParams_, path, bnglFilename, suffix, false, false, false, false, false);
						}
					}

					++numCheckedMessages;
				}
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(FIT_FINISHED);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Master told me to die" << endl;
				}

				exit(0);
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEW_BOOTSTRAP);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Starting a new bootstrapping run" << endl;
				}

				++swarm_->bootstrapCounter;
				currentGeneration_ = 1;

				swarm_->swarmComm->univMessageReceiver.clear();
				return true;
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(BEGIN_NELDER_MEAD);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << "Beginning Nelder-Mead local search" << endl;
				}
				cout << "deserializing simplex" << endl;
				std::stringstream iss;
				iss.str(smhRange.first->second.message[0]); // Simplex is serialized in the first indice of the message vector
				boost::archive::text_iarchive ia(iss);
				map<double, vector<double>> simplex;
				ia >> simplex;

				cout << "running search" << endl;
				runNelderMead(simplex);
				cout << "nelder mead finished. old calc: " << simplex.begin()->first << " new calc: " << fitCalcs[-1] << endl;
				vector<string> paramsStr;
				paramsStr.push_back(toString(fitCalcs[-1]));
				for (auto p = simParams_.begin(); p != simParams_.end(); ++p) {
					paramsStr.push_back(toString(p->second));
					cout << "new param: " << p->second << endl;
				}
				swarm_->swarmComm->sendToSwarm(id_, 0, SIMULATION_END, false, paramsStr);
				swarm_->swarmComm->univMessageReceiver.clear();

				return false;
			}

			smhRange = swarm_->swarmComm->univMessageReceiver.equal_range(NEXT_GENERATION);
			if (smhRange.first != smhRange.second) {

				if (swarm_->options.verbosity >= 3) {
					cout << id_ << ": Master told me to move to the next iteration" << endl;
				}

				swarm_->swarmComm->univMessageReceiver.clear();
				return true;
			}
		}
		swarm_->swarmComm->univMessageReceiver.clear();
	}
}

void Particle::calculateFit(bool local) {
	if (swarm_->options.verbosity >= 3) {
		cout << "Calculating fit" << endl;
	}

	bool usingSD = false;
	bool usingMean = false;

	// Construct pointers to the relevant fitting function.
	// This lets us avoid a switch/case or excess conditionals
	// within the calculation loop
	/*if (swarm_->options.objFunc == 1) {
                objFuncPtr = &Particle::objFunc_sumOfSquares;
       	}
	else if (swarm_->options.objFunc == 2) {
               	objFuncPtr = &Particle::objFunc_chiSquare;
                usingSD = true;
        }
	else if (swarm_->options.objFunc == 3) {
                objFuncPtr = &Particle::objFunc_divByMeasured;
       	}
	else if (swarm_->options.objFunc == 4) {
               	objFuncPtr = &Particle::objFunc_divByMean;
                usingMean = true;
        }
	*/
	switch (swarm_->options.objFunc) {

		case 1: 
			objFuncPtr = &Particle::objFunc_sumOfSquares;
			break;

		case 2:
			objFuncPtr = &Particle::objFunc_chiSquare;
			usingSD = true;
			break;

		case 3:
			objFuncPtr = &Particle::objFunc_divByMeasured;
			break;

		case 4:
			objFuncPtr = &Particle::objFunc_divByMean;
			usingMean = true;
			break;

		default:
			outputError("calculateFit failed. No objFunc");
			break;
	}
	double colSum;
	double setSum;
	double totalSum = 0;
	double divisor = 0;

	// Loop through .exp files. Iterator points to string/dataset pair
	for (auto e = swarm_->options.expFiles.begin(); e != swarm_->options.expFiles.end(); ++e) {
		// Loop through .exp columns. Iterator points to column/map pair
		//cout << "exp loop " << e->first << endl;
		setSum = 0;
		//cout << "first: " << e->second->dataCurrent->begin()->first << endl;
		//cout << "second " <<  e->second->dataCurrent->begin()->second.begin()->first << endl;
		//cout << " third " << e->second->dataCurrent->begin()->second.begin()->second << endl;
		for (std::map<std::string, std::map<double,double> >::iterator exp_col = e->second->dataCurrent->begin(); exp_col != e->second->dataCurrent->end(); ++exp_col) {
			// Loop through timepoints of column.  Iterator points to a timepoint/value pair

			if (usingMean) {
				//cout << "trying to set mean" << endl;
				divisor = e->second->colAverages.at(exp_col->first);
				//cout << "mean set" << endl;
			}

			//cout << "col loop " << exp_col->first << endl;
			colSum = 0;
			for (std::map<double,double>::iterator timepoint = exp_col->second.begin(); timepoint != exp_col->second.end(); ++timepoint) {
				if (std::isnan(timepoint->second)) {
					continue;
				}

				//double exp = timepoint->second;
				//cout << "exp: " << exp << endl;

				if (swarm_->options.smoothing == 1) {
					//cout << "sim: " << dataFiles_.at(e->first).at(swarm_->options.smoothing)->dataCurrent->at(exp_col->first).at(timepoint->first) << endl;
				}
				else {
					//cout << "sim: " << dataFiles_.at(e->first).at(swarm_->options.smoothing + 1)->dataCurrent->at(exp_col->first).at(timepoint->first) << endl;
				}

				if (usingSD) {
					divisor = e->second->standardDeviations.at(exp_col->first).at(timepoint->first);
					//cout << "divisor: " << divisor << endl;
				}

				double sim;
				// TODO: Introduce fudge tolerance to account for precision loss in simulation control column
				if (swarm_->options.smoothing == 1) {
					sim = dataFiles_.at(e->first).at(swarm_->options.smoothing)->dataCurrent->at(exp_col->first).at(timepoint->first);
				}
				else {
					sim = dataFiles_.at(e->first).at(swarm_->options.smoothing+1)->dataCurrent->at(exp_col->first).at(timepoint->first);
				}

				double sum = (this->*objFuncPtr)(sim, timepoint->second, divisor);

				if (swarm_->options.bootstrap) {
					//cout << "all sets: "<< swarm_->bootstrapMaps.size() << endl;
					//cout << "files: " << swarm_->bootstrapMaps[swarm_->bootstrapCounter].size() << endl;
					//cout << "cols: " << (swarm_->bootstrapMaps[swarm_->bootstrapCounter].begin())->first << endl;
					//cout << "tps: " << swarm_->bootstrapMaps[swarm_->bootstrapCounter][exp_col->first][timepoint->first] << endl;
					//cout << "multiplying " << colSum << " by " << swarm_->bootstrapMaps[swarm_->bootstrapCounter - 1][e->first][exp_col->first][timepoint->first] << endl;
					//cout << "m: " << colSum << " by " << swarm_->bootstrapMaps[swarm_->bootstrapCounter].at(e->first).at(exp_col->first).at(timepoint->first) << endl;
					sum *= swarm_->bootstrapMaps[swarm_->bootstrapCounter][e->first][exp_col->first][timepoint->first];
				}
				colSum += sum;
			}
			setSum += colSum;
		}
		totalSum += setSum;
	}

	// Erase our data sets
	dataFiles_.clear();

	// Store our fit calc
	if (!local) {
		fitCalcs[currentGeneration_] = pow(totalSum, 0.5);
	}
	else {
		fitCalcs[-1] = pow(totalSum, 0.5);
	}

	if (swarm_->options.verbosity >= 3) {
		cout << "Fit calculation: " << pow(totalSum, 0.5) << endl;
	}
}


void Particle::finalizeSim() {

	if (swarm_->options.verbosity >= 3) {
		cout << "Finalizing simulation" << endl;
	}

	// Calculate our fit
	calculateFit();

	// Put our fit calc into the message vector
	swarm_->swarmComm->univMessageSender.push_back(toString(fitCalcs.at(currentGeneration_)));
	//cout << "stored fit calc of " << fitCalcs.at(currentGeneration_) << " as " << swarm_->swarmComm->univMessageSender[0] << endl;

	// Put our simulation params into the message vector
	for (map<string,double>::iterator i = simParams_.begin(); i != simParams_.end(); ++i){
		//cout << "stored param of " << i->second << endl;
		swarm_->swarmComm->univMessageSender.push_back(toString(i->second));
	}

	// Tell the swarm master that we're finished
	if (swarm_->options.verbosity >= 3) {
		cout << "Telling master that my simulation is finished" << endl;
	}

	swarm_->swarmComm->sendToSwarm(id_, 0, SIMULATION_END, false, swarm_->swarmComm->univMessageSender);
	// Reset the message vector
	swarm_->swarmComm->univMessageSender.clear();
}

void Particle::smoothRuns() {

	if (swarm_->options.verbosity >= 3) {
		cout << "Smoothing simulation outputs" << endl;
	}

	map<string, map<double, double>> dataSet;
	// For each action/prefix/exp file
	for (auto action = dataFiles_.begin(); action != dataFiles_.end(); ++action) {
		// For each column
		// Insert a new iteration into the prefix set
		//cout << "action loop " << action->first << endl;
		for (auto col = action->second.at(1)->dataCurrent->begin(); col != action->second.at(1)->dataCurrent->end(); ++col) {
			// For each timepoint
			// This map holds time/param value pairs
			//cout << "col loop " << col->first << endl;
			map<double, double> timePairs;
			for (auto time = col->second.begin(); time != col->second.end(); ++time) {
				double sum = 0;
				int i = 0;
				// For each iteration
				for (unsigned int iteration = 1; iteration <= swarm_->options.smoothing; ++iteration) {
					//cout << "it loop " << iteration << endl;
					sum += dataFiles_.at(action->first).at(iteration)->dataCurrent->at(col->first).at(time->first);
					++i;
				}
				double average = sum / (double)i;
				//cout << "average: " << average << endl;
				pair<double, double> timePair;
				timePair = make_pair(time->first, average);
				timePairs.insert(timePair);
			}
			dataSet.insert(pair<string, map<double, double>>(col->first, timePairs));
		}
		action->second.insert(pair<int, Data*>(swarm_->options.smoothing + 1, new Data(dataSet)));
		dataSet.clear();
	}


	/*
	std::cout.precision(18);

	for (auto action = dataFiles_.begin(); action != dataFiles_.end(); ++action) {
		cout << action->first << endl;
		Data * data = action->second.at(swarm_->options.smoothing + 1);
		for (auto col = data->dataCurrent->begin(); col != data->dataCurrent->end(); ++col) {
			cout << col->first << endl;
			for (auto tp = col->second.begin(); tp != col->second.end(); ++tp) {
				cout << tp->second << endl;
			}
		}
	}
	 */
	//cout << id_ << " done smoothing" << endl;
}


#endif











