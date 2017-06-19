/*============================================================================
// Name        : Config.cpp
// Authors     : Brandon Thomas, Abolfazl Razi, Alex Ionkov
// Version     : 2.0
// Last Update : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/


#include <unordered_map>

#include "Config.hh"
#include "Utils.hh"   //razi added to use tolinux function
#include <iostream>
#include <boost/filesystem.hpp>

using namespace std;
namespace fs = boost::filesystem;

Config::Config(string configFile) {
	swarm_ = 0;
	verbose=false; //razi added
	configPath_ = convertToAbsPath(configFile);
	cout << "*** configFile: " << configPath_ << endl;
}



int Config::printDetails(){ //razi added for debugging
	if (swarm_ != NULL){
		swarm_->printDetails();
	}
//	else{
//		cout<<"No Details for config file: "<< configPath_ <<", Swarm object has not created yet !!! "<<endl;
//	}

	return 0;
}



Swarm * Config::createSwarmFromConfig () {

	swarm_ = new Swarm();

	string line;
	ifstream confFile(configPath_);
	boost::smatch matches;
	string name, value;
	unordered_multimap<string,string> pairs;
	unordered_multimap<string,string>::iterator piter; //razi added

	swarm_->setConfigPath(configPath_);
	fs::path cpath(configPath_);



	if(verbose) cout<<"Creating Swarm Object From Config File:"<< configPath_ <<endl;

	// TODO: We really should be assigning iterators inside the conditional
	// statements to avoid the extra map find() each time we find a match
	if (confFile.is_open()) {
		while (getline(confFile, line)) {
			name = "";
			value = "";

			if (line.length() > 0){
				if (line.at(0) == '#')
					continue;

				boost::regex keyVal("^\\s*([a-zA-Z_]+)\\s*=\\s*(.*)\\s*$");
				if (boost::regex_search(line,matches,keyVal)){
					name = matches[1];
					value = matches[2];
					pairs.insert(make_pair(name,value));
				}
				else {
					boost::regex keyVal("^\\s*([a-zA-Z_]+)\\s+(.*)\\s*$");
					if (boost::regex_search(line,matches,keyVal)){
						name = matches[1];
						value = matches[2];
						pairs.insert(make_pair(name,value));
					}
				}
			}
		}
		confFile.close();
	}
	else {
		swarm_->outputError("Error: Couldn't open config file " + configPath_ + " for parsing.");
	}

#ifdef VER2
	if(verbose){
		cout <<endl<<"================================================================================================"<<endl;
		     cout << "=                                  CONFIG FILE CONTENT                                         ="<<endl;
		     cout << "================================================================================================"<<endl;
		for (piter=pairs.begin(); piter!=pairs.end(); piter++){
			cout<<"Item :  " << piter->first <<",   "<< (piter->second) <<""<<endl;
		}
		cout << "================================================================================================"<<endl<<endl;
	}
#endif

	// Set verbosity
	if(pairs.find("verbosity") != pairs.end()) {
		swarm_->options.verbosity = stoi(pairs.find("verbosity")->second);
	}


	// Add our model file to the swarm
	if(pairs.find("model") != pairs.end()) {
		cout << "Processing models include:" << pairs.find("model")->second<< endl;
#ifdef VER2
		

		swarm_->setModels(cpath.parent_path().string(), pairs.find("model")->second, true);
#else
		swarm_->setModel(pairs.find("model")->second);
#endif
	}

	// Update the swarm type
	if(pairs.find("fit_type") != pairs.end()) {
		//cout << "Processing fit type" << endl;
		swarm_->setfitType(pairs.find("fit_type")->second);
	}
	// Update the swarm size
	if(pairs.find("population_size") != pairs.end()) {
		//cout << "Processing swarm size" << endl;
		int swarmSize = stoi(pairs.find("population_size")->second);

		// If the swarm size isn't even, make it even by adding a particle
		//if (swarmSize % 2 != 0 && swarm_->options.fitType == "genetic") {
		//	++swarmSize;
		//}
		swarm_->options.swarmSize = swarmSize;
	}

	// Update the sim path
	if(pairs.find("bng_command") != pairs.end()) {
		//cout << "Processing bng command" << endl;
		if (pairs.find("bng_command")->second.find("random")!=std::string::npos){
			swarm_->options.bngCommand = pairs.find("bng_command")->second;
		}else{//do for all simulators except random
			string bngcmd = pairs.find("bng_command")->second;
			string abspath = convertToAbsPath(bngcmd);
			if (!checkIfFileExists(abspath)) {
				//swarm_->outputError("Error: The specified bng command: " + swarm_->options.bngCommand + " does not exist. Quitting.");
				char *bnfroot = getenv("BNFROOT");
				if (bnfroot != NULL) {
					fs::path bnfr(bnfroot);
					abspath = (bnfr / fs::path(bngcmd)).string();
					swarm_->options.bngCommand = abspath;
				} 
				else {
				swarm_->outputError("Error: The specified bng command: " + swarm_->options.bngCommand + " does not exist. Quitting.");
				}
			}
		}
	}

	// Update the synchronicity
	if(pairs.find("synchronicity") != pairs.end()) {
		//cout << "Processing synchronicity" << endl;
		swarm_->options.synchronicity = (stoi(pairs.find("synchronicity")->second));
	}

	// Update the maximum number of generations
	if(pairs.find("max_generations") != pairs.end()) {
		//cout << "Processing max generations" << endl;
		swarm_->options.maxGenerations = (stoi(pairs.find("max_generations")->second));
	}

	if(pairs.find("output_every") != pairs.end()) {
		//cout << "Processing outputevery" << endl;
		swarm_->options.outputEvery = stoi(pairs.find("output_every")->second.c_str());
	}

	// Tell the swarm if we're using pipes
	if(pairs.find("use_pipes") != pairs.end()) {
		//cout << "Processing usepipes" << endl;
		swarm_->options.usePipes = (stoi(pairs.find("use_pipes")->second) == 1) ? true : false;
	}

	// Tell the swarm if we should delete old files
	if(pairs.find("delete_old_files") != pairs.end()) {
		//cout << "Processing delete old files" << endl;
		swarm_->options.deleteOldFiles = (stoi(pairs.find("delete_old_files")->second) == 1) ? true : false;
	}

	// Tell the swarm if we're using a cluster
	if(pairs.find("use_cluster") != pairs.end()) {
		//cout << "Processing use cluster" << endl;
		if (stoi(pairs.find("use_cluster")->second)) {

			// Set cluster platform if it was specified
			if(pairs.find("cluster_software") != pairs.end()) {
				//cout << "Processing cluster software" << endl;
				swarm_->options.clusterSoftware = pairs.find("cluster_software")->second;
			}

			swarm_->options.useCluster = (stoi(pairs.find("use_cluster")->second) == 1) ? true : false;
			swarm_->getClusterInformation();

			if (pairs.find("email_when_finished") != pairs.end() && stoi(pairs.find("email_when_finished")->second)) {
				swarm_->options.emailWhenFinished = true;
			}

			if (pairs.find("email_address") != pairs.end()) {
				swarm_->options.emailAddress = pairs.find("email_address")->second;
			}
		}
	}

	// Tell the swarm if we should save cluster output
	if(pairs.find("save_cluster_output") != pairs.end()) {
		swarm_->options.saveClusterOutput = (stoi(pairs.find("save_cluster_output")->second) == 1) ? true : false;
	}

	if(pairs.find("cluster_account") != pairs.end()) {
		swarm_->options.clusterAccount = pairs.find("cluster_account")->second;
	}

	if(pairs.find("cluster_queue") != pairs.end()) {
		swarm_->options.clusterQueue = pairs.find("cluster_queue")->second;
	}

	if(pairs.find("hostfile") != pairs.end()) {
		//cout << "Processing use cluster" << endl;
		swarm_->options.hostfile = pairs.find("hostfile")->second;
	}

	// Update seed
	if(pairs.find("seed") != pairs.end()) {
		swarm_->options.seed = stoi(pairs.find("seed")->second);
	}

	// Update swap rate
	if(pairs.find("swap_rate") != pairs.end()) {
		swarm_->options.swapRate = stof(pairs.find("swap_rate")->second);
	}

	// Update number of parents to keep unchanged in breeding
	if(pairs.find("keep_parents") != pairs.end()) {
		swarm_->options.keepParents = stoi(pairs.find("keep_parents")->second);
	}

	// Update extra weight
	if(pairs.find("extra_weight") != pairs.end()) {
		swarm_->options.extraWeight = stoi(pairs.find("extra_weight")->second);
	}

	// Whether or not to force difference parents
	if(pairs.find("force_different_parents") != pairs.end()) {
		swarm_->options.forceDifferentParents = (stoi(pairs.find("force_different_parents")->second) == 1) ? true : false;
	}

	// How many retries when breeding
	if(pairs.find("max_breeding_retries") != pairs.end()) {
		swarm_->options.maxRetryDifferentParents = stoi(pairs.find("max_breeding_retries")->second);
	}

	// Set fit value that will cause fit to end
	if(pairs.find("smoothing") != pairs.end()) {
		swarm_->options.smoothing = stoi(pairs.find("smoothing")->second);

		if (swarm_->options.smoothing == 0) {
			swarm_->options.smoothing = 1;
		}
	}

	// Set output directory
	if(pairs.find("output_dir") != pairs.end()) {
		if (!(checkIfFileExists("output_dir"))) {
                	boost::filesystem::path dir("output_dir");
                	if (boost::filesystem::create_directory("output")) {
                        	std::cout << "created output dir" << "\n";
                	}
       		}
		swarm_->options.outputDir = convertToAbsPath(pairs.find("output_dir")->second);
#if 0
#ifdef VER2
		swarm_->options.outputDir = tolinux(swarm_->options.outputDir);    //razi added for consistency between Linux and Windows
		swarm_->options.outputDir = removeLastSlash(swarm_->options.outputDir);  //razi added later check for linux, works for windows
#endif
#endif
	}

	// Set job name
	if(pairs.find("job_name") != pairs.end()) {
		swarm_->options.jobName = pairs.find("job_name")->second;
	}

	// Set the job output directory
#if(defined (VER2) && defined(PC_VER))
	swarm_->setJobOutputDir(swarm_->options.outputDir + "/" + swarm_->options.jobName);
#else
	swarm_->setJobOutputDir(swarm_->options.outputDir + "/" + swarm_->options.jobName + "/");
#endif

	// Set fit value that will cause fit to end
	if(pairs.find("min_fit") != pairs.end()) {
		swarm_->options.minFit = stod(pairs.find("min_fit")->second);
	}

	// Set maximum number of simulations in an asynchronous genetic fit
	if(pairs.find("max_num_simulations") != pairs.end()) {
		swarm_->options.maxNumSimulations = stol(pairs.find("max_num_simulations")->second);
	}

	// Set maximum fitting time
	// TODO: Implement this in synchronous genetic fitting
	if(pairs.find("max_fit_time") != pairs.end()) {
		vector<string> timeElements;
		split(pairs.find("max_fit_time")->second, timeElements, ":");
		long timeInSeconds = (stol(timeElements[0]) * 3600) + (stol(timeElements[1]) * 60) + stol(timeElements[2]);
		swarm_->options.maxFitTime = timeInSeconds;
	}

	// Update the maximum number of parallel threads (non-cluster only)
	if(pairs.find("parallel_count") != pairs.end()) {
		if (stoi(pairs.find("parallel_count")->second) > swarm_->options.swarmSize) {
			outputError("Error: Parallel_count is set higher than population_size.");
		}
		if (swarm_->options.useCluster) {
			swarm_->options.parallelCount = swarm_->options.swarmSize;
		}
		else {
			swarm_->options.parallelCount = stoi(pairs.find("parallel_count")->second);
		}

	}




	// Whether or not to divide by value at t=0
	if(pairs.find("divide_by_init") != pairs.end()) {
		swarm_->options.divideByInit = (stoi(pairs.find("divide_by_init")->second.c_str()) == 1) ? true : false;
	}

	// Whether or not to log transform simulation output
	if(pairs.find("log_transform_sim_data") != pairs.end()) {
		swarm_->options.logTransformSimData = stoi(pairs.find("log_transform_sim_data")->second);
	}

	// Whether or not to standardize simulation output
	if(pairs.find("standardize_sim_data") != pairs.end()) {
		swarm_->options.standardizeSimData = (stoi(pairs.find("standardize_sim_data")->second) == 1) ? true : false;
	}

	// Whether or not to standardize exp data
	if(pairs.find("standardize_exp_data") != pairs.end()) {
		swarm_->options.standardizeExpData = (stoi(pairs.find("standardize_exp_data")->second) == 1) ? true : false;
	}

	// Update fit calculation method
	if(pairs.find("objfunc") != pairs.end()) {
		swarm_->options.objFunc = stoi(pairs.find("objfunc")->second);
	}


	//////////////////////////////////////////////////////////////////
	// PSO STUFF
	//////////////////////////////////////////////////////////////////

	if (pairs.find("inertia") != pairs.end()) {
		swarm_->options.inertia = stof(pairs.find("inertia")->second);
	}

	if (pairs.find("cognitive") != pairs.end()) {
		swarm_->options.cognitive = stof(pairs.find("cognitive")->second);
	}

	if (pairs.find("social") != pairs.end()) {
		swarm_->options.social = stof(pairs.find("social")->second);
	}

	if (pairs.find("inertia") != pairs.end()) {
		swarm_->options.inertia = stof(pairs.find("inertia")->second);
	}

	if (pairs.find("nmin") != pairs.end()) {
		swarm_->options.nmin = stoi(pairs.find("nmin")->second);
	}

	if (pairs.find("nmax") != pairs.end()) {
		swarm_->options.nmax = stoi(pairs.find("nmax")->second);
	}

	if (pairs.find("inertia_init") != pairs.end()) {
		swarm_->options.inertiaInit = stof(pairs.find("inertia_init")->second);
	}

	if (pairs.find("inertia_final") != pairs.end()) {
		swarm_->options.inertiaFinal = stof(pairs.find("inertia_final")->second);
	}

	if (pairs.find("abs_tolerance") != pairs.end()) {
		swarm_->options.absTolerance = stof(pairs.find("abs_tolerance")->second);
	}

	if (pairs.find("rel_tolerance") != pairs.end()) {
		swarm_->options.relTolerance = stof(pairs.find("rel_tolerance")->second);
	}

	if (pairs.find("mutate_qpso") != pairs.end()) {
		swarm_->options.mutateQPSO = (stoi(pairs.find("mutate_qpso")->second) == 1) ? true : false;
	}

	if (pairs.find("beta_min") != pairs.end()) {
		swarm_->options.betaMin = stof(pairs.find("beta_min")->second);
	}

	if (pairs.find("beta_max") != pairs.end()) {
		swarm_->options.betaMax = stof(pairs.find("beta_max")->second);
	}

	if (pairs.find("topology") != pairs.end()) {
		swarm_->options.topology = pairs.find("topology")->second;
	}

	if (pairs.find("pso_type") != pairs.end()) {
		swarm_->options.psoType = pairs.find("pso_type")->second;
	}

	if(pairs.find("enhanced_stop") != pairs.end()) {
		swarm_->options.enhancedStop = (stoi(pairs.find("enhanced_stop")->second) == 1) ? true : false;
	}

	if(pairs.find("enhanced_inertia") != pairs.end()) {
		swarm_->options.enhancedInertia = (stoi(pairs.find("enhanced_inertia")->second) == 1) ? true : false;

	}

	if(pairs.find("num_islands") != pairs.end()) {
		swarm_->options.numIslands = stoi(pairs.find("num_islands")->second);

		if (swarm_->options.swarmSize % swarm_->options.numIslands != 0 && swarm_->options.fitType == "de") {
			outputError("Error: The number of islands must divide evenly into the population size. Quitting.");
		}
	}

	if(pairs.find("mutate_type") != pairs.end()) {
		swarm_->options.mutateType = stoi(pairs.find("mutate_type")->second);
	}

	if(pairs.find("crossover_rate") != pairs.end()) {
		swarm_->options.cr = stof(pairs.find("crossover_rate")->second);
	}

	if(pairs.find("migration_frequency") != pairs.end()) {
		swarm_->options.migrationFrequency = stoi(pairs.find("migration_frequency")->second);
	}

	if(pairs.find("num_to_migrate") != pairs.end()) {
		swarm_->options.numToMigrate = stoi(pairs.find("num_to_migrate")->second);
	}

	if(pairs.find("min_temp") != pairs.end()) {
		swarm_->options.minTemp = stof(pairs.find("min_temp")->second);
	}

	if(pairs.find("min_radius") != pairs.end()) {
		swarm_->options.minRadius = stof(pairs.find("min_radius")->second);
	}

	if(pairs.find("local_search_probability") != pairs.end()) {
		swarm_->options.localSearchProbability = stof(pairs.find("local_search_probability")->second);
	}

	if(pairs.find("rand_params_probability") != pairs.end()) {
		swarm_->options.randParamsProbability = stof(pairs.find("rand_params_probability")->second);
	}

	/////////////////////////RAQUEL adding weight parameter for ranking models by constraint/////////////////////

	if(pairs.find("constraint_weight") != pairs.end()) {
		swarm_->options.constraintWeight = stof(pairs.find("constraint_weight")->second);
	}


#ifdef VER2
	swarm_->consolidate_model_params();
#endif

	int constraintIndex = 0; //Raquel added a constraint index integer
	// Add any init param generation options
	for (unordered_multimap<string, string>::iterator pair = pairs.begin(); pair != pairs.end(); ++pair) {
		// We should use a map.equalrange to speed this part up
		if (pair->first == "random_var" || pair->first == "lognormrandom_var" || pair->first == "loguniform_var" || pair->first == "constraint" || pair->first == "constraints") {
			vector<string> paramComponents;
			split(pair->second,paramComponents);
			//string genString = pair->first + " " + paramComponents[1] + " " + paramComponents[2];

			// Make sure we have three components to work with
			if (paramComponents.size() == 3) {
				// Make sure first parameter name exists as a free parameter
				if (swarm_->options.freeParams_.count(paramComponents[0]) > 0) {  // it exist in the complete list (so belongs to some of the models)
					// Make sure 2nd and 3rd components are numeric
					if (isFloat(paramComponents[1]) && isFloat(paramComponents[2])) {
						//apply to the common list
						swarm_->options.freeParams_.at(paramComponents[0])->setGenerationMethod(pair->first);
						swarm_->options.freeParams_.at(paramComponents[0])->setParameterName(paramComponents[0]);
						swarm_->options.freeParams_.at(paramComponents[0])->setGenMin(stof(paramComponents[1]));
						swarm_->options.freeParams_.at(paramComponents[0])->setGenMax(stof(paramComponents[2]));
						if (pair->first == "lognormrandom_var" || pair->first == "loguniform_var") {
							swarm_->options.freeParams_.at(paramComponents[0])->setIsLog(true);
						}
						else if (pair->first == "random_var") {
							swarm_->options.freeParams_.at(paramComponents[0])->setParameterName(paramComponents[0]);
						}
						if (verbose)  cout << "setting " << paramComponents[0] << " to " << pair->first << ":" << paramComponents[1] << ":" << paramComponents[2] << endl;
					}
					else {
						swarm_->outputError("Error: Problem parsing your free parameter generation option in your .conf file. The min and/or max values were non-numeric.");
					}
				}
				else {
					cout << "Warning: We found a parameter '" << paramComponents[0] << "' in your .conf file, but don't see a matching free parameter specification in your model file. We will ignore this parameter." << endl;
				}
			}
			else {
				//Raquel: added constraint section to the config file
				if(paramComponents.size() == 1){
					if (pair->first == "constraint" || pair->first == "constraints") {
								//swarm_->options.constraints_.at(constraintIndex) = paramComponents[0];
								swarm_->options.constraints_[constraintIndex]=paramComponents[0];
								constraintIndex++;
								 cout << "setting " << paramComponents[0] << endl;
					}else{
						swarm_->outputError("Error: Problem parsing your free parameter generation option in your .conf file. Each parameter generation option requires three components: the parameter name, minimum, and maximum.");
					}
				}else{
					swarm_->outputError("Error: Problem parsing your free parameter generation option in your .conf file. Each parameter generation option requires three components: the parameter name, minimum, and maximum.");
				}
				//Raquel: end, added constraint section to the config file
				//swarm_->outputError("Error: Problem parsing your free parameter generation option in your .conf file. Each parameter generation option requires three components: the parameter name, minimum, and maximum.");
			}
		}
	}

	for (map<string, FreeParam*>::iterator i = swarm_->options.freeParams_.begin(); i != swarm_->options.freeParams_.end(); ++i) {
		if (!i->second) {
			outputError("Error: We found a free parameter '" + i->first + "' specified in your model file but can't find a matching parameter generator in your .conf file.");
		}
	}
std::pair <std::unordered_multimap<string,string>::iterator, std::unordered_multimap<string,string>::iterator> it;

#ifdef VER2   //razi added to include multiple exp and bngl files in the swarm object, [one exp file per model to avoid ambiguity] 2017-1-7
	if(pairs.find("exp_file") != pairs.end()){
		if(verbose) cout<<"Config::createSwarmFromConfig: EXPs files: "<<pairs.find("exp_file")->second<<".\n";
		swarm_->setExpPath(cpath.parent_path().string(),pairs.find("exp_file")->second,-1);
	}
#else

	// Add any .exp files to the swarm
	it = pairs.equal_range("exp_file");
	for (unordered_multimap<string,string>::iterator exp = it.first; exp != it.second; ++exp) {
		swarm_->addExp(exp->second);
	}
#endif

	if(pairs.find("bootstrap") != pairs.end()) {
		swarm_->options.bootstrap = stoi(pairs.find("bootstrap")->second);
		if (swarm_->options.bootstrap) {
			vector<map<string, map<string, map<double, unsigned int>>>> bootstrapMaps;
			swarm_->generateBootstrapMaps(bootstrapMaps);
			swarm_->bootstrapMaps = bootstrapMaps;
			if (checkIfFileExists(swarm_->options.outputDir + "/" + swarm_->options.jobName + "_bootstrap")) {
				runCommand("rm " + swarm_->options.outputDir + "/" + swarm_->options.jobName + "_bootstrap/*");
			}
		}
	}

	// Add any .mutation rates to the swarm
	it = pairs.equal_range("mutate");

	for (unordered_multimap<string,string>::iterator exp = it.first; exp != it.second; ++exp) {
		swarm_->addMutate(exp->second);
		swarm_->hasMutate = true;

		if (boost::regex_search(exp->second, boost::regex("^default\\s"))) {
			//cout << "found default mutation rate. breaking" << endl;
			break;
		}
	}

	// Link all model actions with their corresponding data sets (.exp files)
	// Remove any .1exp files that aren't pointed to in the model file
	// And Remove any Model::actions without corresponding .exp files specified in the .conf file

#ifdef VER2
	vector<string> prefixedActions;
	vector<vector<string> > prefixedActionsMatrix;
	Model * curmodel;


	vector<Model::action> toDeleteActs;
	//cout<<"EXP files are:"; for (map<string, Data*>::iterator i = swarm_->options.expFiles.begin(); i != swarm_->options.expFiles.end(); ++i) 		cout<<""<< i->first <<", ";	cout<<endl;
/* just for debguging
	for (int mid = 0; mid < swarm_->options.models.size(); ++mid){ //do for all modls independently
		//cout<<"checking model["<<cnt<<"]: "<< swarm_->options.models.at(mid)->getName() <<"="<< swarm_->getModel(mid) <<"., press a key to continue ....";
		curmodel = swarm_->options.models.at(mid);
		for (map<string, Model::action>::iterator i = curmodel->actions.begin(); i != curmodel->actions.end();++i) {
			cout << "The action item["<<i->first<<"] for model:"<<  swarm_->getModel(mid) <<"exists. press a key to continue ...\n"; cin>>input;
		}
	}

*/


	unsigned int mid;

	for (mid = 0; mid < swarm_->options.models.size(); ++mid){ //do for all modls independently
		curmodel = swarm_->options.models.at(mid);
		prefixedActions.clear();
		for (map<string, Model::action>::iterator i = curmodel->actions.begin(); i != curmodel->actions.end();) {
			//cout << "The action item["<<i->first<<"] for model:"<<  swarm_->getModel(mid) <<"exists.\n";
			prefixedActions.push_back(i->first);

			if (swarm_->options.expFiles.count(i->first) == 1) {
				if (swarm_->options.verbosity >=3 ) {
					cout << "Linking action " << i->first << " with exp file: " << swarm_->options.expFiles[i->first]->getPath() << endl;
				}
				i->second.dataSet = swarm_->options.expFiles[i->first];
				++i;
			}
			else { // Have a prefix but no .exp file
				cout << "Warning: The model file specifies an action with the prefix '" << i->first << "' but there isn't a matching .exp file specified in your .conf file. We will ignore this action command." << endl;
				curmodel->actions.erase(++i);
			}
		}
		//cout<<"working on model:"<<mid<<" is finished..\n";
		prefixedActionsMatrix.push_back(prefixedActions);
	}

/*
	if (verbose){  //later delete, just for debugging
		cout<<"The list of prefixes are:";
		for (int mid = 0; mid < swarm_->options.models.size(); ++mid)
			for (int ii=0; ii<prefixedActionsMatrix.at(mid).size(); ++ii) 	cout<< prefixedActionsMatrix.at(mid).at(ii) <<", "; cout<<endl;
	}
*/
	for (map<string, Data*>::iterator i = swarm_->options.expFiles.begin(); i != swarm_->options.expFiles.end(); ++i) {
		bool exp_linked = false;
		for (int mid = 0; mid < swarm_->options.models.size(); ++mid){
			curmodel = swarm_->options.models.at(mid);
			prefixedActions = prefixedActionsMatrix.at(mid);
			// Have exp but no prefix
			if(std::find(prefixedActions.begin(), prefixedActions.end(), i->first) != prefixedActions.end() )
				exp_linked = true; //found
		}
		if (!exp_linked){
			cout << "Warning: The .conf file specifies an .exp file '" << i->first << "' but there isn't a matching action command in your model file specified with the prefix=> argument." << endl;
			swarm_->options.expFiles.erase(++i);
		}
	}

#else
	vector<string> prefixedActions;
	vector<Model::action> toDeleteActs;

	for (map<string, Model::action>::iterator i = swarm_->options.model->actions.begin(); i != swarm_->options.model->actions.end();) {
		prefixedActions.push_back(i->first);

		if (swarm_->options.expFiles.count(i->first) == 1) {
			if (swarm_->options.verbosity >=3 ) {
				cout << "Linking action " << i->first << " with exp file: " << swarm_->options.expFiles[i->first]->getPath() << endl;
			}
			i->second.dataSet = swarm_->options.expFiles[i->first];
			++i;
		}
		else { // Have a prefix but no .exp file
			cout << "Warning: The model file specifies an action with the prefix '" << i->first << "' but there isn't a matching .exp file specified in your .conf file. We will ignore this action command." << endl;
			swarm_->options.model->actions.erase(++i);
		}
	}

	for (map<string, Data*>::iterator i = swarm_->options.expFiles.begin(); i != swarm_->options.expFiles.end(); ++i) {

		// Have exp but no prefix
		if(std::find(prefixedActions.begin(), prefixedActions.end(), i->first) == prefixedActions.end() ) {
			cout << "Warning: The .conf file specifies an .exp file '" << i->first << "' but there isn't a matching action command in your model file specified with the prefix=> argument." << endl;
			swarm_->options.expFiles.erase(++i);
		}
	}
#endif

	/*
	// TODO: Make sure we have either min fit, max time, or max sims when doing an asynchronous fit
	// This whole checker needs to be much more complex
	if (!swarm_->options.synchronicity) {
		if (swarm_->options.minFit == -1 && swarm_->options.maxFitTime == MAX_LONG && swarm_->options.maxNumSimulations == MAX_LONG && !swarm_->options.enhancedStop) {
			outputError("Error: You are running an asynchronous fit, but failed to specify either a minimum fitting value, maximum fit time, or maximum number of simulations in the .conf file.");
		}
	}
	 */

	/*
	for (auto i : swarm_->options.expFiles){
		cout << "Exp: " << i.second->getPath() << endl;
	}

	for (auto i : swarm_->options.modelFile->actions_){
		cout << "Prefix: " << i.prefix << endl;
	}

	for (auto i : swarm_->options.modelFile->actions_){
		cout << "DataSets: " << i.dataSet->getPath() << endl;
	}
	 */
	checkConsistency();

	return swarm_;
}

string Config::getLocation () {
	return configPath_;
}

void Config::checkConsistency() {
	/*
	 * Make sure we have SOME sort of stop criteria
	 */

	if (swarm_->options.fitType.empty()) {
		swarm_->outputError("Error: You didn't specify a fit_type in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.jobName.empty()) {
		swarm_->outputError("Error: You didn't specify a job_name. Quitting.");
	}

	if (swarm_->options.fitType != "ga" && swarm_->options.fitType != "sa" && swarm_->options.fitType != "de" && swarm_->options.fitType != "pso") {
		swarm_->outputError("Error: Unrecognized fit_type. Quitting.");
	}
#ifdef VER2
	if (swarm_->options.models.size() < 1) {
		swarm_->outputError("Error: You didn't specify a model in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}
#else
	if (swarm_->options.model == 0) {
		swarm_->outputError("Error: You didn't specify a model in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}
#endif
	if (swarm_->options.expFiles.empty()) {
		swarm_->outputError("Error: You didn't specify any .exp files in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.outputDir.empty()) {
		swarm_->outputError("Error: You didn't specify an output_dir in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.bngCommand.empty()) {
		swarm_->outputError("Error: You didn't specify a bng_command in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.swarmSize == 0) {
		swarm_->outputError("Error: You didn't specify a population_size in your .conf file. Refer to the documentation with help on setting this option. Quitting.");
	}

	if (swarm_->options.enhancedInertia && swarm_->options.fitType == "pso") {
		if (!swarm_->options.enhancedStop) {
			swarm_->outputError("Error: You set enhanced_inertia to true, but enhanced_stop to false. Enhanced inertia requires enhanced_stop to be set to true. Quitting.");
		}
	}

	if (swarm_->options.enhancedStop  && swarm_->options.fitType == "pso") {
		if (swarm_->options.nmax == 0) {
			swarm_->outputError("Error: You set enhanced_stop to true, but did not set a maximum permance value (nmax) in your .conf file. Quitting.");
		}
	}

	if (swarm_->options.useCluster && swarm_->options.emailWhenFinished && swarm_->options.emailAddress.empty()) {
		swarm_->outputError("Error: You set email_when_finished to true, but did not specify an email address. Quitting.");
	}

	if (swarm_->options.useCluster && !swarm_->options.clusterSoftware.empty() && (swarm_->options.clusterSoftware != "mpi" || swarm_->options.clusterSoftware != "slurm" || swarm_->options.clusterSoftware != "torque")) {
		swarm_->outputError("Error: You set cluster_software to an invalid option. Quitting.");
	}

	if (swarm_->options.fitType == "de") {
		if (swarm_->options.mutateType == 1 && (swarm_->options.swarmSize / swarm_->options.numIslands) <= 2 ) {
			swarm_->outputError("Error: mutate_type=1 requires at least 2 members per island. Quitting.");
		}

		if (swarm_->options.mutateType == 2 && (swarm_->options.swarmSize / swarm_->options.numIslands) <= 4 ) {
			swarm_->outputError("Error: mutate_type=2 requires at least 4 members per island. Quitting.");
		}

		if (swarm_->options.mutateType == 3 && (swarm_->options.swarmSize / swarm_->options.numIslands) <= 3 ) {
			swarm_->outputError("Error: mutate_type=3 requires at least 3 members per island. Quitting.");
		}

		if (swarm_->options.mutateType == 4 && (swarm_->options.swarmSize / swarm_->options.numIslands) <= 5 ) {
			swarm_->outputError("Error: mutate_type=4 requires at least 5 members per island. Quitting.");
		}

		if (swarm_->options.numIslands == 0) {
			swarm_->outputError("Error: You are running a differential evolution fit, but didn't specify the number of islands to use. Quitting.");
		}
	}

	if (swarm_->options.fitType == "pso") {
		if (swarm_->options.fitType == "qpso") {
			if (!swarm_->options.maxNumSimulations) {
				swarm_->outputError("Error: You are using qpso, but didn't set max_num_simulations. Quitting.");
			}
		}
	}


	if (swarm_->options.fitType == "sa") {
#ifdef VER2
		if (swarm_->options.swarmSize < (swarm_->options.models[0]->getNumFreeParams() + 1)) {
			swarm_->outputError("Error: You are using simulated annealing but your swarm size is less than N+1, where N is the number of free parameters in your model. Quitting.");
		}
#else
		if (swarm_->options.swarmSize < (swarm_->options.model->getNumFreeParams() + 1)) {
			swarm_->outputError("Error: You are using simulated annealing but your swarm size is less than N+1, where N is the number of free parameters in your model. Quitting.");
		}
#endif
		if (swarm_->options.synchronicity) {
			swarm_->outputError("Error: Simulated annealing does not support synchronous fitting. Quitting.");
		}
	}

	if (swarm_->options.fitType == "ga") {
		if (swarm_->options.synchronicity && swarm_->options.maxGenerations == 0) {
			swarm_->outputError("Error: You are running a synchronous ga fit, but didn't specify max_generations. Quitting.");
		}
	}

#ifdef VER2
	//razi: check consistency between the model and exp files.
	//TODO more checks may be needed here
	if (!swarm_->check_model_exp_consistency()){
		swarm_->outputError("Error: Models are not consistent with exp files.");
	}
#endif

}
