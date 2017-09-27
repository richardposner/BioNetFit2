/*============================================================================
// Name        : GenFit2.cpp
// Authors     : Brandon Thomas, Abolfazl Razi, Alex Ionkov, Raquel Dias
// Version     : 2.0
// Last Update : 2017-06-22
// Copyright   :
// Description :
//============================================================================*/

#include "GenFit2.hh"

using namespace std;

#ifdef VER2
int main(int argc, char *argv[]) {
	srand(clock());
	int generation;
	string action;
	string configFile;
	string type;
	int pID, ret;//Raquel removed unused nModels,
	bool verbose = false;
	string execPath;

	execPath = boost::filesystem::canonical(boost::filesystem::system_complete(argv[0])).string();

//	SetExecPath(argv[0]);
	vector<string> configFiles; //razi : we could have more than one config file for joint optimization. But we chose to use one config file with multiple models,	However,this is still a vector

	string expFile; //razi : These files are not required for ordinary runs. Included as a reference for artificial simulation i order to generate gdat files from exp file by adding noise !!!!
	//string bnglFile;

	if (argc < 2){
		cout <<"Please provide enough information, type BioNetFit -h for help. I am quitting .... \n";
		return 0;
	}
    //cout << "1" << endl;
	map <string, vector <string>> cmdLine;
	map <string, vector <string>>::iterator it;

	//cout << "2" << endl;

	if (readCommandLine(argc, const_cast<const char **> (argv), cmdLine)==0){
		cout<<"Invalid command line format, quitting ..."<<endl;
		return 0;
	}
	if(verbose){
		cout << "Parsed command line" << endl;
	}
/*
	for(it=cmdLine.begin(); it!=cmdLine.end(); it++){
		cout<<"Found Switch: "<< it->first<<endl;
		if (it->second.size()>0){
			cout<< "     Options: ";
			for (unsigned int ii=0; ii< it->second.size(); ii++){
				cout<< it->second[ii]<<", ";
		./	}
			cout<<endl;
		}else
			cout<< "     No Options for this switch"<<endl;
	}
*/

	//default values
	action="run";
	type = "master";
	generation =0;
	pID = 0;
	int subParID=0; int mid=-1;
	ret = 0;
	verbose=false;


	configFiles.clear();
	for(it=cmdLine.begin(); it!=cmdLine.end(); it++){
		if (it->first=="a"){
				if (it->second.size()==0){
					cout<<"Specify only one action"<<endl;
					return 0;
				}
				action = it->second.at(0);
				if (action != "run" && action != "cluster" && action != "load" && action != "results" && action != "resume") {
					cout << "Error: Couldn't find a valid 'action' in your arguments, add switch -h for more info...";
					cout << "Default action is" << action;
				}
		}

		if (it->first=="t"){
				if (it->second.size()==0){
					cout<<"Specify only one type, check help for more info..."<<endl;
					cout << "type " << it->second.at(0) << endl;
					return 0;
				}
				type = it->second.at(0);
				if (type != "master" && type != "particle") {
					cout << "Error: Couldn't find a valid 'type' in your arguments., add switch -h for more info...";
					cout << "Default type is" << type;
				}
		}

		if (it->first=="c"){
				if (it->second.size()<1){
					cout<<"Define at least one config file..."<<endl;
					return 0;
				}
				for (unsigned int ii=0; ii< it->second.size(); ii++){
					configFiles.push_back(it->second[ii]);
				}
		}
		if (it->first=="e"){
			cout<<"-e  #of files :"<<it->second.size();int ik; cin>>ik;
			if (it->second.size()==1){
				expFile = it->second.at(0);   if (verbose) cout<<"exp file found: -e "<< expFile <<endl;
			}else{
				cout<<"Define only one exp file..., #of files :"<<it->second.size();return 0;
			}
		}


		if (it->first=="p"){
				if (it->second.size()!=1){
					cout<<"Specify the pID [an integer value]..."<<endl;
					return 0;
				}
				subParID = stoi(it->second.at(0));
		}



		if (it->first=="n"){
				if (it->second.size()!=1){
					cout<<"Specify the number of models [an integer value]..."<<endl;
					return 0;
				}
				//nModels = stoi(it->second.at(0));
		}


		if (it->first=="g"){
				if (it->second.size()!=1){
					cout<<"Specify the generation [an integer value]..."<<endl;
					return 0;
				}
				generation = stoi(it->second.at(0));
		}

		if (it->first=="h"){
			outputHelp(); return 0;
		}
		if (it->first=="v"){
			verbose=true;
		}
	}

	if(verbose){
		cout << "Finished command line argument parsing" << endl;
	}

	if ((type != "particle") && (pID!=0)){
		cout<<"particle id entered, but type is not set to particle, I am quitting ..."<<endl;	return 0;
	}

	if(verbose){
		cout << "Finished testing for errors" << endl;
	}

	if (verbose){
		cout <<endl<< "************************************************************************************************"<<endl;
		cout << "*                                   GENERAL PARAMETERS                                         *"<<endl;
		cout << "************************************************************************************************"<<endl;
		cout << "Type: " << type << endl;
		cout << "Action: " << action << endl;
		cout << "Config files: ";


		for (unsigned int ii=0; ii<configFiles.size(); ii++)
			cout<<configFiles.at(ii)<<", ";
		cout<<endl;
		if (type == "particle") cout << "subParID:" << subParID << endl;
		cout << "Generation: " << generation << endl;
		cout << "Show debugging messages: " << verbose<< endl;
		#ifdef WIN_VER
			cout<<"OS: windows"<<endl;
		#else
			#ifdef CYGWIN
		 	 	 cout<<"OS: CYGWIN"<<endl;
			#else
		 	 	 cout<<"OS: Linux/MAC"<<endl;
			#endif
		#endif
 		cout << "************************************************************************************************"<<endl;
	}

	if(verbose){
		cout << "going through config files " << endl;
	}
	//if(configFiles.size()!=1){
	//	cout<<"Error: More than one config file is not supported ...."<<endl;
	//	return 0;
	//}else{
		configFile = configFiles.at(0);  //razi: this line added 2017-1-9, later check for master mode XXX
		if (verbose) cout<<"GenFit: Processing Config File:  "<< configFile<<endl;
		configFile = convertToAbsPath(configFile);

		if (!(checkIfFileExists(configFile))) {
			cout << "Config File: "<<configFile<<" does not exist, I am quitting"<<endl; return 0;
		}
		configFiles.at(0) = configFile;
	//}






	// Regardless of type or action, we need to set up the Swarm
	//Timer tmr;

	//double t = tmr.elapsed();
	//cout << "Adding .conf took " << t << " seconds" << endl;

	//tmr.reset();

	ret = chdir("..");
	if (ret == -1){ //program runs in mainpath/bin
		cout << "Failed moving to the main directory....\n";
	}

	Swarm *s;
	if (type == "master") {
		Config myconfig(configFile);
		if (verbose){
			myconfig.verbose_on();
		}
		if (action != "load" && action != "resume") { //run, cluster, results

			if (verbose) cout<<"Initialization Phase for a new Simulations."<<endl;
			s = myconfig.createSwarmFromConfig();
//cout<<"XXX1-AA1\n"; mypause();
			if (verbose) try{myconfig.printDetails();}catch(int e) {cout<<"error occurred wile debugging Err-no:"<<e<<endl;};

			if (s->options.useCluster && s->options.clusterSoftware!="BNF2mpi") {
				action = "cluster";
			}

			//t = tmr.elapsed();
			//cout << "Processing .conf took " << t << " seconds" << endl;

			s->currentGeneration = 1;

			//if(verbose) cout<<"Execution path to be set to: "<< mainpath()<<endl;
			if(verbose){
				cout << "argv[0]: " << argv[0] << ", absolute: " << execPath << endl;
			}
			s->setExePath(execPath); //mainpath()+"/src/BioNetFit2"

			s->isMaster = true;
			if(verbose) cout<<"Swarm Initialization Started.\n";
			s->initComm();
			if(verbose) cout<<"Swarm Initialization Finished.\n";
			s->isMaster = false;
		}

		if (action == "cluster" || action == "run") {

			int randNum;
			for (unsigned int mid = 0; mid < s->getNumModels(); ++ mid){
				randNum = rand();
				string serializedSwarmPath = s->getModelName(mid,false)+"_"+ to_string(static_cast<long long int>(randNum)) + ".sconf";
				if (verbose) cout<<"Start a New Simulation [run/cluster]: sConf File: "<< serializedSwarmPath << endl;

				std::ofstream ofs(serializedSwarmPath);
				if (ofs.is_open()) {
					s->setsConf(convertToAbsPath(serializedSwarmPath), mid);
					boost::archive::binary_oarchive ar(ofs);
					ar & s;
					ofs.close();
				}else cout<<"Filed opening file:"<<serializedSwarmPath<<endl;
			}
			s->isMaster = true;
		}

		if (action == "cluster") {
			//string runCmd = s->getClusterCommand(string(convertToAbsPath(argv[0])));

			string runCmd = s->getClusterCommand(execPath);

			if (verbose) cout<<"Cluster Command is: "<<runCmd<<endl;

			if (s->options.saveClusterOutput) {
				string outputPath = s->options.outputDir + "/" + s->options.jobName + "_cluster_output";
				//cout << "string: " << outputPath << endl;
				if (!checkIfFileExists(outputPath)) {
					string makeClusterOutputDirCmd = "mkdir " + outputPath;
					if (runCommand(makeClusterOutputDirCmd) != 0) {
						cout << "Warning: Couldn't create cluster output directory with command: " << makeClusterOutputDirCmd << ". Turning off save_cluster_output" << endl;
						s->options.saveClusterOutput = false;
						runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));
					}
				}
			}

			cout << "Running BioNetFit on cluster with command: " << runCmd << endl;

			if(runCommand(runCmd) != 0) {
				outputError("Error: Couldn't launch BioNetFit on cluster with command: " + runCmd + ". Quitting.");
			}

			return 0;
		}
		else if (action == "load") {
			std::ifstream ifs(configFile);

			if (ifs.is_open()) {
				boost::archive::binary_iarchive ar(ifs);

				// Load data
				ar & s;
				ifs.close();
			}
			else {
				outputError("Error: Couldn't load config file: " + configFile + ".");
			}

			if (s->options.useCluster) {
				setenv("OMPI_MCA_mpi_warn_on_fork","0",1);
			}

			s->isMaster = true;
//support cygwin
			s->setExePath(execPath);
			//was s->setExePath(convertToAbsPath(argv[0]));
			s->initComm();
		}
		else if (action == "resume") {
			string swarmState = configFile + "/swarmState.sconf";
			std::ifstream ifs(swarmState);

			if (ifs.is_open()) {
				cout << "trying to load swarm..." << endl;
				boost::archive::binary_iarchive ar(ifs);

				// Load data
				ar & s;
				ifs.close();
			}
			else {
				outputError("Error: Couldn't load swarm state: " + swarmState + ".");
			}

			cout << "loaded" << endl;

			s->resumingSavedSwarm = true;

			if (s->options.useCluster) {
				setenv("OMPI_MCA_mpi_warn_on_fork","0", 0);
				int randNum;

				for (unsigned int mid = 0; mid < s->getNumModels(); ++ mid){
					randNum = rand();
					string serializedSwarmPath = s->getModelName(mid,false)+"_"+ to_string(static_cast<long long int>(randNum)) + ".sconf";
					std::ofstream ofs(serializedSwarmPath);
					if (ofs.is_open()) {
						s->setsConf(convertToAbsPath(serializedSwarmPath), mid);
						//cout << "Path is: " << s->getsConf() << endl;
						cout << "serializedSwarmPath " << serializedSwarmPath << endl;
						boost::archive::binary_oarchive ar(ofs);
						ar & s;
						ofs.close();
					}
				}

				s->isMaster = true;
				//string runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));
				string runCmd = s->generateSlurmMultiProgCmd(execPath);

				if (s->options.saveClusterOutput) {
					string outputPath = s->options.outputDir + "/" + s->options.jobName + "_cluster_output";
					if (!checkIfFileExists(outputPath)) {
						string makeClusterOutputDirCmd = "mkdir " + outputPath;
						if (runCommand(makeClusterOutputDirCmd) != 0) {
							cout << "Warning: Couldn't create cluster output directory with command: " << makeClusterOutputDirCmd << ". Turning off save_cluster_output" << endl;
							s->options.saveClusterOutput = false;
							//runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));
							runCmd = s->generateSlurmMultiProgCmd(execPath);

						}
					}
				}

				cout << "Running BioNetFit on cluster with command: " << runCmd << endl;

				if(runCommand(runCmd) != 0) {
					outputError("Error: Couldn't launch BioNetFit on cluster with command: " + runCmd + ". Quitting.");
				}
			}
			else {
				s->initComm();
				s->initRNGS(s->options.seed);
				s->isMaster = true;
				s->doSwarm();
			}

			return 0;
		}

		if (action != "results") {
			s->doSwarm();
		}
		else {

			string messageFilePath = s->options.jobOutputDir + ".req";

			ofstream outFile;
			outFile.open(messageFilePath);

			if (outFile.is_open()) {
				outFile << "output results";
				outFile.close();
			}
			else {
				string errMsg = "Error: Couldn't open file " + messageFilePath + " to request results from the swarm master.";
				outputError(errMsg);
			}
		}
	}

	// We are a particle
	else if (type == "particle"){

		if (verbose) {
			cout<<"1-GenFit:Particle:   started.... Reading sconfig file:"<<configFile <<endl;
		}
		// Try to open the serialized swarm
		while(1) {
			// Create and input archive
			std::ifstream ifs(configFile);

			if (ifs.is_open()) {
				boost::archive::binary_iarchive ar(ifs);

				// Load data
				ar & s;
				//if (verbose && (subParID==4)) s->printDetails();

				ifs.close();
				break;
			}
		}
		if (verbose) {
			cout<<"2-GenFit:Particle:   File opened ...\n"<<endl;
		}
		//in Particle mode, the model is provided as serialized config file

		if (!expFile.empty()){
			expFile = convertToAbsPath(expFile);
			if (!(checkIfFileExists(expFile))) {
				cout << "Exp File: "<<expFile<<" does not exist, I am quitting"<<endl; return 0;
			}
			if (verbose) {cout<<"Particle: Exp file is: " <<expFile <<". Enter a number \n";} //mypause();
		}else{if(verbose){cout<<"No Exp file is provided\n";}}


		//razi: extract particle Id from subparticle id and number of models M,
		// P1=[s1,s2,...sM], P2=[s2M-1 s2M-2    s2M], ....
		//unsigned int parId = pID; unsigned int subParId = 0; //just one subparticle for each instance of running (in particle mode)
		s->initComm();
		if (s->options.clusterSoftware == "BNF2mpi"){
			subParID = s->swarmComm->getRank();
		}
		if(s->options.verbosity >= 3){
			cout << "RRR subParID " << subParID << endl;
		}
		mid = fcalcMID(subParID, s->getNumModels());
		pID = fcalcParID(subParID,s->getNumModels());

		s->setDefaultModelIndex(mid);
		if (s->options.verbosity >= 3) {cout<<"Slave:  Working on Particle:"<< pID<<"   subParticle:" << subParID<< "   Model id:"<< mid << "  Model file:"<< s->getModelName(mid,false)<<endl;}
		if (s->options.verbosity >= 3) {
			cout<<"setting Exp file:"<< expFile <<" for model id:"<<mid <<endl;
		}
		s->setExpPath("", expFile,mid); //the problem is here--> s->addExp() --> new Data(path, this, true) --> Data::parseData --> the exp file in not in the action list
		if(verbose || s->options.verbosity >= 1){
			s->printDetails();
		}
		s->isMaster = false;

		//support cygwin
		s->setExePath(execPath);
		//cout<<"3-GenFit:Particle:   Exe path is set to : "<<(mainpath()+"/bin/BioNetFit")<<endl;
		//was s->setExePath(convertToAbsPath(argv[0]));
		if (verbose) {
			cout<<"4-GenFit:Particle:   init command"<<endl;
		}
		//s->initComm();
		//subParID = s->swarmComm->getRank();
		if (s->options.verbosity >= 3) {
			cout << "RRR supar" << subParID << endl;
		}
		if (pID == 0) {
			pID = s->swarmComm->getRank();
		}
		if (s->options.verbosity >= 3) {
			cout << " pID = s->swarmComm->getRank(); result: " << pID << endl;
		}
		s->initRNGS(s->options.seed + pID);
		if (verbose) {
			cout<<"5-GenFit:Particle:   create particle."<<endl;
		}



//cout<<"5-1-GenFit:Particle..."<<endl;  mypause();
		Particle *p = s->createParticle(pID);
//cout<<"5-2-GenFit:Particle..."<<endl;  mypause();
		if (s->options.verbosity >= 3) {
			cout<<"Particle mode: The model is:"<< s->options.models.at(mid)->getName()<<endl;
		}
//cout<<"5-3-GenFit:Particle..."<<endl;  mypause();


		p->setModel(s->options.models.at(mid), mid);  //razi: this is critical and must be set
		for(unsigned int mm =0; mm < s->getNumModels(); mm++)   p->setModel(s->options.models.at(mm), mm);  //razi: no harm to set all models

//cout<<"5-4-GenFit:Particle..."<<endl;  mypause();
		subParticle * sp = new subParticle(p, mid, subParID); //razi, was wrong subParticle * sp = new subParticle(p, subParId);

		if(verbose){cout << "Launched subPar: " << sp << endl;}
//cout<<"5-5-GenFit:Particle..."<<endl;  mypause();



		if (s->currentGeneration == 1) {
			if (s->options.verbosity >= 1) {
				cout<<"6-GenFit:Particle:   generate params."<<endl;
			}
			p->generateParams();
			if (s->options.verbosity >= 1) {
				cout<<"6-2GenFit:Particle:   generate params."<<endl;
			}
		}

		if (s->options.useCluster) {
			setenv("OMPI_MCA_mpi_warn_on_fork","0",1);
		}
		if (s->options.verbosity >= 1) {
			cout<<"7-GenFit:Particle:   doParticle started."<<endl;
		}
		p->doParticle(mid);
		if (s->options.verbosity >= 1) {
			cout<<"8-GenFit:Particle:   doParticle finished."<<endl;
		}
	}
	if (s->options.verbosity >= 1) {
		cout<<"GenFit: Fitting Finished Kill Pheromones."<<endl;
	}
	s->swarmComm->~Pheromones();
	//if (verbose) s->printDetails();
	return 0;
}


















#else   //old ver
int main(int argc, char *argv[]) {
	srand(clock());
	int generation = 0;
	string action;
	string configFile;
	string type;
	int pID = 0;
	int ret = 0;  //chdir error var

	// GenFit2 [conf_file]
	if (argc == 2) {
		configFile = argv[1];
	}
	// GenFit2 [action] [conf_file]
	else if (argc == 3) {
		action = argv[1];
		configFile = argv[2];
	}
	// GenFit2 [type] [action] [conf_file]
	// Only used internally -- user should never need to specify type
	else if (argc == 4) {
		type = argv[1];
		action = argv[2];
		configFile = argv[3];
	}
	// GenFit2 [type] [pID] [action] [conf_file]
	else if (argc == 5) {
		type = argv[1];
		pID = atoi(argv[2]);
		action = argv[3];
		configFile = argv[4];
	}
	// GenFit2 [type] [pID] [action] [generation] [conf_file]
	else if (argc == 6) {
		type = argv[1];
		pID = atoi(argv[2]);
		action = argv[3];
		generation = atoi(argv[4]);
		configFile = argv[5];
	}
	// Output help if we have too too few or too many arguments
	else {
		outputHelp();
		return 0;
	}

	// Default action is to 'run'
	if (action.empty())
		action = "run";

	// Default type is 'master'
	if (type.empty())
		type = "master";

	// Be sure action and type are valid
	if (action != "run" && action != "cluster" && action != "load" && action != "results" && action != "resume") {
		outputError("Error: Couldn't find a valid 'action' in your arguments.");
	}
	if (type != "master" && type != "particle") {
		outputError("Error: Couldn't find a valid 'type' in your arguments.");
	}

	if (!(checkIfFileExists(configFile))) {
		cout << configFile;
	}
	configFile = convertToAbsPath(configFile);

	ret = chdir("../");
	if (ret == -1)
		cout << "chdir(\"../\") failed\n";


	/*
	cout << "type: " << type << endl;
	cout << "action: " << action << endl;
	cout << "config: " << configFile << endl;
	cout << "pID: " << pID << endl;
	cout << "generation: " << generation << endl;
	*/

	// Regardless of type or action, we need to set up the Swarm
	//Timer tmr;

	//double t = tmr.elapsed();
	//cout << "Adding .conf took " << t << " seconds" << endl;

	//tmr.reset();

	Swarm *s;
	if (type == "master") {
		Config myconfig(configFile);
		if (action != "load" && action != "resume") {
			s = myconfig.createSwarmFromConfig();

			if (s->options.useCluster) {
				action = "cluster";
			}

			//t = tmr.elapsed();
			//cout << "Processing .conf took " << t << " seconds" << endl;

			s->currentGeneration = 1;
			
			ret = chdir("bin/");
			//use return value
			if (ret == -1)
				cout << "chdir(\"bin/\")  failed (line 116)\n";

			s->setExePath(convertToAbsPath(argv[0]));

			s->isMaster = true;
			s->initComm();
			s->isMaster = false;
		}

		if (action == "cluster" || action == "run") {
			int randNum = rand();
			string serializedSwarmPath = to_string(static_cast<long long int>(randNum)) + ".sconf";

			std::ofstream ofs(serializedSwarmPath);
			if (ofs.is_open()) {
				s->setsConf(convertToAbsPath(serializedSwarmPath));

				boost::archive::binary_oarchive ar(ofs);
				ar & s;
				ofs.close();
			}

			s->isMaster = true;
		}

		if (action == "cluster") {
			string runCmd = s->getClusterCommand(string(convertToAbsPath(argv[0])));

			if (s->options.saveClusterOutput) {
				string outputPath = s->options.outputDir + "/" + s->options.jobName + "_cluster_output";
				//cout << "string: " << outputPath << endl;
				if (!checkIfFileExists(outputPath)) {
					string makeClusterOutputDirCmd = "mkdir " + outputPath;
					if (runCommand(makeClusterOutputDirCmd) != 0) {
						cout << "Warning: Couldn't create cluster output directory with command: " << makeClusterOutputDirCmd << ". Turning off save_cluster_output" << endl;
						s->options.saveClusterOutput = false;
						runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));
					}
				}
			}

			cout << "Running BioNetFit on cluster with command: " << runCmd << endl;

			if(runCommand(runCmd) != 0) {
				outputError("Error: Couldn't launch BioNetFit on cluster with command: " + runCmd + ". Quitting.");
			}

			return 0;
		}
		else if (action == "load") {
			std::ifstream ifs(configFile);

			if (ifs.is_open()) {
				boost::archive::binary_iarchive ar(ifs);

				// Load data
				ar & s;
				ifs.close();
			}
			else {
				outputError("Error: Couldn't load config file: " + configFile + ".");
			}

			if (s->options.useCluster) {
				setenv("OMPI_MCA_mpi_warn_on_fork","0",1);
			}

			s->isMaster = true;

			ret = chdir("bin/");
			//use return value
			if (ret == -1)
				cout << "chdir(\"bin/\") failed (line 185)\n";	

#ifdef VER2
			s->setExePath(execPath);
#else
			s->setExePath(convertToAbsPath(argv[0]));
#endif
			s->initComm();
		}
		else if (action == "resume") {
			string swarmState = configFile + "/swarmState.sconf";
			std::ifstream ifs(swarmState);

			if (ifs.is_open()) {
				cout << "trying to load swarm..." << endl;
				boost::archive::binary_iarchive ar(ifs);

				// Load data
				ar & s;
				ifs.close();
			}
			else {
				outputError("Error: Couldn't load swarm state: " + swarmState + ".");
			}

			cout << "loaded" << endl;

			s->resumingSavedSwarm = true;

			if (s->options.useCluster) {
				setenv("OMPI_MCA_mpi_warn_on_fork","0", 0);
				int randNum = rand();
				string serializedSwarmPath = to_string(static_cast<long long int>(randNum)) + ".sconf";

				std::ofstream ofs(serializedSwarmPath);
				if (ofs.is_open()) {
					s->setsConf(convertToAbsPath(serializedSwarmPath));
					cout << "Path is: " << s->getsConf() << endl;

					boost::archive::binary_oarchive ar(ofs);
					ar & s;
					ofs.close();
				}

				s->isMaster = true;
				string runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));

				if (s->options.saveClusterOutput) {
					string outputPath = s->options.outputDir + "/" + s->options.jobName + "_cluster_output";
					if (!checkIfFileExists(outputPath)) {
						string makeClusterOutputDirCmd = "mkdir " + outputPath;
						if (runCommand(makeClusterOutputDirCmd) != 0) {
							cout << "Warning: Couldn't create cluster output directory with command: " << makeClusterOutputDirCmd << ". Turning off save_cluster_output" << endl;
							s->options.saveClusterOutput = false;
							runCmd = s->generateSlurmMultiProgCmd(string(convertToAbsPath(argv[0])));
						}
					}
				}

				cout << "Running BioNetFit on cluster with command: " << runCmd << endl;

				if(runCommand(runCmd) != 0) {
					outputError("Error: Couldn't launch BioNetFit on cluster with command: " + runCmd + ". Quitting.");
				}
			}
			else {
				s->initComm();
				s->initRNGS(s->options.seed);
				s->isMaster = true;
				s->doSwarm();
			}

			return 0;
		}

		if (action != "results") {
			s->doSwarm();
		}
		else {
			string messageFilePath = s->options.jobOutputDir + ".req";

			ofstream outFile;
			outFile.open(messageFilePath);

			if (outFile.is_open()) {
				outFile << "output results";
				outFile.close();
			}
			else {
				string errMsg = "Error: Couldn't open file " + messageFilePath + " to request results from the swarm master.";
				outputError(errMsg);
			}
		}
	}

	// We are a particle
	else if (type == "particle"){
		// Try to open the serialized swarm
		while(1) {
			// Create and input archive
			std::ifstream ifs(configFile);

			if (ifs.is_open()) {
				boost::archive::binary_iarchive ar(ifs);

				// Load data
				ar & s;
				ifs.close();
				break;
			}
		}


		ret = chdir("bin/");
		if (ret == -1)
			cout << "chdir(\"bin/\") failed (line 295)\n"; 


		s->isMaster = false;
		s->setExePath(convertToAbsPath(argv[0]));

		s->initComm();

		if (pID == 0) {
			pID = s->swarmComm->getRank();
		}

		s->initRNGS(s->options.seed + pID);

		Particle *p = s->createParticle(pID);
		p->setModel(s->options.model);

		if (s->currentGeneration == 1) {
			p->generateParams();
		}

		if (s->options.useCluster) {
			setenv("OMPI_MCA_mpi_warn_on_fork","0",1);
		}


		ret = chdir("bin/");
		if (ret == -1)
			cout << "chdir(\"bin/\") failed (line 321)\n";

		p->doParticle();
	}

	s->swarmComm->~Pheromones();

	return 0;
}
#endif

