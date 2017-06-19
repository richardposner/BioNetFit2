/*============================================================================
// Name        : Utils.cpp
// Authors     : Brandon Thomas, Abolfazl Razi, Alex Ionkov
// Version     : 2.0
// Lat Update: : 2017-01-7
// Copyright   :
// Description :
//============================================================================*/

#include "Utils.hh"
#include "Setting.hh"   //razi addded 2016/11/30
#include <ctime>  //razi added
using namespace boost::filesystem;
using namespace std;

string executableName;





#ifdef VER2
void mypause(){
	int ii;
	cin>>ii;
}
/*
void SetExecName(char *argv0) {
        path full_path = system_complete( path(argv0) );

        executableName = full_path.string();
}

string GetExecName() {
	return executableName;
}

string GetExecPath(void) {
        path full_path = system_complete( path(executableName) );

        //Without file name
        return full_path.stem().string();
}

string GetCurrentPath(void) {
        path full_path(current_path());

        //Without file name
        return full_path.stem().string();
}
*/
string convertToAbsPath(string relPath) {
        path fullPath;

        try{
                fullPath = canonical(relPath);
        }
        catch(...){}

        return fullPath.string();
}

/*
string mainpath(){  //razi addede
	string curdirstr;
	#ifdef WIN_VER
		//get the exe file name and extract the current working directory
		curdirstr=string(GetModuleFileName(NULL, curdir, 1000));
		curdirstr=tolinux(curdirstr);
		curdirstr=curdirstr.substr(0, curdirstr.find_last_of( "\\/" ));
	#else  //Linux or cygwin
		char curdir[1000];

		curdirstr = (getcwd(curdir, 1000) ? std::string(curdir):  std::string(""));
		if(!curdirstr.compare("cygdrive")){
			string Dletter= curdirstr.substr(10,1);
			curdirstr = Dletter + ":/" + curdirstr.substr(12,curdirstr.size()-12);
		}
	#endif
	if (curdirstr.find("/bin")!=string::npos || curdirstr.find("/Bin")!=string::npos){
		return(curdirstr.substr(0, curdirstr.size()-4));
	}
	else{
		#ifndef WIN_VER //Raquel: added this block to solve compatibility issues with Linux
			path p = canonical("./"); //Raquel: gets full path with linux sintax
			curdirstr = p.string(); //Raquel: convert it to string
		#endif
		return(curdirstr);
		//cout<<"Error: The executable file is in an improper folder, I am quitting ...."<<endl; exit(1);
	}


}

bool  CheckFullPath(string inputpath){ //razi: check whether or not in canonical format
	if((inputpath.find(":/")!=string::npos) || (inputpath.find("cygdrive")!=string::npos))
		return true;
	else
		return false;
}

string tolinux(string inputpath){   //razi: replace \ with /
	std::replace(inputpath.begin(), inputpath.end(), '\\', '/');
	return inputpath;
}

string removeDoubleSlash(string inputpath){   //razi: remove last \ or /
	//std::replace(inputpath.begin(), inputpath.end(), '\\\\', '\/');
	//std::replace(inputpath.begin(), inputpath.end(), '\/\/', '\/');
	return inputpath;
}


string removeLastSlash(string inputpath){   //razi: remove last \ or /
	if ((inputpath.substr(inputpath.size()-1,1)=="\\") || (inputpath.substr(inputpath.size()-1,1)=="/"))
		return inputpath.substr(0,inputpath.size()-1);
	else
		return inputpath;
}

string addLastSlash(string inputpath){   //razi: remove last \ or /
	if ((inputpath.substr(inputpath.size()-1,1)!="\\") && (inputpath.substr(inputpath.size()-1,1)!="/"))
		return inputpath+"/";
	else
		return inputpath;
}


string removeCygwinAlias(string inputpath){ //convert /cygdrive/X/myfolder to X:/myfolder
	if(inputpath.compare("cygdrive")){
		//string Dletter= boost::to_upper(inputpath.substr(10,1));
		string Dletter= inputpath.substr(10,1); //Drive letter
		return Dletter + ":/" + inputpath.substr(12,inputpath.size()-12);
	}else{
		return inputpath;
	}
}


string convertToAbsPath(string relPath) {
	string fullPath;

//cout<<"1-converting to canonical path format, Input:"<<relPath<<endl;
	relPath=tolinux(relPath);
//cout<<"2-converting to canonical path format, Input [Linux]:"<<relPath<<endl;
	if (CheckFullPath(relPath)){
//cout<<"3-converting to canonical path format, Already in canonical format:"<<relPath<<endl;
		fullPath = relPath;
	}else{
		string basic_path = mainpath();
//cout<<"4-converting to canonical path format, Reference Path:"<<basic_path<<endl;
		fullPath = basic_path + "/" + relPath;

	#ifndef WIN_VER //Raquel: added this block to solve compatibility issues with Linux
		path fullPath2 = canonical(relPath); //Raquel: retrieves full path for relative directory
		fullPath = fullPath2.string(); //Raquel: converts path to string
	#endif
//cout<<"5-converting to canonical path format, Full Path:"<<fullPath<<endl;
	}
	if (!exists(fullPath)){
		string errMsg = "Error in convertToAbsPath: Can't find the file: " + fullPath;

		cout<<errMsg<<"\n Anyway, Enter 1 to continue for debugging purpose... : "; int tempi; cin>>tempi;
		if (tempi==1)
			return(fullPath);
		else
			outputError(errMsg);
		//was outputError(errMsg);
	}else{
		return(fullPath);
	}

	return fullPath;
}
*/
int generate_gdat_file(string exp_file, string gdat_path, int seed){ //razi added to generate artificial gdat files
	//read exp-file
	//copy to gdat file
cout<<"\nUtils::generate_gdat_file: PiD["<<seed <<"] Starts generating file: "<< gdat_path << "   from file :" << exp_file<<"\n.";
try{
	std::ifstream ifs(exp_file);
	std::ofstream ofs(gdat_path);
	std::string line, oline;
	std::string word;
	int wcnt, linecnt;
	int val, rnd1, rnd;
	std::string::size_type pos = 0;
	std::string::size_type prev = 0;
	boost::random::mt19937 generalRand;

//cout<<"Utils:g-1   try to generate file:" <<gdat_path << " from file: " << exp_file<< endl;
	if ((ifs.is_open()) && (ofs.is_open())) {


		linecnt=0;
		while (getline(ifs, line)) {
			if (line.length() > 0){
				linecnt++;
				if (line.at(0) == '#'){
					oline=line;
//cout<<"Utils:g-2   A header line is read:" <<line << endl;
				}else{

//cout<<"Utils:g-3   A valid line is read:" <<line << endl;
					//data line break line into smaller pieces

					pos = 0; prev = 0; wcnt=0; oline=""; line=line+"\t"; //add tab to the end

					srand(time(0)+seed);
					while ((pos = line.find("\t", prev)) != std::string::npos)
					{
						if (pos != prev) //otherwise, most likely the first tab
						{
							wcnt++;
							word = line.substr(prev, pos - prev);


//cout<<"Utils:g-4   A word is read:" <<word << endl;
							if (wcnt == 1){  //the first word is time, don't touch it
								oline=oline + "\t" + word;
							}else{
								//for all pieces, if numerical value , change all of it except the first one
								val = atoi(word.c_str());


								rnd1 = rand()%100;
								val = val + (int) (rnd1);
								if (val < 0){ val=-val;}

//cout<<"Utils:g-6  the new value is :" <<val<< endl;
								oline = oline +"\t" + toString(val);
//cout<<"Utils:g-7  output line is :" <<oline <<endl;

							}
						}
						prev = pos + 1;
					}
					//write the output line into the output file
				}
//cout<<"Utils:g-9  try to dump output line:" <<oline <<endl;
				ofs << oline << std::endl;
			}
		}
		ifs.close();
		ofs.close();
//cout<<"Utils:g-10   finished generating file:" <<gdat_path << " from file: " << exp_file<< " successfully." <<endl;
		return 0;  //successful
	}else{
//cout<<"Utils:g-11   try to generate file:" <<gdat_path << " from file: " << exp_file<< ", Can't open input or output file !!!" <<endl;
		return -1;
	}
}
catch(int err){
	cout<< "Error occurred when generating artificial gdata file....\n\n"; return -1;
}
}

vector<std::string> split_string(std::string str, std::string delimiter){
	vector<std::string> strs;
	size_t pos = 0;
	std::string token;
	while ((pos = str.find(delimiter)) != std::string::npos) {
		token = str.substr(0, pos);
		token = trim(token,true,true);
		strs.push_back(token);
		str.erase(0, pos + delimiter.length());
	}
	if (!str.empty()){
		strs.push_back(trim(str,true,true));  //add the last piece
	}
	return strs;
}

std::string trim(std::string str, bool ltrim, bool rtrim){
	std::string s =str;
	if(ltrim)  //trim from left
		s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	if(rtrim)  //trim from right
		 s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
//cout<<"string: "<<str <<" len["<< str.size() <<"]is trimed to: "<<s<<" len["<< s.size()<<"]\n.";
	return s;
}
unsigned int fcalcsubParID(unsigned int ParID, unsigned int mid, unsigned int nModels){
	return 1+mid+(ParID-1)* nModels;
}
unsigned int fcalcMID(unsigned int subParID, unsigned int nModels){
	return (unsigned int)((subParID-1) % nModels);
}
unsigned int fcalcParID(unsigned int subParID, unsigned int nModels){
	return (unsigned int)((subParID-1)/ nModels)+1;
}





#else
  //previous version
string convertToAbsPath(string relPath) {
	path fullPath;

	try{
		fullPath = canonical(relPath);
	}
	catch(...){}

	if (!exists(fullPath)){
		string errMsg = "Error: Can't find the file: " + relPath;
		outputError(errMsg);
	}

	return fullPath.string();
}
#endif




int checkIfFileExists(string path) {
	if (exists(path))
		return 1;
	else
		return 0;
}




void split(const string& str, vector<string>& tokens, const string& delimiters) {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

string getFilename(string path) {
	return boost::filesystem::path(path).stem().string();
}

#ifdef VER2
void outputHelp() {
	cout << "BioNetFit Usage:" << endl;
	cout << "BioNetFit switch [for example -a]  [information for example: run] ...\n";
	cout << "\nSwitches:\n";
	cout << "   -a   [action]        choose action type" << endl;
	cout << "        Options are:" << endl;
	cout << "              run" << endl;
	cout << "              cluster" << endl;
	cout << "              load" << endl;
	cout << "              results" << endl;
	cout << "              resume" << endl;
	cout << "   -t   [type]          choose type" << endl;
	cout << "        Options are:" << endl;
	cout << "              master" << endl;
	cout << "              particle" << endl;
	cout << "        Options are:" << endl;
	cout << "              PSO(Particle Swarm Optimization)"<<endl;
	cout << "              SA(simulated Annealing)"<<endl;
	cout << "              GA(Genetic Algorithm)"<<endl;
	cout << "              DE(Differential Evolution) " << endl;
	cout << "   -c   [config file]   enter config files" << endl;
	cout << "   -e   [exp file]      enter reference exp files for particle [optional]" << endl;
	cout << "   -b   [bngl file]     enter reference bngl files for particle [optional]" << endl;
	cout << "   -p   [pID]           enter an integer value for pID" << endl;
	cout << "   -g   [Generation]    enter an integer value for Generation" << endl;
	cout << "   -s   [Simulator]     choose simulators [Not developed yet]" << endl;
	cout << "   -h   []              prints helpful information" << endl;
	cout << "   -v   []              activate debugging messages" <<endl;
	cout << "   -o   [Optimization Technique]  choose optimization techniques, [Not developed yet]" << endl;

}


int readCommandLine(int argc, const char *argv[], map<string,vector<string>> &cmdLine){  //razi added
	int i, j=0; string sw; string temp;
	string s;
	for (i=1; i<argc; i++){
		s=argv[i];
		//cout<<"reading command line str["<<i<<"] = "<<s<<endl;
		if (s.substr(0,1).compare("-")==0){  //a new switch is found
			sw=s.substr(1,s.size()-1);
			if (sw.empty()){
				cout<<"Invalid switch:"<< s  <<"  value not specified, quitting .... "<<endl;
				return 0;
			}else{ //now read the options
				//cout<<"switch found:-"<<sw<<endl;

				vector<string> options; options.clear();
				j=1; bool exitflag=false;
				while(!exitflag){
					if ((i+j) < argc){
						temp = argv[i+j];
						//cout<<"checking options["<<j<<"]=" << temp << endl;
						if  (temp.substr(0,1).compare("-")==0){
							//cout<<"This is not an option, it is a new switch, i set to:"<<i+j<< endl;
							i=i+j-1;  //move the pointer to the new switch
							exitflag=true;
						}
						else{  //a valid option
							//cout<<"Option found: "<<argv[i+j]<<"   added to the list of options"<<endl;
							options.push_back(argv[i+j]);
							j++;  //collect other options
						}
					}else{
						//cout<<"Command line reading finished, reached the end while reading the options i+j:"<<i+j<<endl;
						i=i+j-1;  //move the pointer to the new switch
						exitflag=true;
					}
				}
				//if (options.size()>0), cout<<"# "<<options.size() <<"Options were found and added to the list"<<endl;
				cmdLine[sw]=options;
			}
		}else{
			cout<<"invalid format for switch:"<<s<<" Expecting to start with -, quitting .... "<<endl;
			return 0;
		}
	}
	return 1;
}



#else
void outputHelp() {
	cout << "GenFit2 Usage:" << endl;
	cout << "GenFit2 [config_file]" << endl;
}
#endif



bool createParticlePipe(const char * path) {

	if (checkIfFileExists(path)) {
		unlink(path);
	}

	int fifo_status = mkfifo(path, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
	if (fifo_status) {
		cout << "Warning: Couldn't create pipe with path: " << path << endl;
		return false;
	}
	else {
		return true;
	}
}

bool isFloat(string number) {
	istringstream iss(number);
	float f;
	iss >> noskipws >> f; // noskipws considers leading whitespace invalid
	// Check the entire string was consumed and if either failbit or badbit is set
	return iss.eof() && !iss.fail();
}

int runCommand(string cmd, string &result) {

	std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);

	if (!pipe) {
		return 1;
	}

	char buffer[1024];
	result.clear();

	while (!feof(pipe.get())) {
		if (fgets(buffer, 1024, pipe.get()) != NULL)
			result += buffer;
	}

	return 0;
}

int runCommand(string cmd) {

	//cout << "Running command: " << cmd << endl;


	int ret = system(cmd.c_str());
	return ret;



	/*
	FILE *in;

	//cout << "Running command: " << cmd << endl;

	if(!(in = popen(cmd.c_str(), "r"))){
		return 1;
	}
	int status = pclose(in);
	//printf("Exit code: %d\n", WEXITSTATUS(status));

	return status;
	*/

	/*
	char *c_cmd = new char[cmd.length()+1];
	std::strcpy(c_cmd, cmd.c_str());

	pid_t pid;
	char *argv[] = {"bash", "-c", c_cmd, NULL};
	int status;
	printf("Run command: %s\n", c_cmd);
	status = posix_spawn(&pid, "/bin/bash", NULL, NULL, argv, environ);
	if (status == 0) {
		printf("Child pid: %i\n", pid);
		if (waitpid(pid, &status, 0) != -1) {
			printf("Child exited with status %i\n", status);
		} else {
			perror("waitpid");
		}
	} else {
		printf("posix_spawn: %s\n", strerror(status));
	}

	return status;
	*/
}

void outputError(string errorMessage) {
	cout << errorMessage << endl;

	exit (1);
}

string toString(unsigned int theNumber) {
	return to_string(static_cast<long long int>(theNumber));
}

string toString(int theNumber) {
	return to_string(static_cast<long long int>(theNumber));
}

string toString(double theNumber) {
	return to_string(static_cast<long double>(theNumber));
}

string toString(float theNumber) {
	return to_string(static_cast<long double>(theNumber));
}

string toString(unsigned long theNumber) {
	return to_string(static_cast<long double>(theNumber));
}
