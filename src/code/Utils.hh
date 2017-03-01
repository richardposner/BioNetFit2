/*============================================================================
// Name        : Utils.hh
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Lat Update: : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/


#ifndef UTILS_HH_
#define UTILS_HH_

#include <iostream>
#include <map>
#include <vector>  //razi added
#include <string>  //razi added
#include <sstream>
//did not work #include <windows.h> //razi added to use GetModuleFileName function


#include <spawn.h>
#include <sys/wait.h>
#include <boost/filesystem.hpp>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp> //razi added

#include "FreeParam.hh"
#include "Setting.hh"   //razi added

class FreeParam;

void outputError(std::string errorMessage);
std::string convertToAbsPath(std::string relPath);
std::string getFilename(std::string path);
int checkIfFileExists(std::string path);
void split(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ");
void outputHelp();
bool createParticlePipe(const char * path);
//double pickWeighted(double weightSum, std::multimap<double,double>& weights, int extraWeight, boost::random::mt19937 &generalRand);
bool isFloat(std::string number);
double mutateParamGA(FreeParam* fp, double paramValue);
int runCommand(std::string cmd, std::string &result);
int runCommand(std::string cmd);
std::string toString (unsigned int theNumber);
std::string toString (int theNumber);
std::string toString (float theNumber);
std::string toString (double theNumber);
std::string toString (unsigned long theNumber);

#ifdef VER2
	void mypause();
	int readCommandLine(int argc, const char *argv[], std::map<std::string,std::vector<std::string>> &cmdLine);
	std::string mainpath();
	bool CheckFullPath(std::string inputpath);
	std::string tolinux(std::string inputpath);
	std::string removeCygwinAlias(std::string inputpath);
	std::string removeLastSlash(std::string inputpath);
	std::string addLastSlash(std::string inputpath);
	std::string removeDoubleSlash(std::string inputpath);
	int generate_gdat_file(std::string exp_file, std::string gdat_path, int seed); //razi added to generate artificial gdat files
	std::vector<std::string> split_string(std::string str, std::string delimiter);
	std::string trim(std::string str, bool ltrim, bool rtrim);

	unsigned int fcalcsubParID(unsigned int ParID, unsigned int mid, unsigned int nModels);
	unsigned int fcalcMID(unsigned int subParID, unsigned int nModels);
	unsigned int fcalcParID(unsigned int subParID, unsigned int nModels);
#endif //VER2



#endif /* UTILS_HH_ */
