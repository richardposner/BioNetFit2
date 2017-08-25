//Raquel: added new file
//Name: evaluateResults.hh
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <stdlib.h>
#include "Swarm.hh"

using namespace std;

int constraintFunction(double value1, string constraintOperator, double value2);

float evaluateResults(string inputFile1, string inputFile2, map<int,string> constraints,  string outname, float iteration1, float iteration2);

vector<string> splitResults(const string &s, char delim);
