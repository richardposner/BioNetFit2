/*============================================================================
// Name        : FreeParam.cpp
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Last Update : 2017-01-07
// Copyright   :
// Description :
//============================================================================*/

#include "FreeParam.hh"

using namespace std;

FreeParam::FreeParam(string parameterName) {
	parameterName_ = parameterName;
	generationMethod_ = "";
	hasMutation_ = false;
	mutationRate_ = 0;
	mutationFactor_ = 0;
	genMin_ = 0;
	genMax_ = 0;
	isLog_ = false;
}

FreeParam::FreeParam() {
	parameterName_ = "";
	generationMethod_ = "";
	hasMutation_ = false;
	mutationRate_ = 0;
	mutationFactor_ = 0;
	genMin_ = 0;
	genMax_ = 0;
	isLog_ = false;
}
