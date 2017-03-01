/*============================================================================
// Name        : Data.hh
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Last Update : 2017-01-07
// Copyright   :
// Description :
//============================================================================*/


#ifndef DATA_HH_
#define DATA_HH_

#include "Utils.hh"
#include "Swarm.hh"

class Swarm;
class Data {

public:
#ifdef VER2
	Data(std::string path, Swarm * swarm, bool isExp, unsigned int mid);
#else
	Data(std::string path, Swarm * swarm, bool isExp);
#endif
	Data(std::map<std::string, std::map<double, double>> &dataSet);
	Data();

	void parseData();
	std::string getPath();

	void standardizeData();
	void divideByInit();
	void getColumnAverages();
	void logTransformData();

	// Column->timepoint, value
	std::map<std::string, std::map<double, double>> * dataCurrent; // This points to the most recently modified version of the data
	std::map<std::string, std::map<double, double>> standardDeviations;
	std::map<std::string, double> colAverages;

#ifdef VER2
	unsigned int mid_;
#endif

private:
	friend class boost::serialization::access;

	std::string dataPath;

	std::map<std::string,std::map<double,double> > dataOrig_;
	std::map<std::string,std::map<double,double> > dataStandardized_;
	std::map<std::string,std::map<double,double> > dataLogTransformed_;
	std::map<std::string,std::map<double,double> > dataDividedByInit_;

	Swarm * swarm_;
	bool isExp_;

	template<typename Archive>
	void serialize(Archive& ar, const unsigned version) {
		//std::cout << "serializing data" << std::endl;

		ar & standardDeviations;
		ar & colAverages;
		ar & dataPath;
		ar & dataOrig_;
		ar & dataStandardized_;
		ar & dataLogTransformed_;
		ar & dataDividedByInit_;

		ar & swarm_;
		ar & isExp_;
		ar & dataCurrent;
	}
};

#endif /* DATA_HH_ */
