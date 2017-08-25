//Raquel: added new function for comparing results given a list of constraints Frebuary 1st, 2017
#include "Evaluate.hh"



using namespace std;

//Raquel: this function splits the result line into a vector of strings
vector<string> splitResults(const string &s, char delim) {
    stringstream ss;
    ss.str(s);
    string item;
    vector<string> elems;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


//Raquel: this function apply the logical test based on logical operators parsed from the constraint file
int constraintFunction(double value1, string constraintOperator, double value2){

    bool result;

    if(constraintOperator == "=" || constraintOperator == "=="){
      result = (value1 == value2); return result;
    }

    if(constraintOperator == ">"){
      result = (value1 > value2); return result;
    }

    if(constraintOperator == "<"){
       result = (value1 < value2); return result;
    }

    if(constraintOperator == ">="){
      result = (value1 >= value2); return result;
    }

    if(constraintOperator == "<="){
       result = (value1 <= value2); return result;
    }

    if(constraintOperator == "!="){
       result = (value1 != value2); return result;
    }
   //0 means false
   //1 means true
   //-1 means that the logial operator was not found in this function
   return -1;

}

//Raquel: this function will parse the results, apply the logical operations for testing for constraints
//Raquel: this function will write TRUE or FALSE in the output, depending on whether the constraints are fulfilled
float evaluateResults(string inputFile1, string inputFile2, map<int,string> constraints, string outname, float iteration1, float iteration2) {
  ifstream simResult1;
  ifstream simResult2;
  ofstream outFile;
  /*
  if(argc != 4){
     cout << endl << endl;
     cout << "Usage: ./evaluate_results result_file1.txt result_file2.txt constraint_file.cons";
     cout << endl << endl;
     return 1;
  }
*/
  outFile.open(outname, ios_base::app);
  simResult1.open (inputFile1);
  simResult2.open (inputFile2);

  int count = 0;
  vector<string> constraintValue1;
  vector<string> constraintValue2;
  vector<string> constraintOp;
  vector<string> variableNames1;
  vector<string> variableValues1;
  vector<string> variableNames2;
  vector<string> variableValues2;


  string constraintItem;


  map< string, map <int, double> > namesValues1;
  map< string, map <int, double> > namesValues2;

  string possibleConstraints("<>=!");
  size_t pos = 0;
  size_t newpos;
  size_t newpos2;

  string line;

  string tmpStr;

//just to test if the constraing function works
//int test = constraintFunction(0, "=", 0);
//cout << "test result is: " << test << endl; //the result should be 1

  //tests if input files exist
  if (!simResult1.is_open()){cout << "Error! Can't find simulation result file 1.\n"; return 1;}
  if (!simResult2.is_open()){cout << "Error! Can't find simulation result file 2.\n"; return 1;}

  //reads test constraints file
//  cout << "Reading constraints file..." << endl;


  for(auto it = constraints.begin(); it!=constraints.end(); ++it){
      //parse constraints using a logical expression parser
      //continue here
      //constraintValue1[0] = line.substr(0, line.find(constraints));
	  	  line = it->second;
        newpos = line.find_first_of(possibleConstraints, pos);
        constraintValue1.push_back(line.substr(pos, newpos));

        newpos2 = line.find_last_of(possibleConstraints);
        constraintValue2.push_back(line.substr(newpos2+1));

        //cout << "newpos " << newpos << " newpos2 " << newpos2 << endl;

        if(newpos == newpos2){
	constraintOp.push_back(line.substr(newpos,1));

        }else{constraintOp.push_back(line.substr(newpos,newpos2-newpos+1));}



        count++;


  }



  cout << "Done reading constraints file." << endl;







  //get variable/collumn names
  //cout << "Reading Result File 1..." << endl;
  getline (simResult1,line);
  variableNames1 = splitResults(line, ' ');

  if(variableNames1.size() <= 1){
	  variableNames1 = splitResults(line, '\t');
	  if(variableNames1.size() <= 1){
		cout << "Can't split result table, are you using space or tab as delimiter?" << endl;
		return 1;

	  }

  }

  for (int i = 0; i < variableNames1.size(); i++){
      if (variableNames1[i].empty() || variableNames1[i]=="#"){
         //cout << "Found empty string at " << i << endl;
         variableNames1.erase(variableNames1.begin()+i);
         i--;
      }
   }

/*
  //debug
  for(int i = 0; i < variableNames.size(); i++){
      cout << i << ": !!!!" << variableNames[i] << "@@@@" << endl;

  }
*/

  //get variable values
  count = 0;
  while ( getline (simResult1,line) ){

      //cout << line << '\n';

      variableValues1 =  splitResults(line, ' ');

      if(variableValues1.size() <= 1){
    	  variableValues1 = splitResults(line, '\t');
    	  if(variableValues1.size() <= 1){
    		cout << "Can't split result table, are you using space or tab as delimiter?" << endl;
    		return 1;

    	  }

      }

      for(int i = 0; i < variableValues1.size(); i++){
          if (variableValues1[i].empty()){
             //cout << "Found empty string at " << i << endl;
             variableValues1.erase(variableValues1.begin()+i);
             i--;
          }

       }

       //start making the map here
       //make a map with the variable names (string) vs their index (int), pointing to values (double)
       //example map[variableNames[i]][count] = variableValues[i];
       for(int i = 0; i < variableNames1.size(); i++){

          if(variableNames1[i] != "Iteration" || variableNames1[i] != "iteration" || variableNames1[i] != "time" || variableNames1[i] != "Time"){
             tmpStr = static_cast<string>(variableValues1[i]);
             namesValues1[variableNames1[i]][count] = atof(tmpStr.c_str());
             //cout << "tmpStr: " << tmpStr << endl;
             //cout << "name: " << variableNames[i] << "; value at " << count << ": " << namesValues[variableNames[i]][count]  << endl;
          }


       }


  count++;

  }


  //cout << "FOUND LINES IN THE RESULT FILE: " << count << endl;

  simResult1.close();


  //cout << "Done reading Result File 1." << endl;








  //get variable/collumn names
  //cout << "Reading Result File 2..." << endl;
  getline (simResult2,line);
  variableNames2 = splitResults(line, ' ');

  if(variableNames2.size() <= 1){
	  variableNames2 = splitResults(line, '\t');
	  if(variableNames2.size() <= 1){
		cout << "Can't split result table, are you using space or tab as delimiter?" << endl;
		return 1;

	  }

  }

  for (int i = 0; i < variableNames2.size(); i++){
      if (variableNames2[i].empty() || variableNames2[i]=="#"){
         //cout << "Found empty string at " << i << endl;
         variableNames2.erase(variableNames2.begin()+i);
         i--;
      }
   }

/*//debug
  for(int i = 0; i < variableNames.size(); i++){
      cout << i << ": !!!!" << variableNames[i] << "@@@@" << endl;

  }
*/
  //get variable values
  count = 0;
  while ( getline (simResult2,line) ){

      //cout << line << '\n';

      variableValues2 =  splitResults(line, ' ');

      if(variableValues2.size() <= 1){
    	  variableValues2 = splitResults(line, '\t');
    	  if(variableValues2.size() <= 1){
    		cout << "Can't split result table, are you using space or tab as delimiter?" << endl;
    		return 1;

    	  }

      }

      for(int i = 0; i < variableValues2.size(); i++){
          if (variableValues2[i].empty()){
             //cout << "Found empty string at " << i << endl;
             variableValues2.erase(variableValues2.begin()+i);
             i--;
          }

       }

       //start making the map here
       //make a map with the variable names (string) vs their index (int), pointing to values (double)
       //example map[variableNames[i]][count] = variableValues[i];

      for(int i = 0; i < variableNames2.size(); i++){

          if(variableNames2[i] != "Iteration" || variableNames2[i] != "iteration"){

             tmpStr = static_cast<string>(variableValues2[i]);
             namesValues2[variableNames2[i]][count] = atof(tmpStr.c_str());
            //cout << "tmpStr: " << tmpStr << endl;
            // cout << "name: " << variableNames[i] << "; value at " << count << ": " << namesValues[variableNames[i]][count]  << endl;
          }


      }


  count++;

  }


  simResult2.close();


  //cout << "Done reading Result File 2." << endl;

float constraintDist;
float result;
//added new model checking method
float fulfilledConstrants = 0;
int first = 0;
vector<string> consList;
//Now test if the constraints are fulfilled
 //And return a list of iterations or time points and whether they were fulfilled
 //cout << "Time" << "\t" << "Constraint" << "\t" << "Result(0=FALSE;1=TRUE)" << endl;
 //outFile << "Time" << "\t" << "Constraint" << "\t" << "Result(0=FALSE;1=TRUE)" << endl; //old mode
	outFile << "Time" << "\t" << "Constraints_fulfilled" << endl;
 for(int i = 0; i < constraintValue1.size(); i++){
     //constraintFunction(0, "=", 0);
     //cout << "@@@" << constraintValue1[i] << "@@@map size" << namesValues1[constraintValue1[i]].size() << endl;
 //    cout << constraintValue1[i] << "\t" <<  constraintOp[i] << "\t" << constraintValue2[i] << endl;
     cout << "Starting test for constraint: " << constraintValue1[i] << " " <<  constraintOp[i] << " " << constraintValue2[i] << " Iteration1 " << iteration1 << " Iteration2 " << iteration2 << endl;

     for (int j = 0; j < namesValues1[variableNames1[i]].size(); j++){
    	 //Raquel: added this if to check only the last time point
         //if(j==namesValues1[variableNames1[i]].size()-1){
    	 //if found the desired time point

    	 for(int k = 0; k < namesValues2[variableNames2[i]].size(); k++){
    		 //if we reach the correct time point/iteration specified by the user

    		 if(namesValues1[variableNames1[0]][j] == iteration1 && namesValues2[variableNames2[0]][k] == iteration2){
        	 	 cout << "FOUND TIME POINT 1= " << iteration1 << " TIME POINT2= " << iteration2 << endl;
        	 	 cout << "MODEL 1 " << inputFile1 << endl;
        	 	 cout << "MODEL 2 " << inputFile2 << endl;
        	 	cout << "Original constraint: " << constraintValue1[i] << " " << constraintOp[i] << " " << constraintValue2[i] << endl;

        	 	 //cout << "FOUND LAST TIME POINT j= " << j << " value " << namesValues1[variableNames1[0]][j] << endl;
        	 	 //cout << namesValues[constraintValue1[i]][j] << " " << constraintOp[i] << " " << namesValues[constraintValue2[i]][j];

        	 	 //if both constraint values belong to observables (are not explicit numbers)
        	 	 if(constraintValue1[i].find_first_not_of("-01234567890.") != string::npos && constraintValue2[i].find_first_not_of("-01234567890.") != string::npos){

        		 	 result = constraintFunction(namesValues1[constraintValue1[i]][j], constraintOp[i], namesValues2[constraintValue2[i]][j]);
            	 	 cout << "Constraint values: " << namesValues1[constraintValue1[i]][j] << " " << constraintOp[i] << " " << namesValues2[constraintValue2[i]][j] << endl;
            	 	 // (expected_value/(constraint_diff+expected_value))
            	 	 if(result==0){
            	 		 constraintDist = sqrt( pow( ( namesValues2[constraintValue2[i]][j]-namesValues1[constraintValue1[i]][j] ), 2) );
            	 		 result=( namesValues2[constraintValue2[i]][j] / (constraintDist + namesValues2[constraintValue2[i]][j]) );

            	 	 }


        	 	 }else if(constraintValue1[i].find_first_not_of("-01234567890.") != string::npos && constraintValue2[i].find_first_not_of("-01234567890.") == string::npos){

        	 		 result = constraintFunction(namesValues1[constraintValue1[i]][j], constraintOp[i], atof(constraintValue2[i].c_str()));
            	 	 cout << "Constraint values: " << namesValues1[constraintValue1[i]][j] << " " << constraintOp[i] << " " << atof(constraintValue2[i].c_str()) << endl;
            	 	 if(result==0){
            	 		 constraintDist = sqrt( pow( ( atof(constraintValue2[i].c_str())-namesValues1[constraintValue1[i]][j] ), 2) );
            	 		 result=( atof(constraintValue2[i].c_str()) / (constraintDist + atof(constraintValue2[i].c_str()) ) );

            	 	 }
        	 	 }else if(constraintValue1[i].find_first_not_of("-01234567890.") == string::npos && constraintValue2[i].find_first_not_of("-01234567890.") != string::npos){

        		 	 result = constraintFunction(atof(constraintValue1[i].c_str()), constraintOp[i], namesValues2[constraintValue2[i]][j]);
            	 	 cout << "Constraint values: " << atof(constraintValue1[i].c_str()) << " " << constraintOp[i] << " " << namesValues2[constraintValue2[i]][j] << endl;

            	 	 if(result==0){
            	 		 constraintDist = sqrt( pow( ( namesValues2[constraintValue2[i]][j]-atof(constraintValue1[i].c_str()) ), 2) );
            	 		 result=( namesValues2[constraintValue2[i]][j] / (constraintDist + namesValues2[constraintValue2[i]][j]) );

            	 	 }
        	 	 }else if(constraintValue1[i].find_first_not_of("-01234567890.") == string::npos && constraintValue2[i].find_first_not_of("-01234567890.") == string::npos){

        		 	 cout << "ERROR can't find any constraint value in the observable list!!!" << endl;
        		 	 cout << "constraint value 1: " << constraintValue1[i] << "constraint value 2: " << constraintValue2[i] << endl;

        	 	 }

        	 	 cout << "RESULT=" << result << endl;
        	 //cout << " = " << result << "; iteration = " << j << " constraints = " << constraintValue1[i] << constraintOp[i] << constraintValue2[i] << endl;

        	 	 if(result<=1 && result>0){

        	 		 fulfilledConstrants+=result;
        	 		 if(result==1){
            	 		 consList.push_back(constraintValue1[i]+constraintOp[i]+constraintValue2[i]);
        	 		 }
        	 		 result = 0;
        	 		 constraintDist = 0;
        	 	 }else{
        	 		 cout << "ERROR invalid constraint result: " << result << endl;
        	 		 result = 0;
        	 		 constraintDist=0;
        	 	 }
        	 	 if( (variableNames1[0]=="time" || variableNames1[0]=="Time") && first==0){
        	 		 //cout << namesValues1[variableNames1[0]][j] << "\t";
        	 		 outFile << namesValues1[variableNames1[0]][j] << "\t";
        	 		 first = 1;
        	 	 }else{
        	 		 if(first==0){
        	 			 outFile << j+1 << "\t";
        		 	 	 first = 1;
        	 		 }
        	 		 //cout << j+1 << "\t";

        	 	 }
        	 	 //cout << constraintValue1[i] << constraintOp[i] << constraintValue2[i] << "\t";
        	 	 //cout << result << endl; //old mode
        	 	 //outFile << constraintValue1[i] << constraintOp[i] << constraintValue2[i] << "\t"; //old mode
        	 	 //outFile << result << endl; //old mode




         	 }


    	 }

     }




 }

 outFile << fulfilledConstrants << "\t";

 for (int i = 0; i<consList.size(); i++){

	 outFile << consList[i] << "\t";


 }
 	 outFile << endl;

 	consList.clear();


 outFile.close();


//namesValues[variableNames[i]][count]

  return fulfilledConstrants;
}
