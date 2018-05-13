/*============================================================================
// Name        : Setting.hh
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Lat Update: : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/

#ifndef CODE_SETTING_HH_
#define CODE_SETTING_HH_



//razi: choose compile options





//razi: this includes all changes made on the version delivered by Brandon. Main changes include
// fixing compile problems in windows
// fixing run time problems in windows
// modifying swarm object to support multiple models [each model has a different set of free params and exp files]
// developing subparticle object to support multiple connected models
// message passing based on subparticle id (particle id + model index)
// consolidate free parameters among models and check consistency between exp ad model files






#define TEST_SIMULATOR   //uses test simulator for fast run, only for test  // should be commented in the final version
//#define STANDALONE_PARTICLE_TEST_ACTIVE //this is to sets particle stand-alone// should be commented in the final version


//#define CYGWIN   //some special changes to compile in CygWin
#ifndef  CYGWIN
	#if(defined __unix || defined __unix__ || __linux__ || __MACH__) //razi added for test
		#define LIN_VER   //alias for macros
	#endif
	#if(defined _WIN32 || defined  _WIN64) //razi added for test
		#define WIN_VER
	#endif
#endif

#if (defined(CYGWIN) || defined(WIN_VER))
	#define PC_VER
#endif



#endif /* CODE_SETTING_HH_ */
