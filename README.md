# BioNetFit2

# Operating systems supported

	Linux, tested on Ubuntu and Red Hat
	Windows, tested in Windows 7 and Windows 10 with 64bit cygwin
	MacOS, tested with High Sierra 10.13.3
	
# Requirements

	gcc compiler (GCC), v5.2.0 and v5.4.0 tested and working
	mpicc, compiled with the required version of gcc
	boost library v1.65.0, included and pre-compiled for Linux
	libc6 or glibc 2.14, not requered for installation, but requered for running

# Installation at NAU's cluster (Red Hat 7 GNU/Linux with GCC v5.2.0)

	module add gcc/5.2.0
	module add openmpi/1.8.7-gcc-5.2.0
	cd BioNetFit2/boost_1.65.0
	./install_boost.sh
	cd ..
	make clean
	make

	before running BioNetFit
	module load glibc/2.14

*And add option cluster_software=BNF2mpi in your config file

# Local installation in Linux (Ubuntu 14.04 and Ubuntu 16.04 with GCC v5.4.0

	cd BioNetFit2
	make clean
	make

# Installation on MacOS
1. Install Homebrew package manager if it is not already present (https;//brew.sh)
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

2. Install Open-MPI
	brew install open-mpi
	
3. Compile boost and BioNetFit2
	make clean boost macos_install


# New Features in the .conf file

The implementation for multiple models and model checking is complete for the GA, PSO, and DE. So they can be used as templates to update the SA algorithm.
I have also updated the model checking: now the user can provide the exact time point where the logical operation must be done, and the user can also choose which models must be compared.

Now constraints should be provided in .conf files in the following format:

      constraint=Param1[>,>=,==,<=,<]Param2 MODEL1 MODEL2 TIME1 TIME2

For example:

      constraint=RLbonds>=pR 0 1 0 60

where I'm comparing RLbonds from model 0 at time 0 versus pR from model 1 at time 60, and I'm checking if RLbonds is greater then or equal to pR



Another new option was included so the user can specify the constraint weight, which can be any value between 0 and 1.
A very small number (i.e 0.1) = fit values matter the most. 
A very large number (i.e 0.9) = constraints matter the most.
A balanced value (i.e 0.5) = both fit and contrsaint values are equally important.

For example:

      constraint_weight=0.5



# Other differences compared to BioNetFit1

- If using NFsim, whenever you simulate observables that come from the function section of the BNGL file, don't forget to add the "()" suffix to observable labels in the experimental data file. For example, if you have the following function section:

      begin functions

           pre1_dose()=alpha1_pre*Clusters/f # .scan file output (Fig. 2B)
           pre2_time()=alpha2_pre*Clusters/f # .gdat file output (Fig. 3B)
           pre3_dose()=alpha3_pre*pEGFR/f # .scan file output (Fig. 2D)
           pre4_time()=alpha4_pre*pEGFR/f # .gdat file output (Fig. 3D)
           ...

      end functions


- Your experimental data file (*.exp) should contain the function suffix "()" in the obervables that resulted from functions:

      #	time	pre2_time()	pre4_time()	pre2_time()_SD	pre4_time()_SD
      	0	1.6131558314	0.0015707048	0.2004231309	0.0597428571
      	30	0.9593315136	0.8850918303	0.2615617574	0.0667714286
      	60	0.8707928038	1.11127324	0.163062846	0.0632571429
      	120	NaN	1.2208298593	NaN	0.0702857143
      ...


- Now BioNetFit2 supports multiple models or datasets per run (models comparing wild-type versus mutants, for example). To make things clear, you must specify one experimental file per model in the conf file. For example:

      model = parabolaA.bngl, parabolaB.bngl
      exp_file = parabolaA.exp, parabolaB.exp

where parabolaA.exp is the experimental file for parabolaA.bngl, and parabolaB.exp is the experimental file for the model parabolaB.bngl

- Precision: BNF2 allows 5 decimal characters per value in the experimental file (.exp) or time/parameter_scan data points, the obervable values can have more precision), for example: 

        #time/X_value    Y_observable
        1.00372E-04	6.0310121E-02
        3.14492E-04	2.4841445E-01
        9.85384E-04	3.2585751E-01
        3.18698E-03	3.8716871E-01
        9.98564E-03	5.0867287E-01
        3.12876E-02	4.8708839E-01
        1.01192E-01	5.2860091E-01

If you use too many decimal characters the mapping of the simulated data versus experimental data may not match because of precision limitations. 
Make sure that the simulated parameter_scan points or time points have the same number of decimal characters as the experimental data file as well, otherwise the program may return an "map at out_of_range" error.
