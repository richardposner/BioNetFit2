echo "This is a shell script to link the object files"  
cd bin
rm -f BioNetFit.exe
mpic++ -L/usr/lib -o "BioNetFit"  ./src/code/Model.o ./src/code/Data.o ./src/code/FreeParam.o ./src/code/Model.o ./src/code/Parser.o ./src/code/Particle.o ./src/code/Pheromones.o ./src/code/Swarm.o ./src/code/Model.o  ./src/GenFit2.o  ../lib/libboost_iostreams.a ../lib/libboost_regex.a ../lib/libboost_mpi.a ../lib/libboost_program_options.a ../lib/libboost_filesystem.a ../lib/libboost_system.a ../lib/libboost_serialization.a -lrt
./BioNetFit.exe -a run -t master -c ../examples/ex5/ex5.conf -v
cd ..
echo "Linking is finished"  