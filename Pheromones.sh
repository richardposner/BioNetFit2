echo "This is a shell script to compile only one file"  
cd bin
rm -f src/code/Pheromones.o src/code/Pheromones.d
mpic++ -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -I../include -I/usr/lib/openmpi/include -O3 -g3 -Wall -c -fmessage-length=0 -std=gnu++0x -MMD -MP -MF"src/code/Pheromones.d" -MT"src/code/Pheromones.d" -o "src/code/Pheromones.o" "../src/code/Pheromones.cpp"
cd ..
echo "Operation is finished"  