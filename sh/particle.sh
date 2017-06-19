echo "This is a shell script to compile only one file"  
cd bin
rm -f src/code/Particle.o src/code/Particle.d
mpic++ -fpermissive -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -I../include -I/usr/lib/openmpi/include -O3 -g3 -Wall -c -fmessage-length=0 -std=gnu++0x -MMD -MP -MF"src/code/Particle.d" -MT"src/code/Particle.d" -o "src/code/Particle.o" "../src/code/Particle.cpp"
cd ..
echo "Operation is finished"  
