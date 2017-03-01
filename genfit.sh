echo "This is a shell script to compile only one file"  
cd bin
rm -f GenFit2.o GenFit2.d
mpic++ -fpermissive -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -I../include -I/usr/lib/openmpi/include -O3 -g3 -Wall -c -fmessage-length=0 -std=gnu++0x -MMD -MP -MF"src/GenFit2.d" -MT"src/GenFit2.d" -o "src/GenFit2.o" "../src/GenFit2.cpp"
cd ..
echo "Operation is finished"  
