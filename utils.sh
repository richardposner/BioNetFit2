echo "This is a shell script to compile only one file"  
cd bin
rm -f src/code/Utils.o src/code/Utils.d
mpic++ -fpermissive -D__GXX_EXPERIMENTAL_CXX0X__ -D__cplusplus=201103L -I../include -I/usr/lib/openmpi/include -O3 -g3 -Wall -c -fmessage-length=0 -std=gnu++0x -MMD -MP -MF"src/code/Utils.d" -MT"src/code/Utils.d" -o "src/code/Utils.o" "../src/code/Utils.cpp"
cd ..
echo "Operation is finished"  
