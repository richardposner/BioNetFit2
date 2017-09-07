#!/bin/bash

chmod +x ./bootstrap.sh
./bootstrap.sh --prefix=../ --with-libraries=filesystem,iostreams,mpi,regex,system,program_options,serialization
echo "using mpi ;" >> project-config.jam
./b2 install

