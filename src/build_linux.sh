#!/bin/bash

# Make a build directory and copy source code to that directory
mkdir build
cp *.f* build
cd build
# Compiling here, check and traceback are here to assist the user in case something goes wrong
# They do not slow down the code, but the executable is larger than it would be without them. 
ifort -check all -traceback definitions.f90 My_variables.f90 rest_subroutines.f90 ffsqp.f qld.f all_simul.f90 multireservoir.f90 -o multireservoir
# Remove all source code files from the build directory
rm *.f*
# Move the executable and object files to the linux directory
mv * ../../linux
cd ..
# Remove the empty build directory
rmdir build
