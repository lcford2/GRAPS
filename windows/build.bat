@ECHO OFF

cd ../source
mkdir build
copy *.f* build
cd build

ifort /check:all /traceback definitions.f90 My_variables.f90 ffsqp.f qld.f multireservoir.f90 all_simul.f90 rest_subroutines.f90 /exe:multireservoir
del *.f*

move * ../../windows
cd ..
rmdir build
cd ../windows
