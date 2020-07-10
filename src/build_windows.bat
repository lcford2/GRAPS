@ECHO OFF

rem make a build directory and copy source code to that directory
mkdir build
xcopy *.f* build
cd build
rem Compiling here, check and traceback are here to assist the user if something goes wrong.
rem They do not slow down the program, but it is larger than if it is compiled without those flags.
ifort /check:all /traceback definitions.f90 My_variables.f90 ffsqp.f qld.f multireservoir.f90 all_simul.f90 rest_subroutines.f90 /exe:multireservoir
rem Remove all source code files from the build directory
del *.f*
rem Move the executable and object files to the windows directory
move * ../../windows
cd ..
rem Remove the empty build directory
rmdir build
