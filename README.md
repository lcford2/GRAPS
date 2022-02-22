# GRAPS
Generalized Reservoir Analyses using Probabilistic Streamflow
COREGS branch: Model formulation used with the [COREGS model](https://github.com/lcford2/coregs)

## Binaries and Compilation

To compile GRAPS, you will need the intel oneAPI HPC toolkit, which is dependent on the oneAPI base toolkit. 
You can find information on installing these toolkits at the [Intel oneAPI webpage](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.pvef6v). 
You will also need the [GNU Make Utility](https://www.gnu.org/software/make/).
This will allow you to easily compile GRAPS on a Linux machine.
If you are running Mac or Windows, it will be easiest to use the docker image to run this code [ADD DOCKER IMAGE].
If you would rather compile from source on Mac or Windows, you can follow the directions to get the oneAPI kits, and then use the makefile at `src/makefile` as a guide for compiling the model.

After you have installed the toolkits and sourced the setup script (e.g., `source /opt/intel/oneapi/setvars.sh` if you installed the toolkits under the `/opt` directory), you can compile GRAPS by running `make` from the command line in the root directory of this project.. 
This will compile a shared library and store it at `lib/graps.so` and an executable and store it at `bin/graps`. 
COREGS relies on the shared library and does not use the executable but it is generated for your convenience.

## Graphical User Interface

The complementary GRAPS GUI can be found at [github.com/lcford2/graps_gui](https://github.com/lcford2/graps_gui).

In the future, these two projects will be linked and released together but for the moment they exist separately. 

