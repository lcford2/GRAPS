# GRAPS
Generalized Reservoir Analyses using Probabilistic Streamflow

You can find precompiled binaries for 64-bit architechtures for Linux and Windows under their respective directories. 

An example model is presented under the `brazil_example` folder. This example is setup to simulate the Jaguaribe River Basin in Ceara, Brazil for 1998. 
In this case, the releases for each user are prescribed. To run this example using Feasible Sequential Quadratic Programming (FSQP) to optimize the system, change the `0` in `runflag.dat` to a `1`. When optimizing, the contents of `decisionvar_details.dat` are used as an initial solution for FSQP rather than the prescribed releases as in simulation. 