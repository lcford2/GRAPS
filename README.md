# GRAPS
Generalized Reservoir Analyses using Probabilistic Streamflow

## Binaries

You can find precompiled binaries for 64-bit architechtures for Linux and Windows under their respective directories. Included with the reservoir model is the Feasible Sequential Quadratic Programming (FSQP) optimization algorithim.

## Model Example

An example model is presented under the `brazil_example` folder. This example is setup for the Jaguaribe River Basin in Ceara, Brazil for 1997. To run this example you will need to copy all the files from the folder containing the binaries for your operating system (`windows` or `linux`) to the `brazil_example` folder. You can follow the instructions below to decide how you want the model to run. Once you are satisfied with the setup, you should execute `./multireservoir` from within the `brazil_example` folder from the command line. 

### Different model modes

There are four ways to run the model. You can change how the model will run by changing a few lines in the `runflag.dat`, `watershed_details.dat`, and `input.dat` files. Outlined below are instructions on how to run the model in different modes.

1. Single Trace Simulation
    - **This is how the model is currently setup to run**. 
    - The `runflag.dat` contains a single `0` in the first line. Only the first line will be read so anything else in the file does not matter. 
    - The first line in the `input.dat` file is `12  1  1`
    - All of the specified inflow files in `watershed_details.dat` are of the form `<reservoir_name>_97_inflow.dat` where `reservoir_name` is replaced with each reservoir's name. 
2. Inflow Ensemble Simulation
    - The `runflag.dat` file is the same for Single Trace Simulation
    - The first line in the `input.dat` file is `12  1  100`
    - All of the specified inflow files in `watershed_details.dat` are of the form `<reservoir_name>_ensemble.dat` where `reservoir_name` is replaced with each reservoirs name. 
3. Single Trace Optimization (using FSQP)
   - Everything is setup exactly like Single Trace Simulation except the `runflag.dat` file, which will now contain a single `1` on the first line and the name of the objective function on the second line. Possible objective functions are outlined below. 
4. Inflow Ensemble Optimization (using FSQP)
   - Everything is setup exactly like Ensemble Simulation except the `runflag.dat` file, which will now contain a single `1` on the first line and the name of the objective function on the second line. Possible objective functions are outlined below. 

### Objective Function Choices
The following objective functions are currently implemented and their subroutines can be found in the `all_simul.f90` file. The name in parentheses is the subroutine name and is what should be on line 2 of `runflag.dat` when optimizing. You may also add your own objective functions using the currently implemented ones as templates. If you write your own objective function, you will need to update the conditional logic in `multireservoir.f90` to ensure that the program will use your objective function. 
1. Maximize Expected Benefits (`expected_benefits`)
   - This uses the full water contracts as laid out in the input files. It contains disbenefits for not meeting targets. 
2. Maximize Simple Benefits (`simple_benefits`)
   - This simply uses the water tariffs as coefficients for release for each user with out considering disbenefits for not meeting targets.
3. Maximize Release (`max_release`)
   - This maximizes the raw sum of release. Can be considered an unweighted optimization because it does not favor one user over an over. 
4. Minimize Spill (`spill_objective`)
   - This simply tries to minimize the spill across the network without any regard for system benefits. This is mainly implemented to see what the lower bound on spill can be for a particular scenario. 



### Use of `decisionvar_details.dat` and of `release.out` files

When the model is ran in simulation mode, `decisionvar_details.dat` contains the prescribed release for each user in the system and `release.out` will simply contain the same values as `decisionvar_details.dat`. If the system is being optimized, `decisionvar_details.dat` will contain the initial solution or initial guess for FSQP. In this case, `release.out` will contain the optimal release for each user as determined by FSQP. Both files are ordered in ascending order with respect to user ID and time step. So for a model with 4 users and 12 timesteps the length of each file will be 48 (4 x 12) and the first 12 values will correspond to user 1 and the next 12 to user 2 and so on. 
