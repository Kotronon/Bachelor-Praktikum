MolSim
===

The Molecular Dynamics teaching code.

# Group H #
Members:
* Samhitha Girish Jois
* Katharina Miller
* Anna Lena Müller

# Code #
* Link:     https://github.com/Kotronon/Bachelor-Praktikum
* Branch:   master
* Revision: e80650f
* Compiler: gcc 11.4.0

# Run Code #
* packages to install:
  * sudo apt install libspdlog-dev
  * sudo apt-get install libgtest-dev
  * unfortunately only using fetchContent doesn't prevent from installing them
* to run the program:
  * ./MolSim
  * due to the xml input being incomplete it is only possible at the moment to change parameters in the main method manually as well as create cuboids/disks there
  * this should be fixed in the future, log level can not be chosen either currently as this was planned to be set up completely in the xml files
  * 
* to run the tests:
  * ctest
  

# Report #
## Task 1 ##
* Current temperature / kinetic energy
  * To determine current temperature of a system, the given formulas can be combined to
  * T = (sum from i = 1 to #particles (m_i * <v_i,v_i>)) / #dimensions * #particles
* Initialization
  * There is still an option to choose whether to use Brownian Motion for the initialization or not
  * However, it needs to be used when initial velocities are zero otherwise scaling is just applied to the already existing velocities
  * initTemperature and applyBrownianMotion are new parameters
* Direct temperature setting
  * Temperature can be set to an exact temperature in each time step by either setting a specific target temperature or by setting targetTemperatureExists to false when the initial temperature should be used
  * velocity scaling only needs to calculate the current temperature and can then use it to calculate the scaling factor for all velocities
* Gradual temperature setting
  * By using differenceTemperature you can choose the maximal change in temperature after each time step
  * Gradual setting is only enabled if differenceTemperatureExists is set to true
  * It will then gradually try to set the temperature to the target temperature 
 
## Task 2 ##
* Periodic boundary
  * ghost particles to be generated at the opposite side it needed to be mirrored
  * also ghost particles generated in corners
* Gravity force
  * added at the end of each force calculation per particle on the force y-axis
* Individual epsilon and sigma for each particle is possible
  * they are now parameters of an individual particle itself
* Issues when performing the small experiment
  * We encountered some rather weird issue when mixing boundary types
  * For simulations where all (2d) boundaries where of the same type everything worked fine
  * In the case in which the horizontal boundaries were reflective and the vertical boundaries were periodical, some particles started to disappear in the course of the simulation
  * For some reason the amount of particles that were deleted varied heavily between two different computers even though we were on the exact same commit and running the same experiment
  * While only less than 10 particles were deleted on one machine on the other nearly all except about 400 were deleted when running the small experiment for task 2
  * After a certain number of iterations some calculation massively inflates the velocities of some particles, this seems to be most likely caused by some issues with the periodic boundary in combination with the reflective
  * We are not sure if we accidentally created undefined behaviour in any new function to cause this or if there is another explanation
  * Data for the Rayleigh-Taylor instability experiments shown in the videos were calculated on the pc with fewer issues


## Task 3 ##
* Checkpointing
  * individual numbers of checkpoints
  * saving all parameters of all particle in each checkpoint file
  * structure of input file is similar to the input file used in worksheet 1 with added parameters
* Drop
  * used last checkpoint file as input file to equilibrate fluid at the bottom first
  * let drop fall on fluid without temperature regulation
  * simulation is gradually heating up due to gravitational force adding more energy each time step

## Task 4 ##

## Task 5 ##

## MISC ##
