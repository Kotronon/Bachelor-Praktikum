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
* Revision: e7321c1
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
* Thermostat
 
## Task 2 ##
* periodic boundary
  * ghost particles generating at each other side it needed to be mirrored
  * also ghost particles generated in corners
* Gravity force
  * added at the end of each force calculation per particle on the force y-achsis
* Individual epsilon and sigma for each particle is possible
  * they are now parameters of particle itself

## Task 3 ##
* checkpointing
  * individual numbers of checkpoints
  * saving all parameters of all particle in each checkpoint file
* drop
  * used last checkpoint file as input file
  * let drop file on fluid


## MISC ##
* In the falling drop animation we noticed that particles tend to sometimes get stuck at the upper boundary, we did not yet manage to fix this error
* Due to some restructuring of our LinkedCellContainer we might have forgotten to change some functions accordingly