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
  * ./MolSim <end_time> <delta_t> 
  * options:
  * -level <level> (To choose log level)
  * -f <path_to_file> (To add particles manually)
  * -c <numbers_of_cuboids> (for each cuboid) <coordinates_of_left_corner_cuboid_i> <velocity_of_cuboid_i> <dimension_of_cuboid_i> <h_of_cuboid_i> <mass_of_cuboid_i> 
  * for arrays please use the form x,y,z
  * for example: ./MolSim 5 0.0002 -c 2 0,0,0 0,0,0 40,8,1 1.1225 1 15,15,0 0,-10,0 8,8,1 1.1225 1
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