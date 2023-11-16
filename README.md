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
* Revision: 36ff4bf
* Compiler: gcc 11.4.0

# Run Code #
* packages to install:
  * sudo apt install libspdlog-dev
  * sudo apt-get install libgtest-dev
  * unfortunatly only using fetchContent doesn't prevent from installing them
* for the program:
  * ./MolSim <end_time> <delta_t> 
  * options:
  * -level <level> (To choose log level)
  * -f <path_to_file> (To dd particles manually without cuboids)
  * -c <numbers_of_cuboids> (for each cuboid) <coordinates_of_left_corner_cuboid_i> <velocity_of_cuboid_i> <dimension_of_cuboid_i> <h_of_cuboid_i> <mass_of_cuboid_i> 
  * for arrays please use the form x,y,z
  * for example: ./MolSim 1000 0.014 -c 2 0,0,0 0,0,0 40,8,1 1.1225 1 15,15,0 0,-10,0 8,8,1 1.1225 1
* for the tests:
  * ctest
  

# Report #
## Task 1 ##
* Created unit tests with googletest
* for each important function there is one test
* the calculation of velocity with BrownianMotionInitialization couldn't be tested due to random values
* integration of googletest in cmake through fetchContent 

## Task 2 ##
* CI integration in GitHub
* creating yml file and learning about github actions
* Learning YML syntax and adjusting it to our repository structure
* adjust the yml file to trigger a workflow on push and pull requests
* build cmake and run tests via the yml file
* 

## Task 3 ##
* logging
* integration of spdlog in cmake through fetchContent
* log level can be chosen either 
  * via input with -level 
  * via cmake (either via CMAKE_MESSAGE_LOG_LEVEL or SPDLOG_ACTIVE_LEVEL)
  * of via command cmake --log-level=<level> or cmake --loglevel=<level>

## Task 4 ##
* Collision of two bodies
* Calculation of force according to Lennard-Jones-Potential
* Normal calculation with iteration through all particles
* faster calculation with pairwise iteration
* time difference:
* calculation method:
* video:

## MISK ##
* FetchContent is not enough to avoid installing libraries -> for next time better integration is needed
* changing the folder structure was quite difficult
* doxygen was refactored. Now it only starts with make doc_doxygen and will immediately open the documentation in firefox