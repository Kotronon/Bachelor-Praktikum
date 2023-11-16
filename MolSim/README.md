MolSim
===

The Molecular Dynamics teaching code.

# Group H #
Members:
* Samhitha Girish Jois
* Katharina Miller
* Anna Lena Müller

# Code #
* Link:     https://github.com/Kotronon/Bachelor-Praktikum/tree/main/MolSim
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

## Task 3 ##
* logging
* integration of spdlog in cmake through fetchContent
* log level can be chosen either 
  * via input with -level 
  * via cmake (either via CMAKE_MESSAGE_LOG_LEVEL or SPDLOG_ACTIVE_LEVEL)
  * of via command cmake --log-level=<level> or cmake --loglevel=<level>

## Task 4 ##
* Collision of two bodies
* 