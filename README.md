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
* Revision: bc2d47b
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
### Unfinished part of the assignment : 
* github actions did not run due to following errors :
  * Could NOT find Doxygen (missing: DOXYGEN_EXECUTABLE)
  * CMake Error at CMakeLists.txt:113 (find_package):
  * Could not find a package configuration file provided by "spdlog" with any of the following names: spdlogConfig.cmake spdlog-config.cmake
* the following command : cmake -B .github/build -DCMAKE_BUILD_TYPE=${{env.BUILD_TYPE}} was supposed to :
  *  configure cmake in a seperate build file that does not interrupt the project structure
* Earlier versions of the github-actions.yml file gave error messages such as :
  * "Error : cannot load cache"
  * "make: *** No rule to make target"
* every attempt to fix the yml file resulted in the above-mentioned error messages
* Solution attempts : 
  * make github install doxygen and spdlog into the .github/build directory 
  * circumvent the errors by ignoring them
  * attempt to only build src/tests/MolSim directory
  * make github find the CMakeLists and run cmake .. in the cmake directory
* None the of the above-mentioned solution attempts ran, so we are stuck at this point.

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
* Calculation with simple double iteration through all particles
* New input options needed to allow for both file input and multiple cuboids
* Faster calculation with new pairwise iteration method in particle container using Newtons third law to avoid double calculations
* specs: 
  * System:
    Kernel: 5.15.0-88-generic x86_64
    Desktop: Gnome 3.38.4
    Distro: Zorin OS 16.3
    base: Ubuntu 20.04 LTS Focal 
  * CPU:
    Topology: 8-Core model: AMD Ryzen 7 5825U with Radeon Graphics
* conditions: measuring of time with functions in ctime library, setting loglevel to off and disabling/not including time for I/O operations,
  Simulation used is the collision of two bodies as described in the worksheet
* time difference: 18 min 36s (1116s) with simple calculation, 12 min 33 s (753s) with faster calculation
* => very noticeable speed up, faster loop is nearly 1.5 times faster
* video: https://youtu.be/bzQOXPaK2VI

## MISC ##
* FetchContent is not enough to avoid installing libraries -> for next time better integration is needed
* changing the folder structure (root directory in git) was quite difficult
* doxygen was refactored. Now it only starts with make doc_doxygen and will immediately open the documentation in firefox
* Refactored most code from worksheet 1 and implemented feedback
* Added more missing documentation