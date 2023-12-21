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
  

# Report for Samhitha's Branch #
* Link : https://github.com/Kotronon/Bachelor-Praktikum/tree/Aufgabe-1_XML_input_and_XSD
* latest Revision : 
## Task 4 - Profiling ##
* 1: Added following timestamps in MolSim.cpp's main function : setup time, loop time (entire), single iteration time, position+velocity+force calculation time
* Added following measurements : average iterations per second, average molecule/cells update per second
* 2 : Linux cluster : Unfinished due to to login problems
* 3: Profiling : 
* * To be done with AMD uprof (my processor is an AMD ryzen) : https://www.amd.com/en/developer/uprof.html
* * Downloaded the Windows GUI for better statistics
* * download following package : mingw-64 version 8.0.0-1 
* https://www.mingw-w64.org/downloads/ --> https://launchpad.net/ubuntu/+source/mingw-w64
* * use command "x86_64-w64-mingw32-g++ <source file> -o <target executable name>.exe" instead of normal g++ to generate windows-executable .exe
* * paste source path of .exe to AMD uprof main page, and run analysis

* Execution of this roadmap plan failed due to compilation problems

## Task 5 - Optimisation ##
* Biggest potential for optimisation : 
* Following algorithms with O(n^4) runtime :
* --> all functions that plot the cells in the coordinate system
* in calculator and Thermostat classes
* Basic structure of the loops : 
for (x : amount of cells){
for (y : x){
 for (z : y){
 for (p : z(){
 --- initialise particle ---
}
}
}
 }

* Optimisation suggestion (to proceed with caution due to possible integer overflows!): 
* for (int x = 0, int y = 0, int z = 0, x = amount of cells , x++, y++, z++){
 x = vector<vector<vector>>(initialise x);  y = vector<vector>(initialise y), z = <vector>(initialise z);
 initialise particle()
 }
* One loop will take less space, but the repeated initialisations per iteration will take up a lot of overhead
* Therefore : Parallel processing needed

