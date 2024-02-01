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
* Branch:   Aufgabe_1_Samhitha
* Revision: e80650f
* Compiler: gcc 11.4.0

# Run Code #
* to run the program:
  * remmove any existing build folder with rm -rf build
  * mkdir build ; cd build
  * cmake ..
  * make 
  * ./MolSim
* to run the tests:
  * ctest
  

# Report #
## Task 1 ##
* Membrane
* Implemented given formulae for the force calculations
* Implemented particle-specific application and then container-specific force application
* The force parallel to Z-Axis is also implemented 
* Iteration order : "That force" -> Lennard Jones -> Membrane interaction force
* Position and Velocity calculated with Stoermer-Verlet 
* Hardcoding of all the parameters needed (occasionally also in the function calls)
* ghost cell check and calculation (+deletion) needed for all the force applications (extra run time)
* main function iteration took around 46 minutes 
* Simulation with 24 fps took around 3 hrs 20 minutes
* 
## Task 2 pt 2 : Parallelisation strategy ##
* Profiling : ForceCalculation takes the most time 
* strategy : parallelise not the outer 3, but inner 3 iterations and focus on parallelising the inner force calculations and distributions
* OpenMP cannot run on a student account LRZ Rechnerhalle, which is my way of running the code 
* 

