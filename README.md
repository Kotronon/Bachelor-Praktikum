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
* XML file

## Task 2 ##
* Linked-cell algorithm
* Cell wide in each direction is the cutoff value
* Structure as image below:
* We still provide the possability to use simple sum implementation
* Between the two implementations is a big time difference:
![Screenshot](input/both_rpunded_and_scaled.png)

## Task 3 ##
* Boundaries
* Outflow Boundary;
  * checks if particle is out of boundary
  * delets if it is  the case
* Reflection Boundary
  * Mirroring:
    * just reflecting the particles like a mirror 
    * depending on surface and kind of particle physically inaccurate
  * Ghost cells
    * creating an imaginary cell to create a force from the boundary
    * has an individual influence on each particle
    * physically more accurate
    * made as reflection boundary in the end version

## Task 4 ##
* sphere

## MISC ##