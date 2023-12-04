MolSim
===

The Molecular Dynamics teaching code.

# Group H #
Members:
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
  * ./MolSim <xml_file> 
  * example files can be found in input folder
* to run the tests:
  * ctest
  

# Report #
## Task 1 ##
* XML file
* Due to time pressure, we orientated on the xsd file structure of https://github.com/Dominik-Weinzierl/MolSim/tree/main/src/fileReader/XMLReader/template 
and the XMLReader structure of https://github.com/wngTn/MolSim/tree/main/src/inputReader
* We strore everything in an XMLInfo struct and can get access to them in the main file

## Task 2 ##
* Linked-cell algorithm
* Cell wide in each direction is the cutoff value
* Each cell linked to an exact cell coordinate (vector<vector<vector<vector<Particle>>>> in order x, y, z)
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
    * is mirrored behind the boundary (both particles have same distance to boundary)
    * has an individual influence on each particle
    * physically more accurate
    * made as reflection boundary in the end version

## Task 4 ##
* sphere

## MISC ##