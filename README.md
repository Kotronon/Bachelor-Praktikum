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
  * ./MolSim
  * due to the xml input being incomplete it is only possible at the moment to change
  * parameters in the main method manually as well as create cuboids/disks there
  *
  * 
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