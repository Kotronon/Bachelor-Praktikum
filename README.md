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
* XML input with XSD
* Using CodeSynthesis XSD tool, generated files that parse a XML file with the help of a XSD file
* generated using this command : 
* xsd cxx-parser --generate-test-driver --generate-print-impl --force-overwrite --xml-parser expat <xsd input file>.xsd
* adding destructors manually to input-pimpl.cxx and hxx file
* linked inputnew-pimpl file to main running file MolSim.cpp
* Run the parsing with following commands : 
* c++ -std=c++11 -I.../libxsd -c newinput-driver.cxx newinput-pskel.cxx
* c++ -std=c++11 -o driver newinput-driver.o newinput-pskel.o -lexpat
* ./driver newinput.xml
* latest working version can be found on branch : origin/Aufgabe-1_XML_input_and_XSD and commit : f42de68d
* Attempt/ Idea : 
* 1. create class-specific variables in parameter classes specified in .hxx , implement simple getters to those classes in .cxx file
* 2. call getter functions in a function in MolSim.cpp and instantiate all parameters as variables
* 3. use variable definitions instead of hardcode values in function
* Linking of MolSim.cpp (main executable) to newinput-pimpl.cxx failed  due to :
* * linker error : unresolved references to variables of the respective classes
* * possible troubleshooting and debugging very extensive due to the amount of classes
* * Solution : 
* * * check all dependencies once again and enlist help of gdb
* * * read up and learn about scopes of variables
* * * trial and error
* 
## Task 2 ##
* Linked-cell algorithm
* Cell wide in each direction is the cutoff value
* Structure as image below:
* We still provide the possibility to use simple sum implementation
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