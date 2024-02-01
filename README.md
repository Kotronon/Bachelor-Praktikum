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
* Revision: a0ea2fd
* Compiler: gcc 11.4.0

# Run Code #
* packages to install:
  * sudo apt install libspdlog-dev
  * sudo apt-get install libgtest-dev
  * unfortunately only using fetchContent doesn't prevent from installing them
* to run the program:
  * ./MolSim
  * due to the xml input being incomplete it is only possible at the moment to change parameters in the main method manually as well as create cuboids/disks there
  * this should be fixed in the future, log level can not be chosen either currently as this was planned to be set up completely in the xml files
  * 
* to run the tests:
  * ctest
  

# Report #
## Task 1 ##
* Simulation of a Membrane
  * Worked on by Samhitha, report separate
 
## Task 2 ##
* Parallelization
  * Only one strategy mentioned, other one was worked on by Samhitha
* Performance analysis
  * After analyzing the performance of the Rayleigh-Taylor simulation in 3d with Perf it is very noticeable that the vast majority of time is spent in the force calculation
  * The LennardJonesForceCell function takes up nearly 94% of the calculation time
  * Other calculations with a simple iteration over all particles only take up < 1% of the time
  * -> Due to thread creation taking up time as well the first strategy only focuses on optimizing the force calculation
* Parallelization strategy
  * Parallelize first loop in LennardJonesForceCell, each thread gets assigned a group of cells
  * Each thread gets neighbours of this specific cell and calculates forces between particles in the same cell and neighbour cells
  * Due to Newtons third law being used during the calculations the particle updates have to be considered critical zones to avoid race conditions
  * Due to other calculations/functions being simple loops the thread creation usually takes up more time than it speeds the simulation up
  * -> Because of this only the force calculation is parallelized
* Possible improvements:
  * Depending on the kind of simulation it might be beneficial to add information on what scheduler to use
  * If number of particles is uneven across the cells threads might not be balanced, but schedulers to solve this also take more time
* Speedup of Parallelization:
  * Speedup was calculated as average over 1000 iterations of the 3D Rayleigh-Taylor simulation
  * More iterations were too slow to run on 1 thread in a reasonable time
  * Specs:
    * System:
      Kernel: 5.15.0-88-generic x86_64
      Desktop: Gnome 3.38.4
      Distro: Zorin OS 16.3
      base: Ubuntu 20.04 LTS Focal
    * CPU:
      Topology: 8-Core model: AMD Ryzen 7 5825U with Radeon Graphics
  * Conditions: measuring of time duration with functions in chrono library, disabled I/O operations
  * Results:
    * Result as SpeedupParallelization.png in submission
    * Speedup is steadily increasing with number of threads, reaching its peak on 8 threads with a speedup of about 3.64
    * After this speedup is decreasing again due to the pc this was run on only possessing 8 cores
    * Threads can not be run truly in parallel and the overhead of thread creation and management slows down the simulation more

## Task 3 ##
* Rayleigh-Taylor instability in 3d
  * Simulation could not be run completely due to simulation time not being reasonable (>70h)
  * This is most likely due to the optimization task from the last worksheet not being completed
  * Video provided shows the first 20 000 iterations instead of the full 200 000
  * Link to video: (TODO)

## Task 4 ##
* Nano-scale flow simulation (Option A - not chosen)
* We started with this task but decided to switch to task 5
* The walls are working and the new temperature calculations were done as well
* But due to some issues we did not finish all the tasks
* The uncompleted code can be found in branch Task4 (https://github.com/Kotronon/Bachelor-Praktikum/tree/Task4)

## Task 5 ##
* Crystallization of Argon (Option B - chosen)
* The new force calculations had to be added
* While computing the simulations we noticed multiple mistakes in the 3D computations
* To solve this we needed to add more ghost particles and make sure all functions in LinkedCellContainer work properly in a 3d space
* We improved our LinkedCell structure like you suggested by changing the cell sizes accordingly
* Equilibration: gets gaseous after 8 seconds in our video. That equals 30000 iterations
* Cooling: begins to get solid after 25 seconds in our video. That equals  125000 iterations
* Supercooling: gets solid after 6 seconds in our video. That equals 25000 iterations

## Misc ##
* We needed to comment out the thermostat tests because it suddenly would not work anymore, and we could not find the mistake in time
* All other test worked perfectly fine, so we were quite confused since the temperature in our simulations seemed to be calculated correctly