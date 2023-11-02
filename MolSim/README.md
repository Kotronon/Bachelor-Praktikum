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
* Revision: insert last commit number
* Compiler: gcc 11.4.0

# Report #
## Task 1 ##
* Working with C++ works easily with Linux through VM or dual boot. 
* With WSL there were more problems with installing packages.

## Task 2 ##
* Calculation of force, x and v of each particle for planet and comet simulation. 
* Fi = SUM Fij with Fij = (MiMj * (xj-xi))/(||xi-xj||2^3)
* xi(tn+1) = xi(tn) + ∆t * vi(tn) + (∆t)^2 * Fi(tn) / 2mi
* vi (tn+1) = vi(tn) + ∆t * (Fi(tn) + Fi(tn+1)) / 2mi

## Task 3 ##
* Compile with all options enabled.
* run ./MolSim <path_to_file>/eingabe-sonne.txt 
* If you want to change delta_t please use flag -d <value>. 
* If you want to change end_time please use flag -e <value>.
* Visualization of sun, earth, jupiter and Halley's Comet movements.

## Task 4 ##
* Documentation: 
  * Doxygen option (BUILD_DOC) is on as default. But it can be set off in CMakeLists.txt
  * If you want to use doxygen please make sure you have it installed.
  * run make doc_doxygen 
* Refactoring:
  * Encapsulated Particle into ParticleContainer
  * Input immediately through command line from the beginning not while program runs


