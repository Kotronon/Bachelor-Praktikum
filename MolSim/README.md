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

If you do not specify compiler and relevant flags, and we cannot simply compile with whatever compiler we have immediately on hand, then you will lose marks for uncompilable code!

# Report #
## Task 1 ##
* Working with C++ works easily with Linux through VM. With WSL there were more problems.

## Task 2 ##
* Calculation of force, x and v of each particle for planet and comet simulation. 
* Fi = SUM Fij with Fij = (MiMj * (xj-xi))/(||xi-xj||2^3)
* xi(tn+1) = xi(tn) + ∆t * vi(tn) + (∆t)^2 * Fi(tn) / 2mi
* vi (tn+1) = vi(tn) + ∆t * (Fi(tn) + Fi(tn+1)) / 2mi

## Task 3 ##
* Visualization of sun, earth, jupiter and Halley's Comet movements.
* Compile with all options enabled
* run ./MolSim <path_to_file>/eingabe-sonne.txt
* You will get asked for a delta_t and a t_end. You need to write something to the command line, otherwise the program won't go on.

## Task 4 ##
* Refactoring and documentation
* Doxygen option is on as default. But it can be set off in CMakeLists.txt
* If you want to use doxygen please make sure you have it installed.

