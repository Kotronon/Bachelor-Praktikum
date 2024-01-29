//
// Created by Admin on 22.01.2024.
//

#ifndef PSEMOLDYN_GROUPH_MEMBRANE_H
#define PSEMOLDYN_GROUPH_MEMBRANE_H


#include <vector>
#include "Particle.h"

class Membrane {

private :
      static double k_;
   static double r_0_;





public :

    Membrane(double k, double r_0);

    static double euklid_norm
            (std::array<double, 3> x_i,
             std::array<double, 3> x_j);

    double harmonic_potential(std::array<double, 3> x_i, std::array<double, 3> x_j);

    static void force_calculation(Particle *p1, Particle *p2);


    static void
    diagonal_interaction(Particle *p1, Particle *p2);


};


#endif //PSEMOLDYN_GROUPH_MEMBRANE_H
