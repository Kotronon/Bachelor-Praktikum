//
// Created by Admin on 22.01.2024.
//

#ifndef PSEMOLDYN_GROUPH_MEMBRANE_H
#define PSEMOLDYN_GROUPH_MEMBRANE_H


#include <vector>
#include "Particle.h"

class Membrane {

private :
    const double k_ = 0;
    const double r_0_ = 0;
    const Particle grid[148][148][148];


    double euklid_norm
    (std::array<double, 3> x_i,
     std::array<double, 3> x_j);


public :

    std::array<double, 3> force_calculation(std::array<double, 3> x_i, std::array<double, 3> x_j);


    std::array<double, 3>
    diagonal_interaction(std::array<double, 3> x_i, std::array<double, 3> x_j);



};


#endif //PSEMOLDYN_GROUPH_MEMBRANE_H
