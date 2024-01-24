//
// Created by Admin on 22.01.2024.
//

#include "Membrane.h"
#include "Particle.h"
#include "utils/ArrayUtils.h"
#include <cmath>
#include <vector>
#include <array>

double k_ = 300;
double r_0_ = 3;
double sqrt_2 = std::sqrt(2);

double Membrane::euklid_norm(std::array<double, 3> x_i,
                             std::array<double, 3> x_j) {
    return ArrayUtils::L2Norm(x_i - x_j);
}

std::array<double, 3> Membrane::force_calculation(std::array<double, 3> x_i, std::array<double, 3> x_j) {
    std::array<double, 3> result = {0.0,0.0,0.0};
        double norm = euklid_norm(x_i,x_j);
    double teil1 = Membrane::k_ * (norm - r_0_);



    result[0] = teil1 * (x_j[0] - x_i[0] / norm);
    result[0] = teil1 * (x_j[1] - x_i[1] / norm);
    result[0] = teil1 * (x_j[2] - x_i[2] / norm);



    return result;
}

std::array<double, 3> Membrane::diagonal_interaction(std::array<double, 3> x_i, std::array<double, 3> x_j) {
    std::array<double, 3> result = {0.0,0.0,0.0};
    double norm = euklid_norm(x_i,x_j);
    double teil1 = Membrane::k_ * (norm - (sqrt_2*r_0_));

    result[0] = teil1 * (x_j[0] - x_i[0] / norm);
    result[0] = teil1 * (x_j[1] - x_i[1] / norm);
    result[0] = teil1 * (x_j[2] - x_i[2] / norm);

    return result;
}

double Membrane::harmonic_potential(std::array<double, 3> x_i, std::array<double, 3> x_j) {
    return k_/2 * (euklid_norm(x_i,x_j) - r_0_);
}

void Membrane::createGrid(int x, int y, int z) {
       Particle grid [x][y][z];

    for (int i = 0; i < x; ++i) {
        for (int j = 0; j < y; ++j) {
            for (int k = 0; k < z; ++k) {
                grid[i][j][k] = Particle();
            }
        }
    }

}

int main(int argc, char *argsv[]) {

}






