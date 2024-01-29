//
// Created by Admin on 22.01.2024.
//

#include "Membrane.h"
#include "Particle.h"
#include "utils/ArrayUtils.h"
#include <cmath>
#include <vector>
#include <array>


const double sqrt_2 = std::sqrt(2);

/**
 * Calculates the euclidean norm, eventhough it's already given*/
double Membrane::euklid_norm(std::array<double, 3> x_i,
                             std::array<double, 3> x_j) {
    return ArrayUtils::L2Norm(x_i - x_j);
}

/**force calculation, formula 1 on the worksheet*/
void Membrane::force_calculation(Particle *p1, Particle *p2) {
    std::array<double, 3> x_i, x_j;

    x_i = p1->getF();
    x_j = p2->getF();

    std::array<double, 3> result = {0.0,0.0,0.0};
    double norm = euklid_norm(x_i,x_j);
    double teil1 = Membrane::k_ * (norm - r_0_);



    result[0] = teil1 * (x_j[0] - x_i[0] / norm);
    result[1] = teil1 * (x_j[1] - x_i[1] / norm);
    result[2] = teil1 * (x_j[2] - x_i[2] / norm);

    p1->setOldF(p1->getF());
    p1->setF(result);

    p2->setOldF(p2->getF());
    p2->setF(result);


}


void Membrane::diagonal_interaction(Particle *p1, Particle *p2) {
    std::array<double, 3> x_i, x_j;
    std::array<double, 3> result = {0.0,0.0,0.0};
    double norm = euklid_norm(x_i,x_j);
    double teil1 = Membrane::k_ * (norm - (sqrt_2*r_0_));

    result[0] = teil1 * (x_j[0] - x_i[0] / norm);
    result[0] = teil1 * (x_j[1] - x_i[1] / norm);
    result[0] = teil1 * (x_j[2] - x_i[2] / norm);
    p1->setOldF(p1->getF());
    p1->setF(result);

    p2->setOldF(p2->getF());
    p2->setF(result);
}

double Membrane::harmonic_potential(std::array<double, 3> x_i, std::array<double, 3> x_j) {
    return k_/2 * (euklid_norm(x_i,x_j) - r_0_);
}

Membrane::Membrane(double k, double r_0) {
    k_ = k;
    r_0_ = r_0;
}


int main(int argc, char *argsv[]) {

}






