//
// Created by maraconda on 12.11.23.
//

#include "ForceCalculator.h"
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"

void ForceCalculator::SimpleForceCalculation(ParticleContainer &container) {
    std::array<double, 3> force{};
    for (auto &p1: container) {
        force = {0., 0., 0.};
        for (auto &p2: container) {
            if (!(p1 == p2)) {
                force = force + (((p1.getM() * p2.getM()) / pow((ArrayUtils::L2Norm(p1.getX() - p2.getX())), 3)) * (p2.getX() - p1.getX()));
            }
        }
        p1.setOldF(p1.getF());
        p1.setF(force);
    }
}

void ForceCalculator::LennardJonesForce(ParticleContainer &container, double eps, double sig) {
    std::array<double, 3> force{};
    double L2Norm_p1_p2;
    for (auto &p1: container) {
        force = {0., 0., 0.};
        for (auto &p2: container) {
            if (!(p1 == p2)) {
                L2Norm_p1_p2 = ArrayUtils::L2Norm(p1.getX() - p2.getX());
                force = force + ((-24*eps / pow(L2Norm_p1_p2,2)) * (pow(sig/L2Norm_p1_p2,6) - (2 * pow(sig/L2Norm_p1_p2,12))) * (p1.getX() - p2.getX()));
            }
        }
        p1.setOldF(p1.getF());
        p1.setF(force);
    }
}