//
// Created by maraconda on 12.11.23.
//

#include "ForceCalculator.h"
#include "../ParticleContainer.h"
#include "../utils/ArrayUtils.h"
#include <spdlog/spdlog.h>
#include <cfloat>

double ForceCalculator::epsilon = 5;
double ForceCalculator::sigma = 1;
double ForceCalculator::Ggrav = 0;
double ForceCalculator::cutoff = DBL_MAX;


/**
 * Calculates the gravity force of all Particles in given ParticleContainer
 * @param container
 */
void ForceCalculator::GravityForceCalculation(ParticleContainer &container) {
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

/**
 * Calculates the Lennord Jones force of all Particles in given ParticleContainer
 * @param container
 * @param eps
 * @param sig
 */
void ForceCalculator::LennardJonesForce(ParticleContainer &container, double eps, double sig) {
    ForceCalculator::epsilon = eps;
    ForceCalculator::sigma = sig;

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

/**
 * Faster calculation of the the Lennard Jones force of all Particles in given ParticleContainer
 * @param container
 * @param eps
 * @param sig
 */
void ForceCalculator::LennardJonesForceFaster(ParticleContainer &container, double eps, double sig, double grav) {
    ForceCalculator::epsilon = eps;
    ForceCalculator::sigma = sig;
    ForceCalculator::Ggrav = grav;
    std::array<double, 3> zero = {0,0,0};
    for (auto &p: container) {
        p.setOldF(p.getF());
        p.setF(zero);
    }
    container.applyForcePairwise(ForceCalculator::LennardJonesForcePairwise, Ggrav);
}

/**
 * pairwise calculation of the Lennord Jones force of all Particles in given ParticleContainer
 * @param p1
 * @param p2
 */
void ForceCalculator::LennardJonesForcePairwise(Particle *p1, Particle *p2) {
    std::array<double, 3> force = {0,0,0};
    double L2Norm_p1_p2 = ArrayUtils::L2Norm(p1->getX() - p2->getX());
    //make calculation if simple sum or distance between particles is smaller than the cutoff radius
    if(L2Norm_p1_p2 < cutoff) {
        double eps = sqrt(p1->getEps() * p2->getEps());
        double sig = (p1->getSig() + p2->getSig()) / 2;
        force = force +
                ((-24 * eps / pow(L2Norm_p1_p2, 2)) * (pow(sig / L2Norm_p1_p2, 6) - (2 * pow(sig / L2Norm_p1_p2, 12))) *
                 (p1->getX() - p2->getX()));
        if(!p1->getFixed()) p1->setF(p1->getF() + force);
        if(!p2->getFixed()) p2->setF(p2->getF() - force);
    }
}

void ForceCalculator::LennardJonesForceCell(LinkedCellContainer &cells, double grav){
    ForceCalculator::cutoff = cells.getCutoff();
    ForceCalculator::Ggrav = grav;
    cells.setZero();
    cells.applyForcePairwise(ForceCalculator::LennardJonesForcePairwise, Ggrav);
}
