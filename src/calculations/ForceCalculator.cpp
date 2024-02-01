//
// Created by maraconda on 12.11.23.
//

#include "ForceCalculator.h"
#include "../utils/ArrayUtils.h"
#include <cfloat>

double ForceCalculator::gravity = 0;
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
 * Calculates the Lennard-Jones force of all Particles in given ParticleContainer
 * @param container
 * @param eps epsilon
 * @param sig sigma
 */
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

/**
 * Faster calculation of the the Lennard Jones force of all Particles in given ParticleContainer
 * @param container ParticleContainer
 * @param grav gravitational force
 */
void ForceCalculator::LennardJonesForceFaster(ParticleContainer &container, double grav) {
    ForceCalculator::gravity = grav;
    std::array<double, 3> zero = {0,0,0};
    for (auto &p: container) {
        p.setOldF(p.getF());
        p.setF(zero);
    }
    container.applyForcePairwise(ForceCalculator::LennardJonesForcePairwise, gravity);
}

/**
 * calculation of the Lennard-Jones force for a pair of particles
 * @param p1 particle 1
 * @param p2 particle 2
 */
void ForceCalculator::LennardJonesForcePairwise(Particle *p1, Particle *p2) {
    double L2Norm_p1_p2 = ArrayUtils::L2Norm(p1->getX() - p2->getX());

    //make calculation if simple sum or distance between particles is smaller than the cutoff radius
    if(L2Norm_p1_p2 < cutoff) {
        double eps = sqrt(p1->getEps() * p2->getEps());
        double sig = (p1->getSig() + p2->getSig()) / 2;

        std::array<double, 3> force = ((-24 * eps / pow(L2Norm_p1_p2, 2)) * (pow(sig / L2Norm_p1_p2, 6)
                - (2 * pow(sig / L2Norm_p1_p2, 12))) * (p1->getX() - p2->getX()));

        #pragma omp critical
        {
            p1->setF(p1->getF() + force);
            p2->setF(p2->getF() - force);
        }
    }
}

/**
 * calculation of the Lennard-Jones force for all particles in a given LinkedCellContainer
 * @param cells LinkedCellContainer with particles
 * @param grav gravitational force
 */
void ForceCalculator::LennardJonesForceCell(LinkedCellContainer &cells, double grav){
    ForceCalculator::cutoff = cells.getCutoff();
    ForceCalculator::gravity = grav;
    cells.setZero();
    cells.applyForcePairwise(ForceCalculator::LennardJonesForcePairwise, ForceCalculator::smoothedLennardJonesForcePairwise, gravity);
}

/**
 * calculation of the smoothed Lennard-Jones force for a pair of particles
 * @param p1 particle 1
 * @param p2 particle 2
 * @param cutoff cutoff radius for force calculation
 * @param smoothedParameter smoothed parameter
 */
double ForceCalculator::smoothedLennardJonesPotential(Particle *p1, Particle *p2, double cutoff, double smoothedParameter) {
    double eps = sqrt(p1->getEps() * p2->getEps());
    double sig = (p1->getSig() + p2->getSig()) / 2;
    double L2Norm_p1_p2 = ArrayUtils::L2Norm(p1->getX() - p2->getX());

    double potential = 4 * eps * (pow((sig/L2Norm_p1_p2), 12) - pow((sig/L2Norm_p1_p2), 6));
    if (L2Norm_p1_p2 <= smoothedParameter) return potential;
    else if (L2Norm_p1_p2 >= cutoff) return 0;
    else {
        potential *= (1- (pow(L2Norm_p1_p2-smoothedParameter, 2)
                * (3*cutoff-smoothedParameter-2*L2Norm_p1_p2))/pow(cutoff-smoothedParameter, 3));
    }
    return potential;
}

/**
 * calculation of the smoothed Lennard-Jones force for a pair of particles
 * @param p1 particle 1
 * @param p2 particle 2
 * @param cutoff cutoff radius for force calculation
 * @param smoothedParameter smoothed parameter
 */
void ForceCalculator::smoothedLennardJonesForcePairwise(Particle *p1, Particle *p2, double cutoff,
                                                        double smoothedParameter) {
    std::array<double, 3> force = {0,0,0};
    double L2Norm_p1_p2 = ArrayUtils::L2Norm(p1->getX() - p2->getX());

    //make calculation if simple sum or distance between particles is smaller than the cutoff radius
    if(L2Norm_p1_p2 <= cutoff) {
        double eps = sqrt(p1->getEps() * p2->getEps());
        double sig = (p1->getSig() + p2->getSig()) / 2;

        if(L2Norm_p1_p2 <= smoothedParameter) force = ((-24 * eps / pow(L2Norm_p1_p2, 2)) * (pow(sig / L2Norm_p1_p2, 6) - (2 * pow(sig / L2Norm_p1_p2, 12))) *  (p1->getX() - p2->getX()));
        else {
            force  = force + ((-24*pow(sig, 6) * eps) / (pow(L2Norm_p1_p2, 14) * pow(cutoff-smoothedParameter, 3))) *
                    (cutoff - L2Norm_p1_p2) * (pow(cutoff, 2) * (2*pow(sig, 6) - pow(L2Norm_p1_p2, 6)) + cutoff * (3*smoothedParameter - L2Norm_p1_p2) *
                    (pow(L2Norm_p1_p2, 6) - 2*pow(sig, 6)) + L2Norm_p1_p2 * (5*smoothedParameter * pow(sig, 6) - 2*smoothedParameter * pow(L2Norm_p1_p2, 6)-
                    3*pow(sig, 6) * L2Norm_p1_p2 + pow(L2Norm_p1_p2, 7))) * (p2->getX() - p1->getX());
        }

        #pragma omp critical
        {
            p1->setF(p1->getF() + force);
            p2->setF(p2->getF() - force);
        }
    }
}
