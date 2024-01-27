//
// Created by maraconda on 12.11.23.
//

#pragma once

#include "../ParticleContainer.h"
#include "../LinkedCellContainer.h"

class ForceCalculator {
private:
    static double epsilon;
    static double sigma;
    static double gravity;
    static double cutoff;
public:
    static void GravityForceCalculation(ParticleContainer &container);
    static void LennardJonesForce(ParticleContainer &container, double eps, double sig);
    static void LennardJonesForceFaster(ParticleContainer &container, double eps, double sig, double Grav);
    static void LennardJonesForcePairwise(Particle *p1, Particle *p2);
    static void LennardJonesForceCell(LinkedCellContainer &cells, double Grav);
    static double smoothedLennardJonesPotential(Particle *p1, Particle *p2, double cutoff, double smoothedParameter);
    static void smoothedLennardJonesForcePairwise(Particle *p1, Particle *p2, double cutoff, double smoothedParameter);
};


