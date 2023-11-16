//
// Created by maraconda on 12.11.23.
//

#ifndef PSEMOLDYN_GROUPH_FORCECALCULATOR_H
#define PSEMOLDYN_GROUPH_FORCECALCULATOR_H


#include "../ParticleContainer.h"

class ForceCalculator {
private:
    static double epsilon;
    static double sigma;
public:
    static void SimpleForceCalculation(ParticleContainer &container);
    static void LennardJonesForce(ParticleContainer &container, double eps, double sig);
    static void LennardJonesForceFaster(ParticleContainer &container, double eps, double sig);
    static void LennardJonesForcePairwise(Particle *p1, Particle *p2);
};


#endif //PSEMOLDYN_GROUPH_FORCECALCULATOR_H
