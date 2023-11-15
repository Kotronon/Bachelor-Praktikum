//
// Created by maraconda on 12.11.23.
//

#ifndef PSEMOLDYN_GROUPH_FORCECALCULATOR_H
#define PSEMOLDYN_GROUPH_FORCECALCULATOR_H


#include "ParticleContainer.h"

class ForceCalculator {
private:
public:
    static void SimpleForceCalculation(ParticleContainer &container);
    static void LennardJonesForce(ParticleContainer &container, double eps, double sig);
};


#endif //PSEMOLDYN_GROUPH_FORCECALCULATOR_H
