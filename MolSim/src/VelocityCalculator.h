//
// Created by maraconda on 12.11.23.
//

#ifndef PSEMOLDYN_GROUPH_VELOCITYCALCULATOR_H
#define PSEMOLDYN_GROUPH_VELOCITYCALCULATOR_H


#include "ParticleContainer.h"

class VelocityCalculator {
private:
public:
    void BrownianMotionInitialization(ParticleContainer &container);
    static void VelocityStoermerVerlet(ParticleContainer &container, double delta_t);
};


#endif //PSEMOLDYN_GROUPH_VELOCITYCALCULATOR_H
