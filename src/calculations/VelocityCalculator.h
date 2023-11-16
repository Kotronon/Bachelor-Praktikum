//
// Created by maraconda on 12.11.23.
//

#pragma once


#include "../ParticleContainer.h"

class VelocityCalculator {
private:
public:
    static void VelocityStoermerVerlet(ParticleContainer &container, double delta_t);
    static void BrownianMotionInitialization(ParticleContainer &container, double avg_v, int dim);
};


#endif //PSEMOLDYN_GROUPH_VELOCITYCALCULATOR_H
