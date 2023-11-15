//
// Created by maraconda on 12.11.23.
//

#ifndef PSEMOLDYN_GROUPH_POSITIONCALCULATOR_H
#define PSEMOLDYN_GROUPH_POSITIONCALCULATOR_H


#include "../ParticleContainer.h"

class PositionCalculator {
private:
public:
    static void PositionStoermerVerlet(ParticleContainer &container, double delta_t);
};


#endif //PSEMOLDYN_GROUPH_POSITIONCALCULATOR_H
