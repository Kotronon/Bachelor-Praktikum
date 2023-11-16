//
// Created by maraconda on 12.11.23.
//

#pragma once


#include "../ParticleContainer.h"

class PositionCalculator {
private:
public:
    static void PositionStoermerVerlet(ParticleContainer &container, double delta_t);
};


