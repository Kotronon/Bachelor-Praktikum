//
// Created by maraconda on 12.11.23.
//

#pragma once


#include "../ParticleContainer.h"
#include "../LinkedCellContainer.h"

class PositionCalculator {
private:
public:
    static void PositionStoermerVerlet(ParticleContainer &container, double delta_t);
    static void PositionStoermerVerletCell(LinkedCellContainer &cells, double delta_t);
};


