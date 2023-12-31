//
// Created by maraconda on 12.11.23.
//

#pragma once


#include "../ParticleContainer.h"
#include "../LinkedCellContainer.h"

class VelocityCalculator {
private:
public:
    static void VelocityStoermerVerlet(ParticleContainer &container, double delta_t);
    static void BrownianMotionInitialization(ParticleContainer &container, double avg_v, int dim);
    static void VelocityStoermerVerletCell(LinkedCellContainer &cells, double delta_t);
    static void BrownianMotionInitializationCell(LinkedCellContainer &cells, double avg_v, int dim);
};


