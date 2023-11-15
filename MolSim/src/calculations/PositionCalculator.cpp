//
// Created by maraconda on 12.11.23.
//

#include "PositionCalculator.h"
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"

void PositionCalculator::PositionStoermerVerlet(ParticleContainer &container, double delta_t) {
    for (auto &p: container) {
        p.setX(p.getX() + (delta_t * p.getV()) + ((delta_t*delta_t)/(2*p.getM())) * p.getF());
    }
}
