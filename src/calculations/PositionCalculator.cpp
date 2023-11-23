//
// Created by maraconda on 12.11.23.
//

#include "PositionCalculator.h"
#include "../ParticleContainer.h"
#include "../utils/ArrayUtils.h"
#include <spdlog/spdlog.h>
#include <cmath>

/**
 * Calculation of the new position of all molecules in the given ParticleContainer
 * @param container
 * @param delta_t
 */
void PositionCalculator::PositionStoermerVerlet(ParticleContainer &container, double delta_t) {
    for (auto &p: container) {
        p.setX(p.getX() + (delta_t * p.getV()) + ((delta_t*delta_t)/(2*p.getM())) * p.getF());
    }
}

void PositionCalculator::PositionStoermerVerletCell(LinkedCellContainer &grid, double delta_t) {
    for (auto &c: grid) {
        for (auto &p: c) {
            if (p.getType() != 0) {
                p.setX(p.getX() + (delta_t * p.getV()) +
                                              (((delta_t * delta_t) / (2 * p.getM())) *
                                               p.getF()));
            }
        }
    }
    grid.moveToNeighbour();

}