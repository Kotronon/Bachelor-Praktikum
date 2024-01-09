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

/**
 * Calculation of the new position of all molecules in the given LinkedCellContainer
 * @param cells
 * @param delta_t
 */
void PositionCalculator::PositionStoermerVerletCell(LinkedCellContainer &cells, double delta_t) {
    for (auto &x: cells) {
        for (auto &y: x) {
            for (auto &z: y) {
                for (auto &p: z) {
                    p.setOldX(p.getX());
                    std::array<double, 3> x_new = p.getX() + (delta_t * p.getV()) +
                                                  (((delta_t * delta_t) / (2 * p.getM())) *
                                                   p.getF());
                    p.setX(x_new);
                }
            }
        }
    }
    cells.moveToNeighbour();

}