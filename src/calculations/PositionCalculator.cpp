//
// Created by maraconda on 12.11.23.
//

#include "PositionCalculator.h"
#include "../ParticleContainer.h"
#include "../utils/ArrayUtils.h"


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
    for (int i = 0; i < grid.cell_numbers(); i++) {
        for(int j = 0; j < grid.Particles_in_cell(i); j++) {
            grid.cells[i][j].setX(grid.cells[i][j].getX() + (delta_t * grid.cells[i][j].getV()) +
                                  ((delta_t * delta_t) / (2 * grid.cells[i][j].getM())) * grid.cells[i][j].getF());
        }
    }
}