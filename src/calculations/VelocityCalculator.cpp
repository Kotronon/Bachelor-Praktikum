//
// Created by maraconda on 12.11.23.
//

#include "VelocityCalculator.h"
#include "../ParticleContainer.h"
#include "../utils/ArrayUtils.h"
#include "../utils/MaxwellBoltzmannDistribution.h"

/**
 * Calculation of the new velocity of all molecule in the given ParticleContainer according to Brownian Motion Initialization
 * @param container
 * @param avg_v
 */
void VelocityCalculator::BrownianMotionInitialization(ParticleContainer &container, double avg_v, int dim) {
    std::array<double, 3> brownian_motion{};
    for (auto &p: container) {
        brownian_motion = maxwellBoltzmannDistributedVelocity(avg_v, dim);
        p.setV(p.getV() + brownian_motion);
    }
}

/**
 * Calculation of the new velocity of all molecule in the given ParticleContainer according to Strömer Verlet
 * @param container
 * @param delta_t
 */
void VelocityCalculator::VelocityStoermerVerlet(ParticleContainer &container, double delta_t) {
    for (auto &p: container) {
        //vi (tn+1) = vi(tn) + ∆t * Fi(tn) + Fi(tn+1) / 2mi
        p.setV(p.getV() + ((delta_t / (2 * p.getM())) * (p.getOldF() + p.getF())));
    }
}

/**
 * Calculation of the new velocity of all molecule in the given ParticleContainer according to Brownian Motion Initialization
 * @param container
 * @param avg_v
 */
void VelocityCalculator::BrownianMotionInitializationCell(LinkedCellContainer &grid, double avg_v, int dim) {
    std::array<double, 3> brownian_motion{};
    for (int i = 0; i < grid.cell_numbers(); i++) {
        for (int j = 0; j < grid.Particles_in_cell(i); j++) {
            brownian_motion = maxwellBoltzmannDistributedVelocity(avg_v, dim);
            grid.cells[i][j].setV(grid.cells[i][j].getV() + brownian_motion);
        }

    }
}

/**
 * Calculation of the new velocity of all molecule in the given ParticleContainer according to Strömer Verlet
 * @param container
 * @param delta_t
 */
void VelocityCalculator::VelocityStoermerVerletCell(LinkedCellContainer &grid, double delta_t) {
    for (int i = 0; i < grid.cell_numbers(); i++) {
        for (int j = 0; j < grid.Particles_in_cell(i); j++) {
            //vi (tn+1) = vi(tn) + ∆t * Fi(tn) + Fi(tn+1) / 2mi
            grid.cells[i][j].setV(grid.cells[i][j].getV() + ((delta_t / (2 * grid.cells[i][j].getM())) * (grid.cells[i][j].getOldF() + grid.cells[i][j].getF())));
        }
    }
}