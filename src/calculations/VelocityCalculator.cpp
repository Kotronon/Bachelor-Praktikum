//
// Created by maraconda on 12.11.23.
//

#include "VelocityCalculator.h"
#include "../utils/ArrayUtils.h"
#include "../utils/MaxwellBoltzmannDistribution.h"
#include <spdlog/spdlog.h>

/**
 * Calculation of the new velocity of all particles in the given ParticleContainer according to Brownian Motion initialization
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
 * Calculation of the new velocity of all particles in the given ParticleContainer according to Störmer Verlet
 * @param container LinkedCellContainer
 * @param delta_t time difference
 */
void VelocityCalculator::VelocityStoermerVerlet(ParticleContainer &container, double delta_t) {
    for (auto &p: container) {
        //vi (tn+1) = vi(tn) + ∆t * Fi(tn) + Fi(tn+1) / 2mi
        p.setV(p.getV() + ((delta_t / (2 * p.getM())) * (p.getOldF() + p.getF())));
    }
}

/**
 * Calculation of the new velocity of all particles in the given LinkedCellContainer according to Brownian Motion initialization
 * @param cells LinkedCellContainer
 * @param avg_v average velocity
 * @param dim dimension (2 or 3)
 */
void VelocityCalculator::BrownianMotionInitializationCell(LinkedCellContainer &cells, double avg_v, int dim) {
    std::array<double, 3> brownian_motion{};
    //Calculate Brownian motion and add it to velocities of all particles (without ghost particles)
    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); p++) {
                    brownian_motion = maxwellBoltzmannDistributedVelocity(avg_v, dim);
                    p->setV(p->getV() + brownian_motion);
                }
            }
        }

    }
}

/**
 * Calculation of the new velocity of all particles in the given LinkedCellContainer according to Störmer Verlet
 * @param cells LinkedCellContainer
 * @param delta_t time difference
 */
void VelocityCalculator::VelocityStoermerVerletCell(LinkedCellContainer &cells, double delta_t) {
    //Calculate new velocities of all particles (without ghost particles)
    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); p++) {
                    //vi (tn+1) = vi(tn) + ∆t * Fi(tn) + Fi(tn+1) / 2mi
                    p->setV(p->getV() + ((delta_t / (2 * p->getM())) * (p->getOldF() + p->getF())));
                }
            }
        }
    }
}