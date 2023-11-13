//
// Created by maraconda on 12.11.23.
//

#include "VelocityCalculator.h"
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "utils/MaxwellBoltzmannDistribution.h"

void VelocityCalculator::BrownianMotionInitialization(ParticleContainer &container, double avg_v) {
    for (auto &p: container) {
        p.setV(maxwellBoltzmannDistributedVelocity(avg_v, 3));
    }
}

void VelocityCalculator::VelocityStoermerVerlet(ParticleContainer &container, double delta_t) {
    for (auto &p: container) {
        //vi (tn+1) = vi(tn) + ∆t * Fi(tn) + Fi(tn+1) / 2mi
        p.setV(p.getV() + ((delta_t / (2 * p.getM())) * (p.getOldF() + p.getF())));
    }
}