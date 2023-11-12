//
// Created by maraconda on 12.11.23.
//

#include "VelocityCalculator.h"
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"

void VelocityCalculator::BrownianMotionInitialization(ParticleContainer &container) {

}

void VelocityCalculator::VelocityStoermerVerlet(ParticleContainer &container, double delta_t) {
    for (auto &p: container) {
        //vi (tn+1) = vi(tn) + ∆t * Fi(tn) + Fi(tn+1) / 2mi
        p.setV(p.getV() + ((delta_t / (2 * p.getM())) * (p.getOldF() + p.getF())));
    }
}