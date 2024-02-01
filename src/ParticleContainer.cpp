//
// Created by maraconda on 02.11.23.
//

#include <iostream>
#include <vector>
#include <functional>
#include "spdlog/spdlog.h"
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"

/**
 * creates a new empty ParticleContainer
 */
ParticleContainer::ParticleContainer() = default;

/**
 * adds all particles from other ParticleContainer
 * @param other ParticleContainer to copy
 */
ParticleContainer::ParticleContainer(const ParticleContainer &other) {
    containedParticles = other.containedParticles;
}

/**
 * destructs the ParticleContainer
 */
ParticleContainer::~ParticleContainer() = default;

/**
 * adds an existing Particle
 * @param particle particle to add
 */
void ParticleContainer::addParticle(const Particle &particle) {
    const Particle& new_particle(particle);
    containedParticles.emplace_back(new_particle);
}

/**
 * add a new empty particle to the container
 * @param type_arg particle type
 */
void ParticleContainer::addParticle(int type_arg) {
    Particle new_particle(type_arg);
    containedParticles.emplace_back(new_particle);
}


/**
 * adds new Particle
 * @param coordinate position coordinates
 * @param v_arg velocity
 * @param m_arg mass
 * @param type_arg particle type
 */
void ParticleContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg, double sig, double eps) {
    Particle new_particle(x_arg, v_arg, m_arg, sig, eps, type_arg);
    containedParticles.emplace_back(new_particle);
}

/**
 * add another ParticleContainer to the current container
 * @param container
 */
void ParticleContainer::addParticleContainer(ParticleContainer &container) {
    for (auto &p: container) {
        containedParticles.emplace_back(p);
    }
    spdlog::info("Added particles contained in other container to this container!");
}

/**
 * iterate through the Particles pairwise to calculate the force of each particle
 * @param forceCalculation function to calculate force between two particles
 * @param Ggrav gravity force to add
 */
void ParticleContainer::applyForcePairwise(const std::function<void(Particle*, Particle*)>& forceCalculation, double Ggrav){
    auto first = begin();
    auto last = end();
    for (; first != last; ++first) {
        for(auto next = std::next(first); next != last; ++next)
            forceCalculation(&(*first), &(*next));
        std::array<double, 3> gravity = {0, first->getM(), 0};
        first->setF(first->getF() + gravity);
    }
}

/**
 * returns the vector of all particles in the ParticleContainer with Pointer at first element
 * @return
 */
std::vector<Particle>::iterator ParticleContainer::begin() {
    return containedParticles.begin();
}

/**
 * returns the vector of all particles in the ParticleContainer with Pointer at last element
 * @return
 */
std::vector<Particle>::iterator ParticleContainer::end() {
    return containedParticles.end();
}

/**
 * returns number of Particles in ParticleContainer
 * @return
 */
unsigned long ParticleContainer::size() {
    return containedParticles.size();
}

