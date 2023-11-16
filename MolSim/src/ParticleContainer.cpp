//
// Created by maraconda on 02.11.23.
//

#include <iostream>
#include <vector>
#include <functional>
#include "spdlog/spdlog.h"
#include "ParticleContainer.h"


/**
 * creates a new empty ParticleContainer
 */
ParticleContainer::ParticleContainer() {
   spdlog::info("Empty Particle Container generated!");
}

/**
 * creates a new ParticleContainer as a copy from other ParticleContainer
 * @param other
 */
ParticleContainer::ParticleContainer(const ParticleContainer &other) {
    containedParticles = other.containedParticles;
    spdlog::info("Particle Container generated by copy!");
}

ParticleContainer::~ParticleContainer() { spdlog::info("Particle Container destructed!"); }

/**
 * adds an existing Particle
 * @param particle
 */
void ParticleContainer::addParticle(const Particle &particle) {
    const Particle& new_particle(particle);
    containedParticles.emplace_back(new_particle);
    spdlog::info("Added copy of particle to container!");
}

void ParticleContainer::addParticle(int type_arg) {
    Particle new_particle(type_arg);
    containedParticles.emplace_back(new_particle);
    spdlog::info("Added newly created particle to container!");
}


/**
 * adds new Particle
 * @param x_arg
 * @param v_arg
 * @param m_arg
 * @param type_arg
 */
void ParticleContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg) {
    Particle new_particle(x_arg, v_arg, m_arg, type_arg);
    containedParticles.emplace_back(new_particle);
    spdlog::info("Added newly created particle to container!");
    containedParticles.begin();
}

void ParticleContainer::addParticleContainer(ParticleContainer &container) {
    for (auto &p: container) {
        containedParticles.emplace_back(p);
    }
    spdlog::info("Added particles contained in other container to this container!");
}

void ParticleContainer::applyForcePairwise(const std::function<void(Particle*, Particle*)>& forceCalculation){
    auto first = begin();
    auto last = end();
    for (; first != last; ++first) {
        for(auto next = std::next(first); next != last; ++next)
            forceCalculation(&(*first), &(*next));
    }
}

/**
 * returns the vector of all PArticles in the ParticleContainer with Pointer at first element
 * @return
 */
std::vector<Particle>::iterator ParticleContainer::begin() {
    return containedParticles.begin();
}

/**
 * returns the vector of all PArticles in the ParticleContainer with Pointer at last element
 * @return
 */
std::vector<Particle>::iterator ParticleContainer::end() {
    return containedParticles.end();
}

/**
 * returns number of Particles in ParticleContainer
 * @return
 */
int ParticleContainer::size() {
    return containedParticles.size();
}
