//
// Created by maraconda on 02.11.23.
//

#pragma once

#include <vector>
#include <functional>
#include "Particle.h"

class ParticleContainer {
private:
    /**
   * List of all particles inside this particle container
   */
    std::vector<Particle> containedParticles{};

public:
    ParticleContainer();

    ParticleContainer(const ParticleContainer &other);

    virtual ~ParticleContainer();

    int size();

    void addParticle(const Particle &particle);

    void addParticle(int type_arg);

    void addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg, double sig, double eps);

    std::vector<Particle>::iterator begin();

    std::vector<Particle>::iterator end();

    void addParticleContainer(ParticleContainer &container);

    void applyForcePairwise(const std::function<void(Particle *, Particle *)> &forceCalculation, double Ggrav);

    void SetAllNeighbours(Particle &p);

    void removeDuplicates();

};


