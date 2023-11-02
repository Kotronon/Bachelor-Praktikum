//
// Created by maraconda on 02.11.23.
//

#ifndef PSEMOLDYN_GROUPH_PARTICLECONTAINER_H
#define PSEMOLDYN_GROUPH_PARTICLECONTAINER_H


#include <forward_list>
#include "Particle.h"

class ParticleContainer {
private:
    /**
   * List of all particles inside this particle container
   */
    std::forward_list<Particle> containedParticles{};

    /**
   * Number of contained particles
   */
   int numberOfParticles{};

public:
    ParticleContainer();

    ParticleContainer(const ParticleContainer &other);

    virtual ~ParticleContainer();

    void addParticle(const Particle &particle);
};


#endif //PSEMOLDYN_GROUPH_PARTICLECONTAINER_H
