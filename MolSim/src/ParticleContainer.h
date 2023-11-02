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

    int size();

    void addParticle(const Particle &particle);

    void addParticle(int type_arg);

    void addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg);

    std::_Fwd_list_iterator<Particle> begin();

    std::_Fwd_list_iterator<Particle> end();
};


#endif //PSEMOLDYN_GROUPH_PARTICLECONTAINER_H
