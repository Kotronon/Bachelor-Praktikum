//
// Created by kathi on 01.11.23.
//

#ifndef PSEMOLDYN_PARTICLECONTAINER_H
#define PSEMOLDYN_PARTICLECONTAINER_H


#pragma once

#include <array>
#include <string>
#include <list>
#include "Particle.h"

class ParticleContainer: public Particle {
public:
    ParticleContainer(std::list<Particle> particles);
    ParticleContainer(std::pair<int, Particle> particles);

    void iterateThroughList(std::list<Particle> particles);
    std::list<Particle> iterateThroughPairs(std::pair<int, Particle> particles);
    };





#endif //PSEMOLDYN_PARTICLECONTAINER_H
