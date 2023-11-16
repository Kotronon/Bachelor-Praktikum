//
// Created by maraconda on 12.11.23.
//

#pragma once


#include <array>
#include "ParticleContainer.h"

class ParticleGenerator {
private:
public:
    static ParticleContainer
    createCuboid(std::array<double, 3> x, std::array<double, 3> v, std::array<int, 3> N, double h, double m);
};


#endif //PSEMOLDYN_GROUPH_PARTICLEGENERATOR_H
