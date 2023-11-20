//
// Created by maraconda on 12.11.23.
//

#pragma once


#include <array>
#include "ParticleContainer.h"
#include "LinkedCellContainer.h"

class ParticleGenerator {
private:
public:
    static ParticleContainer
    createCuboid(std::array<double, 3> x, std::array<double, 3> v, std::array<int, 3> N, double h, double m);
    void createCuboidInCells(std::array<double, 3> x, std::array<double, 3> v,
                                                std::array<int, 3> N, double h, double m,
                                                LinkedCellContainer &cells, double  cutoff);
};

