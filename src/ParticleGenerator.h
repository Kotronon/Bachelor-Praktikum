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
    createCuboid(std::array<double, 3> x, std::array<double, 3> v, std::array<int, 3> N, double h, double m, double sig, double eps, int type, bool fixed = false);

    static ParticleContainer createSphere(std::array<double, 3> x, std::array<double, 3> v, double m, int r, double h, double sig, double eps, int type, bool fixed = false);

    static ParticleContainer createDisk(std::array<double, 3> x, std::array<double, 3> v, double m, int r, double h, double sig, double eps, int type, bool fixed = false);

    static void createDiskInCells(std::array<double, 3> x, std::array<double, 3> v, double m, int r, double h,
                           LinkedCellContainer &cells, double sig, double eps, int type, bool fixed = false);

    static void createSphereInCells(std::array<double, 3> center, std::array<double, 3> v, double m, int r, double h,
                             LinkedCellContainer cells, double sig, double eps, int type, bool fixed = false);

    static ParticleContainer
    createDiskAlternative(std::array<double, 3> center, std::array<double, 3> v, double m, int r, double h, double sig, double eps, int type, bool fixed = false);

    static void createCuboidInCells(std::array<double, 3> x, std::array<double, 3> v, std::array<int, 3> N, double h, double m,
                                    LinkedCellContainer &cells, double sig, double eps, int type, bool fixed = false);
};

