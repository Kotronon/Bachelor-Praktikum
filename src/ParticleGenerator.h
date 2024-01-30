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
    createCuboid(std::array<double, 3> x, std::array<double, 3> v, std::array<int, 3> N, double h, double m, double sig, double eps, int type);

    static ParticleContainer createSphere(std::array<double, 3> x, std::array<double, 3> v, double m, int r, double h, double sig, double eps, int type);

    static ParticleContainer createDisk(std::array<double, 3> x, std::array<double, 3> v, double m, int r, double h, double sig, double eps, int type);

    static void createDiskInCells(std::array<double, 3> x, std::array<double, 3> v, double m, int r, double h,
                           LinkedCellContainer &cells, double sig, double eps, int type);

    static void createSphereInCells(std::array<double, 3> center, std::array<double, 3> v, double m, int r, double h,
                             LinkedCellContainer cells, double sig, double eps, int type);

    static ParticleContainer createMembrane(std::array<int, 3> n,std::array<double, 3> x,std::array<double,
            3> v,double m, double h,LinkedCellContainer &cells, double eps, double sig,double k, double r_0,int type);

    static ParticleContainer
    createDiskAlternative(std::array<double, 3> center, std::array<double, 3> v, double m, int r, double h, double sig, double eps, int type);

    static void createCuboidInCells(std::array<double, 3> x, std::array<double, 3> v, std::array<int, 3> N, double h, double m,
                                    LinkedCellContainer &cells, double sig, double eps, int type);
};

