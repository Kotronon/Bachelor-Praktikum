//
// Created by kathi on 20.11.23.
//

#pragma once

#include "Particle.h"
#include <vector>


class LinkedCellContainer {
public:
    std::vector<std::vector<Particle>> cells{};

    LinkedCellContainer(std::array<int, 3> N, double cutoff);

    virtual ~LinkedCellContainer();


    int cell_numbers();

    int Particles_in_cell(int cell);


    void addParticle(int cell, std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg);

    void addParticle(int cell, Particle &p);

    void deleteParticle(int cell, Particle &p);

    void moveToNeighbour(int cell_current, int cell_new, Particle &p);

    std::vector<int> get_Particles_from_next_cells(int cell);

    void setZero();


private:
    int x_cells;
    int y_cells;
    int z_cells;
};

