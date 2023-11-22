//
// Created by kathi on 20.11.23.
//

#pragma once

#include "Particle.h"
#include <vector>
#include <functional>


class LinkedCellContainer {
public:

    std::vector<std::vector<Particle>> cells;

    LinkedCellContainer(std::array<int, 3> N, double cutoff,  std::vector<std::string> b);

    virtual ~LinkedCellContainer();


    int cell_numbers() const;

    int Particles_in_cell(int cell);


    void addParticle(int cell, std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg);

    void addParticle(int cell, Particle &p);

    void deleteParticle(int cell, Particle &p);

    void moveToNeighbour();

    [[nodiscard]] std::vector<int> get_Particles_from_next_cells(int cell) const;

    void setZero();

    int getXMax() const;

    int getYMax() const;

    int getZMax() const;

private:
    int x_cells;
    int y_cells;
    int z_cells;
    double c;
    std::vector<std::string> boundary;
};

