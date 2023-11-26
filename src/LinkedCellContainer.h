//
// Created by kathi on 20.11.23.
//

#pragma once

#include "Particle.h"
#include <vector>
#include <functional>


class LinkedCellContainer {
public:



    LinkedCellContainer(std::array<int, 3> N, double cutoff,  std::vector<std::string> b);

    virtual ~LinkedCellContainer();


    [[nodiscard]] int cell_numbers() const;

    int Particles_in_cell(int x, int y, int z);


    void addParticle(int x, int y, int z, std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg);

    void addParticle(int x, int y, int z, Particle &p);

    void addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
                                      int type_arg);

    void addParticle(Particle &p);

    void deleteParticle(int x, int y, int z, Particle &p);

    void moveToNeighbour();

    [[nodiscard]] std::vector<std::array<int, 3>> get_next_cells(int x, int y, int z) const;

    void setZero();

    [[nodiscard]] int getXMax() const;

    [[nodiscard]] int getYMax() const;

    [[nodiscard]] int getZMax() const;

    std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator begin();

    std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator end();

    void applyForcePairwise(const std::function<void(Particle *, Particle *)> &forceCalculation);


private:
    int x_cells;
    int y_cells;
    int z_cells;
    double c;
    std::vector<std::string> boundary;
    std::vector<std::vector<std::vector<std::vector<Particle>>>> cells;

};

