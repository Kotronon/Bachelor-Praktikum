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


    int cell_numbers() const;

    int Particles_in_cell(int x, int y, int z);


    void addParticle(int x, int y, int z, std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg);

    void addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg);

    void addParticle(int x, int y, int z, Particle &p);

    void deleteParticle(int x, int y, int z, Particle &p);

    void moveToNeighbour();

    [[nodiscard]] std::vector<std::array<int, 3>> get_next_cells(int x, int y, int z) const;

    void setZero();

    int getXMax() const;

    int getYMax() const;

    int getZMax() const;

    std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator begin();

    std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator end();

    void applyForcePairwise(const std::function<void(Particle *, Particle *)> &forceCalculation);

    bool applyMirrorBoundary(int particle, int x, int y, int z);

    void generateGhostCell(int index, int x, int y, int z);


private:
    int x_cells;
    int y_cells;
    int z_cells;
    double c;
    std::vector<std::string> boundary;
    std::vector<std::vector<std::vector<std::vector<Particle>>>> cells;
};

