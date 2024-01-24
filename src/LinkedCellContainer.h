//
// Created by kathi on 20.11.23.
//

#pragma once

#include "Particle.h"
#include "ParticleContainer.h"
#include <vector>
#include <functional>


class LinkedCellContainer {
public:



    LinkedCellContainer(std::array<double, 3> N, double cutoff,  std::array<std::string, 6> b, bool smoothed = false, double sLJparameter = -1);

    virtual ~LinkedCellContainer();


    [[nodiscard]] int cell_numbers() const;

    unsigned long Particles_in_cell(int x, int y, int z);


    void addParticle(int x, int y, int z, std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg, double sig, double eps);

    void addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg, double sig, double eps);

    void addParticle(int x, int y, int z, Particle &p);

    void addParticle(Particle &p);

    ParticleContainer toContainer();

    void addContainer(ParticleContainer &container);

    void moveToNeighbour();

    [[nodiscard]] std::vector<std::array<int, 3>> get_next_cells(int x, int y, int z) const;

    void setZero();

    [[nodiscard]] double getXCellSize() const;

    [[nodiscard]] double getYCellSize() const;

    [[nodiscard]] double getZCellSize() const;

    [[nodiscard]] double getCutoff() const;

    std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator begin();

    std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator end();

    void applyForcePairwise(const std::function<void(Particle *, Particle *)> &forceCalculation,
                            const std::function<void(Particle *, Particle *, double, double)> &smoothedforceCalculation, double Grav);

    bool applyMirrorBoundary(int particle, int x, int y, int z);

    void generateGhostCell(int index, int x, int y, int z);

    void deleteGhostCells();

    void moveIfPeriodic(double x_coordinate, double y_coordinate, double z_coordinate, Particle &p);

    double calculateDiffusion();

    std::vector<double> calculateRDF(int intervalBegin, int intervalEnd, double deltaR,  std::vector<int> x_axis_plot, std::ofstream  &filename);


private:
    int x_cells;
    int y_cells;
    int z_cells;
    double cutoff;
    double x_max;
    double y_max;
    double z_max;
    int it = 0;
    double smoothedRadius;
    bool smoothed;
    std::array<std::string, 6> boundary = {"o", "o", "o", "o", "o", "o"};
    std::vector<std::vector<std::vector<std::vector<Particle>>>> cells;
    double x_cell_size;
    double y_cell_size;
    double z_cell_size;

};

