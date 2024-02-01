//
// Created by kathi on 20.11.23.
//


#include "LinkedCellContainer.h"
#include "utils/ArrayUtils.h"
#include <spdlog/spdlog.h>
#include <cmath>
#include <utility>
#include <matplot/matplot.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * create a new linked cell container based on the dimensions and boundary
 * @param N dimensions of the container
 * @param cutoffRadius radius of cutoff
 * @param boundaryConditions boundary types of each side of the container
 * @param smoothed option to enable smoothed Lennard-Jones potential in force calculation
 * @param sLJparameter parameter for smoothed Lennard-Jones potential
 */
LinkedCellContainer::LinkedCellContainer(std::array<double, 3> N, double cutoffRadius,
                                         std::array<std::string, 6> boundaryConditions, bool smoothed,
                                         double sLJparameter) {
    cutoff = cutoffRadius;
    this->smoothed = smoothed;
    smoothedRadius = sLJparameter;

    //calculate number of cells on each axis
    x_cells = (int) round(N[0] / cutoff);
    y_cells = (int) round(N[1] / cutoff);
    z_cells = (int) round(N[2] / cutoff);

    x_max = N[0];
    y_max = N[1];
    z_max = N[2];

    //calculate cell size in each direction
    x_cell_size = x_max / x_cells;
    y_cell_size = y_max / y_cells;
    z_cell_size = z_max / z_cells;

    //create cell structure according to calculated parameters
    std::vector<std::vector<std::vector<std::vector<Particle>>>> x;
    for (int i = 0; i < x_cells + 2; i++) {
        std::vector<std::vector<std::vector<Particle>>> y;
        for (int j = 0; j < y_cells + 2; j++) {
            std::vector<std::vector<Particle>> z(z_cells + 2);
            y.push_back(z);
        }
        x.push_back(y);
    }

    cells = x;
    boundary = std::move(boundaryConditions);
}

/**
 * destructs LinkedCellContainer
 */
LinkedCellContainer::~LinkedCellContainer() = default;

/**
 * returns the number of cells in the grid
 * @return number of cells in entire grid
 */
int LinkedCellContainer::cell_numbers() const {
    return x_cells * y_cells * z_cells;
}

/**
 * returns the number of particles in the given cell
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 * @return number of molecules in cell
 */
unsigned long LinkedCellContainer::Particles_in_cell(int x, int y, int z) {
    return cells[x][y][z].size();
}

/**
 * adds new particle to specific cell
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 * @param x_arg coordinates of new particle
 * @param v_arg velocity of new particle
 * @param m_arg mass of new particle
 * @param type_arg type of new particle
 * @param sig sigma of new particle
 * @param eps epsilon of new particle
 */
void LinkedCellContainer::addParticle(int x, int y, int z, std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                                      double m_arg, int type_arg, double sig, double eps) {
    Particle new_particle = Particle(x_arg, v_arg, m_arg, sig, eps, type_arg);
    cells[x][y][z].emplace_back(new_particle);
}

/**
 * adds new Particle to a cell based on current position
 * @param x_arg coordinates of new particle
 * @param v_arg velocity of new particle
 * @param m_arg mass of new particle
 * @param type_arg type of new particle
 * @param sig sigma of new particle
 * @param eps epsilon of new particle
 */
void
LinkedCellContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg,
                                 double sig, double eps) {
    Particle new_particle = Particle(x_arg, v_arg, m_arg, sig, eps, type_arg);
    int x = (int) floor(x_arg[0] / x_cell_size) + 1;
    if(x_arg[0] == x_max) x = x_cells;
    if(x_arg[0] == 0) x = 1;
    int y = (int) floor(x_arg[1] / y_cell_size) + 1;
    if(x_arg[1] == y_max) y = y_cells;
    if(x_arg[1] == 0) y = 1;
    int z = (int) floor(x_arg[2] / z_cell_size) + 1;
    if(x_arg[2] == z_max) z = z_cells;
    if(x_arg[2] == 0) z = 1;
    cells[x][y][z].emplace_back(new_particle);
}

/**
 * adds existing particle to specific cell
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 * @param p existing particle to add
 */
void LinkedCellContainer::addParticle(int x, int y, int z, Particle &p) {
    cells[x][y][z].emplace_back(p);
}


/**
 * adds existing particle to a cell based on current position
 * @param p existing particle to add
 */
void LinkedCellContainer::addParticle(Particle &p) {
    int x = (int) floor(p.getX()[0] / x_cell_size) + 1;
    if(p.getX()[0] == x_max) x = x_cells;
    if(p.getX()[0] == 0) x = 1;
    int y = (int) floor(p.getX()[1] / y_cell_size) + 1;
    if(p.getX()[1] == y_max) y = y_cells;
    if(p.getX()[1] == 0) y = 1;
    int z = (int) floor(p.getX()[2] / z_cell_size) + 1;
    if(p.getX()[2] == z_max) z = z_cells;
    if(p.getX()[2] == 0) z = 1;
    cells[x][y][z].emplace_back(p);
}

/**
 * create a new ParticleContainer containing the same particles as this LinkedCellContainer
 * @return ParticleContainer with particles
 */
ParticleContainer LinkedCellContainer::toContainer() {
    ParticleContainer container = ParticleContainer();
    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); p++) {
                    container.addParticle(*p);
                }
            }
        }
    }
    return container;
}

/**
 * adds all Particles of given ParticleContainer
 * @param container ParticleContainer
 */
void LinkedCellContainer::addContainer(ParticleContainer &container) {
    for (auto &particle: container) {
        addParticle(particle);
    }
}

/**
 * checks if particle needs to be moved to another cell and moves or deletes them accordingly
 */
void LinkedCellContainer::moveToNeighbour() {
    //begin at 1 and end at x_cells to avoid moving ghost cells
    for (int x = 1; x < x_cells + 1; x++) {
        for (int y = 1; y < y_cells + 1; y++) {
            for (int z = 1; z < z_cells + 1; z++) {
                for (int p = (int) cells[x][y][z].size() - 1; p >= 0; p--) {
                    //calculate new x, y, z position on cell axis
                    int x_now = floor(cells[x][y][z][p].getX()[0] / x_cell_size);
                    int y_now = floor(cells[x][y][z][p].getX()[1] / y_cell_size);
                    int z_now = floor(cells[x][y][z][p].getX()[2] / z_cell_size);
                    //check if particle is in range of inner and boundary cells
                    if ((x_now < x_cells && x_now >= 0 && cells[x][y][z][p].getX()[0] <= x_max) &&
                        (y_now < y_cells && y_now >= 0 && cells[x][y][z][p].getX()[1] <= y_max) &&
                        (z_now < z_cells && z_now >= 0 && cells[x][y][z][p].getX()[2] <= z_max)) {
                        //check if particle needs to be moved to other cell
                        if (x_now + 1 != x || y_now + 1 != y || z_now + 1 != z) {
                            addParticle(x_now + 1, y_now + 1, z_now + 1, cells[x][y][z][p]);
                            generateGhostCell((int) cells[x_now + 1][y_now + 1][z_now + 1].size() - 1, x_now + 1,
                                              y_now + 1,
                                              z_now + 1);
                            cells[x][y][z].erase(cells[x][y][z].begin() + p);
                        } else {
                            generateGhostCell(p, x, y, z);
                        }
                    } else {
                        //check if boundary is periodic and particle needs to be moved to other side
                        moveIfPeriodic(cells[x][y][z][p].getX()[0], cells[x][y][z][p].getX()[1],
                                       cells[x][y][z][p].getX()[2], cells[x][y][z][p]);
                        cells[x][y][z].erase(cells[x][y][z].begin() + p);
                    }
                }
            }
        }
    }
    it++;
}


/**
 * returns the Particles from the next neighbours of the current cell according to N3L
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 * @return
 */
std::vector<std::array<int, 3>> LinkedCellContainer::get_next_cells(int x, int y, int z) const {
    std::vector<std::array<int, 3>> vec = {};

    bool right = x < x_cells || (x == x_cells && boundary[1] == "p");
    bool up = y < y_cells || (y == y_cells && boundary[2] == "p");
    bool left = x > 1 || (x == 1 && boundary[0] == "p");
    bool before = z < z_cells || (z == z_cells && boundary[5] == "p");
    bool down = y > 1 || (y == 1 && boundary[3] == "p");

    //get neighbour cells according to N3L
    if (right) vec.push_back({x + 1, y, z});
    if (up) vec.push_back({x, y + 1, z});
    if (right && up) vec.push_back({x + 1, y + 1, z});
    if (left && up) vec.push_back({x - 1, y + 1, z});
    if (before) vec.push_back({x, y, z + 1});
    if (right && before) vec.push_back({x + 1, y, z + 1});
    if (up && before) vec.push_back({x, y + 1, z + 1});
    if (right && up && before) vec.push_back({x + 1, y + 1, z + 1});
    if (left && up && before) vec.push_back({x - 1, y + 1, z + 1});
    if (right && down && before) vec.push_back({x + 1, y - 1, z + 1});
    if (down && before) vec.push_back({x, y - 1, z + 1});
    if (left && down && before) vec.push_back({x - 1, y - 1, z + 1});
    if (left && before) vec.push_back({x - 1, y, z + 1});

    //left halo cell
    if (x == 1 && boundary[0] != "o") vec.push_back({0, y, z});
    //right halo cell
    if (x == x_cells && boundary[1] == "r") vec.push_back({x + 1, y, z});
    //below halo cell
    if (y == 1 && boundary[3] != "o") vec.push_back({x, y - 1, z});
    //below right halo cell
    if (y == 1 && x == x_cells && boundary[1] == "p" && boundary[3] == "p") vec.push_back({x + 1, y - 1, z});
    //below left halo cell
    if (y == 1 && x == 0 && boundary[0] == "p" && boundary[3] == "p") vec.push_back({x - 1, y - 1, z});
    //up halo cell
    if (y == y_cells && boundary[2] == "r") vec.push_back({x, y + 1, z});
    //before halo and normal cell
    if (z == z_cells && boundary[5] == "r") vec.push_back({x, y, z + 1});
    //behind halo cell
    if (z == 1 && boundary[4] == "r") vec.push_back({x, y, z - 1});
    //all behind halo cells if periodic
    if (z == 1 && boundary[4] == "p") {
        vec.push_back({x, y, z - 1}); //behind
        vec.push_back({x - 1, y, z - 1}); //left behind
        vec.push_back({x + 1, y, z - 1}); //right behind
        vec.push_back({x, y - 1, z - 1}); //down behind
        vec.push_back({x, y + 1, z - 1}); //up behind
        vec.push_back({x - 1, y - 1, z - 1}); //left down behind
        vec.push_back({x - 1, y + 1, z - 1}); //left up behind
        vec.push_back({x + 1, y - 1, z - 1}); //right down behind
        vec.push_back({x + 1, y + 1, z - 1}); //right up behind
    }
    return vec;
}


/**
 * sets old force to current force and current force to zero
 */
void LinkedCellContainer::setZero() {
    for (int x = 1; x <= x_cells; x++) {
        for (int y = 1; y <= y_cells; y++) {
            for (int z = 1; z <= z_cells; z++) {
                for (auto &p: cells[x][y][z]) {
                    p.setOldF(p.getF());
                    p.setF({0, 0, 0});
                }
            }
        }
    }
}

/**
 * applies the force calculation according to N3L
 * @param forceCalculation a function to apply the force calculations pairwise
 * * @param smoothedForceCalculation a function to apply the smoothed force calculations pairwise
 */
void LinkedCellContainer::applyForcePairwise(const std::function<void(Particle *, Particle *)> &forceCalculation,
                                             const std::function<void(Particle *, Particle *, double,
                                                                      double)> &smoothedForceCalculation,
                                             double Grav) {

    //begin at 1 and end at x_cells to avoid calculating the force of ghost cells
    #pragma omp parallel for collapse(3) default(none) shared(Grav, forceCalculation, smoothedForceCalculation)
    for (int x = 1; x <= x_cells; x++) {
        for (int y = 1; y <= y_cells; y++) {
            for (int z = 1; z <= z_cells; z++) {
                //get neighbour cells
                std::vector<std::array<int, 3>> neighbours = get_next_cells(x, y, z);
                for (int j = 0; j < int(cells[x][y][z].size()); j++) {
                    Particle* currentParticle = &(cells[x][y][z][j]);
                    //for all particles in current cell
                    for (int k = j + 1; k < int(cells[x][y][z].size()); k++) {
                        //calculate force with particles in current cell
                        if (!smoothed)
                            forceCalculation(currentParticle, &(cells[x][y][z][k]));
                        else
                            smoothedForceCalculation(currentParticle, &(cells[x][y][z][k]), cutoff,
                                                     smoothedRadius);
                    }
                    for (auto &neighbour: neighbours) {
                        //with neighbour cells
                        for (auto & l : cells[neighbour[0]][neighbour[1]][neighbour[2]]) {
                            //calculate force if neighbour particle is a normal particle or is the specific ghost cell to current particle
                            //if type is positive, it's a normal or a periodic ghost particle
                            //if its negative, it's a reflective ghost particle and then just the one according to the current particle should be used
                            //  -> that's the index of the current particle in the current cell negated and subtracted with one
                            if (l.getType() >= 0 ||
                                l.getType() == -j - 1) {
                                if (!smoothed)
                                    forceCalculation(currentParticle,
                                                     &l);
                                else
                                    smoothedForceCalculation(currentParticle,
                                                             &l,
                                                             cutoff, smoothedRadius);
                            }
                        }
                    }
                    //adds Ggrav force to force at the end
                    std::array<double, 3> grav = {0, currentParticle->getM() * Grav, 0};
                    currentParticle->setF(currentParticle->getF() + grav);
                }
            }
        }
    }
    deleteGhostCells();
}

/**
 * applies the mirroring if boundary is a reflection boundary
 * @param p the  index of the given particle
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 */
bool LinkedCellContainer::applyMirrorBoundary(int p, int x, int y, int z) {
    bool needs_to_be_deleted = true;
    if (cells[x][y][z][p].getX()[0] > x_max && boundary[1] == "r") {
        cells[x][y][z][p].setX(
                {x_max, cells[x][y][z][p].getX()[1], cells[x][y][z][p].getX()[2]});
        cells[x][y][z][p].setV(
                {-cells[x][y][z][p].getV()[0], cells[x][y][z][p].getV()[1], cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    } else if (cells[x][y][z][p].getX()[0] < 0 && boundary[0] == "r") {
        cells[x][y][z][p].setX({0, cells[x][y][z][p].getX()[1], cells[x][y][z][p].getX()[2]});
        cells[x][y][z][p].setV(
                {-cells[x][y][z][p].getV()[0], cells[x][y][z][p].getV()[1], cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    }
    if (cells[x][y][z][p].getX()[1] > y_max && boundary[2] == "r") {
        cells[x][y][z][p].setX({cells[x][y][z][p].getX()[0], y_max, cells[x][y][z][p].getX()[2]});
        cells[x][y][z][p].setV(
                {cells[x][y][z][p].getV()[0], -cells[x][y][z][p].getV()[1], cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    } else if (cells[x][y][z][p].getX()[1] < 0 && boundary[3] == "r") {
        cells[x][y][z][p].setX({cells[x][y][z][p].getX()[0], 0, cells[x][y][z][p].getX()[2]});
        cells[x][y][z][p].setV(
                {cells[x][y][z][p].getV()[0], -cells[x][y][z][p].getV()[1], cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    }
    if (cells[x][y][z][p].getX()[2] > z_max && boundary[5] == "r") {
        cells[x][y][z][p].setX({cells[x][y][z][p].getX()[0], cells[x][y][z][p].getX()[1], z_max});
        cells[x][y][z][p].setV(
                {cells[x][y][z][p].getV()[0], cells[x][y][z][p].getV()[1], -cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    } else if (cells[x][y][z][p].getX()[2] < 0 && boundary[4] == "r") {
        cells[x][y][z][p].setX({cells[x][y][z][p].getX()[0], cells[x][y][z][p].getX()[1], 0});
        cells[x][y][z][p].setV(
                {cells[x][y][z][p].getV()[0], cells[x][y][z][p].getV()[1], -cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    }
    return needs_to_be_deleted;
}

/**
 * generates ghost particles for given particle
 * @param index the index of the particle
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 */
void LinkedCellContainer::generateGhostCell(int index, int x, int y, int z) {
    int periodic = 0;
    double x_coordinate = cells[x][y][z][index].getX()[0];
    double y_coordinate = cells[x][y][z][index].getX()[1];
    double z_coordinate = cells[x][y][z][index].getX()[2];
    int x_new = x;
    int y_new = y;
    int z_new = z;
    bool threedim = false;
    double boundary_check = pow(2.0, 1.0 / 6.0) * cells[x][y][z][index].getSig();
    //if boundary is reflective, check if particle is nearer than 2^(1/6)*sig to boundary; one ghost particle per boundary
    //if periodic mirror particle to other sides (can be multiple ghost particles)
    if (x == 1) {
        if (boundary[0] == "r" && x_coordinate < boundary_check) {
            std::array<double, 3> ghost_x = {-cells[x][y][z][index].getX()[0],
                                             cells[x][y][z][index].getX()[1], cells[x][y][z][index].getX()[2]};
            std::array<double, 3> ghost_v = {0, 0, 0};
            addParticle(x - 1, y, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), -index - 1,
                        cells[x][y][z][index].getSig(), cells[x][y][z][index].getEps());
        }
        if (boundary[0] == "p") {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0] + x_max,
                                             cells[x][y][z][index].getX()[1], cells[x][y][z][index].getX()[2]};
            addParticle(x_cells + 1, y, z, ghost_x, cells[x][y][z][index].getV(), cells[x][y][z][index].getM(),
                        cells[x][y][z][index].getType(), cells[x][y][z][index].getSig(),
                        cells[x][y][z][index].getEps());
            x_coordinate += x_max;
            periodic++;
            x_new = x_cells + 1;
        }
    }
    if (x == x_cells) {
        if (boundary[1] == "r" && x_coordinate > (x_max - boundary_check)) {
            std::array<double, 3> ghost_x = {x_max + x_max - x_coordinate,
                                             cells[x][y][z][index].getX()[1], cells[x][y][z][index].getX()[2]};
            std::array<double, 3> ghost_v = {0, 0, 0};
            addParticle(x + 1, y, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), -index - 1,
                        cells[x][y][z][index].getSig(), cells[x][y][z][index].getEps());
        }
        if (boundary[1] == "p") {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0] - x_max,
                                             cells[x][y][z][index].getX()[1], cells[x][y][z][index].getX()[2]};
            addParticle(0, y, z, ghost_x, cells[x][y][z][index].getV(), cells[x][y][z][index].getM(),
                        cells[x][y][z][index].getType(),
                        cells[x][y][z][index].getSig(), cells[x][y][z][index].getEps());
            x_coordinate -= x_max;
            periodic++;
            x_new = 0;
        }
    }
    if (y == 1) {
        if (boundary[3] == "r" && y_coordinate < boundary_check) {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0],
                                             -cells[x][y][z][index].getX()[1],
                                             cells[x][y][z][index].getX()[2]};
            std::array<double, 3> ghost_v = {0, 0, 0};
            addParticle(x, y - 1, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), -index - 1,
                        cells[x][y][z][index].getSig(), cells[x][y][z][index].getEps());
        }
        if (boundary[3] == "p") {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0],
                                             cells[x][y][z][index].getX()[1] + y_max,
                                             cells[x][y][z][index].getX()[2]};
            addParticle(x, y_cells + 1, z, ghost_x, cells[x][y][z][index].getV(), cells[x][y][z][index].getM(),
                        cells[x][y][z][index].getType(), cells[x][y][z][index].getSig(),
                        cells[x][y][z][index].getEps());
            y_coordinate += y_max ;
            periodic++;
            y_new = y_cells + 1;
        }
    }
    if (y == y_cells) {
        if (boundary[2] == "r" && y_coordinate > (y_max - boundary_check)) {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0],
                                             y_max + y_max - y_coordinate,
                                             cells[x][y][z][index].getX()[2]};
            std::array<double, 3> ghost_v = {0, 0, 0};
            addParticle(x, y + 1, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), -index - 1,
                        cells[x][y][z][index].getSig(), cells[x][y][z][index].getEps());
        }
        if (boundary[2] == "p") {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0],
                                             cells[x][y][z][index].getX()[1] - y_max,
                                             cells[x][y][z][index].getX()[2]};
            addParticle(x, 0, z, ghost_x, cells[x][y][z][index].getV(), cells[x][y][z][index].getM(),
                        cells[x][y][z][index].getType(),
                        cells[x][y][z][index].getSig(), cells[x][y][z][index].getEps());
            y_coordinate -= y_max;
            periodic++;
            y_new = 0;
        }
    }
    if (z == 1) {
        if (boundary[4] == "r" && z_coordinate < boundary_check) {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0], cells[x][y][z][index].getX()[1],
                                             -cells[x][y][z][index].getX()[2]};
            std::array<double, 3> ghost_v = {0, 0, 0};
            addParticle(x, y, z - 1, ghost_x, ghost_v, cells[x][y][z][index].getM(), -index - 1,
                        cells[x][y][z][index].getSig(), cells[x][y][z][index].getEps());
        }
        if (boundary[4] == "p") {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0],
                                             cells[x][y][z][index].getX()[1],
                                             cells[x][y][z][index].getX()[2] + z_max};
            addParticle(x, y, z_cells + 1, ghost_x, cells[x][y][z][index].getV(), cells[x][y][z][index].getM(),
                        cells[x][y][z][index].getType(), cells[x][y][z][index].getSig(),
                        cells[x][y][z][index].getEps());
            z_coordinate += z_max;
            threedim = true;
            periodic++;
            z_new = z_cells + 1;
        }
    }
    if (z == z_cells) {
        if (boundary[5] == "r" && z_coordinate > z_max - boundary_check) {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0], cells[x][y][z][index].getX()[1],
                                             z_max + z_max - z_coordinate};
            std::array<double, 3> ghost_v = {0, 0, 0};
            addParticle(x, y, z + 1, ghost_x, ghost_v, cells[x][y][z][index].getM(), -index - 1,
                        cells[x][y][z][index].getSig(), cells[x][y][z][index].getEps());
        }
        if (boundary[5] == "p") {
            std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0],
                                             cells[x][y][z][index].getX()[1],
                                             cells[x][y][z][index].getX()[2] - z_max};
            addParticle(x, y, 0, ghost_x, cells[x][y][z][index].getV(), cells[x][y][z][index].getM(),
                        cells[x][y][z][index].getType(), cells[x][y][z][index].getSig(),
                        cells[x][y][z][index].getEps());
            z_coordinate -= z_max;
            threedim = true;
            periodic++;
            z_new = 0;
        }
    }
    if (periodic > 1) {
        //multiple generated periodic ghost particles -> particle is in corner -> mirrored in corner
        if(threedim) {
            if (periodic == 3) {
                addParticle(x_new, y_new, z, {x_coordinate, y_coordinate, cells[x][y][z][index].getX()[2]},
                            cells[x][y][z][index].getV(),
                            cells[x][y][z][index].getM(), cells[x][y][z][index].getType(), cells[x][y][z][index].getSig(),
                            cells[x][y][z][index].getEps());
                addParticle(x_new, y, z_new, {x_coordinate, cells[x][y][z][index].getX()[1], z_coordinate},
                            cells[x][y][z][index].getV(),
                            cells[x][y][z][index].getM(), cells[x][y][z][index].getType(), cells[x][y][z][index].getSig(),
                            cells[x][y][z][index].getEps());
                addParticle(x, y_new, z_new, {cells[x][y][z][index].getX()[0], y_coordinate, z_coordinate},
                            cells[x][y][z][index].getV(),
                            cells[x][y][z][index].getM(), cells[x][y][z][index].getType(), cells[x][y][z][index].getSig(),
                            cells[x][y][z][index].getEps());
            }
        }
        addParticle(x_new, y_new, z_new,{x_coordinate, y_coordinate, z_coordinate}, cells[x][y][z][index].getV(),
                    cells[x][y][z][index].getM(), cells[x][y][z][index].getType(), cells[x][y][z][index].getSig(),
                    cells[x][y][z][index].getEps());
    }
}

/**
 * deletes all ghost particles
 */
void LinkedCellContainer::deleteGhostCells() {
    for (int y = 0; y <= y_cells + 1; y++) {
        for (int z = 0; z <= z_cells + 1; z++) {
            cells[0][y][z].clear();
            cells[x_cells + 1][y][z].clear();
        }
    }
    for (int x = 0; x <= x_cells + 1; x++) {
        for (int z = 0; z <= z_cells + 1; z++) {
            cells[x][y_cells + 1][z].clear();
            cells[x][0][z].clear();
        }
    }
    for (int x = 0; x <= x_cells + 1; x++) {
        for (int y = 0; y <= y_cells + 1; y++) {
            cells[x][y][0].clear();
            cells[x][y][z_cells + 1].clear();
        }
    }
}

/**
 * moves the particle if it is behind a periodic boundary
 * @param x_coordinate x_coordinate of current particle
 * @param y_coordinate y_coordinate of current particle
 * @param z_coordinate z_coordinate of current particle
 * @param p Particle to be moved
 */
void LinkedCellContainer::moveIfPeriodic(double x_coordinate, double y_coordinate, double z_coordinate, Particle &p) {
    bool periodic = false;

    double oldX = p.getOldX()[0];
    double oldY = p.getOldX()[1];
    double oldZ = p.getOldX()[2];

    if (x_coordinate > x_max && boundary[1] == "p") {
        periodic = true;
        oldX -= x_max;
        x_coordinate -= x_max;
        if(x_coordinate > x_max) x_coordinate = std::fmod(x_coordinate, x_max); //0;
    }
    else if (x_coordinate < 0 && boundary[0] == "p") {
        periodic = true;
        oldX += x_max;
        x_coordinate += x_max;
        if(x_coordinate < 0) x_coordinate = x_max - std::fmod(-x_coordinate, x_max);//x_max;
    }
    if (y_coordinate > y_max && boundary[2] == "p") {
        periodic = true;
        oldY -= y_max;
        y_coordinate -= y_max;
        if(y_coordinate > y_max ) y_coordinate = std::fmod(y_coordinate, y_max); // 0;
    }
    else if (y_coordinate < 0 && boundary[3] == "p") {
        periodic = true;
        y_coordinate += y_max;
        oldY += y_max;
        if(y_coordinate < 0 ) y_coordinate = y_max - std::fmod(-y_coordinate, y_max); //y_max;
    }
    if (z_coordinate > z_max && boundary[5] == "p") {
        periodic = true;
        z_coordinate -= z_max;
        oldZ -= z_max;
        if(z_coordinate > z_max) z_coordinate = std::fmod(z_coordinate, z_max); // 0;
    }
    else if (z_coordinate < 0 && boundary[4] == "p") {
        periodic = true;
        z_coordinate += z_max;
        oldZ += z_max;
        if(z_coordinate < 0) z_coordinate = z_max - std::fmod(-z_coordinate, z_max); //z_max;
    }
    if (x_coordinate > x_max || x_coordinate < 0 || y_coordinate > y_max || y_coordinate < 0 || z_coordinate > z_max ||
        z_coordinate < 0) {
        spdlog::info("needs to be deleted because x {} y {} z {}", x_coordinate, y_coordinate, z_coordinate);
        periodic = false;
        return;
    }

    if (periodic) {
        Particle new_particle({x_coordinate, y_coordinate, z_coordinate}, p.getV(), p.getM(), p.getSig(),
                              p.getEps(), p.getType());

        new_particle.setOldX({oldX, oldY, oldZ});

        int x = (int) floor(x_coordinate / x_cell_size) + 1;
        if(x_coordinate == x_max) x = x_cells;
        if(x_coordinate == 0) x = 1;
        int y = (int) floor(y_coordinate / y_cell_size) + 1;
        if(y_coordinate == y_max) y = y_cells;
        if(y_coordinate == 0) y = 1;
        int z = (int) floor(z_coordinate / z_cell_size) + 1;
        if(z_coordinate == z_max) z = z_cells;
        if(z_coordinate == 0) z = 1;

        addParticle(x,y,z,new_particle);
        generateGhostCell((int) cells[x][y][z].size() - 1, x, y, z);
        return;
    }
}

/**
 * calculates the average movement distance (diffusion) of all particles
 * @return the average movement distance
 */
double LinkedCellContainer::calculateDiffusion() {
    double var = 0;
    int particles = 0;
    for (int x = 1; x <= x_cells; x++) {
        for (int y = 1; y <= y_cells; y++) {
            for (int z = 1; z <= z_cells; z++) {
                for (const auto & p : cells[x][y][z]) {
                    var += pow(ArrayUtils::L2Norm(p.getX() - p.getOldX()), 2);
                    particles++;
                }
            }
        }
    }
    var /= particles;
    return var;
}

/**
 * calculates the Radial Distribution Function
 * @param intervalBegin
 * @param intervalEnd
 * @param deltaR
 * @param x_axis_plot plot to add results to
 * @param filename name of file to save results into
 * @return densities
 */
std::vector<double> LinkedCellContainer::calculateRDF(int intervalBegin, int intervalEnd, double deltaR, std::vector<int> x_axis_plot,
                                  std::ofstream &filename) {
    std::vector<double> densities;
    ParticleContainer particles = toContainer();
    auto first = particles.begin();
    auto last = particles.end();
    for (int i = intervalBegin; i <= intervalEnd - deltaR; i++) {
        int num_particles = 0;
        for (auto &p1: particles) {
            for (auto &p2: particles) {
                double distance = ArrayUtils::L2Norm(p1.getX() - p2.getX());
                if (distance >= i && distance <= i + deltaR && !(p1 == p2)) num_particles++;
            }
        }
        double new_density = num_particles / ((4 * M_PI / 3) * (pow(i + deltaR, 3) - pow(i, 3)));
        filename << new_density << ", ";
        densities.emplace_back(new_density);
    }
    filename << "\n";
    matplot::plot(std::move(x_axis_plot), densities);
    matplot::hold(matplot::on);
    return densities;
}

/**
 * return size of cells on x_axis
 * @return size of cells on x_axis
 */
double LinkedCellContainer::getXCellSize() const { return x_cell_size; }

/**
 * return size of cells on y_axis
 * @return size of cells on y_axis
 */
double LinkedCellContainer::getYCellSize() const { return y_cell_size; }

/**
 * return size of cells on z_axis
 * @return size of cells on z_axis
 */
double LinkedCellContainer::getZCellSize() const { return z_cell_size; }

/**
 * return number of cells on x_axis
 * @return number of cells on x_axis
 */
int LinkedCellContainer::getXCells() const {
    return x_cells;
}

/**
 * return number cells on y_axis
 * @return number cells on y_axis
 */
int LinkedCellContainer::getYCells() const {
    return y_cells;
}

/**
 * return number of cells on z_axis
 * @return number of cells on z_axis
 */
int LinkedCellContainer::getZCells() const {
    return z_cells;
}

/**
 * returns a specific cell
 * @param x x-coordinate of cell
 * @param y y-coordinate of cell
 * @param z z-coordinate of cell
 * @return cell
 */
std::vector<Particle> LinkedCellContainer::getCell(int x, int y, int z) {
    return cells[x][y][z];
}

/**
 * returns cutoff radius
 * @return cutoff
 */
double LinkedCellContainer::getCutoff() const { return cutoff; }

/**
 * returns the vector of all particles (including ghost particles) in the ParticleContainer with Pointer at first element
 * @return start pointer
 */
std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator LinkedCellContainer::begin() {
    return cells.begin();
}

/**
 * returns the vector of all particles (including ghost particles) in the LinkedCellContainer with pointer at last element
 * @return end pointer
 */
std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator LinkedCellContainer::end() {
    return cells.end();
}