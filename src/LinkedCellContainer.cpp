//
// Created by kathi on 20.11.23.
//


#include "LinkedCellContainer.h"
#include "calculations/ForceCalculator.h"
#include <spdlog/spdlog.h>
#include <cmath>
#include <utility>
#include "utils/ArrayUtils.h"

/**
 * create a new linked cell container based on the dimensions and boundary
 * @param N dimensions of the container
 * @param cutoff radius of cutoff
 * @param b boundary types of each side of the container
 */
LinkedCellContainer::LinkedCellContainer(std::array<double, 3> N, double cutoff, std::array<std::string, 6> b) {
    //creating list cells[i][j]h length = number of cells
    x_cells = ceil(N[0] / cutoff);
    y_cells = ceil(N[1] / cutoff);
    z_cells = ceil(N[2] / cutoff);
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
    c = cutoff;
    boundary = std::move(b);
    x_max = N[0];
    y_max = N[1];
    z_max = N[2];
}

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
    return cells[x + 1][y + 1][z + 1].size();
}

/**
 * adds new Particle to specific cell
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 * @param x_arg coordinates of new particle
 * @param v_arg velocity of new particle
 * @param m_arg mass of new particle
 * @param type_arg type of new particle
 */
void LinkedCellContainer::addParticle(int x, int y, int z, std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                                      double m_arg,
                                      int type_arg) {
    Particle new_particle = Particle(x_arg, v_arg, m_arg, type_arg);
    cells[x][y][z].emplace_back(new_particle);
}

/**
 * adds new Particle to a cell based on current position
 * @param x_arg coordinates of new particle
 * @param v_arg velocity of new particle
 * @param m_arg mass of new particle
 * @param type_arg type of new particle
 */
void
LinkedCellContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg) {
    Particle new_particle = Particle(x_arg, v_arg, m_arg, type_arg);
    cells[floor(x_arg[0] / c) + 1][floor(x_arg[1] / c) + 1][floor(x_arg[2] / c) + 1].emplace_back(new_particle);
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
    cells[floor(p.getX()[0] / c) +1][floor(p.getX()[1] / c)+1][floor(p.getX()[2] / c)+1].emplace_back(p);
}



/**
 * checks if particle needs to be moved to another cell and moves them accordingly
 */
void LinkedCellContainer::moveToNeighbour() {
    for (int x = 1; x < x_cells + 1; x++) {
        for (int y = 1; y < y_cells + 1; y++) {
            for (int z = 1; z < z_cells + 1; z++) {
                for (int p = cells[x][y][z].size() -1; p >= 0; p--) {
                    int x_now = floor(cells[x][y][z][p].getX()[0] / c);
                    int y_now = floor(cells[x][y][z][p].getX()[1] / c);
                    int z_now = floor(cells[x][y][z][p].getX()[2] / c);
                    if ((x_now < x_cells && x_now >= 0) && (y_now < y_cells && y_now >= 0) &&
                        (z_now < z_cells && z_now >= 0)) {
                        if (x_now + 1 != x || y_now + 1 != y || z_now + 1 != z) {
                            addParticle(x_now+1, y_now+1, z_now+1, cells[x][y][z][p]);
                            cells[x][y][z].erase(cells[x][y][z].begin() + p);
                        }
                        else{
                            generateGhostCell(p, x, y, z);
                        }
                       //
                        //check boundary conditions -> create ghostcell
                    } else {
                        //is outflow boundary
                        //if(applyMirrorBoundary(p, x, y, z))
                       if(needsToBeDeleted(cells[x][y][z][p].getX()[0], cells[x][y][z][p].getX()[1], cells[x][y][z][p].getX()[2]))
                            cells[x][y][z].erase(cells[x][y][z].begin() + p);
                        //}
                    }
                }
            }
        }
    }
}


/**
 * returns the Particles from the next neighbours of the current cell
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 * @return
 */
std::vector<std::array<int, 3>> LinkedCellContainer::get_next_cells(int x, int y, int z) const {
    std::vector<std::array<int, 3>> vec = {};
    bool right = x < x_cells;
    bool up = y < y_cells;
    bool left = x > 1;
    bool before = z < z_cells;
    bool down = y > 1;

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
    if (x == 1) vec.push_back({0, y, z});
    //right halo cell
    if (x == x_cells) vec.push_back({x+1, y, z});
    //below halo cell
    if (y == 1) vec.push_back({x, y - 1, z});
    //up halo cell
    if (y == y_cells) vec.push_back({x, y+1, z});
    //before halo and normal cell
    if (z == z_cells) vec.push_back({x, y, z + 1});
    //behind halo cell
    if (z == 1) vec.push_back({x, y, z - 1});

    return vec;
}
/**
 * sets old force to current force and current force to zero
 */
void LinkedCellContainer::setZero() {
    for (int x = 1; x <= x_cells; x++) {
        for (int y = 1; y <= y_cells; y++) {
            for (int z = 1; z <= z_cells; z++) {
                for (auto & p : cells[x][y][z]) {
                    p.setOldF(p.getF());
                    p.setF({0, 0, 0});
                }
            }
        }
    }
}

int LinkedCellContainer::getXMax() const { return x_cells; }

int LinkedCellContainer::getYMax() const { return y_cells; }

int LinkedCellContainer::getZMax() const { return z_cells; }

/**
 * returns the vector of all particles in the ParticleContainer with Pointer at first element
 * @return
 */
std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator LinkedCellContainer::begin() {
    return cells.begin();
}

/**
 * returns the vector of all particles in the ParticleContainer with Pointer at last element
 * @return
 */
std::vector<std::vector<std::vector<std::vector<Particle>>>>::iterator LinkedCellContainer::end() {
    return cells.end();
}

/**
 * applies the force calculation according to N3L
 * @param forceCalculation a function to apply the force calculations pairwise
 */
void LinkedCellContainer::applyForcePairwise(const std::function<void(Particle *, Particle *)> &forceCalculation) {
    for (int x = 1; x <= x_cells ; x++) {
        for (int y = 1; y <= y_cells; y++) {
            for (int z = 1; z <= z_cells; z++) {
                std::vector<std::array<int, 3>> neighbours = get_next_cells(x, y, z);
                for (int j = 0; j < cells[x][y][z].size(); j++) {
                    //for all particles in current cell
                        for (int k = j + 1; k < cells[x][y][z].size(); k++) {
                            //calculate force with particles in current cell
                            forceCalculation(&(cells[x][y][z][j]), &(cells[x][y][z][k]));
                        }
                        for (auto & neighbour : neighbours) {
                            //with neighbour cells
                            for (int l = 0;
                                 l < cells[neighbour[0]][neighbour[1]][neighbour[2]].size(); l++) {

                                if (cells[neighbour[0]][neighbour[1]][neighbour[2]][l].getType() == 0 ||
                                    cells[neighbour[0]][neighbour[1]][neighbour[2]][l].getType() == j+1) {
                                    forceCalculation(&(cells[x][y][z][j]),
                                                     &(cells[neighbour[0]][neighbour[1]][neighbour[2]][l]));
                                }
                            }
                        }

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
    if (cells[x][y][z][p].getX()[0]  > x_max && boundary[1] == "r") {
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
    } else if (cells[x][y][z][p].getX()[1]  < 0 && boundary[3] == "r") {
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
 * generates ghost cells for given particle
 * @param index the index of the particle
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 */
void LinkedCellContainer::generateGhostCell(int index, int x, int y, int z) {
    if (x == 1 && boundary[0] == "r" && cells[x][y][z][index].getV()[0] < 0) {
        std::array<double, 3> ghost_x = {-cells[x][y][z][index].getX()[0] - 0.0000000001, cells[x][y][z][index].getX()[1], cells[x][y][z][index].getX()[2]};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x - 1, y, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), index+1);
    }
    if (x == x_cells && boundary[1] == "r" && cells[x][y][z][index].getV()[0] > 0) {
        std::array<double, 3> ghost_x = {x_max + c-fmod(cells[x][y][z][index].getX()[0], c) + 0.0000000001, cells[x][y][z][index].getX()[1], cells[x][y][z][index].getX()[2]};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x + 1, y, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), index+1);
    }
    if (y == 1 && boundary[3] == "r" && cells[x][y][z][index].getV()[1] < 0) {
        std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0], -cells[x][y][z][index].getX()[1] - 0.0000000001, cells[x][y][z][index].getX()[2]};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x, y - 1, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), index+1);
    }
    if (y == y_cells && boundary[2] == "r" && cells[x][y][z][index].getV()[1] > 0) {
        std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0],  y_max + c-fmod(cells[x][y][z][index].getX()[1], c) + 0.0000000001, cells[x][y][z][index].getX()[2]};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x, y + 1, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), index+1);
    }
    if (z == 1 && boundary[4] == "r" && cells[x][y][z][index].getV()[2] < 0) {
        std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0], cells[x][y][z][index].getX()[1], -cells[x][y][z][index].getX()[2] - 0.0000000001};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x, y, z - 1, ghost_x, ghost_v, cells[x][y][z][index].getM(), index+1);
    }
    if (z == z_cells && boundary[5] == "r" && cells[x][y][z][index].getV()[2] > 0) {
        std::array<double, 3> ghost_x = {cells[x
        ][y][z][index].getX()[0], cells[x][y][z][index].getX()[1], z_max + c-fmod(cells[x][y][z][index].getX()[2], c) + 0.0000000001};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x, y, z + 1, ghost_x, ghost_v, cells[x][y][z][index].getM(), index+1);
    }
}

/**
 * deletes all ghost cells
 */
void LinkedCellContainer::deleteGhostCells() {
    if (boundary[0] == "r") {
        for (int y = 1; y <= y_cells; y++) {
            for (int z = 1; z <= z_cells; z++) {
                cells[0][y][z].clear();
            }
        }
    }
    if (boundary[1] == "r") {
        for (int y = 1; y <= y_cells; y++) {
            for (int z = 1; z <= z_cells; z++) {
                cells[x_cells + 1][y][z].clear();
            }
        }
    }
    if (boundary[2] == "r") {
        for (int x = 1; x <= x_cells; x++) {
            for (int z = 1; z <= z_cells; z++) {
                cells[x][y_cells + 1][z].clear();
            }
        }
    }
    if (boundary[3] == "r") {
        for (int x = 1; x <= x_cells; x++) {
            for (int z = 1; z <= z_cells; z++) {
                cells[x][0][z].clear();
            }
        }
    }
    if (boundary[4] == "r") {
        for (int x = 1; x <= x_cells; x++) {
            for (int y = 1; y <= y_cells; y++) {
                cells[x][y][0].clear();
            }
        }
    }
    if (boundary[5] == "r") {
        for (int x = 1; x <= x_cells; x++) {
            for (int y = 1; y <= y_cells; y++) {
                cells[x][y][z_cells + 1].clear();
            }
        }
    }
}

/**
 * returns true if particle needs to be deleted
 * @param x_coordinate x_coordinate of current particle
 * @param y_coordinate y_coordinate of current particle
 * @param z_coordinate z_coordinate of current particle
 * @return
 */
bool LinkedCellContainer::needsToBeDeleted(double x_coordinate, double y_coordinate, double z_coordinate) {
    bool needsToBeDeleted = false;
    if ((x_coordinate > x_max) || (x_coordinate < 0)
        || (y_coordinate > y_max) || (y_coordinate < 0)
        || (z_coordinate > z_max) || (z_coordinate < 0 ))
        needsToBeDeleted = true;
    return needsToBeDeleted;
}