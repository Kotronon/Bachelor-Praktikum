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
LinkedCellContainer::LinkedCellContainer(std::array<int, 3> N, double cutoff, std::array<std::string, 6> b) {
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
    three_dim = N[2] > 1;
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
int LinkedCellContainer::Particles_in_cell(int x, int y, int z) {
    return cells[x][y][z].size();
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

void
LinkedCellContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg) {
    Particle new_particle = Particle(x_arg, v_arg, m_arg, type_arg);
    cells[floor(x_arg[0] / c)][floor(x_arg[1] / c)][floor(x_arg[2] / c)].emplace_back(new_particle);
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
 * delete Particle from specific cell
 * @param x index of cell on x axis
 * @param y index of cell on y axis
 * @param z index of cell on z axis
 * @param p particle to deleted
 */
void LinkedCellContainer::deleteParticle(int x, int y, int z, Particle &p) {
    unsigned int pos = 0;
    while (pos < cells[x][y][z].size()) {
        if (cells[x][y][z][pos].operator==(p)) {
            cells[x][y][z].erase(cells[x][y][z].begin() + pos);
            return;
        }
        pos++;
    }
}

/**
 * checks if particle needs to be moved to another cell and moves them accordingly
 */
void LinkedCellContainer::moveToNeighbour() {
    for (int x = 1; x < x_cells + 1; x++) {
        for (int y = 1; y < y_cells + 1; y++) {
            for (int z = 1; z < z_cells + 1; z++) {
                for (int p = 0; p < cells[x][y][z].size(); p++) {
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
                       // if(applyMirrorBoundary(p, x, y, z)) {
                       if(needsToBeDeleted(x_now, y_now, z_now)){
                            cells[x][y][z].erase(cells[x][y][z].begin() + p);
                        }
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
    bool before = z < z_cells && three_dim;
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
    if (x == 1 && boundary[0] == "r") vec.push_back({0, y, z});
    //right halo or normal cell
    if (x == x_cells && boundary[0] != "o") vec.push_back({x+1, y, z});
    //below halo cell
    if (y == 1 && boundary[3] == "r") vec.push_back({x, y - 1, z});
    //before halo and normal cell
    if (z == z_cells && three_dim && boundary[5] != "o" ) vec.push_back({x, y, z + 1});
    //before left normal cells
    if (z == z_cells && three_dim && boundary[5] == "n" && x > 0) vec.push_back({x-1, y, z + 1});
    //before right normal cells
    if (z == z_cells && three_dim && boundary[5] == "n" && x <= x_cells) vec.push_back({x+1, y, z + 1});
    //behind halo cell
    if (z == 1 && three_dim && boundary[4] == "r") vec.push_back({x, y, z - 1});

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
 * @param forceCalculationforceCalculation a function to apply the  force calculations pairwise
 */
void LinkedCellContainer::applyForcePairwise(const std::function<void(Particle *, Particle *)> &forceCalculation) {
    for (int x = 0; x <= x_cells - 1; x++) {
        for (int y = 0; y <= y_cells - 1; y++) {
            for (int z = 0; z <= z_cells-1; z++) {
                std::vector<std::array<int, 3>> neighbours = get_next_cells(x, y, z);
                for (int j = 0; j < cells[x][y][z].size(); j++) {
                    //for all particles in current cell
                    if(cells[x][y][z][j].getType() == 0) {
                        for (int k = j + 1; k < cells[x][y][z].size(); k++) {
                            //calculate force with particles in current cell
                            forceCalculation(&(cells[x][y][z][j]), &(cells[x][y][z][k]));
                        }
                        for (auto & neighbour : neighbours) {
                            //with neighbour cells
                            for (int l = 0;
                                 l < cells[neighbour[0]][neighbour[1]][neighbour[2]].size(); l++) {
                                if (cells[neighbour[0]][neighbour[1]][neighbour[2]][l].getType() == 0 ||
                                    cells[neighbour[0]][neighbour[1]][neighbour[2]][l].getType() == j) {
                                    forceCalculation(&(cells[x][y][z][j]),
                                                     &(cells[neighbour[0]][neighbour[1]][neighbour[2]][l]));
                                }
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
    if (floor(cells[x][y][z][p].getX()[0] / c) > x_cells - 1 && boundary[1] == "r") {
        cells[x][y][z][p].setX(
                {0.0 + c * x_cells - 0.0000001, cells[x][y][z][p].getX()[1], cells[x][y][z][p].getX()[2]});
        cells[x][y][z][p].setV(
                {-cells[x][y][z][p].getV()[0], cells[x][y][z][p].getV()[1], cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    } else if (floor(cells[x][y][z][p].getX()[0] / c) < 0 && boundary[0] == "r") {
        cells[x][y][z][p].setX({0, cells[x][y][z][p].getX()[1], cells[x][y][z][p].getX()[2]});
        cells[x][y][z][p].setV(
                {-cells[x][y][z][p].getV()[0], cells[x][y][z][p].getV()[1], cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    }
    if (floor(cells[x][y][z][p].getX()[1] / c) > y_cells -1 && boundary[2] == "r") {
        cells[x][y][z][p].setX({cells[x][y][z][p].getX()[0], 0.0 + c * y_cells - 0.0000001, cells[x][y][z][p].getX()[2]});
        cells[x][y][z][p].setV(
                {cells[x][y][z][p].getV()[0], -cells[x][y][z][p].getV()[1], cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    } else if (floor(cells[x][y][z][p].getX()[1] / c) < 0 && boundary[3] == "r") {
        cells[x][y][z][p].setX({cells[x][y][z][p].getX()[0], 0, cells[x][y][z][p].getX()[2]});
        cells[x][y][z][p].setV(
                {cells[x][y][z][p].getV()[0], -cells[x][y][z][p].getV()[1], cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    }
    if (floor(cells[x][y][z][p].getX()[2] / c) > z_cells - 1 && boundary[5] == "r" && three_dim) {
        cells[x][y][z][p].setX({cells[x][y][z][p].getX()[0], cells[x][y][z][p].getX()[1], 0.0 + c * z_cells - 0.0000001});
        cells[x][y][z][p].setV(
                {cells[x][y][z][p].getV()[0], cells[x][y][z][p].getV()[1], -cells[x][y][z][p].getV()[2]});
        needs_to_be_deleted = false;
    } else if (floor(cells[x][y][z][p].getX()[2] / c) < 0 && boundary[4] == "r" && three_dim) {
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
        std::array<double, 3> ghost_x = {-cells[x][y][z][index].getX()[0], cells[x][y][z][index].getX()[1],
                                         cells[x][y][z][index].getX()[2]};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x - 1, y, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), index);
    }
    if (x == x_cells && boundary[1] == "r" && cells[x][y][z][index].getV()[0] > 0) {
        std::array<double, 3> ghost_x = {x_cells * c + fmod((cells[x][y][z][index].getX()[0]), c),
                                         cells[x][y][z][index].getX()[1], cells[x][y][z][index].getX()[2]};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x + 1, y, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), index);
    }
    if (y == 1 && boundary[3] == "r" && cells[x][y][z][index].getV()[1] < 0) {
        std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0], -cells[x][y][z][index].getX()[1],
                                         cells[x][y][z][index].getX()[2]};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x, y - 1, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), index);
    }
    if (y == y_cells && boundary[2] == "r" && cells[x][y][z][index].getV()[1] > 0) {
        std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0],
                                         y_cells * c + fmod(cells[x][y][z][index].getX()[0], c),
                                         cells[x][y][z][index].getX()[2]};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x, y + 1, z, ghost_x, ghost_v, cells[x][y][z][index].getM(), index);
    }
    if (z == 1 && boundary[4] == "r" && three_dim && cells[x][y][z][index].getV()[2] < 0) {
        std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0], cells[x][y][z][index].getX()[1],
                                         -cells[x][y][z][index].getX()[2]};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x, y, z - 1, ghost_x, ghost_v, cells[x][y][z][index].getM(), index);
    }
    if (z == z_cells && boundary[5] == "r" && three_dim && cells[x][y][z][index].getV()[2] > 0) {
        std::array<double, 3> ghost_x = {cells[x][y][z][index].getX()[0], cells[x][y][z][index].getX()[1],
                                         z_cells * c + fmod(cells[x][y][z][index].getX()[0], c)};
        std::array<double, 3> ghost_v = {0,0,0};
        addParticle(x, y, z + 1, ghost_x, ghost_v, cells[x][y][z][index].getM(), index);
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
    bool needsToBeDeleted = true;
    if ((x_coordinate >= (x_cells + 1) * c && boundary[1] == "n") || (x_coordinate < 0 && boundary[0] == "n")
        || (y_coordinate >= (y_cells + 1) * c && boundary[2] == "n") || (y_coordinate < 0 && boundary[3] == "n")
        || (z_coordinate >= (z_cells + 1) * c && boundary[5] == "n") || (z_coordinate < 0 && boundary[4] == "n"))
        needsToBeDeleted = false;
    return needsToBeDeleted;
}