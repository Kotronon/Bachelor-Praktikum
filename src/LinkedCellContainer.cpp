//
// Created by kathi on 20.11.23.
//


#include "LinkedCellContainer.h"
#include "calculations/ForceCalculator.h"
#include <spdlog/spdlog.h>
#include <cmath>
#include "utils/ArrayUtils.h"

/**
 * create a new linked cell container based on the dimensions and boundary
 * @param N dimensions of the container
 * @param cutoff radius of cutoff
 * @param b boundary types of each side of the container
 */
LinkedCellContainer::LinkedCellContainer(std::array<int, 3> N, double cutoff, std::vector<std::string> b) {
    //creating list wcells[i][j]h length = number of cells
    x_cells = ceil(N[0] / cutoff);
    y_cells = ceil(N[1] / cutoff);
    z_cells = ceil(N[2] / cutoff);
    std::vector<std::vector<std::vector<std::vector<Particle>>>> x;
    for(int i = 0; i < x_cells; i++){
        std::vector<std::vector<std::vector<Particle>>> y;
        for(int j = 0; j < y_cells; j++){
            std::vector<std::vector<Particle>> z (z_cells);
            y.push_back(z);
        }
        x.push_back(y);
    }
    cells = x;
    c = cutoff;
    boundary = b;
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
 * @param cell
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
void LinkedCellContainer::addParticle(int x, int y, int z, std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
                                      int type_arg) {
    cells[x][y][z].emplace_back(x_arg, v_arg, m_arg, type_arg);
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
 * @param cell
 * @param p
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
    for (int x = 0; x < x_cells; x++) {
        for (int y = 0; y < y_cells; y++) {
            for (int z = 0; z < z_cells; z++) {
                for (int p = 0; p < (int) cells[x][y][z].size(); p++) {
                    int x_now = floor(cells[x][y][z][p].getX()[0] / c);
                    int y_now = floor(cells[x][y][z][p].getX()[1] / c);
                    int z_now = floor(cells[x][y][z][p].getX()[2] / c);
                    if (x_now != x || y_now != y || z_now != z) {
                        if ((x_now < x_cells && x_now >= 0) && (y_now < y_cells && y_now >= 0) && (z_now < z_cells && z_now >= 0)) {
                            addParticle(x_now, y_now, z_now, cells[x][y][z][p]);
                            cells[x][y][z].erase(cells[x][y][z].begin() + p);
                            //check boundary conditions -> create ghostcells

                        } else {
                            //is outflow boundary
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
 * @param cell
 * @return
 */
std::vector<std::array<int, 3>> LinkedCellContainer::get_next_cells(int x, int y, int z) const {
    std::vector<std::array<int, 3>> vec = {};
    bool right = x < x_cells -1;
    bool up = y < y_cells - 1;
    bool left = x > 0;
    bool before = z < z_cells - 1;
    bool after = z > 0;
    //2D
    if (right) vec.push_back({x+1, y, z});
    if (up) vec.push_back({x, y+1, z});
    if (right && up) vec.push_back({x+1, y+1, z});
    if (left && up) vec.push_back({x-1, y+1, z});
    //3D
    if (z_cells > 1) {
        if (before) vec.push_back({x, y, z+1});
        if (right && before) vec.push_back({x+1, y, z+1});
        if (up && before) vec.push_back({x, y+1, z+1});
        if (right && up && before) vec.push_back({x+1, y+1, z+1});
        if (left && up && before) vec.push_back({x-1, y+1, z+1});

        if (after) vec.push_back({x, y, z-1});
        if (right && after) vec.push_back({x+1, y, z-1});
        if (up && after) vec.push_back({x, y+1, z-1});
        if (right && up && after) vec.push_back({x+1, y+1, z-1});
        if (left && up && after) vec.push_back({x-1, y+1, z-1});
    }
    return vec;
}

void LinkedCellContainer::setZero() {
    for (int x = 0; x < x_cells; x++) {
        for (int y = 0; y < y_cells; y++) {
            for(int z = 0; z < z_cells; z++){
                for(auto & p : cells[x][y][z]){
                    p.setOldF(p.getF());
                    p.setF({0,0,0});
                }
            }
        }
    }
}

int LinkedCellContainer::getXMax() const { return x_cells; }

int LinkedCellContainer::getYMax() const { return y_cells; }

int LinkedCellContainer::getZMax() const { return z_cells; }

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

void LinkedCellContainer::applyForcePairwise(const std::function<void(Particle *, Particle *)> &forceCalculation) {
    for (int x = 0; x < x_cells; x++) {
        for (int y = 0; y < y_cells; y++) {
            for (int z = 0; z < z_cells; z++) {
                std::vector<std::array<int, 3>> neighbours = get_next_cells(x, y, z);
                for (int j = 0; j < cells[x][y][z].size(); j++) {
                    //for all particles in current cell
                    for (int k = j + 1; k < cells[x][y][z].size(); k++) {
                        //calculate force with particles in current cell
                        forceCalculation(&(cells[x][y][z][j]), &(cells[x][y][z][k]));
                    }
                    for (auto & neighbour : neighbours) {
                        //with neighbour cells
                        for (int l = 0; l < cells[neighbour[0]][neighbour[1]][neighbour[2]].size(); l++) {
                            forceCalculation(&(cells[x][y][z][j]), &(cells[neighbour[0]][neighbour[1]][neighbour[2]][l]));
                        }
                    }
                }
            }
        }
    }
}


