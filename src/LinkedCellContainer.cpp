//
// Created by kathi on 20.11.23.
//


#include "LinkedCellContainer.h"
#include "calculations/ForceCalculator.h"
#include <spdlog/spdlog.h>
#include <cmath>
#include "utils/ArrayUtils.h"

/**
 * create a cell grid wcells[i][j]h the given numbers o of cells
 * @param number_of_cells
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

LinkedCellContainer::~LinkedCellContainer() {}

/**
 * returns the number of cells in the grid
 * @return
 */
int LinkedCellContainer::cell_numbers() const {
    return x_cells * y_cells * z_cells;
}

/**
 * returns the number od molecules of the given cell
 * @param cell
 * @return
 */
int LinkedCellContainer::Particles_in_cell(int x, int y, int z) {
    return cells[x][y][z].size();
}

/**
 * adds new Particle to specific cell
 * @param cell
 * @param x_arg
 * @param v_arg
 * @param m_arg
 * @param type_arg
 */
void LinkedCellContainer::addParticle(int x, int y, int z, std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
                                      int type_arg) {
    cells[x][y][z].push_back(Particle(x_arg, v_arg, m_arg, type_arg));
}

void
LinkedCellContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg) {
    cells[floor(x_arg[0] / c)][floor(x_arg[1] / c)][floor(x_arg[2] / c)].push_back(
            Particle(x_arg, v_arg, m_arg, type_arg));
}

/**
 * adds existing particle to specific cell
 * @param cell
 * @param p
 */
void LinkedCellContainer::addParticle(int x, int y, int z, Particle &p) {
    cells[x][y][z].push_back(p);
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
 * checks if particle needs to be moved to another cell
 */
void LinkedCellContainer::moveToNeighbour() {
    for (int x = 0; x < x_cells; x++) {
        for (int y = 0; y < y_cells; y++) {
            for (int z = 0; z < z_cells; z++) {
                for (int p = 0; p < cells[x][y][z].size(); p++) {
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
    }
    return vec;
}

void LinkedCellContainer::setZero() {
    for (int x = 0; x < x_cells; x++) {
        for (int y = 0; y < y_cells; y++) {
            for(int z = 0; z < z_cells; z++){
                for(int p = 0; p < cells[x][y][z].size(); p++){
                    cells[x][y][z][p].setOldF(cells[x][y][z][p].getF());
                    cells[x][y][z][p].setF({0,0,0});
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
 * returns the vector of all PArticles in the ParticleContainer with Pointer at last element
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
                    for (int n = 0; n < neighbours.size(); n++) {
                        //with neighbour cells
                        for (int l = 0; l < cells[neighbours[n][0]][neighbours[n][1]][neighbours[n][2]].size(); l++) {
                            forceCalculation(&(cells[x][y][z][j]), &(cells[neighbours[n][0]][neighbours[n][1]][neighbours[n][2]][l]));
                        }
                    }
                }
            }
        }
    }
}


