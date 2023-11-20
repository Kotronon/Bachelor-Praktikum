//
// Created by kathi on 20.11.23.
//


#include "LinkedCellContainer.h"
#include <math.h>
#include "calculations/ForceCalculator.h"
#include <spdlog/spdlog.h>
#include <cmath>
#include "utils/ArrayUtils.h"

/**
 * create a cell grid with the given numbers o of cells
 * @param number_of_cells
 */
LinkedCellContainer::LinkedCellContainer(std::array<int, 3> N, double cutoff) {
    //creating list with length = number of cells
    x_cells = ceil(N[0] / cutoff) + 1;
    y_cells = ceil(N[1] / cutoff) + 1;
    z_cells = ceil(N[2] / cutoff) + 1;
    cells = std::vector<std::vector<Particle>>(x_cells*y_cells*z_cells);
    c = cutoff;
}

LinkedCellContainer::~LinkedCellContainer(){}

/**
 * returns the number of cells in the grid
 * @return
 */
int LinkedCellContainer::cell_numbers() {
    return cells.size();
}

/**
 * returns the number od molecules of the given cell
 * @param cell
 * @return
 */
int LinkedCellContainer::Particles_in_cell(int cell) {
    return cells[cell].size();
}

/**
 * adds new Particle to specific cell
 * @param cell
 * @param x_arg
 * @param v_arg
 * @param m_arg
 * @param type_arg
 */
void LinkedCellContainer::addParticle(int cell, std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg){
    cells[cell].push_back(Particle(x_arg, v_arg, m_arg, type_arg));
}

/**
 * adds existing particle to specific cell
 * @param cell
 * @param p
 */
void LinkedCellContainer::addParticle(int cell, Particle &p){
    cells[cell].push_back(p);
}

/**
 * delets Particle from specific cell
 * @param cell
 * @param p
 */
void LinkedCellContainer::deleteParticle(int cell, Particle &p){
    int pos = 0;
    while(pos < cells[cell].size()){
        if(cells[cell][pos].operator==(p)){
            cells[cell].erase(cells[cell].begin() + pos);
            return;
        }
        pos++;
    }
}

/**
 * moves Particle to neighbour cell
 * @param cell_current
 * @param cell_new
 * @param p
 */
void LinkedCellContainer::moveToNeighbour(int cell_current, int cell_new, Particle &p){
    deleteParticle(cell_current, p);
    addParticle(cell_new, p);
}

/**
 * returns the Particles from the next neighbours of the current cell
 * @param cell
 * @return
 */
std::vector<int> LinkedCellContainer::get_Particles_from_next_cells(int cell){
    std::vector<int> vec;
    bool left = cell%x_cells < x_cells-1;
    bool up = cell < (x_cells*(y_cells-1));
    if(cell == 0){
        if(left) vec.push_back(cell + 1);
        if(up) {
            vec.push_back(cell + x_cells);
            if(left) vec.push_back(cell + 1 + x_cells);
        }
    }
    else if(cell < x_cells -1){
        if(left) vec.push_back(cell + 1);
        if(up) {
            if(left) vec.push_back(cell + 1 + x_cells);
        }
    }
    else if (cell % x_cells == 0){
        if(up) {
            vec.push_back(cell + 1);
            if(left)vec.push_back(cell + 1 + x_cells);
        }
    }
    else {
        if(left && up) vec.push_back(cell + 1 + x_cells);
    }
    return vec;
}

void LinkedCellContainer::setZero() {
    for(int i = 0; i < cells.size(); i++){
        for(int j = 0; j < cells[i].size(); j++){
            cells[i][j].setOldF(cells[i][j].getF());
            cells[i][j].setF({0,0,0});
        }
    }
}

int LinkedCellContainer::getXMax() {return x_cells;}

int LinkedCellContainer::getYMax() {return y_cells;}

int LinkedCellContainer::getZMax() {return z_cells;}