//
// Created by kathi on 20.11.23.
//


#include "LinkedCellContainer.h"
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
    x_cells = ceil(N[0] / cutoff);
    y_cells = ceil(N[1] / cutoff);
    z_cells = ceil(N[2] / cutoff);
    cells = std::vector<std::vector<Particle>>(x_cells* y_cells + x_cells*y_cells*z_cells);
    c = cutoff;
}

LinkedCellContainer::~LinkedCellContainer(){}

/**
 * returns the number of cells in the grid
 * @return
 */
int LinkedCellContainer::cell_numbers() const {
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
    cells[cell].emplace_back((x_arg, v_arg, m_arg, type_arg));
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
 * deletes Particle from specific cell
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
void LinkedCellContainer::moveToNeighbour(){
    for (int i = 0; i < cell_numbers(); i++) {
        for (int j = 0; j < Particles_in_cell(i); j++) {
            int new_cell = i;
            int x_old = i % getXMax();
            int y_old = i / getXMax()*getZMax();
            int z_old = i /(getXMax()*getYMax());
            int x_now =  floor(cells[i][j].getX()[0] / c);
            int y_now = floor(cells[i][j].getX()[1] /c);
            int z_now = floor(cells[i][j].getX()[2] /c);
            if(x_now > x_old){
                if(x_now < getXMax()-1) {
                    if(y_now > y_old && y_now < getYMax()) new_cell+=getXMax();
                    else if(y_now < y_old && y_now >= 0) new_cell-=getXMax();
                    if(z_now > z_old && z_now < getZMax()) new_cell += getXMax()*getYMax();
                    else if(z_now < z_old && z_now >= 0) new_cell -= getXMax()*getYMax();
                    addParticle(new_cell+1, cells[i][j]);
                }
                deleteParticle(i, cells[i][j]);
                j--;
            }
            else if(x_now < x_old){
                if(x_now > 0) {
                    if(y_now > y_old && y_now < getYMax()) new_cell+=getXMax();
                    else if(y_now < y_old && y_now >= 0) new_cell-=getXMax();
                    if(z_now > z_old && z_now < getZMax()) new_cell += getXMax()*getYMax();
                    else if(z_now < z_old && z_now  >= 0) new_cell -= getXMax()*getYMax();
                    addParticle(new_cell-1, cells[i][j]);
                }
                deleteParticle(i, cells[i][j]);
                j--;
            }
            else if(y_now > y_old){
                if(y_now < getYMax()-1) {
                    if(z_now > z_old && z_now < getZMax()) new_cell += getXMax()*getYMax();
                    else if(z_now < z_old && z_now  >= 0) new_cell -= getXMax()*getYMax();
                    addParticle(new_cell+getXMax(), cells[i][j]);
                }
                deleteParticle(i, cells[i][j]);
                j--;
            }
            else if(y_now < y_old){
                if(y_now > 0 ) {
                    if(z_now > z_old && z_now < getZMax()) new_cell += getXMax()*getYMax();
                    else if(z_now < z_old && z_now  >= 0) new_cell -= getXMax()*getYMax();
                    addParticle(new_cell-getXMax(), cells[i][j]);
                }
                deleteParticle(i, cells[i][j]);
                j--;
            }
            else if(z_now > z_old){
                if(z_now < getZMax()-1) {
                    addParticle(i+getXMax()*getYMax(), cells[i][j]);
                }
                deleteParticle(i, cells[i][j]);
                j--;
            }
            else if(z_now < z_old){
                if(z_now > 0 ) {
                    addParticle(i-getXMax()*getYMax(), cells[i][j]);
                }
                deleteParticle(i, cells[i][j]);
                j--;
            }
        }
    }
}

/**
 * returns the Particles from the next neighbours of the current cell
 * @param cell
 * @return
 */
std::vector<int> LinkedCellContainer::get_Particles_from_next_cells(int cell) const{
    std::vector<int> vec;
    bool left = cell%x_cells < x_cells-1;
    bool up = cell / y_cells < y_cells-1;
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
    for(auto & cell : cells){
        for(auto & p : cell){
            p.setOldF(p.getF());
            p.setF({0,0,0});
        }
    }
}

int LinkedCellContainer::getXMax() const {return x_cells;}

int LinkedCellContainer::getYMax() const {return y_cells;}

int LinkedCellContainer::getZMax() const {return z_cells;}