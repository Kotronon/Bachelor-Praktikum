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
    cells = std::vector<std::vector<Particle>>(x_cells*y_cells*z_cells);
    c = cutoff;
    boundary = b;
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
    cells[cell].push_back(Particle(x_arg, v_arg, m_arg, type_arg));
}

void LinkedCellContainer::addParticle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg) {
    cells[floor(x_arg[0]/c)*ceil(x_arg[1]/c)*ceil(x_arg[2]/c)].push_back(Particle(x_arg, v_arg, m_arg, type_arg));
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
 * delete Particle from specific cell
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
 * checks if particle needs to be moved to another cell
 */
void LinkedCellContainer::moveToNeighbour() {
    for (int i = 0; i < cell_numbers(); i++) {
        for (int j = 0; j < Particles_in_cell(i); j++) {
            int x_now = ceil(cells[i][j].getX()[0] / c);
            int y_now = ceil(cells[i][j].getX()[1] / c);
            int z_now = ceil(cells[i][j].getX()[2] / c);
            int new_cell = x_now * y_now * z_now - 1;
            if (new_cell != i) {
                if(new_cell < cell_numbers() && new_cell >= 0){
                    addParticle(new_cell, cells[i][j]);
                    cells[i].erase(cells[i].begin() + j);
                }
                else{
                    bool deleted = false;
                    if(x_now > x_cells){
                        if(boundary[0] == "o") {
                            cells[i].erase(cells[i].begin() + j);
                            deleted = true;
                        }

                        else if(boundary[0] == "r"){
                            cells[i][j].setX({c*x_cells-1, cells[i][j].getX()[1], cells[i][j].getX()[2]});
                            cells[i][j].setV({-cells[i][j].getV()[0], cells[i][j].getV()[1], cells[i][j].getV()[2]});
                        }
                    }
                    else if(x_now < 0){
                        if(boundary[1] == "o"){
                            cells[i].erase(cells[i].begin() + j);
                            deleted = true;
                        }
                        else if(boundary[1] == "r"){
                            cells[i][j].setX({0, cells[i][j].getX()[1], cells[i][j].getX()[2]});
                            cells[i][j].setV({-cells[i][j].getV()[0], cells[i][j].getV()[1], cells[i][j].getV()[2]});
                        }
                    }
                    if(y_now > y_cells){
                        if(boundary[2] == "o" && !deleted){
                            cells[i].erase(cells[i].begin() + j);
                            deleted = true;
                        }
                        else if(boundary[2] == "r" && !deleted){
                            cells[i][j].setX({cells[i][j].getX()[0], c*y_cells-1, cells[i][j].getX()[2]});
                            cells[i][j].setV({cells[i][j].getV()[0], -cells[i][j].getV()[1], cells[i][j].getV()[2]});
                        }
                    }
                    else if(y_now < 0){
                        if (boundary[3] == "o" && !deleted){
                            cells[i].erase(cells[i].begin() + j);
                            deleted = true;
                        }
                        else if (boundary[3] == "r") {
                            cells[i][j].setX({cells[i][j].getX()[0], 0, cells[i][j].getX()[2]});
                            cells[i][j].setV({cells[i][j].getV()[0], -cells[i][j].getV()[1], cells[i][j].getV()[2]});
                        }
                    }
                    if(z_now > z_cells){
                        if(boundary[4] == "o" && !deleted)
                            cells[i].erase(cells[i].begin() + j);
                        else if (boundary[4] == "r" && !deleted) {
                            cells[i][j].setX({cells[i][j].getX()[0], cells[i][j].getX()[1], c*z_cells-1});
                            cells[i][j].setV({cells[i][j].getV()[0], cells[i][j].getV()[1], -cells[i][j].getV()[2]});
                        }
                    }
                    else if(z_now < 0){
                        if(boundary[5] == "o" && !deleted){
                            cells[i].erase(cells[i].begin() + j);
                            deleted = true;
                        }
                        else if (boundary[5] == "r" && !deleted) {
                            cells[i][j].setX({cells[i][j].getX()[0], cells[i][j].getX()[1], 0});
                            cells[i][j].setV({cells[i][j].getV()[0], cells[i][j].getV()[1], -cells[i][j].getV()[2]});
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
std::vector<int> LinkedCellContainer::get_next_cells(int cell) const{
    std::vector<int> vec;
    bool left = cell%x_cells < x_cells-1;
    bool up = cell / x_cells < y_cells-1;
    bool right = cell%x_cells > 0;
    bool down = cell / x_cells > 0;
    bool behind = cell / (x_cells*y_cells) < z_cells-1;
    bool before = cell / (x_cells * y_cells) > 0;
    //2D
    if(left) vec.push_back(cell+1);
    if(up) vec.push_back(cell+x_cells);
    if(left && up) vec.push_back(cell+1+x_cells);
    if(right && up) vec.push_back(cell-1+x_cells);
    //3D
    if(z_cells>1){
    if(before) vec.push_back(cell+x_cells*y_cells);
    if(right && before) vec.push_back(cell-1+x_cells*y_cells);
    if(up && before) vec.push_back(cell+x_cells+x_cells*y_cells);
    if(down && before) vec.push_back(cell-x_cells+x_cells*y_cells);
    if(left && up && before) vec.push_back(cell+1+x_cells+x_cells*y_cells);
    if(right && up && before) vec.push_back(cell-1+x_cells+x_cells*y_cells);
    if(left && down && before) vec.push_back(cell+1-x_cells+x_cells*y_cells);
    if(right && down && before) vec.push_back(cell-1-x_cells+x_cells*y_cells);
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
