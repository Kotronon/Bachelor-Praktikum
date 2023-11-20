//
// Created by maraconda on 12.11.23.
//

#include "PositionCalculator.h"
#include "../ParticleContainer.h"
#include "../utils/ArrayUtils.h"
#include <spdlog/spdlog.h>
#include <cmath>

/**
 * Calculation of the new position of all molecules in the given ParticleContainer
 * @param container
 * @param delta_t
 */
void PositionCalculator::PositionStoermerVerlet(ParticleContainer &container, double delta_t) {
    for (auto &p: container) {
        p.setX(p.getX() + (delta_t * p.getV()) + ((delta_t*delta_t)/(2*p.getM())) * p.getF());
    }
}

void PositionCalculator::PositionStoermerVerletCell(LinkedCellContainer &grid, double delta_t, double cutoff) {
    for (int i = 0; i < grid.cell_numbers(); i++) {
        for(int j = 0; j < grid.Particles_in_cell(i); j++) {
           std::array<double, 3> x_new = grid.cells[i][j].getX() + (delta_t * grid.cells[i][j].getV()) +
                   (((delta_t * delta_t) / (2 * grid.cells[i][j].getM())) * grid.cells[i][j].getF());
           grid.cells[i][j].setX(x_new);
           int new_cell = i;
           int x_old = i% grid.getXMax();
           int y_old = i / grid.getXMax();
           int z_old = i /(grid.getXMax()*grid.getYMax());
           int x_now =  floor(x_new[0] / cutoff);
           int y_now = floor(x_new[1] /cutoff);
           int z_now = floor(x_new[2] /cutoff);
           if(x_now > x_old){
               if(x_now < grid.getXMax()) {
                   if(y_now > y_old) new_cell+=grid.getXMax();
                   else if(y_now < y_old) new_cell-=grid.getXMax();
                   if(z_now > z_old) new_cell += grid.getXMax()*grid.getYMax();
                   else if(z_now < z_old) new_cell -= grid.getXMax()*grid.getYMax();
                   grid.addParticle(new_cell+1, grid.cells[i][j]);
               }
               grid.deleteParticle(i, grid.cells[i][j]);
           }
           else if(x_now < x_old){
               if(x_now > 0) {
                   if(y_now > y_old) new_cell+=grid.getXMax();
                   else if(y_now < y_old) new_cell-=grid.getXMax();
                   if(z_now > z_old) new_cell += grid.getXMax()*grid.getYMax();
                   else if(z_now < z_old) new_cell -= grid.getXMax()*grid.getYMax();
                   grid.addParticle(new_cell-1, grid.cells[i][j]);
               }
               grid.deleteParticle(i, grid.cells[i][j]);
           }
           else if(y_now > y_old){
               if(y_now < grid.getYMax()) {
                   if(z_now > z_old) new_cell += grid.getXMax()*grid.getYMax();
                   else if(z_now < z_old) new_cell -= grid.getXMax()*grid.getYMax();
                   grid.addParticle(new_cell+grid.getXMax(), grid.cells[i][j]);
               }
               grid.deleteParticle(i, grid.cells[i][j]);
            }
           else if(y_now < y_old){
               if(y_now > 0 ) {
                   if(z_now > z_old) new_cell += grid.getXMax()*grid.getYMax();
                   else if(z_now < z_old) new_cell -= grid.getXMax()*grid.getYMax();
                   grid.addParticle(new_cell-grid.getXMax(), grid.cells[i][j]);
               }
               grid.deleteParticle(i, grid.cells[i][j]);
            }
           else if(z_now > z_old){
               if(z_now < grid.getZMax()) {
                   grid.addParticle(i+grid.getXMax()*grid.getYMax(), grid.cells[i][j]);
               }
               grid.deleteParticle(i, grid.cells[i][j]);
            }
           else if(z_now < z_old){
               if(z_now > 0 ) {
                   grid.addParticle(i-grid.getXMax()*grid.getYMax(), grid.cells[i][j]);
               }
               grid.deleteParticle(i, grid.cells[i][j]);
            }
        }
    }
}