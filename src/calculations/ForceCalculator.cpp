//
// Created by maraconda on 12.11.23.
//

#include "ForceCalculator.h"
#include "../ParticleContainer.h"
#include "../utils/ArrayUtils.h"

double ForceCalculator::epsilon = 5;
double ForceCalculator::sigma = 1;

/**
 * Calculates the gravity force of all Particles in given ParticleContainer
 * @param container
 */
void ForceCalculator::SimpleForceCalculation(ParticleContainer &container) {
    std::array<double, 3> force{};
    for (auto &p1: container) {
        force = {0., 0., 0.};
        for (auto &p2: container) {
            if (!(p1 == p2)) {
                force = force + (((p1.getM() * p2.getM()) / pow((ArrayUtils::L2Norm(p1.getX() - p2.getX())), 3)) * (p2.getX() - p1.getX()));
            }
        }
        p1.setOldF(p1.getF());
        p1.setF(force);
    }
}

/**
 * Calculates the Lennord Jones force of all Particles in given ParticleContainer
 * @param container
 * @param eps
 * @param sig
 */
void ForceCalculator::LennardJonesForce(ParticleContainer &container, double eps, double sig) {
    ForceCalculator::epsilon = eps;
    ForceCalculator::sigma = sig;

    std::array<double, 3> force{};
    double L2Norm_p1_p2;
    for (auto &p1: container) {
        force = {0., 0., 0.};
        for (auto &p2: container) {
            if (!(p1 == p2)) {
                L2Norm_p1_p2 = ArrayUtils::L2Norm(p1.getX() - p2.getX());
                force = force + ((-24*eps / pow(L2Norm_p1_p2,2)) * (pow(sig/L2Norm_p1_p2,6) - (2 * pow(sig/L2Norm_p1_p2,12))) * (p1.getX() - p2.getX()));
            }
        }
        p1.setOldF(p1.getF());
        p1.setF(force);
    }
}

/**
 * Faster calculation of the the Lennord Jones force of all Particles in given ParticleContainer
 * @param container
 * @param eps
 * @param sig
 */
void ForceCalculator::LennardJonesForceFaster(ParticleContainer &container, double eps, double sig) {
    ForceCalculator::epsilon = eps;
    ForceCalculator::sigma = sig;
    std::array<double, 3> zero = {0,0,0};
    for (auto &p: container) {
        p.setOldF(p.getF());
        p.setF(zero);
    }
    container.applyForcePairwise(ForceCalculator::LennardJonesForcePairwise);
}

/**
 * pairwise calculation of the Lennord Jones force of all Particles in given ParticleContainer
 * @param p1
 * @param p2
 */
void ForceCalculator::LennardJonesForcePairwise(Particle *p1, Particle *p2) {
    std::array<double, 3> force = {0,0,0};
    double L2Norm_p1_p2 = ArrayUtils::L2Norm(p1->getX() - p2->getX());
    force = force + ((-24*epsilon / pow(L2Norm_p1_p2,2)) * (pow(sigma/L2Norm_p1_p2,6) - (2 * pow(sigma/L2Norm_p1_p2,12))) * (p1->getX() - p2->getX()));
    p1->setF(p1->getF() + force);
    p2->setF(p2->getF() - force);
}

void ForceCalculator::LennardJonesForceCell(LinkedCellContainer &grid, double eps, double sig){
    ForceCalculator::epsilon = eps;
    ForceCalculator::sigma = sig;
    grid.setZero();
    for(int i = 0; i < grid.cell_numbers(); i++){
        std::vector<int> neighbours = grid.get_Particles_from_next_cells(i);
        for(int j = 0; j < grid.Particles_in_cell(i); j++){
            //for all particles in current cell
            for(int k = j+1; k < grid.Particles_in_cell(i); k++){
                //calculate force with particles in current cell
                //LennardJonesForcePairwise(&(grid.cells[i][j]), &(grid.cells[i][k]));
                std::array<double, 3> force = {0,0,0};
                double L2Norm_p1_p2 = ArrayUtils::L2Norm(grid.cells[i][j].getX() - grid.cells[i][k].getX());
                force = force + ((-24*epsilon / pow(L2Norm_p1_p2,2)) * (pow(sigma/L2Norm_p1_p2,6) - (2 * pow(sigma/L2Norm_p1_p2,12))) * (grid.cells[i][j].getX() - grid.cells[i][k].getX()));
                grid.cells[i][j].setF(grid.cells[i][j].getF() + force);
                grid.cells[i][k].setF(grid.cells[i][k].getF() - force);
            }
            for(int neighbour : neighbours){
                //with neighbour cells
                for(int l = 0; l < grid.Particles_in_cell(neighbour); l++){
                   // LennardJonesForcePairwise(&(grid.cells[i][j]), &(grid.cells[neighbour][l]));
                    std::array<double, 3> force = {0,0,0};
                    double L2Norm_p1_p2 = ArrayUtils::L2Norm(grid.cells[i][j].getX() - grid.cells[neighbour][l].getX());
                    force = force + ((-24*epsilon / pow(L2Norm_p1_p2,2)) * (pow(sigma/L2Norm_p1_p2,6) - (2 * pow(sigma/L2Norm_p1_p2,12))) * (grid.cells[i][j].getX() - grid.cells[neighbour][l].getX()));
                    grid.cells[i][j].setF(grid.cells[i][j].getF() + force);
                    grid.cells[neighbour][l].setF(grid.cells[neighbour][l].getF() - force);
                }
            }
        }
    }
}
