//
// Created by maraconda on 12.11.23.
//

#include "ParticleGenerator.h"
#include "ParticleContainer.h"


/**
 * creates a cuboid and stores it in a ParticleContainer
 * @param x
 * @param v
 * @param N
 * @param h
 * @param m
 * @return
 */
ParticleContainer ParticleGenerator::createCuboid(std::array<double, 3> x, std::array<double, 3> v,
                                                  std::array<int, 3> N, double h, double m) {
    ParticleContainer cuboid = ParticleContainer();

    if (N[0] == 0.0 || N[1] == 0.0 || N[2] == 0.0) {
        return cuboid;
    }

    std::array<double, 3> coordinate = x;

    for (int z_i = 0; z_i < N[2]; z_i++) {
        for (int y_i = 0; y_i < N[1]; y_i++) {
            for (int x_i = 0; x_i < N[0]; x_i++) {
                cuboid.addParticle(coordinate, v, m, 0);
                coordinate[0] += h;
            }
            coordinate[0] = x[0];
            coordinate[1] += h;
        }
        coordinate[1] = x[1];
        coordinate[2] += h;
    }

    return cuboid;
}



/**
 * create a cuboid and stores it in the given cell grid
 * @param x
 * @param v
 * @param N
 * @param h
 * @param m
 * @param cells
 * @param cutoff
 */
void ParticleGenerator::createCuboidInCells(std::array<double, 3> x, std::array<double, 3> v,
                                            std::array<int, 3> N, double h, double m,
                                            LinkedCellContainer &cells, double  cutoff){

    if (N[0] == 0.0 || N[1] == 0.0 || N[2] == 0.0) {
        return;
    }

    std::array<double, 3> coordinate = x;
    int x_achsis = N[0] / cutoff;
    int y_achsis = N[1] / cutoff;
    int z_achsis = N[2] / cutoff;
    //to move to next cell
    int x_achsis_tmp = x_achsis;
    int y_achsis_tmp = y_achsis;
    int z_achsis_tmp = z_achsis;

    for (int z_i = 0; z_i < N[2]; z_i++) {
        for (int y_i = 0; y_i < N[1]; y_i++) {
            for (int x_i = 0; x_i < N[0]; x_i++) {
                cells.addParticle((x_achsis_tmp*y_achsis_tmp*z_achsis_tmp), coordinate, v, m, 0);
                coordinate[0] += h;
                if(x_achsis_tmp < coordinate[0] / cutoff) x_achsis_tmp ++;
            }
            coordinate[0] = x[0];
            x_achsis_tmp = x_achsis;
            coordinate[1] += h;
            if(y_achsis_tmp < coordinate[1] / cutoff) y_achsis_tmp ++;
        }
        coordinate[1] = x[1];
        y_achsis_tmp = y_achsis;
        coordinate[2] += h;
        if(z_achsis_tmp < coordinate[2] / cutoff) z_achsis_tmp ++;
    }
}
