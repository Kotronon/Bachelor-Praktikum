//
// Created by maraconda on 12.11.23.
//

#define _USE_MATH_DEFINES
#include "utils/ArrayUtils.h"
#include "ParticleGenerator.h"
#include "ParticleContainer.h"
#include <cmath>


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
 * @param x position of left front corner
 * @param v initial velocity
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
    int x_axis = floor(x[0] / cutoff);
    int y_axis = floor(x[1] / cutoff);
    int z_axis = floor(x[2] / cutoff);
    //to move to next cell
    int x_axis_tmp = x_axis;
    int y_axis_tmp = y_axis;
    int z_axis_tmp = z_axis;

    for (int z_i = 0; z_i < N[2]; z_i++) {
        for (int y_i = 0; y_i < N[1]; y_i++) {
            for (int x_i = 0; x_i < N[0]; x_i++) {
                cells.addParticle(x_axis_tmp, y_axis_tmp, z_axis_tmp, coordinate, v, m, 0);
                coordinate[0] += h;
                if(x_axis_tmp < floor(coordinate[0] / cutoff)) x_axis_tmp ++;
            }
            coordinate[0] = x[0];
            x_axis_tmp = x_axis;
            coordinate[1] += h;
            if(y_axis_tmp < floor(coordinate[1] / cutoff)) y_axis_tmp ++;
        }
        coordinate[1] = x[1];
        y_axis_tmp = y_axis;
        coordinate[2] += h;
        if(z_axis_tmp < ceil(coordinate[2] / cutoff)) z_axis_tmp ++;
    }
}

/**
 * creates a 2-dimensional sphere/disk and stores it in a ParticleContainer
 * @param center position of center
 * @param v initial velocity
 * @param m mass of each particle
 * @param r number of molecules along the radius
 * @param h distance between molecules
 * @return
 */
ParticleContainer ParticleGenerator::createDisk(std::array<double, 3> center, std::array<double, 3> v, double m,
                                                  int r, double h) {
    //based on Bresenham algorithm for circles
    ParticleContainer disk;
    double radius = r*h;
    double d = (5 - (radius * 4)) / 4;
    double x = 0;
    double y = radius;

    do {
        disk.addParticle({center[0] + x, center[1] + y, center[2]},v,m,0);
        disk.addParticle({center[0] + x, center[1] - y, center[2]},v,m,0);
        disk.addParticle({center[0] - x, center[1] + y, center[2]},v,m,0);
        disk.addParticle({center[0] - x, center[1] - y, center[2]},v,m,0);

        disk.addParticle({center[0] + y, center[1] + x, center[2]},v,m,0);
        disk.addParticle({center[0] + y, center[1] - x, center[2]},v,m,0);
        disk.addParticle({center[0] - y, center[1] + x, center[2]},v,m,0);
        disk.addParticle({center[0] - y, center[1] - x, center[2]},v,m,0);

        if (d < 0) {
            d += 2 * x + 1;
        }

        else {
            d += 2 * (x-y) + 1;
            y -= h;
        }

        x += h;

    } while (x <= y);
    return disk;
}

/**
 * creates a 2-dimensional sphere/disk and stores it in a LinkedCellContainer
 * @param x position of center
 * @param v initial velocity
 * @param m mass of each particle
 * @param r number of molecules along the radius
 * @param h distance between molecules
 * @return
 */
void ParticleGenerator::createDiskInCells(std::array<double, 3> x, std::array<double, 3> v, double m,
                                                int r, double h, LinkedCellContainer &cells) {

    ParticleContainer container = ParticleGenerator::createDisk(x,v,m,r,h);
    for (auto p = container.begin(); p <= container.end(); p++) {
        cells.addParticle(*p);
    }
}

/**
 * creates a 3-dimensional sphere and stores it in a ParticleContainer
 * @param x position of center
 * @param v initial velocity
 * @param m mass of each particle
 * @param r number of molecules along the radius
 * @param h distance between molecules
 * @return
 */
ParticleContainer ParticleGenerator::createSphere(std::array<double, 3> x, std::array<double, 3> v, double m,
                                                  int r, double h) {

    ParticleContainer sphere = ParticleContainer();
    //Add center
    //sphere.addParticle(x, v, m, 0);
    //Iterate over each layer pair
    //for (double z_dif; z_dif < r; z_dif += h) {

    //}
    //Create circle for each layer

    return sphere;
}
