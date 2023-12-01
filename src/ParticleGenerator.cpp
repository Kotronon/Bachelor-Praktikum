//
// Created by maraconda on 12.11.23.
//

#define _USE_MATH_DEFINES
#include "utils/ArrayUtils.h"
#include "ParticleGenerator.h"
#include "ParticleContainer.h"
#include <cmath>
#include <iostream>


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
 * @return particle container containing the newly created particles
 */
ParticleContainer ParticleGenerator::createDisk(std::array<double, 3> center, std::array<double, 3> v, double m, int r, double h) {

    auto disk = new ParticleContainer();
    double radius = r * h;
    double d = 3 - (2 * radius);
    double x = 0;
    double y = radius;
    bool fill = false;

    //Add center
    disk->addParticle({center[0], center[1] , center[2]},v,m,0);

    if (r == 0) {
        return *disk;
    }

    //Add center line
    for (double x_minus = center[0] - h, x_plus = center[0] + h; x_plus < center[0] + radius - (h/4); x_minus -= h, x_plus += h) {
        disk->addParticle({x_minus, center[1], center[2]},v,m,0);
        disk->addParticle({x_plus, center[1], center[2]},v,m,0);
    }

    //Add particles for top, bottom, left and right
    disk->addParticle({center[0] + x, center[1] + y, center[2]},v,m,0);
    disk->addParticle({center[0] + x, center[1] - y, center[2]},v,m,0);
    disk->addParticle({center[0] + y, center[1] + x, center[2]},v,m,0);
    disk->addParticle({center[0] - y, center[1] - x, center[2]},v,m,0);

    if (d < 0) {
        d = d + (4 * x) + 6;
        fill = false;
    }
    else {
        d = d + 4 * (x - y) + 10;
        y -= h;
        fill = true;
    }

    x += h;

    while (x < y + (h/4)) {

        //starting from top going left and right
        disk->addParticle({center[0] - x, center[1] + y, center[2]},v,m,0);
        for (double x_temp = center[0] - x + h; fill && x_temp < center[0] + x - (h/4); x_temp += h) {
            disk->addParticle({x_temp, center[1] + y, center[2]},v,m,0);
        }
        disk->addParticle({center[0] + x, center[1] + y, center[2]},v,m,0);

        //starting from bottom going left and right
        disk->addParticle({center[0] - x, center[1] - y, center[2]},v,m,0);
        for (double x_temp = center[0] - x + h; fill && x_temp < center[0] + x - (h/4); x_temp += h) {
            disk->addParticle({x_temp, center[1] - y, center[2]},v,m,0);
        }
        disk->addParticle({center[0] + x, center[1] - y, center[2]},v,m,0);

        //starting from left and right going up
        disk->addParticle({center[0] - y, center[1] + x, center[2]},v,m,0);
        for (double x_temp = center[0] - y + h; x_temp < center[0] + y - (h/4); x_temp += h) {
            disk->addParticle({x_temp, center[1] + x, center[2]},v,m,0);
        }
        disk->addParticle({center[0] + y, center[1] + x, center[2]},v,m,0);

        //starting from left and right going down
        disk->addParticle({center[0] - y, center[1] - x, center[2]},v,m,0);
        for (double x_temp = center[0] - y + h; x_temp < center[0] + y - (h/4); x_temp += h) {
            disk->addParticle({x_temp, center[1] - x, center[2]},v,m,0);
        }
        disk->addParticle({center[0] + y, center[1] - x, center[2]},v,m,0);

        if (d < 0) {
            d = d + (4 * x) + 6;
            fill = false;
        }
        else {
            d = d + 4 * (x - y) + 10;
            y -= h;
            fill = true;
        }

        x += h;

    }

    return *disk;
}



/**
 * creates a 2-dimensional sphere/disk and stores it in a LinkedCellContainer
 * @param x position of center
 * @param v initial velocity
 * @param m mass of each particle
 * @param r number of molecules along the radius
 * @param h distance between molecules
 */
void ParticleGenerator::createDiskInCells(std::array<double, 3> center, std::array<double, 3> v, double m,
                                                int r, double h, LinkedCellContainer &cells) {

    ParticleContainer container = ParticleGenerator::createDisk(center,v,m,r,h);
    for (auto p = container.begin(); p < container.end(); p++) {
        cells.addParticle(*p);
    }
}

/**
 * creates a 3-dimensional sphere and stores it in a ParticleContainer (not functional yet)
 * @param x position of center
 * @param v initial velocity
 * @param m mass of each particle
 * @param r number of molecules along the radius
 * @param h distance between molecules
 * @return particle container containing the newly created particles
 */
ParticleContainer ParticleGenerator::createSphere(std::array<double, 3> center, std::array<double, 3> v, double m,
                                                  int r, double h) {

    ParticleContainer sphere = ParticleContainer();
    double radius = r * h;
    //Create cuboid around center
    ParticleContainer cuboid = createCuboid({center[0] - radius, center[1] - radius, center[2] - radius}, v, {2 * r + 1, 2 * r + 1, 2 * r + 1},h,m);
    //Only add particles that are in radius of center
    for (const Particle &p : cuboid) {
        if (ArrayUtils::L2Norm(p.getX() - center) <= radius) {
            sphere.addParticle(p);
        }
    }
    return sphere;
}

/**
 * creates a 3-dimensional sphere and stores it in a LinkedCellContainer (not functional yet)
 * @param x position of center
 * @param v initial velocity
 * @param m mass of each particle
 * @param r number of molecules along the radius
 * @param h distance between molecules
 */
void ParticleGenerator::createSphereInCells(std::array<double, 3> center, std::array<double, 3> v, double m,
                                                  int r, double h, LinkedCellContainer cells) {

    ParticleContainer container = ParticleGenerator::createSphere(center,v,m,r,h);
    for (auto p = container.begin(); p < container.end(); p++) {
        cells.addParticle(*p);
    }
}
