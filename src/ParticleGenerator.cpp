//
// Created by maraconda on 12.11.23.
//

#define _USE_MATH_DEFINES
#include "utils/ArrayUtils.h"
#include "ParticleGenerator.h"
#include "ParticleContainer.h"
#include "Membrane.h"
#include <cmath>
#include <iostream>
#include <spdlog/spdlog.h>


/**
 * creates a cuboid and stores it in a ParticleContainer
 *  @param x position of left front corner
 * @param v initial velocity
 * @param N dimension of the cuboid
 * @param h distance between two particles
 * @param m mass for each particle
 * @param sig sigma for each particle
 * @param eps epsilon for each particle
 * @param type type for each particle
 * @return
 */
ParticleContainer ParticleGenerator::createCuboid(std::array<double, 3> x, std::array<double, 3> v,
                                                  std::array<int, 3> N, double h, double m, double sig, double eps, int type) {
    ParticleContainer cuboid = ParticleContainer();

    if (N[0] == 0.0 || N[1] == 0.0 || N[2] == 0.0) {
        return cuboid;
    }

    std::array<double, 3> coordinate = x;

    for (int z_i = 0; z_i < N[2]; z_i++) {
        for (int y_i = 0; y_i < N[1]; y_i++) {
            for (int x_i = 0; x_i < N[0]; x_i++) {
                cuboid.addParticle(coordinate, v, m, type, sig, eps);
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
 * @param N dimension of the cuboid
 * @param h distance between two particles
 * @param m mass for each particle
 * @param cells LinkedCellContainer
 * @param sig sigma for each particle
 * @param eps epsilon for each particle
 * @param type type for each particle
 */
void ParticleGenerator::createCuboidInCells(std::array<double, 3> x, std::array<double, 3> v,
                                            std::array<int, 3> N, double h, double m,
                                            LinkedCellContainer &cells, double sig, double eps, int type){

    if (N[0] == 0.0 || N[1] == 0.0 || N[2] == 0.0) {
        return;
    }

    std::array<double, 3> coordinate = x;
    int x_axis = (int) floor(x[0] / cells.getCutoff()) + 1;
    int y_axis = (int) floor(x[1] / cells.getCutoff()) + 1;
    int z_axis = (int) floor(x[2] / cells.getCutoff()) + 1;
    //to move to next cell
    int x_axis_tmp = x_axis;
    int y_axis_tmp = y_axis;
    int z_axis_tmp = z_axis;

    for (int z_i = 0; z_i < N[2]; z_i++) {
        for (int y_i = 0; y_i < N[1]; y_i++) {
            for (int x_i = 0; x_i < N[0]; x_i++) {
                cells.addParticle(x_axis_tmp, y_axis_tmp, z_axis_tmp, coordinate, v, m, type, sig, eps);
                //spdlog::info("added x_cell {}", x_axis_tmp);
                coordinate[0] += h;
                if(x_axis_tmp <= (floor(coordinate[0] / cells.getCutoff()))) x_axis_tmp ++;
            }
            coordinate[0] = x[0];
            x_axis_tmp = x_axis;
            coordinate[1] += h;
            if(y_axis_tmp <= floor(coordinate[1] / cells.getCutoff())) y_axis_tmp ++;
        }
        coordinate[1] = x[1];
        y_axis_tmp = y_axis;
        coordinate[2] += h;
        if(z_axis_tmp <= ceil(coordinate[2] / cells.getCutoff())) z_axis_tmp ++;
    }
}

/**
 * creates a 2-dimensional sphere/disk and stores it in a ParticleContainer based on Bresenham
 * @param center position of center
 * @param v initial velocity
 * @param m mass of each particle
 * @param r number of molecules along the radius
 * @param h distance between molecules
 * @param sig sigma of each particle
 * @param eps epsilon of each particle
 * @param type type of each particle
 * @return particle container containing the newly created particles
 */
ParticleContainer ParticleGenerator::createDisk(std::array<double, 3> center, std::array<double, 3> v, double m, int r, double h, double sig, double eps, int type) {

    auto disk = new ParticleContainer();
    double radius = r * h;
    double d = 3 - (2 * radius);
    double x = 0;
    double y = radius;
    bool fill = false;
    bool end = false;

    //Add center
    disk->addParticle({center[0], center[1] , center[2]},v,m,type, sig, eps);

    if (r == 0) {
        return *disk;
    }

    //Add center line
    for (double x_minus = center[0] - h, x_plus = center[0] + h; x_plus < center[0] + radius - (h/4); x_minus -= h, x_plus += h) {
        disk->addParticle({x_minus, center[1], center[2]},v,m,type, sig, eps);
        disk->addParticle({x_plus, center[1], center[2]},v,m,type, sig, eps);
    }

    //Add particles for top, bottom, left and right
    disk->addParticle({center[0] + x, center[1] + y, center[2]},v,m,type, sig, eps);
    disk->addParticle({center[0] + x, center[1] - y, center[2]},v,m,type, sig, eps);
    disk->addParticle({center[0] + y, center[1] + x, center[2]},v,m,type, sig, eps);
    disk->addParticle({center[0] - y, center[1] - x, center[2]},v,m,type, sig, eps);

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

    while (x < y + (h/4.0)) {

        if (std::abs(x - y) < (h/4.0)) {
            fill = false;
            end = true;
        }

        if (!end) {
            //starting from top going left and right
            disk->addParticle({center[0] - x, center[1] + y, center[2]},v,m,type, sig, eps);
            for (double x_temp = center[0] - x + h; fill && x_temp < center[0] + x - (h/4.0); x_temp += h) {
                disk->addParticle({x_temp, center[1] + y, center[2]},v,m,type, sig, eps);
            }
            disk->addParticle({center[0] + x, center[1] + y, center[2]},v,m,type, sig, eps);

            //starting from bottom going left and right
            disk->addParticle({center[0] - x, center[1] - y, center[2]},v,m,type, sig, eps);
            for (double x_temp = center[0] - x + h; fill && x_temp < center[0] + x - (h/4.0); x_temp += h) {
                disk->addParticle({x_temp, center[1] - y, center[2]},v,m,type, sig, eps);
            }
            disk->addParticle({center[0] + x, center[1] - y, center[2]},v,m,type, sig, eps);
        }

        //starting from left and right going up
        disk->addParticle({center[0] - y, center[1] + x, center[2]},v,m,type, sig, eps);
        for (double x_temp = center[0] - y + h; x_temp < center[0] + y - (h/4.0); x_temp += h) {
            disk->addParticle({x_temp, center[1] + x, center[2]},v,m,type, sig, eps);
        }
        disk->addParticle({center[0] + y, center[1] + x, center[2]},v,m,type, sig, eps);

        //starting from left and right going down
        disk->addParticle({center[0] - y, center[1] - x, center[2]},v,m,type, sig, eps);
        for (double x_temp = center[0] - y + h; x_temp < center[0] + y - (h/4.0); x_temp += h) {
            disk->addParticle({x_temp, center[1] - x, center[2]},v,m,type, sig, eps);
        }
        disk->addParticle({center[0] + y, center[1] - x, center[2]},v,m,type, sig, eps);

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
 * creates a 2-dimensional sphere/disk and stores it in a ParticleContainer based on cuboid generation
 * @param center position of center
 * @param v initial velocity
 * @param m mass of each particle
 * @param r number of molecules along the radius
 * @param h distance between molecules
 * @param sig sigma of each particle
 * @param eps epsilon of each particle
 * @param type type of each particle
 * @return particle container containing the newly created particles
 */
ParticleContainer ParticleGenerator::createDiskAlternative(std::array<double, 3> center, std::array<double, 3> v, double m, int r, double h, double sig, double eps, int type) {
    ParticleContainer disk = ParticleContainer();
    double radius = r * h;
    //Create cuboid around center
    ParticleContainer cuboid = createCuboid({center[0] - radius, center[1] - radius, center[2]}, v, {2 * r + 1, 2 * r + 1, 1},h,m, sig, eps, type);
    //Only add particles that are in radius of center
    for (const Particle &p : cuboid) {
        if (ArrayUtils::L2Norm(p.getX() - center) <= radius + (h/4.0)) {
            disk.addParticle(p);
        }
    }
    return disk;
}



/**
 * creates a 2-dimensional sphere/disk and stores it in a LinkedCellContainer
 * @param x position of center
 * @param v initial velocity
 * @param m mass of each particle
 * @param r number of molecules along the radius
 * @param h distance between molecules
 * @param cells LinkedCellContainer
 * @param sig sigma of each particle
 * @param eps epsilon of each particle
 * @param type type of each particle
 */
void ParticleGenerator::createDiskInCells(std::array<double, 3> center, std::array<double, 3> v, double m,
                                                int r, double h, LinkedCellContainer &cells, double sig, double eps, int type) {

    ParticleContainer container = ParticleGenerator::createDiskAlternative(center,v,m,r,h, sig, eps, type);
    for (auto p = container.begin(); p < container.end(); p++) {
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
 * @param sig sigma of each particle
 * @param eps epsilon of each particle
 * @param type type of each particle
 * @return particle container containing the newly created particles
 */
ParticleContainer ParticleGenerator::createSphere(std::array<double, 3> center, std::array<double, 3> v, double m,
                                                  int r, double h, double sig, double eps, int type) {

    ParticleContainer sphere = ParticleContainer();
    double radius = r * h;
    //Create cuboid around center
    ParticleContainer cuboid = createCuboid({center[0] - radius, center[1] - radius, center[2] - radius}, v, {2 * r + 1, 2 * r + 1, 2 * r + 1},h,m, sig, eps, type);
    //Only add particles that are in radius of center
    for (const Particle &p : cuboid) {
        if (ArrayUtils::L2Norm(p.getX() - center) <= radius + (h/4.0)) {
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
 * @param cells LinkedCellContainer
 * @param sig sigma of each particle
 * @param eps epsilon of each particle
 * @param type type of each particle
 */
void ParticleGenerator::createSphereInCells(std::array<double, 3> center, std::array<double, 3> v, double m,
                                                  int r, double h, LinkedCellContainer cells, double sig, double eps, int type) {

    ParticleContainer container = ParticleGenerator::createSphere(center,v,m,r,h, sig, eps, type);
    for (auto p = container.begin(); p < container.end(); p++) {
        cells.addParticle(*p);
    }
}

/**
 * create the membrane (a ParticleContainer) and sets the neighbours of all the particles,
 * @param n dimension
 * @param x was auch immer das ist
 * @param v velocity
 * @param m
 * @returns a PArticleContainer with particles in it
 * */

ParticleContainer ParticleGenerator::createMembrane(std::array<int, 3> n, std::array<double, 3> x, std::array<double, 3> v, double m, double h,
                                       LinkedCellContainer &cells, double eps, double sig,double k, double r_0,int type) {

    ParticleContainer membrane = ParticleContainer();
    spdlog::info("created the membrane");


    if (n[0] == 0.0 || n[1] == 0.0 || n[2] == 0.0) {
        return membrane;
    }

    std::array<double, 3> coordinate = x;

    for (int z_i = 0; z_i < n[2]; z_i++) {
        for (int y_i = 0; y_i < n[1]; y_i++) {
            for (int x_i = 0; x_i < n[0]; x_i++) {

                Particle p = {coordinate,v,m,sig,eps,type};
                membrane.SetAllNeighbours(p);
                spdlog::info("The neighbours are all set ");
                membrane.addParticle(p);

                coordinate[0] += h;
            }
            coordinate[0] = x[0];
            coordinate[1] += h;
        }
        coordinate[1] = x[1];
        coordinate[2] += h;
    }


    /*for (auto p = membrane.begin(); p < membrane.end(); p++) {
        cells.addParticle(*p);
    }*/

    cells.addContainer(membrane);
    spdlog::info("the container added all the particles");



    return membrane;

}




