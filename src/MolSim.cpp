
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "calculations/ForceCalculator.h"
#include "calculations/VelocityCalculator.h"
#include "spdlog/spdlog.h"
#include "calculations/PositionCalculator.h"
#include <string>
#include <chrono>
#include <iostream>


/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

void plotParticlesInCells(int iteration, LinkedCellContainer &cells);

/**
 * get a specific command option
 */
char *getCmdOption(char **begin, char **end, const std::string &option);

/**
 * check if a specific command option exists
 */
bool cmdOptionExists(char **begin, char **end, const std::string &option);

//Hardcoded values for now:
constexpr double start_time = 0;
double end_time = 10;
double delta_t = 0.00005;

double avg_v = 0.1;
int dim = 2;

double eps = 5;
double sig = 1;

std::array<double, 3> domain_size = {120,50,1};
double cutoff = 3.0;

//boundary order (b):  left, right, up, down, behind, before
std::array<std::basic_string<char>, 6> boundary = {"r", "r", "r", "r", "r", "r"};

//Cuboids/Disks have to be created manually in main

//Creation of particle container to be filled with all relevant particles
ParticleContainer container = ParticleContainer();

int main(int argc, char *argsv[]) {

    auto start_time_setup = std::chrono::high_resolution_clock ::now();
    //Creation of cuboids for simulation with simple particle container
    /*
    ParticleContainer cuboid_1 = ParticleGenerator::createCuboid(x_1,v_1,N_1,h,m);
    ParticleContainer cuboid_2 = ParticleGenerator::createCuboid(x_2,v_2,N_2,h,m);
    container.addParticleContainer(cuboid_1);
    container.addParticleContainer(cuboid_2);
    */

    //Creation of linked-cell container to be filled with all relevant particles
    LinkedCellContainer cells = LinkedCellContainer(domain_size, cutoff, boundary);

    //Creation of cuboids/disks for simulation with linked-cell container
    //Use either ParticleGenerator::createCuboidInCells or ParticleGenerator::createDiskInCells

    //ParticleGenerator::createCuboidInCells({20, 20, 0}, {0,0,0}, {100,20,1}, 1.1225, 1, cells, 3.0);
    //ParticleGenerator::createCuboidInCells({70, 60, 0}, {0,-10,0}, {20,20,1}, 1.1225, 1, cells, 3.0);
    ParticleGenerator::createDiskInCells({60, 25, 0}, {0, -10, 0}, 1, 15, 1.225, cells);

    double current_time = start_time;
    int iteration = 0;

    //Pre-calculation of f
    //ForceCalculator::LennardJonesForceFaster(container, eps, sig);

    ForceCalculator::LennardJonesForceCell(cells, eps, sig);

    //Initialization with Brownian Motion
    //VelocityCalculator::BrownianMotionInitialization(container, avg_v, dim);

    VelocityCalculator::BrownianMotionInitializationCell(cells, avg_v, dim);

    auto end_time_setup = std::chrono::high_resolution_clock ::now();
    auto setup_time = end_time_setup - start_time_setup;

    std::cout<<"The setup took "<< setup_time/std::chrono::nanoseconds(2)<<" nanoseconds to run"<<std::endl;

    auto start_time_loop = std::chrono::high_resolution_clock ::now();
    //For this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
        start_time_loop = std::chrono::high_resolution_clock ::now();

        //Calculate new x
        //PositionCalculator::PositionStoermerVerlet(container, delta_t);

        PositionCalculator::PositionStoermerVerletCell(cells, delta_t);

        //Calculate new f
        //ForceCalculator::LennardJonesForceFaster(container, eps, sig);

        ForceCalculator::LennardJonesForceCell(cells, eps, sig);

        //Calculate new v
        //VelocityCalculator::VelocityStoermerVerlet(container, delta_t);

        VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);

        iteration++;
        if (iteration % 10 == 0) {
            //plotParticles(iteration);
            plotParticlesInCells(iteration, cells);
        }
        if (iteration % 100 == 0) {
            spdlog::info("Iteration " + std::to_string(iteration) + " finished.");
        }

        current_time += delta_t;
    }
    auto end_time_loop = std::chrono::high_resolution_clock ::now();
    auto loop_time = end_time_loop - start_time_loop;
    std::cout<<"The loop took "<< loop_time/std::chrono::nanoseconds(2)<<" nanoseconds to run"<<std::endl;



    spdlog::info("Output written. Terminating...");
    return 0;
}

void plotParticlesInCells(int iteration, LinkedCellContainer &grid) {

    std::string out_name("MD_vtk");

    //outputWriter::XYZWriter writer;
    //outputWriter::XYZWriter::plotParticles(container, out_name, iteration);

    outputWriter::VTKWriter writer2;
    int num_of_particles = 0;
    for (auto &x: grid) {
        for (auto &y: x) {
            for (auto &z: y) {
                num_of_particles += z.size();
            }
        }
    }
    writer2.initializeOutput(num_of_particles);
    for (auto &x: grid) {
        for (auto &y: x) {
            for (auto &z: y) {
                for (auto &p: z) {
                    writer2.plotParticle(p);
                }
            }
        }
    }
    writer2.writeFile(out_name, iteration);
}

void plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    //outputWriter::XYZWriter writer;
    //outputWriter::XYZWriter::plotParticles(container, out_name, iteration);

    outputWriter::VTKWriter writer2;
    writer2.initializeOutput(container.size());
    for (auto &p: container) {
        writer2.plotParticle(p);
    }
    writer2.writeFile(out_name, iteration);
}

//Adapted from https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
char *getCmdOption(char **begin, char **end, const std::string &option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        return *itr;
    }
    return nullptr;
}

//Adapted from https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
bool cmdOptionExists(char **begin, char **end, const std::string &option) {
    return std::find(begin, end, option) != end;
}
