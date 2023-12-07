
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


// here we go
#include "our_input/newinput-pimpl.cxx"


/**
 * get over all the inputs from pimpl
 * */
void collectInput();

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
double avg_v = 0.1;

int dim = 2;

double eps = 5;

double sig = 1;


//Creation of particle container to be filled with all relevant particles
ParticleContainer container = ParticleContainer();

int main(int argc, char *argsv[]) {


    /*ParticleContainer cuboid_1 = ParticleGenerator::createCuboid(x_1,v_1,N_1,h,m);
    ParticleContainer cuboid_2 = ParticleGenerator::createCuboid(x_2,v_2,N_2,h,m);
    container.addParticleContainer(cuboid_1);
    container.addParticleContainer(cuboid_2);
*/
    LinkedCellContainer cells = LinkedCellContainer({120, 50, 1}, 3.0, {"r", "r", "r", "r", "r",
                                                                        "r"}); //boundary left, right, up, down, behind, bevor
    //ParticleGenerator::createCuboidInCells({20, 20, 0}, {0,0,0}, {100,20,1}, 1.1225, 1, cells, 3.0);
    //ParticleGenerator::createCuboidInCells({70, 60, 0}, {0,-10,0}, {20,20,1}, 1.1225, 1, cells, 3.0);
    //ParticleGenerator::createCuboidInCells({20, 20, 0}, {0,-10,0}, {20,10,1}, 1.1225, 1, cells, 3.0);

    ParticleGenerator::createDiskInCells({60, 25, 0}, {0, -10, 0}, 1, 15, 1.225, cells);

    double end_time = 10;
    double delta_t = 0.005;
    double current_time = start_time;
    int iteration = 0;

    //DEBUG
    plotParticlesInCells(iteration, cells);
    //DEBUG

    //Pre-calculation of f
    //ForceCalculator::LennardJonesForceFaster(container, eps, sig);

    ForceCalculator::LennardJonesForceCell(cells, eps, sig);

    //Initialization with Brownian Motion
    //VelocityCalculator::BrownianMotionInitialization(container, avg_v, dim);

    VelocityCalculator::BrownianMotionInitializationCell(cells, avg_v, dim);

    //For this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
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

    spdlog::info("Output written. Terminating...");
    return 0;
}






 void collectInput(){
    //LinkedCellContainer cells = LinkedCellContainer({180, 90, 1}, 3.0, {b1_, b2_, b3_, b4_, b5_, b6_}); //boundary left, right, up, down, behind, bevor
    //ParticleGenerator::createCuboidInCells({20, 20, 0}, {0,0,0}, {100,20,1}, 1.1225, 1, cells, 3.0);
    // ParticleGenerator::createCuboidInCells({70, 60, 0}, {0,-1,0}, {20,20,1}, 1.1225, 1, cells, 3.0);
    /*ParticleContainer cuboid_1 = ParticleGenerator::createCuboid(x_1,v_1,N_1,h,m);
    ParticleContainer cuboid_2 = ParticleGenerator::createCuboid(x_2,v_2,N_2,h,m);
    container.addParticleContainer(cuboid_1);
    container.addParticleContainer(cuboid_2);*/
/*constexpr double start_time = 0;
double avg_v = 0.1;
int dim = 2;
double eps = 5;
double sig = 1;*/

    //TODO: hier kommt noch die implementation rein




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
