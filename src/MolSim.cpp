
#include "FileReader.h"
#include "outputWriter/FileWriter.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "calculations/ForceCalculator.h"
#include "calculations/VelocityCalculator.h"
#include "spdlog/spdlog.h"
#include "calculations/PositionCalculator.h"
#include <string>


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
double end_time = 25;
double delta_t = 0.0005;

double avg_v = 0.1;
int dim = 2;

double eps = 5;
double sig = 1;
double Grav = -12.44;
//if you wanna use directSum please use DBL_MAX for each direction
std::array<double, 3> domain_size = {180, 90, 1};
//if you wanna use directSum please use DBL_MAX
double cutoff = 2.5;
//boundary order (b):  left, right, up, down, behind, before
//if you wanna use directSum please use {"o", "o", "o", "o", "o", "o"}
std::array<std::basic_string<char>, 6> boundary = {"p", "p", "r", "r", "o", "o"};
//input file
std::string inputFile;
//checkpoints
bool checkpointing = true;
int num_checkpoints = 1;

//Cuboids/Disks have to be created manually in main
//Creation of particle container to be filled with all relevant particles
ParticleContainer container = ParticleContainer();

int main(int argc, char *argsv[]) {
    //Creation of linked-cell container to be filled with all relevant particles
    LinkedCellContainer cells = LinkedCellContainer(domain_size, cutoff, boundary);
    //add Particles from input file
    if(!inputFile.empty()){
        FileReader::readFile(container, inputFile.data());
        cells.addContainer(container);
    }
    int steps_between_checkpoints = 0;
    int checkpoint = 0;
    if(num_checkpoints > 1){
        steps_between_checkpoints = (end_time/delta_t) / num_checkpoints;
    }
    //Creation of cuboids/disks for simulation with linked-cell container
    //Use either ParticleGenerator::createCuboidInCells or ParticleGenerator::createDiskInCells

    ParticleGenerator::createCuboidInCells({20, 20, 0}, {0,0,0}, {100,20,1}, 1.2, 1, cells, cutoff, 1, 1);
    ParticleGenerator::createCuboidInCells({70, 60, 0}, {0,-10,0}, {20,20,1}, 1.2, 2, cells, cutoff, 1, 0.9412);

    //ParticleGenerator::createDiskInCells({60, 25, 0}, {0, -10, 0}, 1, 15, 1.225, cells, sig, eps);
    double current_time = start_time;
    int iteration = 0;

    //Pre-calculation of f
    ForceCalculator::LennardJonesForceCell(cells, Grav);

    //Initialization with Brownian Motion
    VelocityCalculator::BrownianMotionInitializationCell(cells, avg_v, dim);

    //For this loop, we assume: current x, current f and current v are known
    while (current_time < 3000) {
        //Calculate new x
        PositionCalculator::PositionStoermerVerletCell(cells, delta_t);
        //Calculate new f
        ForceCalculator::LennardJonesForceCell(cells, Grav);

        //Calculate new v
        VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);

        iteration++;
        if (iteration % 10 == 0) {
            plotParticlesInCells(iteration, cells);
        }
        if (iteration % 100 == 0) {
            spdlog::info("Iteration " + std::to_string(iteration) + " finished.");
        }
        if(iteration % steps_between_checkpoints == 0 && checkpointing){
            std::string filename = "checkpoint" + std::to_string(checkpoint) + ".txt";
            checkpoint ++;
            ParticleContainer currentState = cells.toContainer();
            FileWriter::writeFile(currentState, filename);
        }
        current_time += delta_t;
    }

   if(checkpointing){
    std::string filename = "checkpoint" + std::to_string(checkpoint) + ".txt";
            checkpoint ++;
            ParticleContainer currentState = cells.toContainer();
            FileWriter::writeFile(currentState, filename);
            }
    spdlog::info("Output written. Terminating...");
    return 0;
}

void plotParticlesInCells(int iteration, LinkedCellContainer &cells) {

    std::string out_name("MD_vtk");

    //outputWriter::XYZWriter writer;
    //outputWriter::XYZWriter::plotParticles(container, out_name, iteration);

    outputWriter::VTKWriter writer2;
    int num_of_particles = 0;
    for (auto &x: cells) {
        for (auto &y: x) {
            for (auto &z: y) {
                num_of_particles += z.size();
            }
        }
    }
    writer2.initializeOutput(num_of_particles);
    for (auto &x: cells) {
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
