#include "FileReader.h"
#include "outputWriter/FileWriter.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "Thermostat.h"
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

//Hardcoded values for now:
constexpr double start_time = 0;
double end_time = 100;
double delta_t = 0.0005;

int dim = 3;
double Grav = -12.44;

//if you want to use directSum please use DBL_MAX for each direction
std::array<double, 3> domain_size = {60, 60, 60};
//if you want to use directSum please use DBL_MAX
double cutoff = 3.0 * 1.2;

//boundary order:  left, right, up, down, behind, before
//boundary types: "o"(outflow), "r"(reflective), "p"(periodic)
//(if you want use directSum please use {"o", "o", "o", "o", "o", "o"})
std::array<std::basic_string<char>, 6> boundary = {"p", "p", "p", "p", "p", "p"};

//input file (file will be used if valid path is given and file is not empty)
std::string inputFile = "";

//checkpoints
bool checkpointing = false;
int num_checkpoints = 1;
//path to folder to be used for output of checkpoint files
std::string outputDirectory = "../output";

double initTemperature = 40;
int nThermostat = 1000;
bool applyBrownianMotion = true;

//optional:
bool targetTemperatureExists = false;
double targetTemperature = 0;

//optional:
bool differenceTemperatureExists = false;
double differenceTemperature = 10;

//Cuboids/Disks have to be created manually in main

//Creation of particle container to be filled with all relevant particles
ParticleContainer container = ParticleContainer();

int main(int argc, char *argsv[]) {

    //Creation of linked-cell container to be filled with all relevant particles
    LinkedCellContainer cells = LinkedCellContainer(domain_size, cutoff, boundary);

    //Add Particles from input file
    if (!inputFile.empty()) {
        FileReader::readFile(container, inputFile.data());
        cells.addContainer(container);
    }

    int checkpoint = 1;
    if (num_checkpoints < 0) checkpointing = false;
    int steps_between_checkpoints = int(end_time / delta_t) / num_checkpoints;

    //Creation of cuboids/disks for simulation with linked-cell container
    //Use either ParticleGenerator::createCuboidInCells or ParticleGenerator::createDiskInCells

    ParticleGenerator::createCuboidInCells({0.6, 0.6, 0.6}, {0, 0, 0}, {50, 20, 50}, 1.2, 1.0, cells, 1.2, 1, 1);
    ParticleGenerator::createCuboidInCells({0.6, 24.0, 0.6}, {0, 0, 0}, {50, 20, 50}, 1.2, 2.0, cells, 1.1, 1, 2);

    double current_time = start_time;
    int iteration = 0;

    //Pre-calculation of f
    ForceCalculator::LennardJonesForceCell(cells, Grav);

    //Initialization with Brownian Motion (without temperature and with average velocity)
    //VelocityCalculator::BrownianMotionInitializationCell(cells, 1.1, dim);

    //Initialization with Brownian Motion / temperature
    if (applyBrownianMotion) {
        Thermostat::initializeTemperatureWithBrownianMotion(initTemperature, dim, cells);
    } else {
        Thermostat::initializeTemperature(initTemperature, dim, cells);
    }

    if (!targetTemperatureExists) {
        targetTemperature = initTemperature;
    }


    //For this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {

        if (iteration != 0 && iteration % nThermostat == 0) {
            if (differenceTemperatureExists) {
                Thermostat::setTemperatureGradually(targetTemperature, differenceTemperature, dim, cells);
            } else {
                Thermostat::setTemperatureDirectly(targetTemperature, dim, cells);
            }
            spdlog::info("Set temperature to " + std::to_string(Thermostat::calculateCurrentTemperature(2, cells)) +
                         " Kelvin.");
        }

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
        if (iteration % 10 == 0) {
            spdlog::info("Iteration " + std::to_string(iteration) + " finished.");
        }
        if (iteration % steps_between_checkpoints == 0 && checkpointing) {
            std::string filename = outputDirectory + "/checkpoint" + std::to_string(checkpoint) + ".txt";
            checkpoint++;
            ParticleContainer currentState = cells.toContainer();
            FileWriter::writeFile(currentState, filename);
        }
        current_time += delta_t;
    }

    /* if(checkpointing){
      std::string filename = "../input/checkpointNew.txt";
              ParticleContainer currentState = cells.toContainer();
              FileWriter::writeFile(currentState, filename);
              }*/

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
                num_of_particles += (int) z.size();
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
