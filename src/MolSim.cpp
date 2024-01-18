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
#include "gnuplot-iostream.h"
#include <string>
#include <vector>
#include <numeric>
#include <matplot/matplot.h>



/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

void plotParticlesInCells(int iteration, LinkedCellContainer &cells);

//Hardcoded values for now:
constexpr double start_time = 0;
double end_time = 500;
double delta_t = 0.0005;

int dim = 3;
double Grav = -0.8;

//if you want to use directSum please use DBL_MAX for each direction
std::array<double, 3> domain_size = {30, 30, 12};
//if you want to use directSum please use DBL_MAX
double cutoff = 2.5 * 1.1;

//boundary order:  left, right, up, down, behind, before
//boundary types: "o"(outflow), "r"(reflective), "p"(periodic)
//(if you want use directSum please use {"o", "o", "o", "o", "o", "o"})
std::array<std::basic_string<char>, 6> boundary = {"o", "o", "p", "p", "p", "p"};

//input file (file will be used if valid path is given and file is not empty)
std::string inputFile = "";// "../input/checkpoint1.txt";

//checkpoints
bool checkpointing = false;
int num_checkpoints = 1;
//path to folder to be used for output of checkpoint files
std::string outputDirectory = "../input";

double initTemperature = 40;
int nThermostat = 10;
bool applyBrownianMotion = true;

//optional:
bool targetTemperatureExists = false;
double targetTemperature = 0.02;

//optional:
bool differenceTemperatureExists = false;
double differenceTemperature = 2.5 / 1000;

Thermostat thermostat;

std::string filename = "../input/density_velocity_profile.csv";

bool calcDV = true;
int number_of_bins = 50;
double length_bin = domain_size[0] / number_of_bins;
std::ofstream File(filename);

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

    ParticleGenerator::createCuboidInCells({1.0, 0.5, 0.5}, {0,0,0}, {2, 30, 12}, 1, 1, cells, 1.1, 2.0, 1, true);
    ParticleGenerator::createCuboidInCells({27.2, 0.5, 0.5}, {0,0,0}, {2, 30, 12}, 1, 1, cells, 1.1, 2.0, 1, true);
    ParticleGenerator::createCuboidInCells({3.2, 0.6, 0.6}, {0,0,0}, {20, 25, 10}, 1.2, 1, cells, 1.0, 1.0, 2, false);

    double current_time = start_time;
    int iteration = 0;

    for(int i = 0; i < 49; i++){
        File << i << ", ";
    }
    File << 50 << "\n";
    //Pre-calculation of f
    ForceCalculator::LennardJonesForceCell(cells, Grav);

    //Initialization with Brownian Motion (without temperature and with average velocity)
    //VelocityCalculator::BrownianMotionInitializationCell(cells, 1.1, dim);

    //Initialization with Brownian Motion / temperature
    if (applyBrownianMotion) {
        thermostat.initializeTemperatureWithBrownianMotion(initTemperature, dim, cells);
    } else {
        thermostat.initializeTemperature(initTemperature, dim, cells);
    }

    if (!targetTemperatureExists) {
        targetTemperature = initTemperature;
    }

    //plotParticlesInCells(0, cells);
    //For this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {

        if (iteration != 0 && iteration % nThermostat == 0) {
            if (differenceTemperatureExists) {
                initTemperature = thermostat.setTemperatureGradually(targetTemperature, differenceTemperature, dim, cells, initTemperature);
            } else {
                thermostat.setTemperatureDirectly(targetTemperature, dim, cells);
            }
            spdlog::info("temperature with kinetic energy: " + std::to_string(thermostat.calculateCurrentTemperature(3, cells)) +
                         " Kelvin.");
            spdlog::info("new temperature: " + std::to_string(initTemperature) + " Kelvin");
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
        if (iteration % 100 == 0) {
            spdlog::info("Iteration " + std::to_string(iteration) + " finished.");
        }
        if (iteration % steps_between_checkpoints == 0 && checkpointing) {
            std::string filename = outputDirectory + "/checkpoint" + std::to_string(checkpoint) + ".txt";
            checkpoint++;
            ParticleContainer currentState = cells.toContainer();
            FileWriter::writeFile(currentState, filename);
        }

        if(iteration % 10000 == 0 && calcDV){
            cells.calcDVProfile(File, length_bin, number_of_bins);
            spdlog::info("calculation of the density-velocity profile");
        }
        current_time += delta_t;
    }

    File.close();
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
