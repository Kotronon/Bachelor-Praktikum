
#include "FileReader.h"

#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "calculations/ForceCalculator.h"
#include "calculations/VelocityCalculator.h"
#include "spdlog/spdlog.h"
#include "calculations/PositionCalculator.h"
#include "Thermostat.h"
#include "outputWriter/FileWriter.h"
#include <string>
#include <chrono>
#include <cfloat>
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
double end_time = 25;
double delta_t = 0.0005;

double avg_v = 0.1;
int dim = 2;

double eps = 5;
double sig = 1;
double Grav = -12.44;
//if you wanna use directSum please use DBL_MAX for each direction
std::array<double, 3> domain_size = {63,36,1};
//if you wanna use directSum please use DBL_MAX
double cutoff = 2.5 * 1;

//boundary order (b):  left, right, up, down, behind, before
//if you wanna use directSum please use {"o", "o", "o", "o", "o", "o"}
std::array<std::basic_string<char>, 6> boundary = {"p", "p", "r", "r", "o", "o"};
//input file
std::string inputFile = "";
//checkpoints
bool checkpointing = true;
int num_checkpoints = 2;

double initTemperature = 40;
int nThermostat = 1000;
bool applyBrownianMotion = true;

//optional:
bool targetTemperatureExists = true;
double targetTemperature;

//optional:
bool differenceTemperatureExists = true;
double differenceTemperature;


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
    //add Particles from input file
    if(!inputFile.empty()){
        FileReader::readFile(container, inputFile.data());
        cells.addContainer(container);
    }
    int checkpoint = 1;
    if(num_checkpoints < 0) checkpointing = false;
    int steps_between_checkpoints = int(end_time/delta_t) / num_checkpoints;

   //Creation of cuboids/disks for simulation with linked-cell container
    //Use either ParticleGenerator::createCuboidInCells or ParticleGenerator::createDiskInCells
    ParticleGenerator::createCuboidInCells({0.6, 2, 0}, {0,0,0}, {50, 14,1}, 1.2, 1.0, cells, cutoff, 1, 1, 0);
    ParticleGenerator::createCuboidInCells({0.6, 19, 0}, {0,0,0}, {50,14,1}, 1.2, 2, cells, cutoff, 0.9412, 1, 0);
    //ParticleGenerator::createDiskInCells({150, 150, 0}, {0, 0, 0}, 1, 20, 1.2, cells, 1.2, 1, 1);


    double current_time = start_time;
    int iteration = 0;

    //Pre-calculation of f
    ForceCalculator::LennardJonesForceCell(cells, Grav);

    //Initialization with Brownian Motion
    //VelocityCalculator::BrownianMotionInitializationCell(cells, avg_v, dim);
    //Initialization with Brownian Motion / temperature
    if (applyBrownianMotion) {
        Thermostat::initializeTemperatureWithBrownianMotion(initTemperature, dim, cells);
    }
       else {
        Thermostat::initializeTemperature(initTemperature, dim, cells);
    }

    if (!targetTemperatureExists) {
        targetTemperature = initTemperature;
    }

    auto end_time_setup = std::chrono::high_resolution_clock::now();
    auto time_setup = end_time_setup - start_time_setup;

    spdlog::info(&"The setup took this amount of nanoseconds:"[time_setup/std::chrono::nanoseconds(2)]);

    //For this loop, we assume: current x, current f and current v are known
    auto start_time_loop =  std::chrono::high_resolution_clock ::now();
    while (current_time < end_time) {
        auto start_time_loop_inside =  std::chrono::high_resolution_clock ::now();

        if (iteration != 0 && iteration % nThermostat == 0) {
            if (differenceTemperatureExists) {
                Thermostat::setTemperatureGradually(targetTemperature, differenceTemperature, dim, cells);
            }
            else {
                Thermostat::setTemperatureDirectly(targetTemperature, dim, cells);
            }
            spdlog::info("Set temperature to " + std::to_string(Thermostat::calculateCurrentTemperature(2,cells)) + " Kelvin.");
        }

        auto start_time_loop_calc =  std::chrono::high_resolution_clock ::now();
        //Calculate new x
        PositionCalculator::PositionStoermerVerletCell(cells, delta_t);
        //Calculate new f
        ForceCalculator::LennardJonesForceCell(cells, Grav);

        //Calculate new v
        VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);

        auto end_time_loop_calc =  std::chrono::high_resolution_clock ::now();
        auto calculations = end_time_loop_calc - start_time_loop_calc;

        spdlog::info(&"The current calculation for this iteration took this amount of milliseconds: "[calculations/std::chrono::milliseconds (2)]);


        iteration++;
        if (iteration % 10 == 0) {
            plotParticlesInCells(iteration, cells);
        }
        if (iteration % 100 == 0) {
            auto end_time_loop_inside =  std::chrono::high_resolution_clock ::now();
            spdlog::info("Iteration " + std::to_string(iteration) + " finished.");

            auto current_iteration = end_time_loop_inside - start_time_loop_inside;

            spdlog::info(&"The current iteration took this amount of milliseconds: "[current_iteration/std::chrono::milliseconds (2)]);

        }
        if(iteration % steps_between_checkpoints == 0 && checkpointing){
            std::string filename = "../output/checkpoint" + std::to_string(checkpoint) + ".txt";
            checkpoint ++;
            ParticleContainer currentState = cells.toContainer();
            FileWriter::writeFile(currentState, filename);
        }
        current_time += delta_t;

    }
    auto end_time_loop = std::chrono::high_resolution_clock::now();
    auto time_loop = end_time_loop - start_time_loop;


  /* if(checkpointing){
    std::string filename = "../input/checkpointNew.txt";
            ParticleContainer currentState = cells.toContainer();
            FileWriter::writeFile(currentState, filename);
            }*/
    spdlog::info("Output written. Terminating...");
    int64_t  t = time_loop/std::chrono::seconds(2);
    int size = container.size();
    spdlog::info(&"The loop took this amount of milliseconds:"[time_loop/std::chrono::milliseconds(2)]);
    spdlog::info(&"The average amount of iterations per second were:"[ (int64_t) iteration / t]);
    spdlog::info(&"The average amount of molecules updated per second were:"[ (int64_t) size / t]);
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
