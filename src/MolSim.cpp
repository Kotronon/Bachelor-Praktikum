#include "outputWriter/FileWriter.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "FileReader.h"
#include "spdlog/spdlog.h"
#include "ParticleContainer.h"
#include "ParticleGenerator.h"
#include "calculations/ForceCalculator.h"
#include "calculations/VelocityCalculator.h"
#include "calculations/PositionCalculator.h"
#include "Thermostat.h"
#include "gnuplot-iostream.h"
#include <matplot/matplot.h>
#include <string>
#include <vector>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

void plotParticlesInCells(int iteration, LinkedCellContainer &cells);

//Simulation parameters:


//Time:
//----------------------------------------------------------------------------------------------------------------------
constexpr double start_time = 0;
double end_time = 100;
double delta_t = 0.001;

//Number of threads
int num_threads = 8;
//----------------------------------------------------------------------------------------------------------------------


//Checkpointing:
//----------------------------------------------------------------------------------------------------------------------
bool checkpointing = false;
int num_checkpoints = 1;

//Path to folder to be used for output of checkpoint files
std::string outputDirectory = "../input";

//Input file (file will be used if valid path is given and file is not empty)
bool useInput = true;
std::string inputFile = "../input/checkpoint1.txt";
//----------------------------------------------------------------------------------------------------------------------


//Domain with boundaries:
//----------------------------------------------------------------------------------------------------------------------
//Dimensions (Choose between 2 and 3):
int dim = 3;

std::array<double, 3> domain_size = {9.2, 9.2, 9.2};
//If you want to use directSum please use DBL_MAX for each direction

//Boundary types: "o"(outflow), "r"(reflective), "p"(periodic)
//Boundary order:  left, right, up, down, behind, before
std::array<std::basic_string<char>, 6> boundary = {"p", "p", "p", "p", "p", "p"};
//If you want use directSum please use {"o", "o", "o", "o", "o", "o"}

double cutoff = 2.3;
//If you want to use directSum please use DBL_MAX
//----------------------------------------------------------------------------------------------------------------------


//Force Calculation:
//----------------------------------------------------------------------------------------------------------------------
//Addition of gravitational force, set to 0 if it should not be added
double Grav = 0;

//Set to true if you want to use the smoothed Lennard-Jones potential
bool smoothLJ = true;
double sLJRadius = 1.9;
//----------------------------------------------------------------------------------------------------------------------

//Thermostat:
//----------------------------------------------------------------------------------------------------------------------
//Initialization
bool applyBrownianMotion = false;
double initTemperature = 3.0;

//Time steps between applications
int nThermostat = 25;

//Optional target temperature
bool targetTemperatureExists = true;
double targetTemperature = 0.02;

//Optional temperature difference
bool differenceTemperatureExists = true;
double differenceTemperature = 2.5 * pow(10, -3);
//----------------------------------------------------------------------------------------------------------------------

//Diffusion:
//----------------------------------------------------------------------------------------------------------------------
bool calculateDiffusion = true;
int intervalBegin = 0;
int intervalEnd = 80;
double deltaR = 1;
std::string filename = "../input/RDF1.xsl";
//----------------------------------------------------------------------------------------------------------------------


//Cuboids/Disks have to be created manually in main!


//Creation of particle container to be filled with all relevant particles
ParticleContainer container = ParticleContainer();

int main(int argc, char *argsv[]) {

    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    //Creation of linked-cell container to be filled with all relevant particles
    LinkedCellContainer cells = LinkedCellContainer(domain_size, cutoff, boundary, smoothLJ, sLJRadius);

    //Add Particles from input file
    if (useInput && !inputFile.empty()) {
        FileReader::readFile(container, inputFile.data());
        cells.addContainer(container);
    }

    //Creation of cuboids/disks for simulation with linked-cell container
    //------------------------------------------------------------------------------------------------------------------
    //Use either ParticleGenerator::createCuboidInCells or ParticleGenerator::createDiskInCells

    ParticleGenerator::createCuboidInCells({0.575, 0.575, 0.575}, {0, 0, 0}, {8,8,8}, 1.15, 1.0, cells, 1, 1, 1);
    //ParticleGenerator::createCuboidInCells({0.6, 0.6, 0.6}, {0, 0, 0}, {50, 20, 50}, 1.2, 1.0, cells, 1.2, 1, 1);
    //ParticleGenerator::createCuboidInCells({0.6, 24.6, 0.6}, {0, 0, 0}, {50, 20, 50}, 1.2, 2.0, cells, 1.1, 1, 2);
    //------------------------------------------------------------------------------------------------------------------

    double current_time = start_time;
    int iteration = 0;

    int checkpoint = 1;
    if (num_checkpoints < 0) checkpointing = false;
    int steps_between_checkpoints = int(end_time / delta_t) / num_checkpoints;

    std::ofstream RDFFile (filename);
    std::vector<int> x_axis_plot;

    if (calculateDiffusion) {
        std::iota(std::begin(x_axis_plot), std::end(x_axis_plot), intervalBegin);
    }

    //Pre-calculation of f
    ForceCalculator::LennardJonesForceCell(cells, Grav);

    /*
    //Initialization with Brownian Motion / temperature
    if (applyBrownianMotion) {
        Thermostat::initializeTemperatureWithBrownianMotion(initTemperature, dim, cells);
    } else {
        Thermostat::initializeTemperature(initTemperature, dim, cells);
    }
    */

    if (!targetTemperatureExists) {
        targetTemperature = initTemperature;
    }

    while (current_time < end_time) {

        //Apply thermostat every nThermostat iterations
        if (iteration != 0 && iteration % nThermostat == 0) {
            if (differenceTemperatureExists) {
                initTemperature = Thermostat::setTemperatureGradually(targetTemperature, differenceTemperature, dim, cells, initTemperature);
            } else {
                Thermostat::setTemperatureDirectly(targetTemperature, dim, cells);
            }
            spdlog::info("Temperature with kinetic energy: " + std::to_string(Thermostat::calculateCurrentTemperature(3, cells)));
        }

        //Calculate new positions
        PositionCalculator::PositionStoermerVerletCell(cells, delta_t);
        //Calculate new forces
        ForceCalculator::LennardJonesForceCell(cells, Grav);
        //Calculate new velocities
        VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);

        iteration++;

        if (iteration % 10 == 0) {
            plotParticlesInCells(iteration, cells);
        }
        if (iteration % 100 == 0) {
            spdlog::info("Iteration " + std::to_string(iteration) + " finished.");
        }


        if (iteration % steps_between_checkpoints == 0 && checkpointing) {
            std::string file = outputDirectory + "/checkpoint" + std::to_string(checkpoint) + ".txt";
            checkpoint++;
            ParticleContainer currentState = cells.toContainer();
            FileWriter::writeFile(currentState, file);
        }

        if(iteration % 1000 == 0 && calculateDiffusion){
            double diffusion = cells.calculateDiffusion();
            std::vector< double> densities = cells.calculateRDF(intervalBegin, intervalEnd, deltaR, x_axis_plot, RDFFile);
            spdlog::info("Diffusion: {}", diffusion);
        }

        current_time += delta_t;
    }

    if (calculateDiffusion) {
        RDFFile.close();
        matplot::title("RDF");
        matplot::save("../input/plot1.pdf");
        const std::vector<double> leg ({0.5, -0.5});
        matplot::legend();
        matplot::save("../input/plot1_legend.pdf");
        matplot::show();
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
