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
#include <omp.h>
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
double end_time = 150;
double delta_t = 0.001;

int dim = 3;
double Grav = 0;

//if you want to use directSum please use DBL_MAX for each direction
std::array<double, 3> domain_size = {9.2, 9.2, 9.2};
//if you want to use directSum please use DBL_MAX
double cutoff = 2.3;
//if you want to use smoothed Lennard-Jones Potential set smoothLJ to true
double sLJRadius = 1.9;
bool smoothLJ = true;

//boundary order:  left, right, up, down, behind, before
//boundary types: "o"(outflow), "r"(reflective), "p"(periodic)
//(if you want use directSum please use {"o", "o", "o", "o", "o", "o"})
std::array<std::basic_string<char>, 6> boundary = {"p", "p", "p", "p", "p", "p"};

//input file (file will be used if valid path is given and file is not empty)
std::string inputFile = "";// "../input/checkpoint1.txt";

//checkpoints
bool checkpointing = true;
int num_checkpoints = 1;
//path to folder to be used for output of checkpoint files
std::string outputDirectory = "../input";

double initTemperature = 0.01;
int nThermostat = 40;
bool applyBrownianMotion = true;

//optional:
bool targetTemperatureExists = true;
double targetTemperature = 3.0;

//optional:
bool differenceTemperatureExists = true;
double differenceTemperature = 0.001;

Thermostat thermostat;

int intervalBegin = 0;
int intervalEnd = 20;
double deltaR = 1;
std::string filename = "../input/RDF.xsl";

//Cuboids/Disks have to be created manually in main

//Creation of particle container to be filled with all relevant particles
ParticleContainer container = ParticleContainer();

int main(int argc, char *argsv[]) {

    //Creation of linked-cell container to be filled with all relevant particles
    LinkedCellContainer cells = LinkedCellContainer(domain_size, cutoff, boundary, smoothLJ, sLJRadius);

    //Add Particles from input file
    if (!inputFile.empty()) {
        FileReader::readFile(container, inputFile.data());
        cells.addContainer(container);
    }
    int checkpoint = 1;
    if (num_checkpoints < 0) checkpointing = false;
    int steps_between_checkpoints = int(end_time / delta_t) / num_checkpoints;

    std::ofstream RDFFile (filename);

    //Creation of cuboids/disks for simulation with linked-cell container
    //Use either ParticleGenerator::createCuboidInCells or ParticleGenerator::createDiskInCells

    ParticleGenerator::createCuboidInCells({0.575, 0.575, 0.575}, {0, 0, 0}, {8,8,8}, 1.15, 1.0, cells, 1, 1, 1);
    //ParticleGenerator::createCuboidInCells({0.6, 24.6, 0.6}, {0, 0, 0}, {50, 20, 50}, 1.2, 2.0, cells, 1.1, 1, 2);

    double current_time = start_time;
    int iteration = 0;

    //std::ofstream diffusion_file;
    //diffusion_file.open("../input/diffusion.xls");
    //std::ofstream RDF_file("../input/RDF.xls");
    std::vector<int> x_axis_plot;
    std::iota(std::begin(x_axis_plot), std::end(x_axis_plot), intervalBegin);

    /*RDF_file << " ";
    for(int i = intervalBegin; i <= intervalEnd; i++){
        RDF_file << i;
    }
    RDF_file << '\t';*/

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

    //plotParticlesInCells(0, cells);

    //For this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {

        if (iteration != 0 && iteration % nThermostat == 0) {
            if (differenceTemperatureExists) {
                initTemperature = thermostat.setTemperatureGradually(targetTemperature, differenceTemperature, dim, cells, initTemperature);
            } else {
                Thermostat::setTemperatureDirectly(targetTemperature, dim, cells);
            }
            //spdlog::info("Temperature with kinetic energy: " + std::to_string(Thermostat::calculateCurrentTemperature(3, cells)));
            //spdlog::info("New temperature: " + std::to_string(initTemperature));
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

        if(iteration % 1000 == 0){
           // diffusion_file << std::to_string(iteration);
            double diffusion = cells.calculateDiffusion();
            //diffusion_file << std::to_string(diffusion) << '\t';
            //RDF_file << iteration;
            std::vector< double> densities = cells.calculateRDF(intervalBegin, intervalEnd, deltaR, x_axis_plot, RDFFile);
            /*for(int i = 0; i < densities.size(); i++){
                RDF_file << densities[i];
            }
            RDF_file << '\t';*/
            spdlog::info("Diffusion: {}", diffusion);
        }


        current_time += delta_t;
    }

    /* if(checkpointing){
      std::string filename = "../input/checkpointNew.txt";
              ParticleContainer currentState = cells.toContainer();
              FileWriter::writeFile(currentState, filename);
              }*/
    RDFFile.close();
    matplot::title("RDF");
    matplot::save("../input/plot1.pdf");
    //matplot::set_ylabel("Density");
    const std::vector<double> leg ({0.5, -0.5});
    matplot::legend();
    matplot::save("../input/plot1_legend.pdf");
    matplot::show();
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
