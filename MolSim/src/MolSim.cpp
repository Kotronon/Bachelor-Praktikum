
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "ParticleContainer.h"
#include "ForceCalculator.h"
#include "VelocityCalculator.h"
#include "spdlog/spdlog.h"
#include "PositionCalculator.h"
#include <iostream>
#include <string>

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

/**
 * get a specific command option
 */
char* getCmdOption(char ** begin, char ** end, const std::string & option);

/**
 * check if a specific command option exists
 */
bool cmdOptionExists(char** begin, char** end, const std::string& option);

constexpr double start_time = 0;
double end_time = 1000;
double delta_t = 0.014;

ParticleContainer container = ParticleContainer();

int main(int argc, char *argsv[]) {
    spdlog::info("Hello from MolSim for PSE!");
    if (argc <= 2 || argc >= 8) {
        spdlog::info("Erroneous programme call! ");
        spdlog::info("./MolSim <filepath/filename> [options]");
        spdlog::info("Options: ");
        spdlog::info("-e : The end time of the simulation, default value is 1000");
        spdlog::info("-d : ∆time of the simulation, default value is 0.014");
    }

    FileReader fileReader;
    FileReader::readFile(container, argsv[1]);


    //Getting end time and delta t command options if specified

    if(cmdOptionExists(argsv, argsv+argc, "-e"))
    {
        end_time = std::stod(getCmdOption(argsv, argsv + argc, "-e"));
    }

    spdlog::info("end_time: ", end_time);

    if(cmdOptionExists(argsv, argsv+argc, "-d"))
    {
        delta_t = std::stod(getCmdOption(argsv, argsv + argc, "-d"));
    }

    spdlog::info("delta_time: ", delta_t);


    double current_time = start_time;

    int iteration = 0;

    //pre-calculation of f
    ForceCalculator::SimpleForceCalculation(container);

    // for this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
        // calculate new x
        PositionCalculator::PositionStoermerVerlet(container, delta_t);
        // calculate new f
        ForceCalculator::SimpleForceCalculation(container);
        // calculate new v
        VelocityCalculator::VelocityStoermerVerlet(container, delta_t);

        iteration++;
        if (iteration % 10 == 0) {
            plotParticles(iteration);
        }
        spdlog::info("Iteration ", iteration, " finished.");

        current_time += delta_t;
    }


    spdlog::info("output written. Terminating..." );
    return 0;
}

void plotParticles(int iteration) {

    std::string out_name("MD_vtk");

    outputWriter::XYZWriter writer;
    outputWriter::XYZWriter::plotParticles(container, out_name, iteration);

    outputWriter::VTKWriter writer2;
    writer2.initializeOutput(container.size());
    for (auto &p: container) {
        writer2.plotParticle(p);
    }
    writer2.writeFile(out_name, iteration);
}

//Adapted from https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return nullptr;
}

//Adapted from https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c
bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}
