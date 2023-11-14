
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
#include "spdlog/sinks/stdout_color_sinks.h"


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
    auto console = spdlog::stdout_color_mt("consol_logger");
    console->info("Hello from MolSim for PSE!");
    if (argc <= 2 || argc >= 8) {
        console->info("Erroneous programme call! ");
        console->info("./MolSim <filepath/filename> [options]");
        console->info("Options: ");
        console->info("-e : The end time of the simulation, default value is 1000");
        console->info("-d : ∆time of the simulation, default value is 0.014");
        console->info("-level : ∆the log leve, default is info");
    }

    // FileReader fileReader;
   //  FileReader::readFile(container, argsv[1]);

    //Getting end time and delta t command options if specified

    if(cmdOptionExists(argsv, argsv+argc, "-e"))
    {
        end_time = std::stod(getCmdOption(argsv, argsv + argc, "-e"));
    }

    console->info("end_time: {}" , end_time);

    if(cmdOptionExists(argsv, argsv+argc, "-d"))
    {
        delta_t = std::stod(getCmdOption(argsv, argsv + argc, "-d"));
    }
    if(cmdOptionExists(argsv, argsv+argc, "-d"))
    {
        delta_t = std::stod(getCmdOption(argsv, argsv + argc, "-d"));
    }
    console->info("delta_time: {}", delta_t);
    {
        std::string level = getCmdOption(argsv, argsv + argc, "-level");
        if(level == "info") console->set_level(spdlog::level::info);
        if(level == "debug") console->set_level(spdlog::level::debug);
        if(level == "criticalr") console->set_level(spdlog::level::critical);
        if(level == "err") console->set_level(spdlog::level::err);
        if(level == "n_levels") console->set_level(spdlog::level::n_levels);
        if(level == "off") console->set_level(spdlog::level::off);
        if(level == "trace") console->set_level(spdlog::level::trace);
        if(level == "warn") console->set_level(spdlog::level::warn);
    }
    console->info("log_level: {}", console->level());
    console->dump_backtrace();
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
        console->info("Iteration {} finished", iteration);

        current_time += delta_t;
    }


    console->info("output written. Terminating..." );

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
