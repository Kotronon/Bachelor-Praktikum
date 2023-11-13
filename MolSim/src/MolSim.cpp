
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"
#include "ParticleContainer.h"
#include "spdlog/spdlog.h"

#include <iostream>
#include <string>

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 *  calculate the position for all particles
 */
void calculateV();

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

    // FileReader fileReader;
   //  FileReader::readFile(container, argsv[1]);


    //getting end time and delta t command options if specified

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

    if(cmdOptionExists(argsv, argsv+argc, "-level"))
    {
        std::string level = getCmdOption(argsv, argsv + argc, "-level");
        if(level == "info") spdlog::set_level(spdlog::level::info);
        if(level == "debug") spdlog::set_level(spdlog::level::debug);
        if(level == "criticalr") spdlog::set_level(spdlog::level::critical);
        if(level == "err") spdlog::set_level(spdlog::level::err);
        if(level == "n_levels") spdlog::set_level(spdlog::level::n_levels);
        if(level == "off") spdlog::set_level(spdlog::level::off);
        if(level == "trace") spdlog::set_level(spdlog::level::trace);
        if(level == "warn") spdlog::set_level(spdlog::level::warn);
    }
    spdlog::info("log_level: ", spdlog::get_level());

    double current_time = start_time;

    int iteration = 0;

    // for this loop, we assume: current x, current f and current v are known
    while (current_time < 0) {
        // calculate new x
        calculateX();
        // calculate new f
        calculateF();
        // calculate new v
        calculateV();

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

void calculateF() {
    std::array<double, 3> force{};
    for (auto &p1: container) {
        force = {0., 0., 0.};
        for (auto &p2: container) {
            if (p1 == p2) {}
            else {
                //Fi = SUM Fij
                //Fij = (MiMj * (xj-xi))/(||xi-xj||2^3)
                double mass = p1.getM() * p2.getM();
                std::array<double, 3> diffvec{};
                std::array<double, 3> x1 = p1.getX();
                std::array<double, 3> x2 = p2.getX();
                //calculating xj - xi
                diffvec[0] = x2[0] - x1[0];
                diffvec[1] = x2[1] - x1[1];
                diffvec[2] = x2[2] - x1[2];
                //calculating ||xi - xj||2
                double before_sqrt = diffvec[0] * diffvec[0] + diffvec[1] * diffvec[1] + diffvec[2] * diffvec[2];
                double length_dif = sqrt(before_sqrt);
                double length_pow = pow(length_dif, 3);
                //calculating old force + (MiMj * (xj-xi))/(||xi-xj||2^3)
                force[0] += diffvec[0] * (mass / length_pow);
                force[1] += diffvec[1] * (mass / length_pow);
                force[2] += diffvec[2] * (mass / length_pow);
            }
        }
        p1.setOldF(p1.getF());
        p1.setF(force);

    }
}

void calculateX() {
    for (auto &p: container) {
        //xi(tn+1) = xi(tn) + ∆t * vi(tn) + (∆t)^2 * Fi(tn) /2mi
        std::array<double, 3> x_old = p.getX();
        std::array<double, 3> v = p.getV();
        std::array<double, 3> f = p.getF();
        //calculating (∆t)^2 /2mi
        double t_mul_m = delta_t * delta_t / (2 * p.getM());
        std::array<double, 3> x_new{};
        //calculating xi(tn) + ∆t * vi(tn) + (∆t)^2 * Fi(tn) /2mi
        x_new[0] = x_old[0] + delta_t * v[0] + t_mul_m * f[0];
        x_new[1] = x_old[1] + delta_t * v[1] + t_mul_m * f[1];
        x_new[2] = x_old[2] + delta_t * v[2] + t_mul_m * f[2];
        p.setX(x_new);
    }
}

void calculateV() {
    for (auto &p: container) {
        //vi (tn+1) = vi(tn) + ∆t * Fi(tn) + Fi(tn+1) / 2mi
        std::array<double, 3> v_old = p.getV();
        std::array<double, 3> f = p.getF();
        std::array<double, 3> f_old = p.getOldF();
        std::array<double, 3> v_new{};
        //calculating vi(tn) + ∆t * Fi(tn) + Fi(tn+1) / 2mi
        v_new[0] = v_old[0] + (delta_t / (2 * p.getM())) * (f_old[0] + f[0]);
        v_new[1] = v_old[1] + (delta_t / (2 * p.getM())) * (f_old[1] + f[1]);
        v_new[2] = v_old[2] + (delta_t / (2 * p.getM())) * (f_old[2] + f[2]);
        p.setV(v_new);

    }
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
