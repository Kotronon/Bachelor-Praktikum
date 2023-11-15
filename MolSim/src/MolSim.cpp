
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

//Hardcoded values for now:
constexpr double start_time = 0;
double avg_v = 0.1;
int dim = 3;
int eps = 5;
int sig = 1;
//Creation of particle container to be filled with all relevant particles
ParticleContainer container = ParticleContainer();

int main(int argc, char *argsv[]) {
    spdlog::info("Hello from MolSim for PSE!");
    if (argc <= 3 || ((argc - 3)%6 != 0 && (argc - 4)%6 != 0)) {
        spdlog::info("Erroneous programme call! ");
        spdlog::info("./MolSim end_time delta_time [options]");

        spdlog::info("Parameters: ");
        spdlog::info("end_time: The end time of the simulation");
        spdlog::info("delta_time: ∆time of the simulation");

        spdlog::info("Options: ");
        spdlog::info("-level : The log level. The default log level is info");
        spdlog::info("-f : The filepath and filename of the file to be used in the simulation in the format filepath/filename");
        spdlog::info("-c : Create and add a cuboid to the simulation. Has to be followed by 5 specific values describing the cuboid:");
        spdlog::info(" -> the number of cuboids you want to create");
        spdlog::info(" -> the coordinate of the lower left front-side corner of the cuboid with each values separated by a comma (x,y,z)");
        spdlog::info(" -> the initial velocity of all particles in the cuboid, each value separated by a comma (x,y,z)");
        spdlog::info(" -> the number of particles per dimension in the cuboid, each value separated by a comma (N_x,N_y,N_z)");
        spdlog::info(" -> the distance h of the particles");
        spdlog::info(" -> the mass of one particle in the cuboid");
        spdlog::info(" -> everything from the left_corner of the next cuboid");
    }

    //TODO: Add input for file and/or cuboid

    //Getting parameters end_time and delta_t
    double end_time = std::stod(argsv[1]);
    spdlog::info("end_time: ", end_time);

    double delta_t = std::stod(argsv[2]);
    spdlog::info("delta_time: ", delta_t);

    //Check for optional file input
    if(cmdOptionExists(argsv, argsv+argc, "-f"))
    {
        FileReader fileReader;
        FileReader::readFile(container, getCmdOption(argsv, argsv + argc, "-f"));
        spdlog::info("Read given file");
    }


    //Check for additional cuboids to be created
    if(cmdOptionExists(argsv + 3, argsv+argc, "-c"))
    {
        int i = 3;
        while(argsv[i] != "-c"){
            i++;
        }
        i++;
        int numbers_cuboid = std::stod(argsv[i]);
        std::list<std::string> cuboids_left_corner;
        std::list<std::string> cuboids_velocity;
        std::list<std::string> cuboids_dimension;
        std::list<std::string> cuboids_h;
        std::list<std::string> cuboids_mass;
        for(int j = 0; j < numbers_cuboid; j++){
            cuboids_left_corner.push_back(argsv[i]);
            i++;
            cuboids_velocity.push_back(argsv[i]);
            i++;
            cuboids_dimension.push_back(argsv[i]);
            i++;
            cuboids_h.push_back(argsv[i]);
            i++;
            cuboids_mass.push_back(argsv[i]);

        }
        //ParticleContainer cuboid = ParticleGenerator::createCuboid();
        for(int j = 0; j < numbers_cuboid; j++){
            std::string first_particle = cuboids_left_corner.front();
            cuboids_left_corner.pop_front();
            first_particle.erase(0,1);
            first_particle.pop_back();
            std::string s;
            std::array<double, 3> coordinates_left_corner;
            int k = 0;
            i = 0;
            while (first_particle[i] != '\0') {
                if (first_particle[i] != ',') {
                    // Append the char to the temp string.
                    s += first_particle[i];
                } else {
                    coordinates_left_corner[k] = std::stod(s);
                    k++;
                    s.clear();
                }
                i++;

            }
            std::string velocity_particle = cuboids_velocity.front();
            cuboids_velocity.pop_front();
            velocity_particle.erase(0,1);
            velocity_particle.pop_back();
            std::array<double, 3> coordinates_velocity;
            k = 0;
            i = 0;
            while (velocity_particle[i] != '\0') {
                if (velocity_particle[i] != ',') {
                    // Append the char to the temp string.
                    s += velocity_particle[i];
                } else {
                    coordinates_velocity[k] = std::stod(s);
                    k++;
                    s.clear();
                }
                i++;

            }
            std::string dimension = cuboids_dimension.front();
            cuboids_dimension.pop_front();
            dimension.erase(0,1);
            dimension.pop_back();
            std::array<int, 3> dimensions;
            k = 0;
            i = 0;
            while (dimension[i] != '\0') {
                if (dimension[i] != ',') {
                    // Append the char to the temp string.
                    s += dimension[i];
                } else {
                    dimensions[k] = std::stod(s);
                    k++;
                    s.clear();
                }
                i++;

            }
            double h = std::stod(cuboids_h.front());
            cuboids_h.pop_front();

            double m = std::stod(cuboids_mass.front());
            cuboids_mass.pop_front();

            ParticleContainer new_cuboid = ParticleGenerator::createCuboid(coordinates_left_corner, coordinates_velocity, dimensions, h, m);
            container.addParticleContainer(new_cuboid);
        }
    }



    double current_time = start_time;
    int iteration = 0;

    //Pre-calculation of f
    ForceCalculator::LennardJonesForce(container, eps, sig);
    //Initialization with Brownian Motion
    VelocityCalculator::BrownianMotionInitialization(container, avg_v, dim);
    //For this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
        // calculate new x
        PositionCalculator::PositionStoermerVerlet(container, delta_t);
        // calculate new f
        ForceCalculator::LennardJonesForce(container, eps, sig);
        // calculate new v
        VelocityCalculator::VelocityStoermerVerlet(container, delta_t);

        iteration++;
        if (iteration % 10 == 0) {
            plotParticles(iteration);
        }
        spdlog::info("Iteration " + std::to_string(iteration) + " finished.");

        current_time += delta_t;
    }


    spdlog::info("Output written. Terminating..." );
    return 0;
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
