
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

void plotParticlesInCells(int iteration, LinkedCellContainer &cells);

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
double eps = 5;
double sig = 1;
//Creation of particle container to be filled with all relevant particles
ParticleContainer container = ParticleContainer();

int main(int argc, char *argsv[]) {
    spdlog::info("Hello from MolSim for PSE!");
   /* if (argc <= 3) {
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
        spdlog::info(" -> the coordinate of the lower left front-side corner of the cuboid with each values separated by a comma x,y,z");
        spdlog::info(" -> the initial velocity of all particles in the cuboid, each value separated by a comma x,y,z");
        spdlog::info(" -> the number of particles per dimension in the cuboid, each value separated by a comma N_x,N_y,N_z");
        spdlog::info(" -> the distance h of the particles");
        spdlog::info(" -> the mass of one particle in the cuboid");
        spdlog::info(" -> everything from the left_corner of the next cuboid");
    }

    //TODO: Add input for file and/or cuboid

    //Getting parameters end_time and delta_t
   /* double end_time = std::stod(argsv[1]);
    spdlog::info("end_time: {}", end_time);

    double delta_t = std::stod(argsv[2]);
    spdlog::info("delta_time: {}", delta_t);

    //Check for optional file input
    if(cmdOptionExists(argsv, argsv+argc, "-f"))
    {
        FileReader fileReader;
        FileReader::readFile(container, getCmdOption(argsv, argsv + argc, "-f"));
        spdlog::info("Read given file");
    }

    if(cmdOptionExists(argsv, argsv+argc, "-level")){
        std::string level = getCmdOption(argsv, argsv + argc, "-level");
        if(level == "off") spdlog::set_level(spdlog::level::off);
        else if(level == "info") spdlog::set_level(spdlog::level::info);
        else if(level == "warn") spdlog::set_level(spdlog::level::warn);
        else if(level == "trace") spdlog::set_level(spdlog::level::trace);
        else if(level == "n_levels") spdlog::set_level(spdlog::level::n_levels);
        else if(level == "err") spdlog::set_level(spdlog::level::err);
        else if(level == "critical") spdlog::set_level(spdlog::level::critical);
        else if(level == "debug") spdlog::set_level(spdlog::level::debug);
        else spdlog::error("Level does not exist");
    }
    spdlog::info("log level: {}", spdlog::get_level());


    //Check for additional cuboids to be created
    if(cmdOptionExists(argsv, argsv+argc, "-c"))
    {
        int i = 3;
        std::string check = argsv[i];
        while(check != "-c"){
            i++;
            check = argsv[i];
        }
        i++;
        int numbers_cuboid = std::stoi(argsv[i]);

        i++;
        //iterating through each cuboid thet needs to be generated
        for(int j = 0; j < numbers_cuboid; j++){
            std::string first_particle = argsv[i];
            std::string s;
            //get x coordinates of the particle in the left corner of the cuboid
            std::array<double, 3> coordinates_left_corner = {0,0,0};
            int k = 0;
            int it = 0;
            while (first_particle[it] != '\0') {
                if (first_particle[it] != ',') {
                    // Append the char to the temp string.
                    s += first_particle[it];
                } else {
                    coordinates_left_corner[k] = std::stod(s);
                    k++;
                    s.clear();
                }
                it++;

            }
            coordinates_left_corner[k] = std::stod(s);
            i++;
            //get velocity of the cuboid
            std::string velocity_particle = argsv[i];
            std::array<double, 3> coordinates_velocity = {0,0,0};
            k = 0;
            it = 0;
            while (velocity_particle[it] != '\0') {
                if (velocity_particle[it] != ',') {
                    // Append the char to the temp string.
                    s += velocity_particle[it];
                } else {
                    coordinates_velocity[k] = std::stod(s);
                    k++;
                    s.clear();
                }
                it++;

            }
            coordinates_velocity[k] = std::stod(s);
            i++;
            //get the dimension of the cuboid
            std::string dimension = argsv[i];
            //dimension.erase(0,1);
            //dimension.pop_back();
            std::array<int, 3> dimensions = {0,0,0};
            k = 0;
            it = 0;
            while (dimension[it] != '\0') {
                if (dimension[it] != ',') {
                    // Append the char to the temp string.
                    s += dimension[it];
                } else {
                    dimensions[k] = std::stoi(s);
                    k++;
                    s.clear();
                }
                it++;

            }
            dimensions[k] = std::stoi(s);
            i++;
            //get the distance between the particles in the cuboid
            double h = std::stod(argsv[i]);
            i++;
            //get the mass of the particles in the cuboid
            double m = std::stod(argsv[i]);
            i++;
            //generate the cuboid and add it to the ParticleContainer
            ParticleContainer new_cuboid = ParticleGenerator::createCuboid(coordinates_left_corner, coordinates_velocity, dimensions, h, m);
            container.addParticleContainer(new_cuboid);
        }
    }

/*

    ParticleContainer cuboid_1 = ParticleGenerator::createCuboid(x_1,v_1,N_1,h,m);
    ParticleContainer cuboid_2 = ParticleGenerator::createCuboid(x_2,v_2,N_2,h,m);
    container.addParticleContainer(cuboid_1);
    container.addParticleContainer(cuboid_2);
*/
   LinkedCellContainer cells = LinkedCellContainer({180, 90, 1}, 3.0);
   ParticleGenerator::createCuboidInCells({20,20,0}, {0,0,0}, {100,20,1}, 1.1225, 1, cells, 3.0);
    ParticleGenerator::createCuboidInCells({70,60,0}, {0,-10,0}, {20,20,1}, 1.1225, 1, cells, 3.0);
   double end_time = 20;
   double delta_t = 0.0005;
    double current_time = start_time;
    int iteration = 0;
      //Pre-calculation of f
    //ForceCalculator::LennardJonesForceFaster(container, eps, sig);
   ForceCalculator::LennardJonesForceCell(cells, eps, sig);
    //Initialization with Brownian Motion
    //VelocityCalculator::BrownianMotionInitialization(container, avg_v, dim);
   /* VelocityCalculator::BrownianMotionInitializationCell(cells, avg_v, dim);
    //For this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
        // calculate new x
        //PositionCalculator::PositionStoermerVerlet(container, delta_t);
        PositionCalculator::PositionStoermerVerletCell(cells, delta_t);
        // calculate new f
        //ForceCalculator::LennardJonesForceFaster(container, eps, sig);
        ForceCalculator::LennardJonesForceCell(cells, eps, sig);
        // calculate new v
        //VelocityCalculator::VelocityStoermerVerlet(container, delta_t);
        VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);

        iteration++;
        if (iteration % 10 == 0) {
            //plotParticles(iteration);
            plotParticlesInCells(iteration, cells);
        }
        if (iteration % 100 == 0) {
            spdlog::info("Iteration " + std::to_string(iteration) + " finished.");
        }

        current_time += delta_t;
    }

*/
    spdlog::info("Output written. Terminating..." );
    return 0;
}

void plotParticlesInCells(int iteration, LinkedCellContainer &grid) {

    std::string out_name("MD_vtk");

    //outputWriter::XYZWriter writer;
    //outputWriter::XYZWriter::plotParticles(container, out_name, iteration);

    outputWriter::VTKWriter writer2;
    writer2.initializeOutput(grid.cell_numbers());
    for (int i = 0; i < grid.cell_numbers(); i++) {
        for(int j = 0; j < grid.Particles_in_cell(i); j++){
            writer2.plotParticle(grid.cells[i][j]);
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
