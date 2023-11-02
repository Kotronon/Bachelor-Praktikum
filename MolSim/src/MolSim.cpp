
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>
#include <list>
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

std::list<Particle> particles;

int main(int argc, char *argsv[]) {

    std::cout << "Hello from MolSim for PSE!" << std::endl;
    if (argc <= 2 || argc >= 8) {
        std::cout << "Erroneous programme call! " << std::endl;
        std::cout << "./MolSim <filepath/filename> [options]"<< std::endl;
        std::cout << "Options: "<< std::endl;
        std::cout << "-e : The end time of the simulation, default value is 1000"<< std::endl;
        std::cout << "-d : ∆time of the simulation, default value is 0.014"<< std::endl;
    }

    FileReader fileReader;
    FileReader::readFile(particles, argsv[1]);


    //getting end time and delta t command options if specified

    if(cmdOptionExists(argsv, argsv+argc, "-e"))
    {
        end_time = std::stod(getCmdOption(argsv, argsv + argc, "-e"));
    }

    if(cmdOptionExists(argsv, argsv+argc, "-d"))
    {
        end_time = std::stod(getCmdOption(argsv, argsv + argc, "-d"));
    }

    /*
    std::string input;
    std::cout << "Please enter the end time. If you want to use the default value 1000 please enter x." << std::endl;
    bool valid_input_received = false;

    while (!valid_input_received){
        std::cin >> input;
        if (input != "x"){
            try{
                end_time = std::stod(input);
                valid_input_received = true;
            }
            catch(std::invalid_argument& e){
                std::cout << "Invalid input. Please enter a valid value or x for the default value."<< std::endl;
            }
        }
        else{
            valid_input_received = true;
        }
    }


    std::cout << "Please enter the delta time. If you want to use the default value 0.014 please enter x." << std::endl;
    valid_input_received = false;

    while (!valid_input_received){
        std::cin >> input;
        if (input != "x"){
            try{
                delta_t = stod(input);
                valid_input_received = true;
            }
            catch(std::invalid_argument& e){
                std::cout << "Invalid input. Please enter a valid value or x for the default value."<< std::endl;
            }
        }
        else{
            valid_input_received = true;
        }
    }
    */

    double current_time = start_time;

    int iteration = 0;

    // for this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
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
        std::cout << "Iteration " << iteration << " finished." << std::endl;

        current_time += delta_t;
    }


    std::cout << "output written. Terminating..." << std::endl;
    return 0;
}

void calculateF() {
    std::array<double, 3> force{};
    for (auto &p1: particles) {
        force = {0., 0., 0.};
        for (auto &p2: particles) {
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
    for (auto &p: particles) {
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
    for (auto &p: particles) {
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
    outputWriter::XYZWriter::plotParticles(particles, out_name, iteration);

    outputWriter::VTKWriter writer2;
    writer2.initializeOutput(particles.size());
    for (auto &p: particles) {
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
