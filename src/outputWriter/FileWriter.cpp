#include "FileWriter.h"
#include <iostream>
#include <fstream>

void FileWriter::writeFile(ParticleContainer &container, const std::string &filename) {
    //create and open file
     std::ofstream MyFile(filename);
     //amound of particles
     MyFile << container.size() << "\n";
     for(auto &particle : container){
         //write to file
         std::string input;
         input += std::to_string(particle.getX()[0]) + " " + std::to_string(particle.getX()[1]) + " " + std::to_string(particle.getX()[2])
         + " " + std::to_string(particle.getV()[0]) + " " + std::to_string(particle.getV()[1]) + " " + std::to_string(particle.getV()[2])
         + " " + std::to_string(particle.getF()[0]) + " " + std::to_string(particle.getF()[1]) + " " + std::to_string(particle.getF()[2])
         + " " + std::to_string(particle.getOldF()[0]) + " " + std::to_string(particle.getOldF()[1]) + " " + std::to_string(particle.getOldF()[2])
         + " " + std::to_string(particle.getM()) + " " + std::to_string(particle.getSig()) + " " + std::to_string(particle.getEps())
         + " " + std::to_string(particle.getType()) + "\n";
         //std::cout << input << std::endl;
         MyFile << input << std::endl;

     }
     //close file for good practice
     MyFile.close();
}