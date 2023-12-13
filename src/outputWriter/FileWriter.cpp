#include "FileWriter.h"
#include <iostream>
#include <fstream>

void FileWriter::writeFile(ParticleContainer &container, const std::string &filename) {
    //create and open file
     std::ofstream MyFile(filename);
     MyFile << container.size() << "\n";
     for(auto &particle : container){
         //write to file
         MyFile << particle.getX()[0] << " " << particle.getX()[1] << " " << particle.getX()[2]
         << " " << particle.getV()[0] << " " << particle.getV()[1] << " " << particle.getV()[2]
         << " " << particle.getF()[0] << " " << particle.getF()[1] << " " << particle.getF()[2]
         << " " << particle.getOldF()[0] << " " << particle.getOldF()[1] << " " << particle.getOldF()[2]
         << " " << particle.getM() << " " << particle.getSig() << " " << particle.getEps()
         << " " << particle.getType() << "\n";
     }
     MyFile.close();
}