/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "ParticleContainer.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

FileReader::FileReader() = default;

FileReader::~FileReader() = default;

/**
 * reads file and stores particles in given ParticleContainer
 * @param container
 * @param filename
 */

void FileReader::readFile(ParticleContainer &container, char *filename) {
  std::array<double, 3> x{};
  std::array<double, 3> v{};
  std::array<double, 3> f{};
  std::array<double, 3> oldF{};
  double m;
  double sig;
  double eps;
  int num_particles = 0;
  int type;

  std::ifstream input_file(filename);
  std::string tmp_string;

  if (input_file.is_open()) {

    getline(input_file, tmp_string);
    std::cout << "Read line: " << tmp_string << std::endl;

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << std::endl;
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    std::cout << "Reading " << num_particles << "." << std::endl;
    getline(input_file, tmp_string);
    std::cout << "Read line: " << tmp_string << std::endl;

    for (int i = 0; i < num_particles; i++) {
      std::istringstream datastream(tmp_string);

      for (auto &xj : x) {
        datastream >> xj;
      }
      for (auto &vj : v) {
        datastream >> vj;
      }
      for(auto &fj : f){
          datastream >> fj;
      }
      for(auto &oldFj : oldF){
          datastream >> oldFj;
      }
      if (datastream.eof()) {
        std::cout
            << "Error reading file: eof reached unexpectedly reading from line "
            << i << std::endl;
        exit(-1);
      }
      datastream >> m;
      datastream >> sig;
      datastream >> eps;
      datastream >> type;
      Particle new_particle = Particle(x, v, m, sig, eps, type);
      new_particle.setF(f);
      new_particle.setOldF(oldF);
      container.addParticle(new_particle);
      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << std::endl;
    }
  } else {
    std::cout << "Error: could not open file " << filename << std::endl;
    exit(-1);
  }
}
