
#pragma once

#include "Particle.h"
#include "ParticleContainer.h"

#include <list>
#include <forward_list>

class FileWriter {

public:
  FileWriter() = default;
  virtual ~FileWriter() = default;

  static void writeFile(ParticleContainer &container, const std::string &filename);

};
