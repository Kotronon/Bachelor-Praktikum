/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include "C:\Users\katha\Bachelor-Praktikum\MolSim\src\Particle.h"

#include <fstream>
#include <list>

namespace outputWriter {

class XYZWriter {

public:
  XYZWriter();

  virtual ~XYZWriter();

  void plotParticles(std::list<Particle> particles, const std::string &filename,
                     int iteration);
};

} // namespace outputWriter
