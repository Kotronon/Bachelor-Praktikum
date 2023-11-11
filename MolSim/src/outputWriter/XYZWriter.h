/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include  "Particle.h"
#include "ParticleContainer.h"

#include <fstream>
#include <list>
#include <forward_list>

namespace outputWriter {

class XYZWriter {

public:
  XYZWriter();

  virtual ~XYZWriter();

    static void plotParticles(ParticleContainer &container, const std::string &filename, int iteration);
};

} // namespace outputWriter
