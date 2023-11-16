/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include "Particle.h"
#include "ParticleContainer.h"

#include <list>
#include <forward_list>

class FileReader {

public:
  FileReader();
  virtual ~FileReader();

  static void readFile(ParticleContainer &container, char *filename);

};
