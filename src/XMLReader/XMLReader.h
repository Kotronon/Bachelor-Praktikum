#pragma once

#include "../ParticleGenerator.h"
#include "../LinkedCellContainer.h"
#include "../ParticleContainer.h"
#include "input.hxx"

#include <vector>
#include <string>
#include  <spdlog/spdlog.h>

class XMLReader {

public:


	struct XMLInfo {
		// general parameters of the simulation
		double epsilon;
		double sigma;
		double t_end;
		double delta_t;
		// all inputFiles (format of week 1)
		std::vector<std::string> inputFiles;
        double cutoff;
		double delta_r;
		double maxDistance;
        std::array<double, 3> dim;
        bool linkedCells;
        std::array<std::string, 6> boundary;
        LinkedCellContainer cells;
        ParticleContainer container;
        spdlog::level::level_enum level = spdlog::level::off;
	};

	static XMLInfo readFile(const std::string &s);

};