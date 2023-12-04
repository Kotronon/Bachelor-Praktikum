#include <iostream>

#include "XMLReader.h"

XMLReader::XMLInfo XMLReader::readFile(const std::string &s) {
    auto sim = Simulation(s, xml_schema::flags::dont_validate);

	XMLInfo  info = {};
    double endTime = sim->endTime();
    double deltaT = sim->deltaT();
    if (sim->Strategy().LinkedCell().present()) {
      info.linkedCells = true;
      info.cutoff = sim->Strategy().LinkedCell().get().cutoff();
		auto dimension = sim->Strategy().LinkedCell().get().Domain();
        info.dim = {dimension.x(), dimension.y(), dimension.z()};
		auto borderType = sim->Strategy().LinkedCell().get().Boundary();
		std::array<std::string, 6> boundaryConds = {"o", "o", "o", "o", "o", "o"};
		boundaryConds[0] = borderType.boundary_left().get();
		boundaryConds[1] = borderType.boundary_right().get();
		boundaryConds[2] = borderType.boundary_top().get();
		boundaryConds[3] = borderType.boundary_bottom().get();
        boundaryConds[4] = borderType.boundary_back().get();
		boundaryConds[5] = borderType.boundary_front().get();

		info.boundary = boundaryConds;
        info.cells = LinkedCellContainer(info.dim, info.cutoff, boundaryConds);
	} else {
		info.linkedCells = false;
        info.container = ParticleContainer();
	}
     for (auto &it: sim->Shapes()) {
      for (auto &sphere: it.Sphere()) {
          if(info.linkedCells)
            ParticleGenerator::createDiskInCells({sphere.Center().x(), sphere.Center().y(), sphere.Center().z()}, {sphere.Velocity().x(),sphere.Velocity().y(), sphere.Velocity().z()},
                                                 sphere.mass(), sphere.radius(), sphere.distance(), info.cells);
          else {
              ParticleContainer helper = ParticleGenerator::createDisk({sphere.Center().x(), sphere.Center().y(), sphere.Center().z()}, {sphere.Velocity().x(),sphere.Velocity().y(), sphere.Velocity().z()},
                                                                         sphere.mass(), sphere.radius(), sphere.distance());
              info.container.addParticleContainer(helper);
          }
      }
    }
     for (auto &it: sim->Shapes()) {
      for (auto &cuboid: it.Cuboid()) {
          if (info.linkedCells)
              ParticleGenerator::createCuboidInCells(
                      {cuboid.Position().x(), cuboid.Position().y(), cuboid.Position().z()},
                      {cuboid.Velocity().x(), cuboid.Velocity().y(), cuboid.Velocity().z()},
                      {cuboid.Dimension().x(), cuboid.Dimension().y(), cuboid.Dimension().z()}, cuboid.distance(),
                      cuboid.mass(),
                      info.cells, info.cutoff);
          else {
             ParticleContainer new_container = ParticleGenerator::createCuboid(
                      {cuboid.Position().x(), cuboid.Position().y(), cuboid.Position().z()},
                      {cuboid.Velocity().x(), cuboid.Velocity().y(), cuboid.Velocity().z()},
                      {cuboid.Dimension().x(), cuboid.Dimension().y(), cuboid.Dimension().z()}, cuboid.distance(),
                      cuboid.mass());
             info.container.addParticleContainer(new_container);
          }
      }
    }
    if(sim->logLevel().present()){
        auto level = sim->logLevel().get();
        if(level == "off"){
            info.level= spdlog::level::off;
        }
        if(level == "trace"){
            info.level= spdlog::level::trace;
        }
        if(level == "debug"){
            info.level= spdlog::level::debug;
        }
        if(level == "info"){
            info.level= spdlog::level::info;
        }
        if(level == "warn"){
            info.level= spdlog::level::warn;
        }
        if(level == "error"){
            info.level= spdlog::level::err;
        }
        if(level == "critical"){
            info.level= spdlog::level::critical;
        }
    }
    return info;
}