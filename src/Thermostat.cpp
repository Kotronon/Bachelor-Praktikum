//
// Created by maraconda on 10.12.23.
//

#include "Thermostat.h"
#include "LinkedCellContainer.h"
#include "cmath"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/ArrayUtils.h"

void Thermostat::initializeTemperature(double initialTemperature, int dimension, LinkedCellContainer &cells) {
    double factor;
    std::array<double, 3> brownian_motion{};

    for (auto &x: cells) {
        for (auto &y: x) {
            for (auto &z: y) {
                for (auto &p: z) {

                    factor = std::sqrt(initialTemperature / p.getM());
                    brownian_motion = maxwellBoltzmannDistributedVelocity(factor, dimension);
                    p.setV(p.getV() + brownian_motion);
                }
            }
        }
    }
}

void Thermostat::setTemperatureDirectly(double temperature, LinkedCellContainer &cells) {}

void Thermostat::setTemperatureGradually(double targetTemperature, int timeSteps, LinkedCellContainer &cells){}

void Thermostat::setTemperatureGradually(double targetTemperature, int timeSteps, double temperatureDifference, LinkedCellContainer &cells){}
