//
// Created by maraconda on 10.12.23.
//

#include "Thermostat.h"
#include "LinkedCellContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/ArrayUtils.h"
#include <cmath>

/**
 * initialize the temperature in a LinkedCellContainer (velocities of particles should not be (0,0,0))
 * @param initialTemperature temperature to initialize
 * @param dimension dimension of the simulation (possible values: 2 or 3)
 * @param cells LinkedCellContainer
 */
void Thermostat::initializeTemperature(double initialTemperature, int dimension, LinkedCellContainer &cells) {
   setTemperatureDirectly(initialTemperature, dimension,cells);
}

/**
 * initialize the temperature in a LinkedCellContainer with Brownian Motion (velocities of particles should be (0,0,0))
 * @param initialTemperature temperature to initialize
 * @param dimension dimension of the simulation (possible values: 2 or 3)
 * @param averageVelocity average velocity of the Brownian Motion
 * @param cells LinkedCellContainer
 */
void Thermostat::initializeTemperatureWithBrownianMotion(double initialTemperature, int dimension, LinkedCellContainer &cells) {
    double factor;
    std::array<double, 3> brownian_motion{};

    //iterate over all particles (without ghost particles)
    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); p++) {

                    factor = std::sqrt(initialTemperature / p->getM());
                    brownian_motion = maxwellBoltzmannDistributedVelocity(factor, dimension);
                    p->setV(p->getV() + brownian_motion);
                }
            }
        }
    }
}

/**
 * set the temperature in a LinkedCellContainer directly to a certain value with velocity scaling (velocities should already be initialized)
 * @param newTemperature temperature
 * @param dimension dimension of the simulation (possible values: 2 or 3)
 * @param cells LinkedCellContainer
 */
void Thermostat::setTemperatureDirectly(double newTemperature, int dimension, LinkedCellContainer &cells) {

    //Calculate current temperature and velocity scaling factor
    double currentTemperature = calculateCurrentTemperature(dimension, cells);
    double factor = std::sqrt( newTemperature/currentTemperature);

    //Scale velocities of all particles (without ghost particles)
    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); p++) {
                    p->setV(factor * (p->getV()));
                }
            }
        }
    }
}

/**
 * set the temperature in a LinkedCellContainer gradually to a certain value with velocity scaling (velocities should already be initialized)
 * @param targetTemperature temperature to be reached eventually
 * @param temperatureDifference maximal absolute temperature change allowed for one application of the thermostat
 * @param dimension dimension of the simulation (possible values: 2 or 3)
 * @param cells LinkedCellContainer
 * @param newTemperature the last applied temperature
 * @returns new applied temperature
 */
double Thermostat::setTemperatureGradually(double targetTemperature, double temperatureDifference, int dimension, LinkedCellContainer &cells, double newTemperature) {

    //Calculate current temperature
    double currentTemperature = calculateCurrentTemperature(dimension, cells);

    //Calculate the new temperature to set based on the allowed difference
    if (std::abs(targetTemperature - currentTemperature) <= temperatureDifference) {
        newTemperature = targetTemperature;
    }
    else {
        if (targetTemperature > newTemperature) {
            newTemperature += temperatureDifference;
        }
        else if (targetTemperature < newTemperature){
            newTemperature -= temperatureDifference;
            }
    }
    if (newTemperature < 0) {
        newTemperature = 0;
    }

    //Calculate the velocity scaling factor
    double factor = std::sqrt(std::abs(newTemperature)/currentTemperature);

    //Scale velocities of all particles (without ghost particles)
    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); p++) {

                    p->setV(factor * p->getV());
                }
            }
        }
    }

    return newTemperature;
}

/**
 * Calculate the current temperature in a LinkedCellContainer
 * @param dimension dimension of the simulation (possible values: 2 or 3)
 * @param cells LinkedCellContainer
 * @return current temperature in the LinkedCellContainer
 */
double Thermostat::calculateCurrentTemperature(int dimension, LinkedCellContainer cells) {
    double kineticEnergy = 0;
    int numberOfParticles = 0;
    std::array<double, 3> v_multiplication{};

    //Calculate current temperature with
    //T = (sum from i = 1 to #particles (m_i * <v_i,v_i>)) / #dimensions * #particles

    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); p++) {
                    v_multiplication = p->getV() * p->getV();
                    kineticEnergy += p->getM() * (v_multiplication[0] + v_multiplication[1] + v_multiplication[2]);
                    numberOfParticles++;
                }
            }
        }
    }

    //Calculate current temperature
    double currentTemperature = kineticEnergy / (numberOfParticles * dimension);
    return currentTemperature;
}