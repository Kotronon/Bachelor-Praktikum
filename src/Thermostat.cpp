//
// Created by maraconda on 10.12.23.
//

#include "Thermostat.h"
#include "LinkedCellContainer.h"
#include <cmath>
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/ArrayUtils.h"

/**
 * initialize the temperature in a LinkedCellContainer (velocities of particles should not be (0,0,0))
 * @param initialTemperature temperature to initialize in Kelvin
 * @param dimension dimension of the simulation (possible values: 2 or 3)
 * @param cells LinkedCellContainer
 */
void Thermostat::initializeTemperature(double initialTemperature, int dimension, LinkedCellContainer &cells) {
    setTemperatureDirectly(initialTemperature, dimension,cells);
}

/**
 * initialize the temperature in a LinkedCellContainer with Brownian Motion (velocities of particles should be (0,0,0))
 * @param initialTemperature temperature to initialize in Kelvin
 * @param dimension dimension of the simulation (possible values: 2 or 3)
 * @param averageVelocity average velocity of the Brownian Motion
 * @param cells LinkedCellContainer
 */
void Thermostat::initializeTemperatureWithBrownianMotion(double initialTemperature, int dimension, double averageVelocity, LinkedCellContainer &cells) {
    double factor;
    std::array<double, 3> brownian_motion{};

    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); z++) {
                    factor = std::sqrt(initialTemperature / p->getM());
                    brownian_motion = maxwellBoltzmannDistributedVelocity(factor * averageVelocity, dimension);
                    p->setV(p->getV() + brownian_motion);
                }
            }
        }
    }
}

/**
 * set the temperature in a LinkedCellContainer directly to a certain value with velocity scaling (velocities should already be initialized)
 * @param newTemperature temperature in Kelvin
 * @param dimension dimension of the simulation (possible values: 2 or 3)
 * @param cells LinkedCellContainer
 */
void Thermostat::setTemperatureDirectly(double newTemperature, int dimension, LinkedCellContainer &cells) {
    double sum = 0;
    double factor;
    int numberOfParticles = 0;
    std::array<double, 3> v_multiplication{};

    //Calculate current temperature with
    //T = (sum from i = 1 to #particles (m_i * <v_i,v_i>)) / #dimensions * #particles

    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); z++) {
                    v_multiplication = p->getV() * p->getV();
                    sum += p->getM() * (v_multiplication[0] + v_multiplication[1] + v_multiplication[2]);
                    numberOfParticles++;
                }
            }
        }
    }

    //Calculate current temperature and velocity scaling factor
    double currentTemperature = sum / (numberOfParticles * dimension);
    factor = std::sqrt(newTemperature/currentTemperature);

    //Scale velocities of all particles
    for (auto x = cells.begin() + 1; x < cells.end() - 1; x++) {
        for (auto y = x->begin() + 1; y < x->end() - 1; y++) {
            for (auto z = y->begin() + 1; z < y->end() - 1; z++) {
                for (auto p = z->begin(); p < z->end(); z++) {
                    p->setV(factor * (p->getV()));
                }
            }
        }
    }

}

void Thermostat::setTemperatureGradually(double targetTemperature, int timeSteps, LinkedCellContainer &cells){}

void Thermostat::setTemperatureGradually(double targetTemperature, int timeSteps, double temperatureDifference, LinkedCellContainer &cells){}
