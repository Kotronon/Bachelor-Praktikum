//
// Created by maraconda on 10.12.23.
//

#ifndef PSEMOLDYN_GROUPH_THERMOSTAT_H
#define PSEMOLDYN_GROUPH_THERMOSTAT_H


#include "LinkedCellContainer.h"

class Thermostat {
public:

    void setTemperatureGradually(double targetTemperature, int timeSteps, LinkedCellContainer &cells);

    void setTemperatureGradually(double targetTemperature, int timeSteps, double temperatureDifference,
                                 LinkedCellContainer &cells);

    void initializeTemperature(double initialTemperature, int dimension, LinkedCellContainer &cells);


    void initializeTemperatureWithBrownianMotion(double initialTemperature, int dimension, double averageVelocity,
                                                 LinkedCellContainer &cells);

    void setTemperatureDirectly(double temperature, int dimension, LinkedCellContainer &cells);
};


#endif //PSEMOLDYN_GROUPH_THERMOSTAT_H
