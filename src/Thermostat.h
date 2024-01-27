//
// Created by maraconda on 10.12.23.
//

#ifndef PSEMOLDYN_GROUPH_THERMOSTAT_H
#define PSEMOLDYN_GROUPH_THERMOSTAT_H


#include "LinkedCellContainer.h"

class Thermostat {
public:

    Thermostat() = default;

    static void initializeTemperature(double initialTemperature, int dimension, LinkedCellContainer &cells);

    static void setTemperatureDirectly(double temperature, int dimension, LinkedCellContainer &cells);

    static double
    setTemperatureGradually(double targetTemperature, double temperatureDifference, int dim,
                            LinkedCellContainer &cells, double deltaT);

    static double calculateCurrentTemperature(int dimension, LinkedCellContainer cells);

    static void initializeTemperatureWithBrownianMotion(double initialTemperature, int dimension, LinkedCellContainer &cells);

private:

};


#endif //PSEMOLDYN_GROUPH_THERMOSTAT_H
