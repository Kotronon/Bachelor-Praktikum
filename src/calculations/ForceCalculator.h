//
// Created by maraconda on 12.11.23.
//

#pragma once

#include "../ParticleContainer.h"
#include "../LinkedCellContainer.h"
#include "../Membrane.h"

class ForceCalculator {
private:
    static double epsilon;
    static double sigma;
    static double Ggrav;
    static double cutoff;
    static double k;
    static double r_0;
    static double f_z;
public:
    static void LateralForceCalculation(Particle *p1, Particle *p2);
    static void DiagonalForceCalculation(Particle *p1, Particle *p2);
    static void ThatOneForceCalculation(Particle *p1, Particle *p2, double fZ);
    static void GravityForceCalculation(ParticleContainer &container);
    static void LennardJonesForce(ParticleContainer &container, double eps, double sig);
    static void LennardJonesForceFaster(ParticleContainer &container, double eps, double sig, double Grav);
    static void LennardJonesForcePairwise(Particle *p1, Particle *p2);
    static void LennardJonesForceCell(LinkedCellContainer &cells, double Grav);
    static void LennardJonesForceMembrane(LinkedCellContainer &cells, double Grav);
    static void MembraneForceCalculation(LinkedCellContainer &cells, double Grav, double h);
    static void ThatOneMembraneForceCalculation(LinkedCellContainer &cells, double Grav, double fZ);



};


