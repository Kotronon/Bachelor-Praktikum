//
// Created by kathi on 14.11.23.
//
#include "../src/ParticleContainer.h"
#include "../src/PositionCalculator.h"
#include "../src/VelocityCalculator.h"
#include "../src/ForceCalculator.h"
#include "gmock/gmock.h"


class MockParticleContainer: public ParticleContainer {
public:
    MOCK_METHOD(void , addPartocle, (const Particle &particle));
    MOCK_METHOD((std::vector<Particle>::iterator) , begin, ());
    MOCK_METHOD((std::vector<Particle>::iterator) , end, ());
    MOCK_METHOD(int , size, ());
};

class MockForceCalculator: public  ForceCalculator {
public:
    MOCK_METHOD(void, SimpleForceCalculation, (ParticleContainer &container));
    MOCK_METHOD(void, LennardJonesForce, (ParticleContainer &container, double eps, double sig));
};

class MockPositionCalculator: public  PositionCalculator {
public:
    MOCK_METHOD(void, PositionStoermerVerlet, (ParticleContainer &container, double delta_t));
};

class MockVelocityCalculator: public VelocityCalculator {
public:
    MOCK_METHOD(void, BrownianMotionInitialization, (ParticleContainer &container, double avg_v));
    MOCK_METHOD(void, VelocityStoermerVerlet, (ParticleContainer &container, double delta_t));
};