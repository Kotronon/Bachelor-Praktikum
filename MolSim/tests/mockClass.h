//
// Created by kathi on 14.11.23.
//
#include "../src/ParticleContainer.h"
#include "../src/PositionCalculator.h"
#include "../src/VelocityCalculator.h"
#include "../src/ForceCalculator.h"
#include "../src/Particle.h"
#include "gmock/gmock.h"

namespace {

    class MockParticle : public Particle{
    public:
        MOCK_METHOD(void, Particle, ((std::array<double, 3>) x_arg, (std::array<double, 3>) v_arg, double m_arg, int type_arg), ());
        MOCK_METHOD((std::array<double, 3>),getX, (), (const));
        MOCK_METHOD((std::array<double, 3>),getV, (), (const));
        MOCK_METHOD((std::array<double, 3>),getF, (), (const));
        //MOCK_METHOD(bool, (operator==), (Particle & other), ());
    };

    class MockParticleContainer : public ParticleContainer {
    public:
        MOCK_METHOD(void, addPartocle,((std::array<double, 3>) x_arg, (std::array<double, 3>) v_arg, double m_arg, int type_arg), ());
        MOCK_METHOD((std::vector<Particle>::iterator), begin,(), ());
        MOCK_METHOD((std::vector<Particle>::iterator), end, (), ());
        MOCK_METHOD(int, size,(), ());
    };

    class MockForceCalculator : public ForceCalculator {
    public:
        MOCK_METHOD(void, SimpleForceCalculation,(ParticleContainer &container), ());

        MOCK_METHOD(void, LennardJonesForce,(ParticleContainer &container, double eps, double sig), ());
    };

    class MockPositionCalculator : public PositionCalculator {
    public:
        MOCK_METHOD(void, PositionStoermerVerlet,(ParticleContainer &container, double delta_t),());
    };

    class MockVelocityCalculator : public VelocityCalculator {
    public:
        MOCK_METHOD(void, BrownianMotionInitialization,(ParticleContainer &container, double avg_v),());

        MOCK_METHOD(void, VelocityStoermerVerlet,(ParticleContainer &container, double delta_t),());
    };

}