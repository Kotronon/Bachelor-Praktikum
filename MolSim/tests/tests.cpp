
#include <gtest/gtest.h>
#include "gmock/gmock.h"
#include "../src/ParticleContainer.h"
#include "../src/PositionCalculator.h"
#include "../src/VelocityCalculator.h"
#include "../src/ForceCalculator.h"
#include "../src/Particle.h"



//@TODO write tests with TEST()
//To compile tests write cmake --build . in terminal and afterwarts ctest works


TEST(ParticleContainerTest, ParticleContainer){
    ParticleContainer particles = ParticleContainer();
    //test addParticle
    Particle particle({0, 0, 0}, {0, 3, 0}, 50, 0);
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    //test begin
    ASSERT_TRUE(particles.begin().base()->operator==(particle));
    //test end
    ASSERT_FALSE(particles.end().base()->operator==(particle));
    //test size
    EXPECT_EQ(1, particles.size());
}

TEST(PositionTest, stroemerVelvet){
    //Checking the correctness of the stoemer velvet position calculation
    //for force = 0
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    std::array<double, 3> res = {0, 0.042, 0};
    PositionCalculator calculator;
    calculator.PositionStoermerVerlet(particles, 0.014);
    std::vector<Particle>::iterator particleVector =  particles.begin();
    EXPECT_EQ(res, particleVector.base()->getX());

}

TEST(VelocityTest, BrownianMotionInitialization){
    //Checking the correctness of the Brownian Motion Initialization velocity calculation
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
}

TEST(VelocityTest, stroemerVelvet) {
    //Checking the correctness of the stoemer velvet velocity calculation
    //for force = 0
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    std::array<double, 3> res = {0, 3, 0};
    VelocityCalculator calculator;
    calculator.VelocityStoermerVerlet(particles, 0.014);
    std::vector<Particle>::iterator particleVector = particles.begin();
    EXPECT_EQ(res, particleVector.base()->getV());
}

TEST(ForceTest, SimpleForceCalculation){
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    particles.addParticle({0, 3, 0}, {4, 3, 5}, 20, 0);
    std::array<double, 3> res1 = {0, 111.111, 0};
    std::array<double, 3> res2 = {0, -111.111, 0};
    ForceCalculator calculator;
    calculator.SimpleForceCalculation(particles);
    std::vector<Particle>::iterator particleVector = particles.begin();
    double x1;
    double x2;
    for (size_t i=0; i<3; i++) {
        x1 = (res1[i] - particleVector.base()->getF()[i]);
        x2 = (res1[i] - std::next(particleVector.base())->getF()[i]);
    }
    EXPECT_EQ(x1, 0);
    EXPECT_EQ(x2, 0);
}

TEST(ForceTest, LennardJonesForce){
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    particles.addParticle({0, 3, 0}, {4, 3, 5}, 20, 0);
    std::array<double, 3> res1 = {0, 0.0109438, 0};
    std::array<double, 3> res2 = {0, -0.0109438, 0};
    ForceCalculator calculator;
    calculator.LennardJonesForce(particles, 1, 1);
    std::vector<Particle>::iterator particleVector = particles.begin();
    double x1;
    double x2;
    for (size_t i=0; i<3; i++) {
        x1 = (res1[i] - particleVector.base()->getF()[i]);
        x2 = (res1[i] - std::next(particleVector.base())->getF()[i]);
    }
    EXPECT_EQ(x1, 0);
    EXPECT_EQ(x2, 0);
}

int main(){
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}