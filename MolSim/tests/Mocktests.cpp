
#include <gtest/gtest.h>
#include "gmock/gmock.h"
#include "mockClass.h"




//@TODO write tests with TEST()
//To compile tests write cmake --build . in terminal and afterwarts ctest works
using ::testing::Mock;

TEST(ParticleContainerTest, ParticleContainer){
    MockParticleContainer particles;
    //test addParticle
    MockParticle particle;
    particle.Particle({0, 0, 0}, {0, 3, 0}, 50, 0);
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    //test begin
    ASSERT_TRUE(particles.begin().base()->operator==(particle));
    //test end
    ASSERT_TRUE(particles.end().base()->operator==(particle));
    //test size
    EXPECT_EQ(1, particles.size());
}

/*TEST(PositionTest, stroemerVelvet){
    //Checking the correctness of the stoemer velvet position calculation
    //for force = 0
    ::testing::NiceMock<MockParticleContainer> particles;
    //MockParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    std::array<double, 3> res = {0, 0.042, 0};
    MockPositionCalculator calculator;
    calculator.PositionStoermerVerlet(particles, 0.014);
    std::vector<Particle>::iterator particleVector =  particles.begin();
    EXPECT_EQ(res, particleVector.base()->getX());

}

TEST(VelocityTest, BrownianMotionInitialization){
    //Checking the correctness of the Brownian Motion Initialization velocity calculation
    MockParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
}

TEST(VelocityTest, stroemerVelvet) {
    //Checking the correctness of the stoemer velvet velocity calculation
    //for force = 0
    MockParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    std::array<double, 3> res = {0, 3, 0};
    MockVelocityCalculator calculator;
    calculator.VelocityStoermerVerlet(particles, 0.014);
    std::vector<Particle>::iterator particleVector = particles.begin();
    EXPECT_EQ(res, particleVector.base()->getV());
}

TEST(ForceTest, SimpleForceCalculation){
    MockParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    particles.addParticle({0, 3, 0}, {4, 3, 5}, 20, 0);
    std::array<double, 3> res1 = {0, 1000, 0};
    std::array<double, 3> res2 = {0, -1000, 0};
    MockForceCalculator calculator;
    calculator.SimpleForceCalculation(particles);
    std::vector<Particle>::iterator particleVector = particles.begin();
    EXPECT_EQ(res1, particleVector.base()->getF());
    EXPECT_EQ(res2, std::next(particleVector.base())->getF());
}

TEST(ForceTest, LennardJonesForce){
    MockParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    particles.addParticle({0, 3, 0}, {4, 3, 5}, 20, 0);
    std::array<double, 3> res1 = {0, -8/531441, 0};
    std::array<double, 3> res2 = {0, 8/531441, 0};
    MockForceCalculator calculator;
    calculator.LennardJonesForce(particles, 1, 1);
    std::vector<Particle>::iterator particleVector = particles.begin();
    EXPECT_EQ(res1, particleVector.base()->getF());
    EXPECT_EQ(res2, std::next(particleVector.base())->getF());
}*/

int main(){
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}