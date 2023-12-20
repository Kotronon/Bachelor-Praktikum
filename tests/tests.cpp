
#include "gtest/gtest.h"
#include "../src/ParticleContainer.h"
#include "../src/calculations/PositionCalculator.h"
#include "../src/calculations/VelocityCalculator.h"
#include "../src/calculations/ForceCalculator.h"
#include "../src/ParticleGenerator.h"
#include "../src/LinkedCellContainer.h"
#include "src/Thermostat.h"
#include <math.h>
//@TODO write tests with TEST()
//To compile tests write cmake --build . in terminal and afterwarts ctest works

/**
 * Test for plain Particle Container for addParticle, begin(), end() and size() functions
 */
TEST(ParticleContainerTest, ParticleContainer){
    ParticleContainer particles = ParticleContainer();
    //test addParticle
    Particle particle({0, 0, 0}, {0, 3, 0}, 50, 0, 1, 5);
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0, 1, 5);
    //test begin
    //ASSERT_TRUE(particles.begin().base()->operator==(particle)) << "didn't point on only element in container";
    //test end
    ASSERT_FALSE(particles.end().base()->operator==(particle)) << "didn't point on only element in container";
    //test size
    EXPECT_EQ(1, particles.size()) << "wrong size";
}

/**
 * Test for linked cell container for addParticle, begin(), end() and size() functions as well as to calculate position and move particle to next cell
 */
TEST(cellTest, LinkedCellContainer){
    LinkedCellContainer cells = LinkedCellContainer({180, 90, 1}, 3.0, {"r", "r", "r", "o", "r", "r"});
    //test cells size
    EXPECT_EQ(1800, cells.cell_numbers()) << "wrong size";
    //test addParticle
    Particle particle({0, 0, 0}, {0, 3, 0}, 50, 0, 1, 5);
    cells.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0, 1, 5);
    EXPECT_EQ(1, cells.Particles_in_cell(1,1,1));
    //test calculate position and move particle to next cell
    Particle particle2({0, 3, 0}, {0, 3, 0}, 50, 0, 1, 5);
    PositionCalculator::PositionStoermerVerletCell(cells, 1);
    EXPECT_EQ(0, cells.Particles_in_cell(1,1,1));
    EXPECT_EQ(1, cells.Particles_in_cell(1,2,1));
    EXPECT_EQ(0, cells.Particles_in_cell(1,3,1));

}
/**
 * Test for outflow, reflective and periodic boundary
 */
TEST(Boundary, Boundry){
    LinkedCellContainer cells = LinkedCellContainer({180, 90, 1}, 3.0, {"r", "p", "p", "o", "r", "r"});
    //test outflow
    cells.addParticle({0, 0, 0}, {0, -1, 0}, 50, 0, 1, 5);
    EXPECT_EQ(1, cells.Particles_in_cell(1,1,1));
    PositionCalculator::PositionStoermerVerletCell(cells, 1);
    EXPECT_EQ(0, cells.Particles_in_cell(1,1,1));
    //test reflective boundary
    cells.addParticle({0, 0, 0}, {-1, 3, 0}, 50, 0, 1, 5);
    ASSERT_EQ(1, cells.Particles_in_cell(1,1,1));
    cells.generateGhostCell(0, 1, 1, 1);
    Particle ghost ({-0.0000000001, 0, 0}, {0, 0, 0}, 50, 1, 5, -1);
    EXPECT_EQ(1, cells.Particles_in_cell(0,1,1));
    //test periodic
    cells.addParticle({179, 89, 0}, {1,1,0}, 50, 0, 1, 5);
    cells.generateGhostCell(0, 60, 30, 1);
    ASSERT_EQ(1, cells.Particles_in_cell(60, 0, 1));
    ASSERT_EQ(1, cells.Particles_in_cell(0, 30, 1));
    ASSERT_EQ(1, cells.Particles_in_cell(0, 0, 1));

}

/**
 * test for creating a cuboid and adding it to a plain container
 * test for creating a cuboid and adding it to a linked cell container
 * test for creating a sphere and adding it to a linked cell container
 */
TEST(ParticleGeneratorTest, ParticleGenerator){
    //Test of adding the right numbers of particles when generating a cuboid to plain container
    ParticleContainer particles = ParticleGenerator::createCuboid({0,0,0}, {0,0,0}, {40,8,1}, 2, 3, 1, 5, 0);
    //test size
    EXPECT_EQ(320, particles.size()) << "wrong number of particles generated in plain container";
    //Test of adding the right numbers of particles when generating a cuboid to linked cell container
    LinkedCellContainer cells = LinkedCellContainer({180, 90, 1}, 3.0, {"r", "r", "r", "r", "r", "r"});
    ParticleGenerator::createCuboidInCells({0,0,0}, {0,0,0}, {40,8,1}, 2, 3, cells, 3.0, 1, 5, 0);
    int particles_num = 0;
    for (auto &x: cells) {
        for (auto &y: x) {
            for (auto &z: y) {
                for (auto &p: z) {
                    particles_num++;
                }
            }
        }
    }
    EXPECT_EQ(320, particles_num) << "wrong number of particles generated in linked cell container";
}
/**
 * Test for calculation the new Position of a particle
 */
TEST(PositionTest, stroemerVelvet){
    //Checking the correctness of the stoemer velvet position calculation
    //for force = 0
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0, 1, 5);
    std::array<double, 3> res = {0, 0.042, 0};
    PositionCalculator::PositionStoermerVerlet(particles, 0.014);
    auto particleVector =  particles.begin();
    EXPECT_EQ(res, particleVector.base()->getX()) << "wrong position calculated";

}
/**
 * Test for the calculation of the velocity of a particle
 */
TEST(VelocityTest, stroemerVelvet) {
    //Checking the correctness of the stoemer velvet velocity calculation
    //for force = 0
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0, 1, 5);
    std::array<double, 3> res = {0, 3, 0};
    VelocityCalculator::VelocityStoermerVerlet(particles, 0.014);
    auto particleVector = particles.begin();
    EXPECT_EQ(res, particleVector.base()->getV()) << "wrong velocity calculated";
}

/**
 * Test for the calculation of gravity force
 */
TEST(ForceTest, SimpleForceCalculation){
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0, 1, 5);
    particles.addParticle({0, 3, 0}, {4, 3, 5}, 20, 0, 1, 5);
    double res = (50.0*20.0*3.0)/27.0;
    std::array<double, 3> res1 = {0, res, 0};
    std::array<double, 3> res2 = {0, -res, 0};
    ForceCalculator::GravityForceCalculation(particles);
    auto particleVector = particles.begin();
    EXPECT_EQ(res1, particleVector.base()->getF());
    EXPECT_EQ(res2, std::next(particleVector.base())->getF()) << "wrong force calculated";
}

/**
 * Test for the Lennard Jones force calculation
 */
TEST(ForceTest, LennardJonesForce){
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0, 1, 5);
    particles.addParticle({0, 3, 0}, {4, 3, 5}, 20, 0, 1, 5);
    double res = (24.0/9.0)*(pow(1.0/3.0, 6.0)-2*pow(1.0/3.0, 12.0)) * 3;
    std::array<double, 3> res1 = {0, res, 0};
    std::array<double, 3> res2 = {0, -res, 0};
    ForceCalculator::LennardJonesForce(particles, 1, 1);
    auto particleVector = particles.begin();
    EXPECT_EQ(res1, particleVector.base()->getF());
    EXPECT_EQ(res2, std::next(particleVector.base())->getF()) << "wrong force calculated";
}

/**
 * tests for the thermostat including heating up, cooling down, holding a temperature
 */
TEST(ThermostatTest, Thermostat){

    //Setting up LinkedCellContainer and some particles
    LinkedCellContainer cells = LinkedCellContainer({30,30,1},2.5,{"r,r,r,r,o,o"});
    ParticleGenerator::createDiskInCells({15, 15, 1}, {0, 0, 0}, 1.0, 3, 1.2, cells, 1.0, 1.0, 1);

    //Setting up parameters
    double grav = -12.44;
    double delta_t = 0.0005;
    double init_T = 20;
    double higher_T = 30;
    double lower_T = 10;

    //Setup before first iteration and initialization with Brownian Motion to a temperature of 20 Kelvin
    ForceCalculator::LennardJonesForceCell(cells, grav);
    Thermostat::initializeTemperatureWithBrownianMotion(init_T, 2, cells);

    //Check initial setting of temperature
    ASSERT_TRUE(std::abs(Thermostat::calculateCurrentTemperature(2, cells) - 20) < 5) << "wrong temperature after temperature initialization";

    //Iteration 1 (Holding the temperature)
    PositionCalculator::PositionStoermerVerletCell(cells, delta_t);
    ForceCalculator::LennardJonesForceCell(cells, grav);
    VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);
    Thermostat::setTemperatureDirectly(init_T,2,cells);

    //Check temperature after first iteration
    ASSERT_TRUE(std::abs(Thermostat::calculateCurrentTemperature(2, cells) - 20) < 0.05) << "wrong temperature while holding temperature";

    //Iteration 2 (Holding the temperature)
    PositionCalculator::PositionStoermerVerletCell(cells, delta_t);
    ForceCalculator::LennardJonesForceCell(cells, grav);
    VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);
    Thermostat::setTemperatureDirectly(init_T,2,cells);

    //Check temperature after second iteration
    ASSERT_TRUE(std::abs(Thermostat::calculateCurrentTemperature(2, cells) - 20) < 0.05) << "wrong temperature while holding temperature";

    //Iteration 3 (heating up the simulation)
    PositionCalculator::PositionStoermerVerletCell(cells, delta_t);
    ForceCalculator::LennardJonesForceCell(cells, grav);
    VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);
    Thermostat::setTemperatureDirectly(higher_T,2,cells);

    //Check temperature after third iteration
    ASSERT_TRUE(std::abs(Thermostat::calculateCurrentTemperature(2, cells) - 30) < 0.05) << "wrong temperature while heating up";

    //Iteration 4 (cooling down the simulation)
    PositionCalculator::PositionStoermerVerletCell(cells, delta_t);
    ForceCalculator::LennardJonesForceCell(cells, grav);
    VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);
    Thermostat::setTemperatureDirectly(lower_T,2,cells);

    //Check temperature after fourth iteration
    ASSERT_TRUE(std::abs(Thermostat::calculateCurrentTemperature(2, cells) - 10) < 0.05) << "wrong temperature while cooling down";
}

int main(){
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
   // return 0;
}