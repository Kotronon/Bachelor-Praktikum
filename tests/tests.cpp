
#include "gtest/gtest.h"
#include "../src/ParticleContainer.h"
#include "../src/calculations/PositionCalculator.h"
#include "../src/calculations/VelocityCalculator.h"
#include "../src/calculations/ForceCalculator.h"
#include "../src/ParticleGenerator.h"
#include "../src/LinkedCellContainer.h"
#include <math.h>
//@TODO write tests with TEST()
//To compile tests write cmake --build . in terminal and afterwarts ctest works

/**
 * Test for plain Particle Container for addParticle, begin(), end() and size() functions
 */
TEST(ParticleContainerTest, ParticleContainer){
    ParticleContainer particles = ParticleContainer();
    //test addParticle
    Particle particle({0, 0, 0}, {0, 3, 0}, 50, 0);
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    //test begin
    ASSERT_TRUE(particles.begin().base()->operator==(particle)) << "didn't point on only element in container";
    //test end
    ASSERT_FALSE(particles.end().base()->operator==(particle)) << "didn't point on only element in container";
    //test size
    EXPECT_EQ(1, particles.size()) << "wrong size";
}

/**
 * Test for linked cell container for addParticle, begin(), end() and size() functions as well as to calculate position and move particle to next cell and generating ghostcell and outflow boundary
 */
TEST(cellTest, LinkedCellContainer){
    LinkedCellContainer cells = LinkedCellContainer({180, 90, 1}, 3.0, {"r", "r", "r", "o", "r", "r"});
    //test cells size
    EXPECT_EQ(1800, cells.cell_numbers()) << "wrong size";
    //test addParticle
    Particle particle({0, 0, 0}, {0, 3, 0}, 50, 0);
    cells.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    //test begin
    ASSERT_TRUE(std::next(std::next(std::next(cells.begin())->begin())->begin())->begin().base()->operator==(particle)) << "didn't point on only element in container";
    //test end
    ASSERT_FALSE(std::next(std::next(std::next(cells.begin())->begin())->begin())->end().base()->operator==(particle)) << "didn't point on only element in container";
    //test calculate position and move particle to next cell
    Particle particle2({0, 3, 0}, {0, 3, 0}, 50, 0);
    PositionCalculator::PositionStoermerVerletCell(cells, 1);
    auto it = std::next(std::next(cells.begin())->begin());
    ASSERT_TRUE(std::next(std::next(it)->begin())->begin().base()->operator==(particle2)) << "particle wasn't moved to next cell";
    cells.addParticle({1, 0, 0}, {-1, 3, 0}, 50, 0);
    ASSERT_EQ(1, cells.Particles_in_cell(1,1,1));
    cells.generateGhostCell(0, 1, 1, 1);
    Particle ghost ({-1.0000000001, 0, 0}, {0, 0, 0}, 50, 1);
    ASSERT_TRUE(std::next(std::next(cells.begin()->begin())->begin())->begin().base()->operator==(ghost));
    //test outflow
    cells.addParticle({0, 0, 0}, {0, -1, 0}, 50, 0);
    EXPECT_EQ(2, cells.Particles_in_cell(1,1,1));
    PositionCalculator::PositionStoermerVerletCell(cells, 1);
    EXPECT_EQ(0, cells.Particles_in_cell(1,1,1));
}

/**
 * test for creating a cuboid and adding it to a plain container
 * test for creating a cuboid and adding it to a linked cell container
 * test for creating a sphere and adding it to a linked cell container
 */
TEST(ParticleGeneratorTest, ParticleGenerator){
    //Test of adding the right numbers of particles when generating a cuboid to plain container
    ParticleContainer particles = ParticleGenerator::createCuboid({0,0,0}, {0,0,0}, {40,8,1}, 2, 3);
    //test size
    EXPECT_EQ(320, particles.size()) << "wrong number of particles generated in plain container";
    //Test of adding the right numbers of particles when generating a cuboid to linked cell container
    LinkedCellContainer cells = LinkedCellContainer({180, 90, 1}, 3.0, {"r", "r", "r", "r", "r", "r"});
    ParticleGenerator::createCuboidInCells({0,0,0}, {0,0,0}, {40,8,1}, 2, 3, cells, 3.0);
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
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
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
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
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
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    particles.addParticle({0, 3, 0}, {4, 3, 5}, 20, 0);
    double res = (50.0*20.0*3.0)/27.0;
    std::array<double, 3> res1 = {0, res, 0};
    std::array<double, 3> res2 = {0, -res, 0};
    ForceCalculator::SimpleForceCalculation(particles);
    auto particleVector = particles.begin();
    EXPECT_EQ(res1, particleVector.base()->getF());
    EXPECT_EQ(res2, std::next(particleVector.base())->getF()) << "wrong force calculated";
}

/**
 * Test for the Lennard Jones force calculation
 */
TEST(ForceTest, LennardJonesForce){
    ParticleContainer particles;
    particles.addParticle({0, 0, 0}, {0, 3, 0}, 50, 0);
    particles.addParticle({0, 3, 0}, {4, 3, 5}, 20, 0);
    double res = (24.0/9.0)*(pow(1.0/3.0, 6.0)-2*pow(1.0/3.0, 12.0)) * 3;
    std::array<double, 3> res1 = {0, res, 0};
    std::array<double, 3> res2 = {0, -res, 0};
    ForceCalculator::LennardJonesForce(particles, 1, 1);
    auto particleVector = particles.begin();
    EXPECT_EQ(res1, particleVector.base()->getF());
    EXPECT_EQ(res2, std::next(particleVector.base())->getF()) << "wrong force calculated";
}

int main(){
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}