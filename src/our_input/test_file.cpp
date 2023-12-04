//
// Created by Samsamsamsamsamsamsamsamsa on 03.12.2023.
//

#include "test_file.h"
#include "newinput-pimpl.hxx"
#include "newinput-pskel.hxx"
#include "newinput-pimpl.cxx"
#include <string>
#include <iostream>

//alle parameter!
std::string inp_p_algorithm_option_, inp_p_output_file_param,bound_p_b1_,bound_p_b2_,bound_p_b3_,bound_p_b4_,bound_p_b5_,bound_p_b6_;
float inp_p_write_frequency_,inp_p_log_level_,
        sim_p_dimension_,
        sim_p_avg_velocity_,sim_p_epsilon_,sim_p_delta_t_,sim_p_t_end_,sim_p_sigma_,
        sim_p_r_cutoff_,
        cub_p_h_, cub_p_m_,
        sph_p_h, sph_p_dimension, sph_p_m;
short sim_p_domain_size_l_x_;

signed char sim_p_domain_size_l_y_, sim_p_domain_size_l_z_,
        cub_p_x1_x_,cub_p_x1_y_,cub_p_x1_z_,
        cub_p_x2_x_,cub_p_x2_y_,cub_p_x2_z_,
        cub_p_v1_x_,cub_p_v1_y_,cub_p_v1_z_,
        cub_p_v2_x_,cub_p_v2_y_,cub_p_v2_z_,
        cub_p_N1_x_,cub_p_N1_y_,cub_p_N1_z_,
        cub_p_N2_x_,cub_p_N2_y_,cub_p_N2_z_,
        sph_p_x_center_x_,sph_p_x_center_y_,sph_p_x_center_z_,
        sph_p_v_x_,sph_p_v_y_,sph_p_v_z_
;

void setup(){
    //TODO: copy this to void collect_input();
    input_parameters_pimpl inputParametersPimpl;
    inp_p_algorithm_option_ = inputParametersPimpl.getAlgorithmOption();
    inp_p_log_level_ = inputParametersPimpl.getLogLevel();
    inp_p_output_file_param = inputParametersPimpl.getOutputFileName();
    inp_p_write_frequency_ = inputParametersPimpl.getWriteFrequency();

    //simulation params
    simulation_input_parameters_pimpl simulationInputParametersPimpl;
    sim_p_avg_velocity_ = simulationInputParametersPimpl.getAvgVelocity();
    sim_p_delta_t_ = simulationInputParametersPimpl.getDeltaT();
    sim_p_dimension_ = simulationInputParametersPimpl.getDimension();
    sim_p_epsilon_ = simulationInputParametersPimpl.getEpsilon();
    sim_p_t_end_ = simulationInputParametersPimpl.getTEnd();
    sim_p_sigma_ = simulationInputParametersPimpl.getSigma();
    sim_p_r_cutoff_ = simulationInputParametersPimpl.getRCutoff();
    sim_p_domain_size_l_x_ = simulationInputParametersPimpl.getDomainSizeLX();
    sim_p_domain_size_l_y_ =simulationInputParametersPimpl.getDomainSizeLY();
    sim_p_domain_size_l_z_ = simulationInputParametersPimpl.getDomainSizeLZ();
    sim_p_dimension_ = simulationInputParametersPimpl.getDimension();
    sim_p_delta_t_ = simulationInputParametersPimpl.getDeltaT();

    //innput boundaries
    input_boundary_options_pimpl inputBoundaryOptionsPimpl;
    bound_p_b1_ = inputBoundaryOptionsPimpl.getB1();
    bound_p_b2_ = inputBoundaryOptionsPimpl.getB2();
    bound_p_b3_ = inputBoundaryOptionsPimpl.getB3();
    bound_p_b4_ = inputBoundaryOptionsPimpl.getB4();
    bound_p_b5_ = inputBoundaryOptionsPimpl.getB5();
    bound_p_b6_ = inputBoundaryOptionsPimpl.getB6();

    //cuboid parameters
    cuboid_input_parameters_pimpl cuboidInputParametersPimpl;

    cub_p_N1_x_ = cuboidInputParametersPimpl.getN1X();
    cub_p_N1_y_ = cuboidInputParametersPimpl.getN1Y();
    cub_p_N1_z_ = cuboidInputParametersPimpl.getN1Z();

    cub_p_N2_x_ = cuboidInputParametersPimpl.getN2X();
    cub_p_N2_x_ = cuboidInputParametersPimpl.getN2Y();
    cub_p_N2_x_ = cuboidInputParametersPimpl.getN2Z();

    cub_p_v1_x_ = cuboidInputParametersPimpl.getV1X();
    cub_p_v1_y_ = cuboidInputParametersPimpl.getV1Y();
    cub_p_v1_z_ = cuboidInputParametersPimpl.getV1Z();

    cub_p_v2_x_ = cuboidInputParametersPimpl.getV2X();
    cub_p_v2_y_ = cuboidInputParametersPimpl.getV2Y();
    cub_p_v2_z_ = cuboidInputParametersPimpl.getV2Z();

    cub_p_x1_x_ = cuboidInputParametersPimpl.getX1X();
    cub_p_x1_y_ = cuboidInputParametersPimpl.getX1Y();
    cub_p_x1_z_ = cuboidInputParametersPimpl.getX1Z();

    cub_p_x2_x_ = cuboidInputParametersPimpl.getX2X();
    cub_p_x2_y_ = cuboidInputParametersPimpl.getX2Y();
    cub_p_x2_z_ = cuboidInputParametersPimpl.getX2Z();

    cub_p_h_ = cuboidInputParametersPimpl.getH();
    cub_p_m_ = cuboidInputParametersPimpl.getM();



    //sphere parameters
    sphere_input_parameters_pimpl sphereInputParametersPimpl;

    sph_p_v_x_ = sphereInputParametersPimpl.getVX();
    sph_p_v_y_ = sphereInputParametersPimpl.getVY();
    sph_p_v_z_ = sphereInputParametersPimpl.getVZ();

    sph_p_x_center_x_ = sphereInputParametersPimpl.getXCenterX();
    sph_p_x_center_y_ = sphereInputParametersPimpl.getXCenterY();
    sph_p_x_center_z_ = sphereInputParametersPimpl.getXCenterZ();

    sph_p_dimension = sphereInputParametersPimpl.getDimension();
    sph_p_h = sphereInputParametersPimpl.getH();
    sph_p_m = sphereInputParametersPimpl.getM();




}


int main(int argc, char *argsv[]){
setup();

//copy and paste this into main and pray for the best

/*spdlog::info("Hello from MolSim for PSE!");


    ParticleContainer cuboid_1 = ParticleGenerator::createCuboid(cub_p_,v_1,N_1,h,m);
    ParticleContainer cuboid_2 = ParticleGenerator::createCuboid(x_2,v_2,N_2,h,m);
    container.addParticleContainer(cuboid_1);
    container.addParticleContainer(cuboid_2);

    LinkedCellContainer cells = LinkedCellContainer({sim_p_domain_size_l_x_, sim_p_domain_size_l_y_, sim_p_domain_size_l_z_}, sim_p_r_cutoff_, {bound_p_b1_,bound_p_b2_,bound_p_b3_,bound_p_b4_,bound_p_b5_,bound_p_b6_}); //boundary left, right, up, down, behind, bevor
    ParticleGenerator::createCuboidInCells({cub_p_x1_x_,cub_p_x1_y_,cub_p_x1_z_}, {cub_p_v1_x_,cub_p_v1_y_,cub_p_v1_z_}, {cub_p_N1_x_,cub_p_N1_y_,cub_p_N1_z_}, cub_p_h_, cub_p_m_, cells, sim_p_r_cutoff_);
    ParticleGenerator::createCuboidInCells({70, 60, 0}, {0,-1,0}, {20,20,1}, 1.1225, 1, cells, sim_p_r_cutoff_);
    double end_time = sim_p_t_end_;
    double delta_t = sim_p_delta_t_;
    double current_time = start_time;
    int iteration = inp_p_write_frequency_;
    //Pre-calculation of f
    //ForceCalculator::LennardJonesForceFaster(container, sim_p_epsilon_, sim_p_sigma_);
    ForceCalculator::LennardJonesForceCell(cells, sim_p_epsilon_, sim_p_sigma_);
    //Initialization with Brownian Motion
    //VelocityCalculator::BrownianMotionInitialization(container, sim_p_avg_velocity_, sim_p_dimension_);
    VelocityCalculator::BrownianMotionInitializationCell(cells, sim_p_avg_velocity_, sim_p_dimension_);
    //For this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
        if (iteration == 30)
            spdlog::info("get in losers, we're going shopping");
        //Calculate new x
        //PositionCalculator::PositionStoermerVerlet(container, delta_t);
        PositionCalculator::PositionStoermerVerletCell(cells, delta_t);
        //Calculate new f
        //ForceCalculator::LennardJonesForceFaster(container, eps, sig);
        ForceCalculator::LennardJonesForceCell(cells, eps, sig);
        //Calculate new v
        //VelocityCalculator::VelocityStoermerVerlet(container, delta_t);
        VelocityCalculator::VelocityStoermerVerletCell(cells, delta_t);

        iteration++;
        if (iteration % 10 == 0) {
            //plotParticles(iteration);
            plotParticlesInCells(iteration, cells);
        }
        if (iteration % 100 == 0) {
            spdlog::info("Iteration " + std::to_string(iteration) + " finished.");
        }

        current_time += delta_t;
    }

    spdlog::info("Output written. Terminating..." );
    return 0;*/


}
