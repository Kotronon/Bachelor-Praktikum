//
// Created by Admin on 11.12.2023.
//

#include "implementation_input.h"
#include "newinput-pimpl.hxx"

dimension_pimpl Dimension = {dimension_pimpl::get_dimension()};
avg_velocity_pimpl avgVelocity = {avg_velocity_pimpl::get_avg_velocity()};
r_cutoff_pimpl rCutoff = {r_cutoff_pimpl::get_r_cutoff()};
epsilon_pimpl epsilon = {epsilon_pimpl :: get_epsilon()};
delta_t_pimpl deltaT_ = {delta_t_pimpl::get_delta_t()};
t_end_pimpl tEndPimpl = {t_end_pimpl::get_t_end()};
sigma_pimpl sigma = {sigma_pimpl::get_sigma()};
domain_size_l_pimpl domainSize = {domain_size_l_pimpl::get_domain_size_x(),domain_size_l_pimpl::get_domain_size_y(),domain_size_l_pimpl::get_domain_size_z()};

h_pimpl hParams = {h_pimpl::get_h()};
m_pimpl mParams = {m_pimpl::get_m()};

x1_pimpl x1Params = {x1_pimpl::get_x1_x(), x1_pimpl::get_x1_y(), x1_pimpl::get_x1_z()};
x2_pimpl x2Params = {x2_pimpl::get_x2_x(), x2_pimpl::get_x2_y(), x2_pimpl::get_x2_z()};
v1_pimpl v1Params = {v1_pimpl::get_v1_x(), v1_pimpl::get_v1_y(), v1_pimpl::get_v1_z()};
v2_pimpl v2Params = {v2_pimpl::get_v2_x(), v2_pimpl::get_v2_y(), v2_pimpl::get_v2_z()};
N1_pimpl n1Params = {N1_pimpl::get_n1_x(), N1_pimpl::get_n1_y(), N1_pimpl::get_n1_z()};
N2_pimpl n2Params = {N2_pimpl::get_n2_x(), N2_pimpl::get_n2_y(), N2_pimpl::get_n2_z()};

h_pimpl sphHparams =  {h_pimpl::get_h()};
m_pimpl sphMparams = {m_pimpl::get_m()};
dimension_pimpl sphDim = {dimension_pimpl::get_dimension()};
x_center_pimpl xCenterParams ={x_center_pimpl::get_x_center_x(), x_center_pimpl::get_x_center_y(), x_center_pimpl::get_x_center_z()};
v_pimpl vParams = {v_pimpl::get_v_x(), v_pimpl::get_v_y(), v_pimpl::get_v_z()};

input_boundary_options_pimpl boundaryParams ={input_boundary_options_pimpl::getB1(),
                                              input_boundary_options_pimpl::getB2(),input_boundary_options_pimpl::getB3(),
                                              input_boundary_options_pimpl::getB4(),input_boundary_options_pimpl::getB5(),
                                              input_boundary_options_pimpl::getB6()};

simulation_input_parameters_pimpl simulationParams =
        {Dimension,
         avgVelocity,
         epsilon,
         deltaT_,
         tEndPimpl,
         sigma,
         rCutoff,
         domainSize};


cuboid_input_parameters_pimpl cuboidParams = {hParams,mParams,x1Params,x2Params,
                                              v1Params,v2Params,n1Params,n2Params};
sphere_input_parameters_pimpl sphereParams = {sphHparams,sphDim, sphMparams,xCenterParams,vParams};

input_parameters_pimpl input = {input_parameters_pimpl::get_algorithm_option(),
                                input_parameters_pimpl::get_write_frequency(),input_parameters_pimpl::get_output_file_name(), input_parameters_pimpl::get_log_level(),
                                simulationParams,boundaryParams,cuboidParams,sphereParams};



input_parameters_pimpl::input_parameters_pimpl(std::string algo, float write_f, std::string outputfilename,
                                               float loglevel, simulation_input_parameters_pskel sim,
                                               input_boundary_options_pskel boundaries,
                                               cuboid_input_parameters_pskel cuboid,
                                               sphere_input_parameters_pskel sphere) {

}



simulation_input_parameters_pimpl::simulation_input_parameters_pimpl(dimension_pskel dim, avg_velocity_pskel avgV,
                                                                     epsilon_pskel eps, delta_t_pskel deltaT,
                                                                     t_end_pskel tEnd, sigma_pskel sigma,
                                                                     r_cutoff_pskel rCut,
                                                                     domain_size_l_pskel domainSizeL) {

}

input_boundary_options_pimpl::input_boundary_options_pimpl(std::string b1, std::string b2, std::string b3,
                                                           std::string b4, std::string b5, std::string b6) {

}

cuboid_input_parameters_pimpl::cuboid_input_parameters_pimpl(h_pskel, m_pskel, x1_pskel, x2_pskel, v1_pskel, v2_pskel,
                                                             N1_pskel, N2_pskel) {

}

sphere_input_parameters_pimpl::sphere_input_parameters_pimpl(h_pskel, dimension_pskel, m_pskel, x_center_pskel,
                                                             v_pskel) {

}

dimension_pimpl::dimension_pimpl(float value) {

}

avg_velocity_pimpl::avg_velocity_pimpl(float value) {

}

epsilon_pimpl::epsilon_pimpl(float value) {

}

delta_t_pimpl::delta_t_pimpl(float value) {

}

t_end_pimpl::t_end_pimpl(float value) {

}

sigma_pimpl::sigma_pimpl(float value) {

}

r_cutoff_pimpl::r_cutoff_pimpl(float value) {

}

domain_size_l_pimpl::domain_size_l_pimpl(short x, signed char y, signed char z) {

}


h_pimpl::h_pimpl(float value) {

}

m_pimpl::m_pimpl(float value) {

}

x1_pimpl::x1_pimpl(signed char x, signed char y, signed char z) {

}

x2_pimpl::x2_pimpl(signed char x, signed char y, signed char z) {

}

v1_pimpl::v1_pimpl(signed char x, signed char y, signed char z) {

}