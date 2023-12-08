// Not copyrighted - public domain.
//
// This sample parser implementation was generated by CodeSynthesis XSD,
// an XML Schema to C++ data binding compiler. You may use it in your
// programs without any restrictions.
//

#include "newinput-pimpl.hxx"

#include <iostream>

// input_parameters_pimpl
//


//variables :

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

void input_parameters_pimpl::
pre ()
{
}

void input_parameters_pimpl::
algorithm_option (const ::std::string& algorithm_option)
{
    std::cout << "algorithm_option: " << algorithm_option << std::endl;
    inp_p_algorithm_option_ = algorithm_option;
}

std::string input_parameters_pimpl::
get_algorithm_option(){
    return inp_p_algorithm_option_;
}

void input_parameters_pimpl::
write_frequency (float write_frequency)
{
    std::cout << "write_frequency: " << write_frequency << std::endl;
    inp_p_write_frequency_ = write_frequency;
}

float input_parameters_pimpl::get_write_frequency() {
    return inp_p_write_frequency_;
}



void input_parameters_pimpl::
output_file_param (const ::std::string& output_file_param)
{
    std::cout << "output_file_param: " << output_file_param << std::endl;
    inp_p_output_file_param = output_file_param;
}

std::string input_parameters_pimpl::get_output_file_name() {
    return inp_p_output_file_param;
}

void input_parameters_pimpl::
log_level (float log_level)
{
    std::cout << "log_level: " << log_level << std::endl;
    inp_p_log_level_ = log_level;
}

float input_parameters_pimpl::get_log_level() {
    return inp_p_log_level_;
}

void input_parameters_pimpl::
simulation_input_parameters ()
{
}

void input_parameters_pimpl::
input_boundary_options ()
{
}

void input_parameters_pimpl::
cuboid_input_parameters ()
{
}

void input_parameters_pimpl::
sphere_input_parameters ()
{
}

void input_parameters_pimpl::
post_input_parameters ()
{
}

input_parameters_pimpl::~input_parameters_pimpl() {

}






// simulation_input_parameters_pimpl
//

void simulation_input_parameters_pimpl::
pre ()
{
}

void simulation_input_parameters_pimpl::
dimension ()
{
}

void simulation_input_parameters_pimpl::
avg_velocity ()
{
}

void simulation_input_parameters_pimpl::
epsilon ()
{
}

void simulation_input_parameters_pimpl::
delta_t ()
{
}

void simulation_input_parameters_pimpl::
t_end ()
{
}

void simulation_input_parameters_pimpl::
sigma ()
{
}

void simulation_input_parameters_pimpl::
r_cutoff ()
{
}

void simulation_input_parameters_pimpl::
domain_size_l ()
{
}

void simulation_input_parameters_pimpl::
post_simulation_input_parameters ()
{
}

simulation_input_parameters_pimpl::~simulation_input_parameters_pimpl() {

}

// input_boundary_options_pimpl
//

void input_boundary_options_pimpl::
pre ()
{
}

void input_boundary_options_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void input_boundary_options_pimpl::
b1 (const ::std::string& b1)
{
    std::cout << "b1: " << b1 << std::endl;
    bound_p_b1_ = b1;
}

std::string input_boundary_options_pimpl::getB1() {
    return bound_p_b1_;
}

void input_boundary_options_pimpl::
b2 (const ::std::string& b2)
{
    std::cout << "b2: " << b2 << std::endl;
    bound_p_b2_ = b2;
}

std::string input_boundary_options_pimpl::getB2() {
    return bound_p_b2_;
}

void input_boundary_options_pimpl::
b3 (const ::std::string& b3)
{
    std::cout << "b3: " << b3 << std::endl;
    bound_p_b3_ = b3;
}

std::string input_boundary_options_pimpl::getB3() {
    return bound_p_b3_;
}

void input_boundary_options_pimpl::
b4 (const ::std::string& b4)
{
    std::cout << "b4: " << b4 << std::endl;
    bound_p_b4_ = b4;
}

std::string input_boundary_options_pimpl::getB4() {
    return bound_p_b4_;
}

void input_boundary_options_pimpl::
b5 (const ::std::string& b5)
{
    std::cout << "b5: " << b5 << std::endl;
    bound_p_b5_ = b5;
}

std::string input_boundary_options_pimpl::getB5() {
    return bound_p_b5_;
}

void input_boundary_options_pimpl::
b6 (const ::std::string& b6)
{
    std::cout << "b6: " << b6 << std::endl;
    bound_p_b6_ = b6;
}

std::string input_boundary_options_pimpl::getB6() {
    return bound_p_b6_;
}

void input_boundary_options_pimpl::
post_input_boundary_options ()
{
}

input_boundary_options_pimpl::~input_boundary_options_pimpl() {

}



// cuboid_input_parameters_pimpl
//

void cuboid_input_parameters_pimpl::
pre ()
{
}

void cuboid_input_parameters_pimpl::
h ()
{
}

void cuboid_input_parameters_pimpl::
m ()
{
}

void cuboid_input_parameters_pimpl::
x1 ()
{
}

void cuboid_input_parameters_pimpl::
x2 ()
{
}

void cuboid_input_parameters_pimpl::
v1 ()
{
}

void cuboid_input_parameters_pimpl::
v2 ()
{
}

void cuboid_input_parameters_pimpl::
N1 ()
{
}

void cuboid_input_parameters_pimpl::
N2 ()
{
}

void cuboid_input_parameters_pimpl::
post_cuboid_input_parameters ()
{
}

cuboid_input_parameters_pimpl::~cuboid_input_parameters_pimpl() {

}

// sphere_input_parameters_pimpl
//

void sphere_input_parameters_pimpl::
pre ()
{
}

void sphere_input_parameters_pimpl::
h ()
{
}

void sphere_input_parameters_pimpl::
dimension ()
{
}

void sphere_input_parameters_pimpl::
m ()
{
}

void sphere_input_parameters_pimpl::
x_center ()
{
}

void sphere_input_parameters_pimpl::
v ()
{
}

void sphere_input_parameters_pimpl::
post_sphere_input_parameters ()
{
}

sphere_input_parameters_pimpl::~sphere_input_parameters_pimpl() {

}

// dimension_pimpl
//

void dimension_pimpl::
pre ()
{
}

void dimension_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void dimension_pimpl::
value (float value)
{
    std::cout << "value: " << value << std::endl;
    sim_p_dimension_ = value;
}

float dimension_pimpl::get_dimension() {
    return sim_p_dimension_;
}


void dimension_pimpl::
post_dimension ()
{
}

dimension_pimpl::~dimension_pimpl() {

}


// avg_velocity_pimpl
//

void avg_velocity_pimpl::
pre ()
{
}

void avg_velocity_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void avg_velocity_pimpl::
value (float value)
{
    std::cout << "value: " << value << std::endl;
    sim_p_avg_velocity_ = value;
}

float avg_velocity_pimpl::get_avg_velocity() {
    return sim_p_avg_velocity_;
}

void avg_velocity_pimpl::
post_avg_velocity ()
{
}

avg_velocity_pimpl::~avg_velocity_pimpl() {

}



// epsilon_pimpl
//

void epsilon_pimpl::
pre ()
{
}

void epsilon_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void epsilon_pimpl::
value (float value)
{
    std::cout << "value: " << value << std::endl;
    sim_p_epsilon_ = value;
}

float epsilon_pimpl::get_epsilon() {
    return sim_p_epsilon_;
}

void epsilon_pimpl::
post_epsilon ()
{
}

epsilon_pimpl::~epsilon_pimpl() {

}



// delta_t_pimpl
//

void delta_t_pimpl::
pre ()
{
}

void delta_t_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void delta_t_pimpl::
value (float value)
{
    std::cout << "value: " << value << std::endl;
    sim_p_delta_t_ = value;
}

float delta_t_pimpl::get_delta_t() {
    return sim_p_delta_t_;
}

void delta_t_pimpl::
post_delta_t ()
{
}

delta_t_pimpl::~delta_t_pimpl() {

}



// t_end_pimpl
//

void t_end_pimpl::
pre ()
{
}

void t_end_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void t_end_pimpl::
value (float value)
{
    std::cout << "value: " << value << std::endl;
    sim_p_t_end_ = value;
}

float t_end_pimpl::get_t_end() {
    return sim_p_t_end_;
}

void t_end_pimpl::
post_t_end ()
{
}

t_end_pimpl::~t_end_pimpl() {

}



// sigma_pimpl
//

void sigma_pimpl::
pre ()
{
}

void sigma_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void sigma_pimpl::
value (float value)
{
    std::cout << "value: " << value << std::endl;
    sim_p_sigma_ = value;
}

float sigma_pimpl::get_sigma() {
    return sim_p_sigma_;
}

void sigma_pimpl::
post_sigma ()
{
}

sigma_pimpl::~sigma_pimpl() {

}



// r_cutoff_pimpl
//

void r_cutoff_pimpl::
pre ()
{
}

void r_cutoff_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void r_cutoff_pimpl::
value (float value)
{
    std::cout << "value: " << value << std::endl;
    sim_p_r_cutoff_ = value;
}

float r_cutoff_pimpl::get_r_cutoff() {
    return sim_p_r_cutoff_;
}

void r_cutoff_pimpl::
post_r_cutoff ()
{
}

r_cutoff_pimpl::~r_cutoff_pimpl() {

}



// domain_size_l_pimpl
//

void domain_size_l_pimpl::
pre ()
{
}

void domain_size_l_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void domain_size_l_pimpl::
x (short x)
{
    std::cout << "x: " << static_cast<short> (x) << std::endl;
    sim_p_domain_size_l_x_ = x;
}

short domain_size_l_pimpl::get_domain_size_x() {
    return sim_p_domain_size_l_x_;
}

void domain_size_l_pimpl::
y (signed char y)
{
    std::cout << "y: " << static_cast<short> (y) << std::endl;
    sim_p_domain_size_l_y_ = y;
}

signed char domain_size_l_pimpl::get_domain_size_y() {
    return sim_p_domain_size_l_y_;
}

void domain_size_l_pimpl::
z (signed char z)
{
    std::cout << "z: " << static_cast<short> (z) << std::endl;
    sim_p_domain_size_l_z_ = z;
}

signed char domain_size_l_pimpl::get_domain_size_z() {
    return sim_p_domain_size_l_z_;
}

void domain_size_l_pimpl::
post_domain_size_l ()
{
}

domain_size_l_pimpl::~domain_size_l_pimpl() {

}







// h_pimpl
//

void h_pimpl::
pre ()
{
}

void h_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void h_pimpl::
value (float value)
{
    std::cout << "value: " << value << std::endl;
    cub_p_h_ = value;
}

float h_pimpl::get_h() {
    return cub_p_h_;
}

void h_pimpl::
post_h ()
{
}

h_pimpl::~h_pimpl() {

}



// m_pimpl
//

void m_pimpl::
pre ()
{
}

void m_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void m_pimpl::
value (float value)
{
    std::cout << "value: " << value << std::endl;
    cub_p_m_ = value;
}

float m_pimpl::get_m() {
    return cub_p_m_;
}

void m_pimpl::
post_m ()
{
}

m_pimpl::~m_pimpl() {

}



// x1_pimpl
//

void x1_pimpl::
pre ()
{
}

void x1_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void x1_pimpl::
x (signed char x)
{
    std::cout << "x: " << static_cast<short> (x) << std::endl;
    cub_p_x1_x_ = x;
}

signed char x1_pimpl::get_x1_x() {
    return cub_p_x1_x_;
}

signed char x1_pimpl::get_x1_y() {
    return cub_p_x1_y_;
}

signed char x1_pimpl::get_x1_z() {
    return cub_p_x1_z_;
}

void x1_pimpl::
y (signed char y)
{
    std::cout << "y: " << static_cast<short> (y) << std::endl;
    cub_p_x1_y_ = y;
}

void x1_pimpl::
z (signed char z)
{
    std::cout << "z: " << static_cast<short> (z) << std::endl;
    cub_p_x1_z_ = z;
}

void x1_pimpl::
post_x1 ()
{
}

x1_pimpl::~x1_pimpl() {

}



// x2_pimpl
//

void x2_pimpl::
pre ()
{
}

void x2_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void x2_pimpl::
x (signed char x)
{
    std::cout << "x: " << static_cast<short> (x) << std::endl;
    cub_p_x2_x_ = x;
}

signed char x2_pimpl::get_x2_x() {
    return cub_p_x2_x_;
}

void x2_pimpl::
y (signed char y)
{
    std::cout << "y: " << static_cast<short> (y) << std::endl;
    cub_p_x2_y_ = y;
}

signed char x2_pimpl::get_x2_y() {
    return cub_p_x2_y_;
}

void x2_pimpl::
z (signed char z)
{
    std::cout << "z: " << static_cast<short> (z) << std::endl;
    cub_p_x2_z_ = z;
}

signed char x2_pimpl::get_x2_z() {
    return cub_p_x2_z_;
}

void x2_pimpl::
post_x2 ()
{
}




// v1_pimpl
//

void v1_pimpl::
pre ()
{
}

void v1_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void v1_pimpl::
x (signed char x)
{
    std::cout << "x: " << static_cast<short> (x) << std::endl;
    cub_p_v1_x_ = x;
}

signed char v1_pimpl::get_v1_x() {
    return cub_p_v1_x_;
}

void v1_pimpl::
y (signed char y)
{
    std::cout << "y: " << static_cast<short> (y) << std::endl;
    cub_p_v1_y_ = y;
}

signed char v1_pimpl::get_v1_y() {
    return cub_p_v1_y_;
}

void v1_pimpl::
z (signed char z)
{
    std::cout << "z: " << static_cast<short> (z) << std::endl;
    cub_p_v1_z_ = z;
}

signed char v1_pimpl::get_v1_z() {
    return cub_p_v1_z_;
}

void v1_pimpl::
post_v1 ()
{
}


// v2_pimpl
//

void v2_pimpl::
pre ()
{
}

void v2_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void v2_pimpl::
x (signed char x)
{
    std::cout << "x: " << static_cast<short> (x) << std::endl;
    cub_p_v2_x_ = x;
}

signed char v2_pimpl::get_v2_x() {
    return cub_p_v2_x_;
}

void v2_pimpl::
y (signed char y)
{
    std::cout << "y: " << static_cast<short> (y) << std::endl;
    cub_p_v2_y_= y;
}

signed char v2_pimpl::get_v2_y() {
    return cub_p_v2_y_;
}

void v2_pimpl::
z (signed char z)
{
    std::cout << "z: " << static_cast<short> (z) << std::endl;
    cub_p_v2_z_ = z;
}

signed char v2_pimpl::get_v2_z() {
    return cub_p_v2_z_;
}

void v2_pimpl::
post_v2 ()
{
}

// N1_pimpl
//

void N1_pimpl::
pre ()
{
}

void N1_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void N1_pimpl::
x (signed char x)
{
    std::cout << "x: " << static_cast<short> (x) << std::endl;
    cub_p_N1_x_ = x;
}
signed char N1_pimpl::get_n1_x() {
    return cub_p_N1_x_;
}

void N1_pimpl::
y (signed char y)
{
    std::cout << "y: " << static_cast<short> (y) << std::endl;
    cub_p_N1_y_ = y;
}

signed char N1_pimpl::get_n1_y() {
    return cub_p_N1_y_;
}

void N1_pimpl::
z (signed char z)
{
    std::cout << "z: " << static_cast<short> (z) << std::endl;
    cub_p_N1_z_ = z;
}

signed char N1_pimpl::get_n1_z() {
    return cub_p_N1_z_;
}

void N1_pimpl::
post_N1 ()
{
}


// N2_pimpl
//

void N2_pimpl::
pre ()
{
}

void N2_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void N2_pimpl::
x (signed char x)
{
    std::cout << "x: " << static_cast<short> (x) << std::endl;
    cub_p_N2_x_ = x;
}

signed char N2_pimpl::get_n2_x() {
    return cub_p_N2_x_;
}

void N2_pimpl::
y (signed char y)
{
    std::cout << "y: " << static_cast<short> (y) << std::endl;
    cub_p_N2_y_ = y;
}

signed char N2_pimpl::get_n2_y() {
    return cub_p_N2_y_;
}

void N2_pimpl::
z (signed char z)
{
    std::cout << "z: " << static_cast<short> (z) << std::endl;
    cub_p_N2_z_ = z;
}

signed char N2_pimpl::get_n3_z() {
    return cub_p_N2_z_;
}

void N2_pimpl::
post_N2 ()
{
}







// x_center_pimpl
//

void x_center_pimpl::
pre ()
{
}

void x_center_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void x_center_pimpl::
x (signed char x)
{
    std::cout << "x: " << static_cast<short> (x) << std::endl;
    sph_p_x_center_x_ = x;
}

signed char x_center_pimpl::get_x_center_x() {
    return sph_p_x_center_x_;
}

void x_center_pimpl::
y (signed char y)
{
    std::cout << "y: " << static_cast<short> (y) << std::endl;
    sph_p_x_center_y_ = y;
}

signed char x_center_pimpl::get_x_center_y() {
    return sph_p_x_center_y_;
}


void x_center_pimpl::
z (signed char z)
{
    std::cout << "z: " << static_cast<short> (z) << std::endl;
    sph_p_x_center_z_ = z;
}

signed char x_center_pimpl::get_x_center_z() {
    return sph_p_x_center_z_;
}


void x_center_pimpl::
post_x_center ()
{
}


// v_pimpl
//

void v_pimpl::
pre ()
{
}

void v_pimpl::
name (const ::std::string& name)
{
    std::cout << "name: " << name << std::endl;
}

void v_pimpl::
x (signed char x)
{
    std::cout << "x: " << static_cast<short> (x) << std::endl;
    sph_p_v_x_ = x;
}

signed char v_pimpl::get_v_x() {
    return sph_p_v_x_;
}

void v_pimpl::
y (signed char y)
{
    std::cout << "y: " << static_cast<short> (y) << std::endl;
    sph_p_v_y_ = y;
}

signed char v_pimpl::get_v_y() {
    return sph_p_v_y_;
}

void v_pimpl::
z (signed char z)
{
    std::cout << "z: " << static_cast<short> (z) << std::endl;
    sph_p_v_z_ = z;
}

signed char v_pimpl::get_v_z() {
    return sph_p_v_z_;
}

void v_pimpl::
post_v ()
{
}







