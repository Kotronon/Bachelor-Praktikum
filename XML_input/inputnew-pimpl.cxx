// Not copyrighted - public domain.
//
// This sample parser implementation was generated by CodeSynthesis XSD,
// an XML Schema to C++ data binding compiler. You may use it in your
// programs without any restrictions.
//

#include "inputnew-pimpl.hxx"

#include <iostream>

// parameters_pimpl
//

void parameters_pimpl::
pre ()
{
}

void parameters_pimpl::
algorithm_option (const ::std::string& algorithm_option)
{
  std::cout << "algorithm_option: " << algorithm_option << std::endl;
}

void parameters_pimpl::
write_frequency (unsigned char write_frequency)
{
  std::cout << "write_frequency: " << static_cast<unsigned short> (write_frequency) << std::endl;
}

void parameters_pimpl::
output_file_name (const ::std::string& output_file_name)
{
  std::cout << "output_file_name: " << output_file_name << std::endl;
}

void parameters_pimpl::
log_level (unsigned char log_level)
{
  std::cout << "log_level: " << static_cast<unsigned short> (log_level) << std::endl;
}

void parameters_pimpl::
simulation_parameters ()
{
}

void parameters_pimpl::
boundaries ()
{
}

void parameters_pimpl::
cuboid_parameters ()
{
}

void parameters_pimpl::
sphere_parameters ()
{
}

void parameters_pimpl::
post_parameters ()
{
}

parameters_pimpl::~parameters_pimpl() {

}

// simulation_parameters_pimpl
//

void simulation_parameters_pimpl::
pre ()
{
}

void simulation_parameters_pimpl::
name ()
{
}

void simulation_parameters_pimpl::
post_simulation_parameters ()
{
}

simulation_parameters_pimpl::~simulation_parameters_pimpl() {

}

// boundaries_pimpl
//

void boundaries_pimpl::
pre ()
{
}

void boundaries_pimpl::
b1 (const ::std::string& b1)
{
  std::cout << "b1: " << b1 << std::endl;
}

void boundaries_pimpl::
b2 (const ::std::string& b2)
{
  std::cout << "b2: " << b2 << std::endl;
}

void boundaries_pimpl::
b3 (const ::std::string& b3)
{
  std::cout << "b3: " << b3 << std::endl;
}

void boundaries_pimpl::
b4 (const ::std::string& b4)
{
  std::cout << "b4: " << b4 << std::endl;
}

void boundaries_pimpl::
b5 (const ::std::string& b5)
{
  std::cout << "b5: " << b5 << std::endl;
}

void boundaries_pimpl::
b6 (const ::std::string& b6)
{
  std::cout << "b6: " << b6 << std::endl;
}

void boundaries_pimpl::
post_boundaries ()
{
}

boundaries_pimpl::~boundaries_pimpl() {

}

// cuboid_parameters_pimpl
//

void cuboid_parameters_pimpl::
pre ()
{
}

void cuboid_parameters_pimpl::
name ()
{
}

void cuboid_parameters_pimpl::
post_cuboid_parameters ()
{
}

cuboid_parameters_pimpl::~cuboid_parameters_pimpl() {

}

// sphere_parameters_pimpl
//

void sphere_parameters_pimpl::
pre ()
{
}

void sphere_parameters_pimpl::
name ()
{
}

void sphere_parameters_pimpl::
post_sphere_parameters ()
{
}

sphere_parameters_pimpl::~sphere_parameters_pimpl() {

}

// name_pimpl
//

void name_pimpl::
pre ()
{
}

void name_pimpl::
x (unsigned char x)
{
  std::cout << "x: " << static_cast<unsigned short> (x) << std::endl;
}

void name_pimpl::
y (unsigned char y)
{
  std::cout << "y: " << static_cast<unsigned short> (y) << std::endl;
}

void name_pimpl::
z (unsigned char z)
{
  std::cout << "z: " << static_cast<unsigned short> (z) << std::endl;
}

void name_pimpl::
value (double value)
{
  std::cout << "value: " << value << std::endl;
}

void name_pimpl::
post_name ()
{
}

name_pimpl::~name_pimpl() {

}

