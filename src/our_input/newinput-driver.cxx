// Not copyrighted - public domain.
//
// This sample parser implementation was generated by CodeSynthesis XSD,
// an XML Schema to C++ data binding compiler. You may use it in your
// programs without any restrictions.
//

#include "newinput-pimpl.hxx"
#include "newinput-pimpl.cxx"

#include <iostream>

int
main (int argc, char* argv[])
{
  if (argc != 2)
  {
    std::cerr << "usage: " << argv[0] << " file.xml" << std::endl;
    return 1;
  }

  try
  {
    // Instantiate individual parsers.
    //
    ::input_parameters_pimpl input_parameters_p;
    ::xml_schema::string_pimpl string_p;
    ::xml_schema::float_pimpl float_p;
    ::simulation_input_parameters_pimpl simulation_input_parameters_p;
    ::dimension_pimpl dimension_p;
    ::avg_velocity_pimpl avg_velocity_p;
    ::epsilon_pimpl epsilon_p;
    ::delta_t_pimpl delta_t_p;
    ::t_end_pimpl t_end_p;
    ::sigma_pimpl sigma_p;
    ::r_cutoff_pimpl r_cutoff_p;
    ::domain_size_l_pimpl domain_size_l_p;
    ::xml_schema::byte_pimpl byte_p;
    ::xml_schema::short_pimpl short_p;
    ::input_boundary_options_pimpl input_boundary_options_p;
    ::cuboid_input_parameters_pimpl cuboid_input_parameters_p;
    ::h_pimpl h_p;
    ::m_pimpl m_p;
    ::x1_pimpl x1_p;
    ::x2_pimpl x2_p;
    ::v1_pimpl v1_p;
    ::v2_pimpl v2_p;
    ::N1_pimpl N1_p;
    ::N2_pimpl N2_p;
    ::sphere_input_parameters_pimpl sphere_input_parameters_p;
    ::x_center_pimpl x_center_p;
    ::v_pimpl v_p;

    // Connect the parsers together.
    //
    input_parameters_p.parsers (string_p,
                                float_p,
                                string_p,
                                float_p,
                                simulation_input_parameters_p,
                                input_boundary_options_p,
                                cuboid_input_parameters_p,
                                sphere_input_parameters_p);

    simulation_input_parameters_p.parsers (dimension_p,
                                           avg_velocity_p,
                                           epsilon_p,
                                           delta_t_p,
                                           t_end_p,
                                           sigma_p,
                                           r_cutoff_p,
                                           domain_size_l_p);

    dimension_p.parsers (string_p,
                         float_p);

    avg_velocity_p.parsers (string_p,
                            float_p);

    epsilon_p.parsers (string_p,
                       float_p);

    delta_t_p.parsers (string_p,
                       float_p);

    t_end_p.parsers (string_p,
                     float_p);

    sigma_p.parsers (string_p,
                     float_p);

    r_cutoff_p.parsers (string_p,
                        float_p);

    domain_size_l_p.parsers (string_p,
                             short_p,
                             byte_p,
                             byte_p);

    input_boundary_options_p.parsers (string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p,
                                      string_p);

    cuboid_input_parameters_p.parsers (h_p,
                                       m_p,
                                       x1_p,
                                       x2_p,
                                       v1_p,
                                       v2_p,
                                       N1_p,
                                       N2_p);

    h_p.parsers (string_p,
                 float_p);

    m_p.parsers (string_p,
                 float_p);

    x1_p.parsers (string_p,
                  byte_p,
                  byte_p,
                  byte_p);

    x2_p.parsers (string_p,
                  byte_p,
                  byte_p,
                  byte_p);

    v1_p.parsers (string_p,
                  byte_p,
                  byte_p,
                  byte_p);

    v2_p.parsers (string_p,
                  byte_p,
                  byte_p,
                  byte_p);

    N1_p.parsers (string_p,
                  byte_p,
                  byte_p,
                  byte_p);

    N2_p.parsers (string_p,
                  byte_p,
                  byte_p,
                  byte_p);

    sphere_input_parameters_p.parsers (h_p,
                                       dimension_p,
                                       m_p,
                                       x_center_p,
                                       v_p);

    x_center_p.parsers (string_p,
                        byte_p,
                        byte_p,
                        byte_p);

    v_p.parsers (string_p,
                 byte_p,
                 byte_p,
                 byte_p);

    // Parse the XML document.
    //
    ::xml_schema::document doc_p (input_parameters_p, "input_parameters");

    input_parameters_p.pre ();
    doc_p.parse (argv[1]);
    input_parameters_p.post_input_parameters ();
  }
  catch (const ::xml_schema::exception& e)
  {
    std::cerr << e << std::endl;
    return 1;
  }
  catch (const std::ios_base::failure&)
  {
    std::cerr << argv[1] << ": error: io failure" << std::endl;
    return 1;
  }
}

