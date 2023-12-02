// Not copyrighted - public domain.
//
// This sample parser implementation was generated by CodeSynthesis XSD,
// an XML Schema to C++ data binding compiler. You may use it in your
// programs without any restrictions.
//

#ifndef DEBUG
#  define DEBUG 1
#endif


#include "inputnew-pimpl.hxx"
#include "inputnew-pimpl.cxx"



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
    ::parameters_pimpl parameters_p;
    ::xml_schema::string_pimpl string_p;
    ::xml_schema::unsigned_byte_pimpl unsigned_byte_p;
    ::simulation_parameters_pimpl simulation_parameters_p;
    ::name_pimpl name_p;
    ::xml_schema::decimal_pimpl decimal_p;
    ::boundaries_pimpl boundaries_p;
    ::cuboid_parameters_pimpl cuboid_parameters_p;
    ::sphere_parameters_pimpl sphere_parameters_p;

    // Connect the parsers together.
    //
    parameters_p.parsers (string_p,
                          unsigned_byte_p,
                          string_p,
                          unsigned_byte_p,
                          simulation_parameters_p,
                          boundaries_p,
                          cuboid_parameters_p,
                          sphere_parameters_p);

    simulation_parameters_p.parsers (name_p);

    name_p.parsers (unsigned_byte_p,
                    unsigned_byte_p,
                    unsigned_byte_p,
                    decimal_p);

    boundaries_p.parsers (string_p,
                          string_p,
                          string_p,
                          string_p,
                          string_p,
                          string_p);

    cuboid_parameters_p.parsers (name_p);

    sphere_parameters_p.parsers (name_p);

    // Parse the XML document.
    //
    ::xml_schema::document doc_p (parameters_p, "parameters");

    parameters_p.pre ();
    doc_p.parse (argv[1]);
    parameters_p.post_parameters ();
  }
  catch (const ::xml_schema::exception& e)
  {
    std::cerr << e << std::endl;
    return 1;
  }
  catch (const std::ios_base::failure&)
  {
     throw xml_schema::parsing();
    std::cerr << argv[1] << ": error: io failure" << std::endl;
    return 1;
  }
}

