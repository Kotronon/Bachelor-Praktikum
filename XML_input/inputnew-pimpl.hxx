// Not copyrighted - public domain.
//
// This sample parser implementation was generated by CodeSynthesis XSD,
// an XML Schema to C++ data binding compiler. You may use it in your
// programs without any restrictions.
//

#ifndef INPUTNEW_PIMPL_HXX
#define INPUTNEW_PIMPL_HXX

#include "inputnew-pskel.hxx"

class parameters_pimpl: public virtual parameters_pskel
{
public:

    ~parameters_pimpl();
    virtual void
    pre ();

    virtual void
    algorithm_option (const ::std::string&);

    virtual void
    write_frequency (unsigned char);

    virtual void
    output_file_name (const ::std::string&);

    virtual void
    log_level (unsigned char);

    virtual void
    simulation_parameters ();

    virtual void
    boundaries ();

    virtual void
    cuboid_parameters ();

    virtual void
    sphere_parameters ();

    virtual void
    post_parameters ();
};

class simulation_parameters_pimpl: public virtual simulation_parameters_pskel
{
public:
    ~simulation_parameters_pimpl();
    virtual void
    pre ();

    virtual void
    name ();

    virtual void
    post_simulation_parameters ();
};

class boundaries_pimpl: public virtual boundaries_pskel
{
public:
    ~boundaries_pimpl();
    virtual void
    pre ();

    virtual void
    b1 (const ::std::string&);

    virtual void
    b2 (const ::std::string&);

    virtual void
    b3 (const ::std::string&);

    virtual void
    b4 (const ::std::string&);

    virtual void
    b5 (const ::std::string&);

    virtual void
    b6 (const ::std::string&);

    virtual void
    post_boundaries ();
};

class cuboid_parameters_pimpl: public virtual cuboid_parameters_pskel
{
public:
    ~cuboid_parameters_pimpl();
    virtual void
    pre ();

    virtual void
    name ();

    virtual void
    post_cuboid_parameters ();
};

class sphere_parameters_pimpl: public virtual sphere_parameters_pskel
{
public:
    ~sphere_parameters_pimpl ();
    virtual void
    pre ();

    virtual void
    name ();

    virtual void
    post_sphere_parameters ();
};

class name_pimpl: public virtual name_pskel
{
public:
public:
    ~name_pimpl();

    virtual void
    pre ();

    virtual void
    x (unsigned char);

    virtual void
    y (char);

    virtual void
    z (unsigned char);

    virtual void
    value (double);

    virtual void
    post_name ();
};

#endif // INPUTNEW_PIMPL_HXX
