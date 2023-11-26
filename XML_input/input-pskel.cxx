// Copyright (c) 2005-2023 Code Synthesis.
//
// This program was generated by CodeSynthesis XSD, an XML Schema to
// C++ data binding compiler.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// In addition, as a special exception, Code Synthesis gives permission
// to link this program with the Xerces-C++ library (or with modified
// versions of Xerces-C++ that use the same license as Xerces-C++), and
// distribute linked combinations including the two. You must obey the GNU
// General Public License version 2 in all respects for all of the code
// used other than Xerces-C++. If you modify this copy of the program, you
// may extend this exception to your version of the program, but you are
// not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// Furthermore, Code Synthesis makes a special exception for the Free/Libre
// and Open Source Software (FLOSS) which is described in the accompanying
// FLOSSE file.
//

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/pre.hxx>

#include "input-pskel.hxx"

// input_pskel
//

void input_pskel::
vector_parameters_parser (::vector_parameters_pskel& p)
{
  this->vector_parameters_parser_ = &p;
}

void input_pskel::
value_parameters_parser (::value_parameters_pskel& p)
{
  this->value_parameters_parser_ = &p;
}

void input_pskel::
parsers (::vector_parameters_pskel& vector_parameters,
         ::value_parameters_pskel& value_parameters)
{
  this->vector_parameters_parser_ = &vector_parameters;
  this->value_parameters_parser_ = &value_parameters;
}

input_pskel::
input_pskel ()
: vector_parameters_parser_ (0),
  value_parameters_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// vector_parameters_pskel
//

void vector_parameters_pskel::
param_parser (::param_pskel& p)
{
  this->param_parser_ = &p;
}

void vector_parameters_pskel::
parsers (::param_pskel& param)
{
  this->param_parser_ = &param;
}

vector_parameters_pskel::
vector_parameters_pskel ()
: param_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// value_parameters_pskel
//

void value_parameters_pskel::
param_parser (::param1_pskel& p)
{
  this->param_parser_ = &p;
}

void value_parameters_pskel::
parsers (::param1_pskel& param)
{
  this->param_parser_ = &param;
}

value_parameters_pskel::
value_parameters_pskel ()
: param_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// param_pskel
//

void param_pskel::
value_parser (::xml_schema::string_pskel& p)
{
  this->value_parser_ = &p;
}

void param_pskel::
type_parser (::xml_schema::string_pskel& p)
{
  this->type_parser_ = &p;
}

void param_pskel::
parsers (::xml_schema::string_pskel& value,
         ::xml_schema::string_pskel& type)
{
  this->value_parser_ = &value;
  this->type_parser_ = &type;
}

param_pskel::
param_pskel ()
: value_parser_ (0),
  type_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// param1_pskel
//

void param1_pskel::
value_parser (::xml_schema::float_pskel& p)
{
  this->value_parser_ = &p;
}

void param1_pskel::
type_parser (::xml_schema::string_pskel& p)
{
  this->type_parser_ = &p;
}

void param1_pskel::
parsers (::xml_schema::float_pskel& value,
         ::xml_schema::string_pskel& type)
{
  this->value_parser_ = &value;
  this->type_parser_ = &type;
}

param1_pskel::
param1_pskel ()
: value_parser_ (0),
  type_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// input_pskel
//

void input_pskel::
vector_parameters (const std::array<int,3>&)
{
}

void input_pskel::
value_parameters (const std::int8_t&)
{
}

void input_pskel::
post_input ()
{
}

// vector_parameters_pskel
//

void vector_parameters_pskel::
param ()
{
}

void vector_parameters_pskel::
post_vector_parameters ()
{
}

// value_parameters_pskel
//

void value_parameters_pskel::
param ()
{
}

void value_parameters_pskel::
post_value_parameters ()
{
}

// param_pskel
//

void param_pskel::
value (const ::std::string&)
{
}

void param_pskel::
type (const ::std::string&)
{
}

void param_pskel::
post_param ()
{
}

// param1_pskel
//

void param1_pskel::
value (float)
{
}

void param1_pskel::
type (const ::std::string&)
{
}

void param1_pskel::
post_param1 ()
{
}

#include <cassert>

// Element validation and dispatch functions for input_pskel.
//
bool input_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  if (vd->func == 0 && vd->state == 0)
  {
    if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
      return true;
    else
      vd->state = 1;
  }

  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, ns, n, t, true);

    vd = vs.data + (vs.size - 1);

    if (vd->state == ~0UL)
      vd = vs.data + (--vs.size - 1);
    else
      break;
  }

  if (vd->func == 0)
  {
    if (vd->state != ~0UL)
    {
      unsigned long s = ~0UL;

      if (n == "vector_parameters" && ns.empty ())
        s = 0UL;

      if (s != ~0UL)
      {
        vd->count++;
        vd->state = ~0UL;

        vd = vs.data + vs.size++;
        vd->func = &input_pskel::sequence_0;
        vd->state = s;
        vd->count = 0;

        this->sequence_0 (vd->state, vd->count, ns, n, t, true);
      }
      else
      {
        if (vd->count < 1UL)
          this->_expected_element (
            "", "vector_parameters",
            ns, n);
        return false;
      }
    }
    else
      return false;
  }

  return true;
}

bool input_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size - 1];

  if (vd.func == 0 && vd.state == 0)
  {
    if (!::xml_schema::complex_content::_end_element_impl (ns, n))
      assert (false);
    return true;
  }

  assert (vd.func != 0);
  (this->*vd.func) (vd.state, vd.count, ns, n, 0, false);

  if (vd.state == ~0UL)
    vs.size--;

  return true;
}

void input_pskel::
_pre_e_validate ()
{
  this->v_state_stack_.push ();
  static_cast< v_state_* > (this->v_state_stack_.top ())->size = 0;

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size++];

  vd.func = 0;
  vd.state = 0;
  vd.count = 0;
}

void input_pskel::
_post_e_validate ()
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  ::xml_schema::ro_string empty;
  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, empty, empty, 0, true);
    assert (vd->state == ~0UL);
    vd = vs.data + (--vs.size - 1);
  }

  if (vd->count < 1UL)
    this->_expected_element (
      "", "vector_parameters");

  this->v_state_stack_.pop ();
}

void input_pskel::
sequence_0 (unsigned long& state,
            unsigned long& count,
            const ::xml_schema::ro_string& ns,
            const ::xml_schema::ro_string& n,
            const ::xml_schema::ro_string* t,
            bool start)
{
  XSD_UNUSED (t);

  switch (state)
  {
    case 0UL:
    {
      if (n == "vector_parameters" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->vector_parameters_parser_;

          if (this->vector_parameters_parser_)
            this->vector_parameters_parser_->pre ();
        }
        else
        {
          if (this->vector_parameters_parser_)
          {
            this->vector_parameters_parser_->post_vector_parameters ();
            this->vector_parameters (std::array<int,3>() );
          }

          count = 0;
          state = 1UL;
        }

        break;
      }
      else
      {
        assert (start);
        if (count < 1UL)
          this->_expected_element (
            "", "vector_parameters",
            ns, n);
        count = 0;
        state = 1UL;
      }
    }
    // Fall through.
    case 1UL:
    {
      if (n == "value_parameters" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->value_parameters_parser_;

          if (this->value_parameters_parser_)
            this->value_parameters_parser_->pre ();
        }
        else
        {
          if (this->value_parameters_parser_)
          {
            this->value_parameters_parser_->post_value_parameters ();
            this->value_parameters (std::int8_t());
          }

          count = 0;
          state = ~0UL;
        }

        break;
      }
      else
      {
        assert (start);
        if (count < 1UL)
          this->_expected_element (
            "", "value_parameters",
            ns, n);
        count = 0;
        state = ~0UL;
      }
    }
    // Fall through.
    case ~0UL:
      break;
  }
}

// Element validation and dispatch functions for vector_parameters_pskel.
//
bool vector_parameters_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  if (vd->func == 0 && vd->state == 0)
  {
    if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
      return true;
    else
      vd->state = 1;
  }

  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, ns, n, t, true);

    vd = vs.data + (vs.size - 1);

    if (vd->state == ~0UL)
      vd = vs.data + (--vs.size - 1);
    else
      break;
  }

  if (vd->func == 0)
  {
    if (vd->state != ~0UL)
    {
      unsigned long s = ~0UL;

      if (n == "param" && ns.empty ())
        s = 0UL;

      if (s != ~0UL)
      {
        vd->count++;
        vd->state = ~0UL;

        vd = vs.data + vs.size++;
        vd->func = &vector_parameters_pskel::sequence_0;
        vd->state = s;
        vd->count = 0;

        this->sequence_0 (vd->state, vd->count, ns, n, t, true);
      }
      else
      {
        return false;
      }
    }
    else
      return false;
  }

  return true;
}

bool vector_parameters_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size - 1];

  if (vd.func == 0 && vd.state == 0)
  {
    if (!::xml_schema::complex_content::_end_element_impl (ns, n))
      assert (false);
    return true;
  }

  assert (vd.func != 0);
  (this->*vd.func) (vd.state, vd.count, ns, n, 0, false);

  if (vd.state == ~0UL)
    vs.size--;

  return true;
}

void vector_parameters_pskel::
_pre_e_validate ()
{
  this->v_state_stack_.push ();
  static_cast< v_state_* > (this->v_state_stack_.top ())->size = 0;

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size++];

  vd.func = 0;
  vd.state = 0;
  vd.count = 0;
}

void vector_parameters_pskel::
_post_e_validate ()
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  ::xml_schema::ro_string empty;
  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, empty, empty, 0, true);
    assert (vd->state == ~0UL);
    vd = vs.data + (--vs.size - 1);
  }


  this->v_state_stack_.pop ();
}

void vector_parameters_pskel::
sequence_0 (unsigned long& state,
            unsigned long& count,
            const ::xml_schema::ro_string& ns,
            const ::xml_schema::ro_string& n,
            const ::xml_schema::ro_string* t,
            bool start)
{
  XSD_UNUSED (t);

  switch (state)
  {
    case 0UL:
    {
      if (n == "param" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->param_parser_;

          if (this->param_parser_)
            this->param_parser_->pre ();
        }
        else
        {
          if (this->param_parser_)
          {
            this->param_parser_->post_param ();
            this->param ();
          }

          count++;
        }

        break;
      }
      else
      {
        assert (start);
        count = 0;
        state = ~0UL;
      }
    }
    // Fall through.
    case ~0UL:
      break;
  }
}

// Element validation and dispatch functions for value_parameters_pskel.
//
bool value_parameters_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  if (vd->func == 0 && vd->state == 0)
  {
    if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
      return true;
    else
      vd->state = 1;
  }

  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, ns, n, t, true);

    vd = vs.data + (vs.size - 1);

    if (vd->state == ~0UL)
      vd = vs.data + (--vs.size - 1);
    else
      break;
  }

  if (vd->func == 0)
  {
    if (vd->state != ~0UL)
    {
      unsigned long s = ~0UL;

      if (n == "param" && ns.empty ())
        s = 0UL;

      if (s != ~0UL)
      {
        vd->count++;
        vd->state = ~0UL;

        vd = vs.data + vs.size++;
        vd->func = &value_parameters_pskel::sequence_0;
        vd->state = s;
        vd->count = 0;

        this->sequence_0 (vd->state, vd->count, ns, n, t, true);
      }
      else
      {
        return false;
      }
    }
    else
      return false;
  }

  return true;
}

bool value_parameters_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size - 1];

  if (vd.func == 0 && vd.state == 0)
  {
    if (!::xml_schema::complex_content::_end_element_impl (ns, n))
      assert (false);
    return true;
  }

  assert (vd.func != 0);
  (this->*vd.func) (vd.state, vd.count, ns, n, 0, false);

  if (vd.state == ~0UL)
    vs.size--;

  return true;
}

void value_parameters_pskel::
_pre_e_validate ()
{
  this->v_state_stack_.push ();
  static_cast< v_state_* > (this->v_state_stack_.top ())->size = 0;

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size++];

  vd.func = 0;
  vd.state = 0;
  vd.count = 0;
}

void value_parameters_pskel::
_post_e_validate ()
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  ::xml_schema::ro_string empty;
  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, empty, empty, 0, true);
    assert (vd->state == ~0UL);
    vd = vs.data + (--vs.size - 1);
  }


  this->v_state_stack_.pop ();
}

void value_parameters_pskel::
sequence_0 (unsigned long& state,
            unsigned long& count,
            const ::xml_schema::ro_string& ns,
            const ::xml_schema::ro_string& n,
            const ::xml_schema::ro_string* t,
            bool start)
{
  XSD_UNUSED (t);

  switch (state)
  {
    case 0UL:
    {
      if (n == "param" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->param_parser_;

          if (this->param_parser_)
            this->param_parser_->pre ();
        }
        else
        {
          if (this->param_parser_)
          {
            this->param_parser_->post_param1 ();
            this->param ();
          }

          count++;
        }

        break;
      }
      else
      {
        assert (start);
        count = 0;
        state = ~0UL;
      }
    }
    // Fall through.
    case ~0UL:
      break;
  }
}

// Element validation and dispatch functions for param_pskel.
//
bool param_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  if (vd->func == 0 && vd->state == 0)
  {
    if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
      return true;
    else
      vd->state = 1;
  }

  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, ns, n, t, true);

    vd = vs.data + (vs.size - 1);

    if (vd->state == ~0UL)
      vd = vs.data + (--vs.size - 1);
    else
      break;
  }

  if (vd->func == 0)
  {
    if (vd->state != ~0UL)
    {
      unsigned long s = ~0UL;

      if (n == "value" && ns.empty ())
        s = 0UL;

      if (s != ~0UL)
      {
        vd->count++;
        vd->state = ~0UL;

        vd = vs.data + vs.size++;
        vd->func = &param_pskel::sequence_0;
        vd->state = s;
        vd->count = 0;

        this->sequence_0 (vd->state, vd->count, ns, n, t, true);
      }
      else
      {
        return false;
      }
    }
    else
      return false;
  }

  return true;
}

bool param_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size - 1];

  if (vd.func == 0 && vd.state == 0)
  {
    if (!::xml_schema::complex_content::_end_element_impl (ns, n))
      assert (false);
    return true;
  }

  assert (vd.func != 0);
  (this->*vd.func) (vd.state, vd.count, ns, n, 0, false);

  if (vd.state == ~0UL)
    vs.size--;

  return true;
}

void param_pskel::
_pre_e_validate ()
{
  this->v_state_stack_.push ();
  static_cast< v_state_* > (this->v_state_stack_.top ())->size = 0;

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size++];

  vd.func = 0;
  vd.state = 0;
  vd.count = 0;
}

void param_pskel::
_post_e_validate ()
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  ::xml_schema::ro_string empty;
  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, empty, empty, 0, true);
    assert (vd->state == ~0UL);
    vd = vs.data + (--vs.size - 1);
  }


  this->v_state_stack_.pop ();
}

void param_pskel::
sequence_0 (unsigned long& state,
            unsigned long& count,
            const ::xml_schema::ro_string& ns,
            const ::xml_schema::ro_string& n,
            const ::xml_schema::ro_string* t,
            bool start)
{
  XSD_UNUSED (t);

  switch (state)
  {
    case 0UL:
    {
      if (n == "value" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->value_parser_;

          if (this->value_parser_)
            this->value_parser_->pre ();
        }
        else
        {
          if (this->value_parser_)
          {
            this->value (this->value_parser_->post_string ());
          }

          count++;
        }

        break;
      }
      else
      {
        assert (start);
        count = 0;
        state = ~0UL;
      }
    }
    // Fall through.
    case ~0UL:
      break;
  }
}

// Element validation and dispatch functions for param1_pskel.
//
bool param1_pskel::
_start_element_impl (const ::xml_schema::ro_string& ns,
                     const ::xml_schema::ro_string& n,
                     const ::xml_schema::ro_string* t)
{
  XSD_UNUSED (t);

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  if (vd->func == 0 && vd->state == 0)
  {
    if (this->::xml_schema::complex_content::_start_element_impl (ns, n, t))
      return true;
    else
      vd->state = 1;
  }

  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, ns, n, t, true);

    vd = vs.data + (vs.size - 1);

    if (vd->state == ~0UL)
      vd = vs.data + (--vs.size - 1);
    else
      break;
  }

  if (vd->func == 0)
  {
    if (vd->state != ~0UL)
    {
      unsigned long s = ~0UL;

      if (n == "value" && ns.empty ())
        s = 0UL;

      if (s != ~0UL)
      {
        vd->count++;
        vd->state = ~0UL;

        vd = vs.data + vs.size++;
        vd->func = &param1_pskel::sequence_0;
        vd->state = s;
        vd->count = 0;

        this->sequence_0 (vd->state, vd->count, ns, n, t, true);
      }
      else
      {
        if (vd->count < 1UL)
          this->_expected_element (
            "", "value",
            ns, n);
        return false;
      }
    }
    else
      return false;
  }

  return true;
}

bool param1_pskel::
_end_element_impl (const ::xml_schema::ro_string& ns,
                   const ::xml_schema::ro_string& n)
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size - 1];

  if (vd.func == 0 && vd.state == 0)
  {
    if (!::xml_schema::complex_content::_end_element_impl (ns, n))
      assert (false);
    return true;
  }

  assert (vd.func != 0);
  (this->*vd.func) (vd.state, vd.count, ns, n, 0, false);

  if (vd.state == ~0UL)
    vs.size--;

  return true;
}

void param1_pskel::
_pre_e_validate ()
{
  this->v_state_stack_.push ();
  static_cast< v_state_* > (this->v_state_stack_.top ())->size = 0;

  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_& vd = vs.data[vs.size++];

  vd.func = 0;
  vd.state = 0;
  vd.count = 0;
}

void param1_pskel::
_post_e_validate ()
{
  v_state_& vs = *static_cast< v_state_* > (this->v_state_stack_.top ());
  v_state_descr_* vd = vs.data + (vs.size - 1);

  ::xml_schema::ro_string empty;
  while (vd->func != 0)
  {
    (this->*vd->func) (vd->state, vd->count, empty, empty, 0, true);
    assert (vd->state == ~0UL);
    vd = vs.data + (--vs.size - 1);
  }

  if (vd->count < 1UL)
    this->_expected_element (
      "", "value");

  this->v_state_stack_.pop ();
}

void param1_pskel::
sequence_0 (unsigned long& state,
            unsigned long& count,
            const ::xml_schema::ro_string& ns,
            const ::xml_schema::ro_string& n,
            const ::xml_schema::ro_string* t,
            bool start)
{
  XSD_UNUSED (t);

  switch (state)
  {
    case 0UL:
    {
      if (n == "value" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->value_parser_;

          if (this->value_parser_)
            this->value_parser_->pre ();
        }
        else
        {
          if (this->value_parser_)
          {
            this->value (this->value_parser_->post_float ());
          }

          count = 0;
          state = ~0UL;
        }

        break;
      }
      else
      {
        assert (start);
        if (count < 1UL)
          this->_expected_element (
            "", "value",
            ns, n);
        count = 0;
        state = ~0UL;
      }
    }
    // Fall through.
    case ~0UL:
      break;
  }
}

// Attribute validation and dispatch functions for param_pskel.
//
bool param_pskel::
_attribute_impl_phase_one (const ::xml_schema::ro_string& ns,
                           const ::xml_schema::ro_string& n,
                           const ::xml_schema::ro_string& s)
{
  if (n == "type" && ns.empty ())
  {
    if (this->type_parser_)
    {
      this->type_parser_->pre ();
      this->type_parser_->_pre_impl ();
      this->type_parser_->_characters (s);
      this->type_parser_->_post_impl ();
      this->type (this->type_parser_->post_string ());
    }

    return true;
  }

  return false;
}

// Attribute validation and dispatch functions for param1_pskel.
//
bool param1_pskel::
_attribute_impl_phase_one (const ::xml_schema::ro_string& ns,
                           const ::xml_schema::ro_string& n,
                           const ::xml_schema::ro_string& s)
{
  if (n == "type" && ns.empty ())
  {
    if (this->type_parser_)
    {
      this->type_parser_->pre ();
      this->type_parser_->_pre_impl ();
      this->type_parser_->_characters (s);
      this->type_parser_->_post_impl ();
      this->type (this->type_parser_->post_string ());
    }

    return true;
  }

  return false;
}

// Character validation functions for vector_parameters_pskel.
//
bool vector_parameters_pskel::
_characters_impl (const ::xml_schema::ro_string& s)
{
  this->_any_characters (s);
  return true;
}

// Character validation functions for value_parameters_pskel.
//
bool value_parameters_pskel::
_characters_impl (const ::xml_schema::ro_string& s)
{
  this->_any_characters (s);
  return true;
}

// Character validation functions for param_pskel.
//
bool param_pskel::
_characters_impl (const ::xml_schema::ro_string& s)
{
  this->_any_characters (s);
  return true;
}

// Character validation functions for param1_pskel.
//
bool param1_pskel::
_characters_impl (const ::xml_schema::ro_string& s)
{
  this->_any_characters (s);
  return true;
}

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

