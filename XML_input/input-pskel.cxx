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

// parameters_pskel
//

void parameters_pskel::
coordinates_parser (::coordinates_pskel& p)
{
  this->coordinates_parser_ = &p;
}

void parameters_pskel::
values_parser (::values_pskel& p)
{
  this->values_parser_ = &p;
}

void parameters_pskel::
parsers (::coordinates_pskel& coordinates,
         ::values_pskel& values)
{
  this->coordinates_parser_ = &coordinates;
  this->values_parser_ = &values;
}

parameters_pskel::
parameters_pskel ()
: coordinates_parser_ (0),
  values_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// coordinates_pskel
//

void coordinates_pskel::
name_parser (::name_pskel& p)
{
  this->name_parser_ = &p;
}

void coordinates_pskel::
parsers (::name_pskel& name)
{
  this->name_parser_ = &name;
}

coordinates_pskel::
coordinates_pskel ()
: name_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// values_pskel
//

void values_pskel::
name_parser (::name1_pskel& p)
{
  this->name_parser_ = &p;
}

void values_pskel::
parsers (::name1_pskel& name)
{
  this->name_parser_ = &name;
}

values_pskel::
values_pskel ()
: name_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// name_pskel
//

void name_pskel::
x_parser (::xml_schema::short_pskel& p)
{
  this->x_parser_ = &p;
}

void name_pskel::
y_parser (::xml_schema::short_pskel& p)
{
  this->y_parser_ = &p;
}

void name_pskel::
z_parser (::xml_schema::short_pskel& p)
{
  this->z_parser_ = &p;
}

void name_pskel::
parsers (::xml_schema::short_pskel& x,
         ::xml_schema::short_pskel& y,
         ::xml_schema::short_pskel& z)
{
  this->x_parser_ = &x;
  this->y_parser_ = &y;
  this->z_parser_ = &z;
}

name_pskel::
name_pskel ()
: x_parser_ (0),
  y_parser_ (0),
  z_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// name1_pskel
//

void name1_pskel::
value_parser (::xml_schema::float_pskel& p)
{
  this->value_parser_ = &p;
}

void name1_pskel::
parsers (::xml_schema::float_pskel& value)
{
  this->value_parser_ = &value;
}

name1_pskel::
name1_pskel ()
: value_parser_ (0),
  v_state_stack_ (sizeof (v_state_), &v_state_first_)
{
}

// parameters_pskel
//

void parameters_pskel::
coordinates ()
{
}

void parameters_pskel::
values ()
{
}

void parameters_pskel::
post_parameters ()
{
}

// coordinates_pskel
//

void coordinates_pskel::
name ()
{
}

void coordinates_pskel::
post_coordinates ()
{
}

// values_pskel
//

void values_pskel::
name ()
{
}

void values_pskel::
post_values ()
{
}

// name_pskel
//

void name_pskel::
x (short)
{
}

void name_pskel::
y (short)
{
}

void name_pskel::
z (short)
{
}

void name_pskel::
post_name ()
{
}

// name1_pskel
//

void name1_pskel::
value (float)
{
}

void name1_pskel::
post_name1 ()
{
}

#include <cassert>

// Element validation and dispatch functions for parameters_pskel.
//
bool parameters_pskel::
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

      if (n == "coordinates" && ns.empty ())
        s = 0UL;

      if (s != ~0UL)
      {
        vd->count++;
        vd->state = ~0UL;

        vd = vs.data + vs.size++;
        vd->func = &parameters_pskel::sequence_0;
        vd->state = s;
        vd->count = 0;

        this->sequence_0 (vd->state, vd->count, ns, n, t, true);
      }
      else
      {
        if (vd->count < 1UL)
          this->_expected_element (
            "", "coordinates",
            ns, n);
        return false;
      }
    }
    else
      return false;
  }

  return true;
}

bool parameters_pskel::
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

void parameters_pskel::
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

void parameters_pskel::
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
      "", "coordinates");

  this->v_state_stack_.pop ();
}

void parameters_pskel::
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
      if (n == "coordinates" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->coordinates_parser_;

          if (this->coordinates_parser_)
            this->coordinates_parser_->pre ();
        }
        else
        {
          if (this->coordinates_parser_)
          {
            this->coordinates_parser_->post_coordinates ();
            this->coordinates ();
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
            "", "coordinates",
            ns, n);
        count = 0;
        state = 1UL;
      }
    }
    // Fall through.
    case 1UL:
    {
      if (n == "values" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->values_parser_;

          if (this->values_parser_)
            this->values_parser_->pre ();
        }
        else
        {
          if (this->values_parser_)
          {
            this->values_parser_->post_values ();
            this->values ();
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
            "", "values",
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

// Element validation and dispatch functions for coordinates_pskel.
//
bool coordinates_pskel::
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

      if (n == "name" && ns.empty ())
        s = 0UL;

      if (s != ~0UL)
      {
        vd->count++;
        vd->state = ~0UL;

        vd = vs.data + vs.size++;
        vd->func = &coordinates_pskel::sequence_0;
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

bool coordinates_pskel::
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

void coordinates_pskel::
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

void coordinates_pskel::
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

void coordinates_pskel::
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
      if (n == "name" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->name_parser_;

          if (this->name_parser_)
            this->name_parser_->pre ();
        }
        else
        {
          if (this->name_parser_)
          {
            this->name_parser_->post_name ();
            this->name ();
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

// Element validation and dispatch functions for values_pskel.
//
bool values_pskel::
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

      if (n == "name" && ns.empty ())
        s = 0UL;

      if (s != ~0UL)
      {
        vd->count++;
        vd->state = ~0UL;

        vd = vs.data + vs.size++;
        vd->func = &values_pskel::sequence_0;
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

bool values_pskel::
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

void values_pskel::
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

void values_pskel::
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

void values_pskel::
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
      if (n == "name" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->name_parser_;

          if (this->name_parser_)
            this->name_parser_->pre ();
        }
        else
        {
          if (this->name_parser_)
          {
            this->name_parser_->post_name1 ();
            this->name ();
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

// Element validation and dispatch functions for name_pskel.
//
bool name_pskel::
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

      if (n == "x" && ns.empty ())
        s = 0UL;

      if (s != ~0UL)
      {
        vd->count++;
        vd->state = ~0UL;

        vd = vs.data + vs.size++;
        vd->func = &name_pskel::sequence_0;
        vd->state = s;
        vd->count = 0;

        this->sequence_0 (vd->state, vd->count, ns, n, t, true);
      }
      else
      {
        if (vd->count < 1UL)
          this->_expected_element (
            "", "x",
            ns, n);
        return false;
      }
    }
    else
      return false;
  }

  return true;
}

bool name_pskel::
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

void name_pskel::
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

void name_pskel::
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
      "", "x");

  this->v_state_stack_.pop ();
}

void name_pskel::
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
      if (n == "x" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->x_parser_;

          if (this->x_parser_)
            this->x_parser_->pre ();
        }
        else
        {
          if (this->x_parser_)
          {
            this->x (this->x_parser_->post_short ());
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
            "", "x",
            ns, n);
        count = 0;
        state = 1UL;
      }
    }
    // Fall through.
    case 1UL:
    {
      if (n == "y" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->y_parser_;

          if (this->y_parser_)
            this->y_parser_->pre ();
        }
        else
        {
          if (this->y_parser_)
          {
            this->y (this->y_parser_->post_short ());
          }

          count = 0;
          state = 2UL;
        }

        break;
      }
      else
      {
        assert (start);
        if (count < 1UL)
          this->_expected_element (
            "", "y",
            ns, n);
        count = 0;
        state = 2UL;
      }
    }
    // Fall through.
    case 2UL:
    {
      if (n == "z" && ns.empty ())
      {
        if (start)
        {
          this->::xml_schema::complex_content::context_.top ().parser_ = this->z_parser_;

          if (this->z_parser_)
            this->z_parser_->pre ();
        }
        else
        {
          if (this->z_parser_)
          {
            this->z (this->z_parser_->post_short ());
          }

          count = 0;
          state = ~0UL;
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

// Element validation and dispatch functions for name1_pskel.
//
bool name1_pskel::
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
        vd->func = &name1_pskel::sequence_0;
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

bool name1_pskel::
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

void name1_pskel::
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

void name1_pskel::
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

void name1_pskel::
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

// Character validation functions for name_pskel.
//
bool name_pskel::
_characters_impl (const ::xml_schema::ro_string& s)
{
  this->_any_characters (s);
  return true;
}

// Character validation functions for name1_pskel.
//
bool name1_pskel::
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

