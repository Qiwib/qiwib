#include "grid1d.hh"

using namespace std;
#include <stdio.h>

#define gridfunction_member(T) template <class space> T gridfunction<space>::
#define grid_member(T) template <typename value_t> T Grid1D<value_t>::
#define gridfunction_t typename Grid1D<value_t>::function_t

gridfunction_member(gridfunction<space>&) operator+=(const gridfunction<space>& g)
{
  const_iterator y(g.begin());
  for(iterator x(this->begin()); x!=this->end();x++,y++) *x += *y;

  return *this;
}
gridfunction_member(gridfunction<space>&) operator-=(const gridfunction& g)
{
  const_iterator y(g.begin());
  for(iterator x(this->begin()); x!=this->end();x++,y++) *x -= *y;

  return *this;
}

gridfunction_member(gridfunction<space>&) operator*=(const gridfunction& g)
{
  const_iterator y(g.begin());
  for(iterator x(this->begin()); x!=this->end();x++,y++) *x *= *y;

  return *this;
}

gridfunction_member(gridfunction<space>&) operator*=(const gridfunction::scalar_t& s)
{
  for(iterator x(this->begin()); x!=this->end();x++) *x *= s;

  return *this;
}

gridfunction_member(gridfunction<space>) operator+(const gridfunction& g) const 
{
  gridfunction z(*this);
  return (z += g);
}

gridfunction_member(gridfunction<space>) operator-(const gridfunction& g) const 
{
  gridfunction z(*this);
  return (z -= g);
}

gridfunction_member(gridfunction<space>) operator*(const gridfunction& g) const 
{
  gridfunction z(*this);
  return (z *= g);
}

gridfunction_member(gridfunction<space>) operator*(const scalar_t& s) const 
{
  gridfunction z(*this);
  return (z *= s);
}

grid_member(value_t) integrate(const function_t& f) const 
{
  value_t sum(0);
  for(typename function_t::const_iterator x(f.begin()); x != f.end(); x++)
    sum += *x;
  return sum*dx;
}


grid_member(value_t) inner(const function_t& f, const function_t& g) const 
{
  value_t sum(0);
  for(typename function_t::const_iterator x(f.begin()), y(g.begin()); x != f.end(); x++, y++)
    sum += (*x) * (*y);
  return sum*dx;
}

grid_member(value_t) inner(const function_t& f, const function_t& g, const function_t& h) const 
{
  value_t sum(0);
  for(typename function_t::const_iterator x(f.begin()), y(g.begin()), z(h.begin()); x != f.end(); x++, y++, z++)
    sum += (*x) * (*y) * (*z);

  return sum*dx;
}

grid_member(gridfunction_t) derivative(const function_t& f) const
{
  function_t df(*this);
  typename function_t::const_iterator fx(f.begin());
  typename function_t::iterator dfx(df.begin());

  *dfx++ = ((*(fx+1)) -  (*fx))/dx;
  for(;dfx+1 != df.end();fx++,dfx++)
    *dfx = ((*(fx+2)) - (*fx))/(2*dx);
  *dfx = ((*(fx+1)) - (*fx))/dx;

  return df;
}
