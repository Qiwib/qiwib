#include "grid1d.hh"
#include "pointwise_operator.hh"
#include <stdlib.h>

using namespace std;
#include <stdio.h>

#define gridfunction_member(T) template <class space> T gridfunction<space>::
#define grid_member(T) template <typename value_t> T Grid1D<value_t>::
#define gridfunction_t typename Grid1D<value_t>::function_t

// Define point-wise binary operations using the macros in pointwise_operator.hh.
pointwise_operator(+)
pointwise_operator(-)
pointwise_operator(*)
scalar_operator(*)  

gridfunction_member(gridfunction<space>) conj() const 
{
  gridfunction z(*this);

  for(size_t i=0;i<this->size();i++) z[i] = space::conj(z[i]);

  return z;
}

grid_member(value_t) integrate(const function_t& f) const 
{
  value_t sum(0);
  for(unsigned int i=0;i<f.size();i++) sum += f[i];

  return sum*dx;
}


grid_member(value_t) inner(const function_t& f, const function_t& g) const 
{
  value_t sum(0);
  for(unsigned int i=0;i<f.size();i++) sum += conj(f[i])*g[i];

  return sum*dx;
}

grid_member(value_t) inner(const function_t& f, const function_t& g, const function_t& h) const 
{
  value_t sum(0);
  for(unsigned int i=0;i<f.size();i++) sum += conj(f[i])*g[i]*h[i];

  return sum*dx;
}

// Definition of conjugate depends on the scalar field.
template<> double Grid1D<double>::conj(const double& v) { return v; }
template<> complex<double> Grid1D< complex<double> >::conj(const complex<double>& v) { return complex<double>(v.real(),-v.imag()); }

grid_member(double) h2_derivative1_stencil[3] = {-1/2.,0,1/2.};
grid_member(double) h2_derivative2_stencil[3] = {1,-2,1};
grid_member(double) h2_derivative3_stencil[5] = {-1/2.,1,0,-1,1/2.};

grid_member(double) h4_derivative1_stencil[5] = { 1/12.,-2/3.,  0,  2/3.,-1/12.};
grid_member(double) h4_derivative2_stencil[5] = {-1/12., 4/3.,-5/2.,4/3.,-1/12.};
grid_member(double) h4_derivative3_stencil[7] = {1/8.,-1,13/8.,0,-13/8.,1,-1/8.};


grid_member(gridfunction_t) derivative(const function_t& f, bool periodic) const
{
#ifdef DERIVATIVE_FIVE_POINT
  return stencil_operator(f,h4_derivative1_stencil,5,1.0/dx,periodic);
#else
  return stencil_operator(f,h2_derivative1_stencil,3,1.0/dx,periodic);
#endif
}

grid_member(gridfunction_t) second_derivative(const function_t& f, bool periodic) const
{
#ifdef DERIVATIVE_FIVE_POINT
  return stencil_operator(f,h4_derivative2_stencil,5,1.0/(dx*dx),periodic);
#else
  return stencil_operator(f,h2_derivative2_stencil,3,1.0/(dx*dx),periodic);
#endif
}

grid_member(gridfunction_t) third_derivative(const function_t& f, bool periodic) const
{
#ifdef DERIVATIVE_FIVE_POINT
  return stencil_operator(f,h4_derivative3_stencil,7,1.0/(dx*dx*dx),periodic);
#else
  return stencil_operator(f,h2_derivative3_stencil,5,1.0/(dx*dx*dx),periodic);
#endif
}

// C/C++ modulus is defined wrongly for all negative integers. This small modification works
// for -n < i < \infty.
#define mod(i,n) (i+n)%n

grid_member(gridfunction_t) stencil_operator(const function_t& f, const double *stencil, size_t stencil_length,
					     double delta, bool periodic) const
{
  function_t df(*this);
  unsigned int n = f.size();
  int max = stencil_length/2;
  
  if(periodic){
    n = n-1;			// Assume that end point is copy of start point
    for(size_t i=0;i<n;i++){
      value_t sum(0);
      for(size_t j=0;j<stencil_length;j++)
	sum += stencil[j]*f[mod(i+j-max,n)];
      df[i] = sum*delta;
    }
    df[n] = df[0];
  } else {
    for(size_t i=max;i<n-max;i++){
      value_t sum(0);
      for(size_t j=0;j<stencil_length;j++){
	sum += stencil[j]*f[i+j-max];
      }
      df[i] = sum*delta;
      for(size_t i=0;i<max;i++){
	df[i]     = df[max];
	df[n-i-1] = df[n-max-1];
      }
    }
  }

  return df;
}
