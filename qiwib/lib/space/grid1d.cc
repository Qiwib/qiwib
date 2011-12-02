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
scalar_operator(/)
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

grid_member(value_t) inner(const function_t& f, const function_t& g, const function_t& h, const function_t& m) const 
{
  value_t sum(0);
  for(unsigned int i=0;i<f.size();i++) sum += conj(f[i])*conj(g[i])*h[i]*m[i];

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


grid_member(gridfunction_t) derivative(const function_t& f) const
{
#ifdef DERIVATIVE_FIVE_POINT
  return stencil_operator(f,h4_derivative1_stencil,5,1.0/dx);
#else
  return stencil_operator(f,h2_derivative1_stencil,3,1.0/dx);
#endif
}

grid_member(gridfunction_t) second_derivative(const function_t& f) const
{
#ifdef DERIVATIVE_FIVE_POINT
  return stencil_operator(f,h4_derivative2_stencil,5,1.0/(dx*dx));
#else
  return stencil_operator(f,h2_derivative2_stencil,3,1.0/(dx*dx));
#endif
}

grid_member(gridfunction_t) third_derivative(const function_t& f) const
{
#ifdef DERIVATIVE_FIVE_POINT
  return stencil_operator(f,h4_derivative3_stencil,7,1.0/(dx*dx*dx));
#else
  return stencil_operator(f,h2_derivative3_stencil,5,1.0/(dx*dx*dx));
#endif
}

// C/C++ modulus is defined wrongly for all negative integers. This small modification works
// for -n < i < \infty.
#define mod(i,n) (i+n)%n

grid_member(gridfunction_t) stencil_operator(const function_t& f, const double *stencil, size_t stencil_length, double delta) const
{
  function_t df(*this);
  unsigned int n = f.size();
  unsigned int max = stencil_length/2;

  switch(boundary_condition){	
  case OPEN:                            // Open boundary condition 
                                        // (assume function continues outside [a;b]; C1-continuous on boundary)
    for(size_t i=max;i<n-max;i++){	// For this boundary type, numerical derivative is only
      value_t sum(0);			// defined on inner points
      for(size_t j=0;j<stencil_length;j++)
	sum += stencil[j]*f[i+j-max];	

      df[i] = sum*delta;	
    }
    for(size_t i=0;i<max;i++){        // Boundary point derivatives are set to same as last inner point
      df[i]     = df[max];	      // (TODO: C2-continuity instead of C1)
      df[n-i-1] = df[n-max-1];		
    }
    break;

  case BOX:		        // Homogeneous Dirichlet boundary condition (f(x) = 0 outside [a;b])
    for(size_t i=0;i<n;i++){	// End point is included in interval [a;b]
      value_t sum(0);
      for(size_t j=0;j<stencil_length;j++)
	if(i+j-max>=0 && i+j-max<n)
	  sum += stencil[j]*f[i+j-max,n];
      df[i] = sum*delta;
    }
    break;
  case PERIODIC:	        
    for(size_t i=0;i<n;i++){	// End point is not included in interval [a;b[
      value_t sum(0);
      for(size_t j=0;j<stencil_length;j++)
	sum += stencil[j]*f[mod(i+j-max,n)];
      df[i] = sum*delta;
    }
    break;

  default:
    fprintf(stderr,"Invalid input: boundary type %d\n",boundary_condition);
    abort();
  }

  return df;
}
