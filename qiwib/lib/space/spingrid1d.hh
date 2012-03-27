#ifndef SPINGRID1D_HH
# define SPINGRID1D_HH

#include "grid1d.hh"

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

template <typename scalar_t, unsigned int n> class SpinValue {
public:
  scalar_t x[n];
  unsigned int N;

  SpinValue(const scalar_t v=0) : N(n)  {  for(unsigned int i=0;i<N;i++) x[i] = v; }//fill(&x[0],&x[n-1],v); }
  SpinValue(const scalar_t *v, unsigned int m) : N(n)  { assert(m==n); memcpy(x,v,n*sizeof(scalar_t)); }
  SpinValue(const SpinValue& v) : N(n) { copy(v.x,v.x+n,x); }
  SpinValue(const std::vector<scalar_t>& v) : N(n) { copy(v.begin(),v.end(),x); }
  
  inline SpinValue& operator+=(const SpinValue& y) { for(unsigned int i=0;i<n;i++) x[i] += y.x[i]; return *this; }
  inline SpinValue& operator-=(const SpinValue& y) { for(unsigned int i=0;i<n;i++) x[i] -= y.x[i]; return *this; }
  inline SpinValue& operator*=(const SpinValue& y) { for(unsigned int i=0;i<n;i++) x[i] *= y.x[i]; return *this; }
  inline SpinValue& operator/=(const SpinValue& y) { for(unsigned int i=0;i<n;i++) x[i] /= y.x[i]; return *this; }
  inline SpinValue& operator*=(const scalar_t& s)  { for(unsigned int i=0;i<n;i++) x[i] *= s;      return *this; }
  inline SpinValue operator+(const SpinValue& y) const { SpinValue result(*this); return result += y; }
  inline SpinValue operator-(const SpinValue& y) const { SpinValue result(*this); return result -= y; }
  inline SpinValue operator*(const SpinValue& y) const { SpinValue result(*this); return result *= y; }
  inline SpinValue operator/(const SpinValue& y) const { SpinValue result(*this); return result /= y; }
  inline SpinValue operator*(const scalar_t& s)  const { SpinValue result(*this); return result *= s; }
  inline friend SpinValue operator*(const scalar_t& s, const SpinValue& v){ return v*s; }
  
  inline SpinValue conj() const {
    SpinValue result(*this);
    for(unsigned int i=0;i<n;i++) result.x[i] = FieldTraits<scalar_t>::conj(result.x[i]);
    return result;
  }

  inline scalar_t contract() const {
    scalar_t sum = 0;
    for(unsigned int i=0;i<n;i++) sum += x[i];
    return sum;
  }

  inline double abs() const {
    double sum = 0;
    for(unsigned int i=0;i<n;i++) sum += FieldTraits<scalar_t>::abs(x[i]);
    return sum;
  }

  inline scalar_t& operator[](const unsigned int i){ return x[i]; }
  inline scalar_t  operator[](const unsigned int i) const { return x[i]; }

  friend ostream& operator<<(ostream& stream, const SpinValue& v){
    stream << "("; 
    for(unsigned int i=0;i<n;i++) stream << v.x[i] << (i+1<n?", ":")");
    return stream;
  }
};

template <typename scalar_t, unsigned int n> class SpinGrid1D : public Grid1D< scalar_t, SpinValue<scalar_t,n> > {
public:
  typedef Grid1D< scalar_t, SpinValue<scalar_t,n> > basetype;
  typedef Grid1D< scalar_t > space_component_t;
  typedef gridfunction<space_component_t> function_component_t;
  typedef typename basetype::value_t    value_t;
  typedef typename basetype::function_t function_t;

//   class Term {
//   public:
//     scalar_t v;
//     function_t f;
//     Array2D<scalar_t> A;    
//     Array2D<function_component_t> Af;
//     int type;
//   
//     Term(const scalar_t& v=0) : v(v), type(0) {}
//     Term(const function_t& f) : f(f), type(1) {}
//     Term(const Array2D<scalar_t>& A) : A(A), type(2) {}
//     Term(const Array2D<function_component_t>& Af) : Af(Af), type(3) {}
//     
//     function_t operator *(const function_t& g) const {
//       switch(type){
//       case 0: return g*v;
//       case 1: return g*f;
//       case 2: 
//       case 3: {
// 	function_t h(g.size());
// 	
// 	for(int ix=0;ix<g.size();ix++){
// 	  value_t hx(0);
// 
// 	  if(type == 2){
// 	    for(int i=0;i<A.rows();i++)
// 	      for(int j=0;j<A.columns();j++)
// 		hx[i] += A(i,j) *g[ix][j];
// 	  } else {
// 	    for(int i=0;i<Af.rows();i++)
// 	      for(int j=0;j<Af.columns();j++)
// 		hx[i] += Af(i,j)[ix] * g[ix][j];
// 	  }
// 	  h[ix] = hx;
// 	}
// 	return h;
//       }
//       default:
// 	fprintf(stderr,"gridterm type %d not implemented.\n",type);
// 	abort();
//       }
//     }
//   };


  SpinGrid1D(double xmin=0,double xmax=1, unsigned int Nx=1, boundary_t boundary = PERIODIC_BOUNDARY) : basetype(xmin,xmax,Nx,boundary)
  {}
};

template <typename scalar_t, unsigned int n> class spingridfunction : public gridfunction<SpinGrid1D<scalar_t,n> > {
public:
  typedef Grid1D< scalar_t > space_component_t;
  typedef gridfunction<space_component_t> function_component_t;
  typedef SpinGrid1D<scalar_t,n> Space;
  typedef gridfunction<Space> basetype;
  typedef typename Space::value_t value_t;

  spingridfunction() {  }
  spingridfunction(const Space& space) : basetype(space) {}
  spingridfunction(const std::vector<value_t>& v) : basetype(v) { }

  spingridfunction(const Array2D<scalar_t>& vs) {
    std::vector<value_t> v(vs.rows());
    for(unsigned int i=0;i<vs.rows();i++)
      v[i] = value_t(&vs[i*n],n);
    *this = v;
  }

  function_component_t component(unsigned int j) const { 
    function_component_t f(this->size());
    for(unsigned int i=0;i<f.size();i++)
      f[i] = (*this)[i][j];
    return f;
  }

};




#endif
