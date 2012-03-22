#ifndef GRID1D_HH
# define GRID1D_HH

#include <vector>
#include <list>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <complex>

typedef enum { OPEN_BOUNDARY, BOX_BOUNDARY, PERIODIC_BOUNDARY }  boundary_t;
typedef std::complex<double> complex_t; // TODO: organize

template <typename scalar_t, typename value_t = scalar_t> struct FieldTraits {
  static scalar_t contract(const value_t& x){ return x.contract(); }
  static value_t  conj(const value_t& x)    { return x.conj(); }
  static double   abs(const value_t& x)     { return x.abs(); }
};

template <> struct FieldTraits<double,double> {
  static double contract(const double& x){ return x; }
  static double conj(const double& x){ return x; }
  static double abs(const double& x) { return fabs(x); }
};

template <> struct FieldTraits<complex_t,complex_t> {
  static complex_t contract(const complex_t& x){ return x; }
  static complex_t conj(const complex_t& x){ return complex_t(x.real(),-x.imag()); }
  static double    abs(const complex_t& x) { return std::abs(x); }
};



template <typename Space>
class gridfunction : public std::vector<typename Space::value_t> {
public:
  typedef Space space_t;
  typedef typename Space::value_t  value_t;
  typedef typename Space::scalar_t scalar_t;
  typedef std::vector<value_t> basetype;
  typedef FieldTraits<scalar_t,value_t> Traits;

  gridfunction& operator += (const gridfunction& g);
  gridfunction& operator -= (const gridfunction& g);
  gridfunction& operator *= (const gridfunction& g);
  gridfunction& operator /= (const scalar_t& g);
  gridfunction& operator *= (const scalar_t& g);

  gridfunction operator + (const gridfunction& g) const;
  gridfunction operator - (const gridfunction& g) const;
  gridfunction operator * (const gridfunction& g) const;
  gridfunction operator / (const scalar_t& g) const;
  gridfunction operator * (const scalar_t& g) const;

  gridfunction conj () const;

  gridfunction(const unsigned int Nx=0) : basetype(Nx) {
    //    fprintf(stderr,"default constructor gridfunction()\n");
  }

  gridfunction(const Space& space) : std::vector<value_t>(space.Nx) {
    //    fprintf(stderr,"gridfunction(space)\n");
  }
  gridfunction(const std::vector<value_t>& v) : std::vector<value_t>(v.begin(),v.end()) { 
    //    fprintf(stderr,"gridfunction(vectord) : size(%d)\n",this->size());
  }

  std::vector<value_t> get_data() const {
    return std::vector<value_t>(this->begin(),this->end());
  }

  friend std::ostream& operator<<(std::ostream& s, const gridfunction& f) {
    s << "[";
    for(unsigned int i=0;i<f.size();i++) s << f[i] << (i+1<f.size()?", ":"]");
    return s;
  }
  
};


template <typename Scalar = double, typename Value = Scalar> class Grid1D {
public:
  typedef Scalar scalar_t;
  typedef Value value_t;
  typedef gridfunction<Grid1D> function_t;
  typedef Grid1D<Value> self_t;
  typedef FieldTraits<scalar_t,value_t> Traits;
  
  double xmin, xmax, dx;
  unsigned int Nx;
  boundary_t boundary_condition;

  class Term {
  public:
    scalar_t   v;
    function_t f;
    int type;
  
    Term(const scalar_t& v=0) : v(v), type(0) {}
    Term(const function_t& f) : f(f), type(1) {}
    
    function_t operator *(const function_t& g) const {
      if(type==0) return g*v;
      if(type==1) return g*f;
      fprintf(stderr,"gridterm type %d not implemented.\n",type);
      abort();
    }
  };

  
  Grid1D(double xmin=0,double xmax=1, unsigned int Nx=1, boundary_t boundary = PERIODIC_BOUNDARY) : xmin(xmin), xmax(xmax), dx((xmax-xmin)/static_cast<double>(Nx+(boundary==PERIODIC_BOUNDARY)-1)), Nx(Nx), boundary_condition(boundary) {}  

  void set_boundary(const boundary_t boundary){ 
    boundary_condition = boundary; 

    dx = (xmax-xmin)/static_cast<double>(Nx+(boundary_condition==PERIODIC_BOUNDARY)-1);
  }
  void set_boundary(const std::string& boundary){
  set_boundary(boundary == "open"? OPEN_BOUNDARY : (boundary == "box"? BOX_BOUNDARY : PERIODIC_BOUNDARY) );
  }

  value_t  integrate(const function_t& f) const;
  scalar_t inner(const function_t& f, const function_t& g) const;
  scalar_t inner(const function_t& f, const function_t& g, const function_t& h) const;
  scalar_t inner(const function_t& f, const function_t& g, const function_t& h, const function_t& m) const;

  double norm(const function_t& f) const;

  static inline scalar_t contract(const value_t& x){ return Traits::contract(x); }
  static inline value_t  conj(const value_t& x){ return Traits::conj(x); }

  function_t derivative       (const function_t& f) const;
  function_t second_derivative(const function_t& f) const;
  function_t third_derivative (const function_t& f) const;

  function_t stencil_operator(const function_t& f, const double *stencil, size_t stencil_length, 
			      double delta) const;


  std::vector<double> get_xs() const {
    std::vector<double> xs(Nx);
    for(size_t i=0;i<Nx;i++) xs[i] = xmin+i*dx;

    return xs;
  }

private:
  static double h2_derivative1_stencil[3], h2_derivative2_stencil[3], h2_derivative3_stencil[5];
  static double h4_derivative1_stencil[5], h4_derivative2_stencil[5], h4_derivative3_stencil[7];
};

#endif
