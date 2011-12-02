#ifndef GRID1D_HH
# define GRID1D_HH

#include <vector>
#include <list>
#include <string>
#include <stdio.h>

template <typename Space>
class gridfunction : public std::vector<typename Space::value_t> {
public:
  typedef Space space_t;
  typedef typename Space::value_t  value_t;
  typedef typename Space::scalar_t scalar_t;

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

  gridfunction(const Space& space) : std::vector<value_t>(space.Nx) {
    //    fprintf(stderr,"gridfunction(space)\n");
  }
  gridfunction(const std::vector<value_t>& v) : std::vector<value_t>(v.begin(),v.end()) { 
    //    fprintf(stderr,"gridfunction(vectord) : size(%d)\n",this->size());
  }

  gridfunction(const double *v, unsigned int n) : std::vector<value_t>(v,v+n) { 
    //    fprintf(stderr,"gridfunction(pointer) : size(%d)\n",this->size());
  }

  gridfunction() {
    // fprintf(stderr,"default constructor gridfunction()\n");
  }

  std::vector<value_t> get_data() const {
    return std::vector<value_t>(this->begin(),this->end());
  }
};


template <typename Value = double> class Grid1D {
public:
  typedef Value scalar_t;
  typedef Value value_t;
  typedef gridfunction<Grid1D> function_t;
  typedef Grid1D<Value> self_t;
  typedef enum { OPEN, BOX, PERIODIC }  boundary_t;
  
  double xmin, xmax, dx;
  unsigned int Nx;
  boundary_t boundary_condition;
  
  Grid1D(double xmin=0,double xmax=1, unsigned int Nx=1, boundary_t boundary = PERIODIC) : xmin(xmin), xmax(xmax), dx((xmax-xmin)/static_cast<double>(Nx+(boundary==PERIODIC)-1)), Nx(Nx), boundary_condition(boundary) {}  

  void set_boundary(const boundary_t boundary){ 
    boundary_condition = boundary; 

    dx = (xmax-xmin)/static_cast<double>(Nx+(boundary_condition==PERIODIC)-1);
  }
  void set_boundary(const std::string& boundary){
    set_boundary(boundary == "open"? OPEN : (boundary == "box"? BOX : PERIODIC) );
  }

  value_t  integrate(const function_t& f) const;
  scalar_t inner(const function_t& f, const function_t& g) const;
  scalar_t inner(const function_t& f, const function_t& g, const function_t& h) const;
  scalar_t inner(const function_t& f, const function_t& g, const function_t& h, const function_t& m) const;
  static value_t  conj(const value_t& v);

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
