#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/lapack/gbsv.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>

using namespace std;
namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;


class FDInterval {
public:
  const double a, b;
  vector<double> x;
  vector<bool> refine_flags, coarsen_flags;
  size_t next_length;

  class RealFunction {
  public:
    virtual double operator()(double x) const = 0;
  };
  
  class FDFunction {
  public: 
    const FDInterval& space;
    ublas::vector<double> y;

    FDFunction(const FDInterval& space) : space(space), y(space.x.size()) {}
    FDFunction(const FDInterval& space, const ublas::vector<double>& y) : space(space), y(y) {}
    FDFunction(const FDInterval& space, const RealFunction& fun) : space(space) 
    {
      y.resize(space.x.size());

      for(size_t i=0;i<space.x.size();i++)
	y[i] = fun(space.x[i]);
    }

    FDFunction derivative() const {
      FDFunction df(space);


      for(size_t i=1;i<space.x.size()-1;i++)
	df.y[i] = derivative(i);
      
      // zero curvature at endpoints
      df.y[0] = df.y[1]; 	
      df.y[y.size()-1] = df.y[y.size()-2];
      
      return df;
    }

    double derivative(size_t i) const {
      // Central finite difference      
      return (y[i+1]-y[i-1])/(space.x[i+1]-space.x[i-1]);

      // // Forward difference      
      // return (y[i+1]-y[i])/(space.x[i+1]-space.x[i]);
    }

    FDFunction second_derivative() const {
      return derivative().derivative();
    }

    friend ostream& operator << (ostream& s, const FDFunction& f)
    {
      ios::fmtflags stream_flags(s.flags()); 

      s << std::fixed << "{";
      const vector<double>& x(f.space.x);
      for(size_t i=0;i<x.size();i++)
	s<< "{" << x[i] << "," << f.y[i] << "}" << (i+1<x.size()? "," : "}");

      cout.flags( stream_flags );
      return s;
    }
  };

  FDInterval(double a, double b) : a(a), b(b), x(2)
  {
    x[0] = a;
    x[1] = b;
    clear_refine_flags();
  }
  
  void clear_refine_flags()
  {
    refine_flags.resize(x.size());
    coarsen_flags.resize(x.size());

    fill(refine_flags.begin(),  refine_flags.end(), false);
    fill(coarsen_flags.begin(), coarsen_flags.end(), false);
    next_length = x.size();
  }

  void set_refine_flag(size_t i)
  {
    assert(i+1<x.size());

    if(not refine_flags[i]){
      next_length++;
      refine_flags[i] = true;
    }

  }

  void set_coarsen_flag(size_t i)
  {
    assert(i+1<x.size());

    if(not coarsen_flags[i]){
      next_length--;
      coarsen_flags[i] = true;
    }

  }

  void refine_cells()
  {
    vector<double> new_x(next_length);

    for(size_t i=0, j=0;i < x.size();i++){

      if(!(i+1<x.size() && coarsen_flags[i]))
	new_x[j++] = x[i];

      if(i+1<x.size() && refine_flags[i]){
	double cell_midpoint = (x[i]+x[i+1])/2.0;
	new_x[j++] = cell_midpoint;
      }

    }
    
    x = new_x;
    clear_refine_flags();
  }

  void refine_around_points(const vector<double>& xs, double distance, double spacing)
  {
    bool all_done = false;

    while(not all_done){
      all_done = true;
      
      for(size_t i =0;i<x.size()-1;i++)
	for(size_t j=0;j<xs.size();j++)
	  if(
	      ((xs[j] >= x[i] && xs[j] <= x[i+1]) // xs[j] is inside cell
	       || fabs(x[i]-xs[j]) < distance || fabs(x[i+1]-xs[j]) < distance) // or xs[j] is close to cell
	      && (x[i+1]-x[i] > 2.0*spacing) // and cell is too large
	      ){
	    all_done = false;
	    set_refine_flag(i);
	  }

      refine_cells();
    }
  }

  void refine_to_function(const RealFunction& fun, double fun_min, double fun_max)
  {
    bool all_done = false;

    while(not all_done){
      all_done = true;

      for(size_t i=0;i<x.size()-1;i++){
	double cell_size     = x[i+1]-x[i];

	// Approximate cell integral of function
	double cell_integral = cell_size*(fun(x[i])+fun(x[i+1])/2.0);
	//	fprintf(stderr,"cell integral %g: %g*avg(%g,%g) = %g\n",x[i],cell_size,fun(x[i]),fun(x[i+1]),cell_integral);

	if(cell_integral > 2*fun_max){
	  all_done = false;
	  set_refine_flag(i);
	} 

	if(i>=1) {
	  double left_cell_integral = cell_size*(fun(x[i-1])+fun(x[i])/2.0);
	  if(left_cell_integral + cell_integral < fun_min){
	    all_done = false;
	    set_coarsen_flag(i);
	  }
	}	
      }
      refine_cells();
    }
  }

  ublas::banded_matrix<double> laplace_matrix() const 
  {
    ublas::banded_matrix<double> laplace(x.size(),x.size(),1,1);

    // Generate finite difference Laplace matrix:
    for(size_t i=1;i<x.size()-1;i++){
      double hleft = x[i]-x[i-1],
   	    hright = x[i+1]-x[i];

      laplace(i,i)   = 1/hleft + 1/hright;
      laplace(i,i-1) = -1/hleft;
      laplace(i-1,i) = -1/hleft;
    }

    // We treat the endpoints separately.
    laplace(0,0) = 1/(x[1]-x[0]);
    size_t last = x.size()-1;
    laplace(last,last)   =  1/(x[last]-x[last-1]);
    laplace(last,last-1) = -1/(x[last]-x[last-1]);

    return laplace;
  }

  ublas::vector<double> solve_poisson(const FDFunction& density) const 
  {
    // This is an example of how to solve a differential equation
    // on the grid. Your operator is different, but you'll do 
    // essentially the same thing. 
    // 
    // Because of the uniform structure, you don't need to actually
    // construct the matrix, but can use very efficient multigrid
    // solvers. If you were using hundreds of thousands points, you
    // should do it that way, but since you aren't, the extra
    // performance doesn't matter, and we'll just do it straight
    // forwardly.

    ublas::vector<double> solution(x.size());
    ublas::banded_matrix<double> laplace(laplace_matrix());

    lapack::gbtrf(laplace,solution,ublas::lower_tag());

    return solution;
  }
};
