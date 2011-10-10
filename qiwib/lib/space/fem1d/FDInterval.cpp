#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <mkl.h>

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

using namespace std;
namespace ublas = boost::numeric::ublas;


class FDInterval {
  void print_matrix(FILE *stream, const ublas::banded_matrix<double>& A, const string name) const 
  {
    fprintf(stream,"%s = {",name.c_str());
    for(size_t i=0;i<A.size1();i++){
      fprintf(stream,"{");
      for(size_t j=0;j<A.size2();j++)
	fprintf(stream,"%f%s",A(i,j),j+1<A.size2()?",":"}");
      fprintf(stream,"%s",i+1<A.size1()?",\n":"\n};\n");
    }
  }

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
	y[i] = fun(space.x[i]);	// Value in the i'th node.
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
 v	s<< "{" << x[i] << "," << f.y[i] << "}" << (i+1<x.size()? "," : "}");

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

  FDInterval(const vector<double>& points) : a(*points.begin()), b(*points.end()), x(points)
  {
    clear_refine_flags();
  }

  // Approximate overlap of continuous function 'fun' with the
  // i'th shape function. 
  // If exact integrals are required, then use e.g. Gaussian 
  // quadrature instead. This is just a simple, *very* fast
  // and "good enough" method.
  inline double calculate_shape_overlap(size_t i, const RealFunction& fun) const 
  {
    double funval = fun(x[i]);
    double shape_integral = 0;

    if(i>=1)			// Not at left boundary
      shape_integral += 0.5*(x[i]-x[i-1])*funval; // Trapezoidal rule
    if(i<=x.size()-1)
      shape_integral += 0.5*(x[i+1]-x[i])*funval; // Trapezoidal rule

    return shape_integral;
  }

  ublas::vector<double> shape_overlaps(const RealFunction& fun) const 
  {
    ublas::vector<double> overlaps(x.size());
    for(size_t i=0;i<x.size();i++) 
      overlaps[i] = calculate_shape_overlap(i,fun);

    return overlaps;
  }

  // Approximate integral of continuous function on j'th cell.
  // Same considerations apply as for shape_overlap().
  inline double calculate_cell_integral(size_t j, const RealFunction& fun) const 
  {
    double cell_size = x[j+1]-x[j],
           avg_function_value = (fun(x[j])+fun(x[j+1]))/2.0;

    return cell_size*avg_function_value;
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
	double cell_integral = calculate_cell_integral(i,fun);

	if(cell_integral > 2*fun_max){
	  all_done = false;
	  set_refine_flag(i);
	} 

	if(i>=1) {
	  double left_cell_integral = calculate_cell_integral(i-1,fun);
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
    laplace(last-1,last) = -1/(x[last]-x[last-1]);

    return laplace;
  }

  // Solve Ax = b
  // I've implicitly added homogeneous Dirichlet boundary conditions, i.e.
  // the solution always has value zero at endpoints. This is the easiest case.
  // If other boundary conditions are necessary, a little extra work is necessary.
  ublas::vector<double> linear_solve(const ublas::banded_matrix<double>& A, 
			      const ublas::vector<double>& b) const
  {
    const unsigned int n = A.size1();

    const int N(n-2),NRHS(1);
    int INFO(0);
    double diagonal[n-2], subdiagonal[n-3];
    double rhs[n-2];
						
    for(size_t i=1;i<n-1;i++){
      diagonal[i-1] = A(i,i);      
      if(i>=2) subdiagonal[i-2] = A(i,i-1);
    }
    std::copy(b.begin(),b.end(),rhs);

   
    dptsv(&N,&NRHS,diagonal,subdiagonal,rhs,&N,&INFO);

    if(INFO){
      fprintf(stderr,"LAPACK dptsv (symmetric tridiagonal linear solver) failed with INFO = %d\n",INFO);
      abort();
    }

    ublas::vector<double> result(n);
    result[0]   = 0;		// Homogeneous Dirichlet
    result[n-1] = 0;
    copy(rhs,rhs+n-2,result.begin()+1);

    return result;
  };


  FDFunction solve_poisson(const RealFunction& density) const 
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

    print_matrix(stdout,laplace_matrix(),"A");
    solution = linear_solve(laplace_matrix(),(4*M_PI)*shape_overlaps(density));

    return FDFunction(*this,solution);
  }

};
