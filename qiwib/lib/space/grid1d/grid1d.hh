#ifndef GRID1D_HH
# define GRID1D_HH

#include <vector>
#include <string>
#include <stdio.h>

template <typename Space>
class gridfunction : public std::vector<typename Space::value_t> {
public:
  typedef typename Space::value_t  value_t;
  typedef typename Space::scalar_t scalar_t;

  typedef typename std::vector<value_t>::const_iterator const_iterator;
  typedef typename std::vector<value_t>::iterator iterator; 

  gridfunction& operator += (const gridfunction& g);
  gridfunction& operator -= (const gridfunction& g);
  gridfunction& operator *= (const gridfunction& g);
  gridfunction& operator *= (const scalar_t& g);

  gridfunction operator + (const gridfunction& g) const;
  gridfunction operator - (const gridfunction& g) const;
  gridfunction operator * (const gridfunction& g) const;
  gridfunction operator * (const scalar_t& g) const;

  gridfunction(const Space& space) : std::vector<value_t>(space.Nx) {
    //    fprintf(stderr,"gridfunction(space)\n");
  }
  gridfunction(const std::vector<value_t>& values) : std::vector<value_t>(values.begin(),values.end()) { 
    //    fprintf(stderr,"gridfunction(vectord) : size(%d)\n",values.size());
  }

  const std::vector<value_t> get_data() const {
    return std::vector<value_t>(this->begin(), this->end());
  }
};

template <typename Value = double> class Grid1D {
public:
  typedef Value scalar_t;
  typedef Value value_t;
  typedef gridfunction<Grid1D> function_t;
  typedef Grid1D<Value> self_t;

  double xmin, xmax, dx;
  unsigned int Nx;

  Grid1D(double xmin=0,double xmax=1, unsigned int Nx=1) : xmin(xmin), xmax(xmax),  
							   dx((xmax-xmin)/static_cast<double>(Nx)),
							   Nx(Nx) {}
  
  value_t  integrate(const function_t& f) const;
  scalar_t inner(const function_t& f, const function_t& g) const;
  scalar_t inner(const function_t& f, const function_t& g, const function_t& h) const;

  function_t derivative(const function_t& f) const;
  //  function_t second_derivative(const function_t& f) const;
};


#endif
