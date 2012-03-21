#ifndef BASISSET_HH
# define BASISSET_HH

#include <vector>
#include <iostream>
#include <assert.h>
#include "array.hh"

template <typename Function> class basisset : public std::vector<Function> {
public:
  typedef typename Function::space_t Space;
  typedef typename Space::value_t  value_t;
  typedef typename Space::scalar_t scalar_t;
  typedef typename std::vector<Function>::const_iterator const_iterator;
  typedef typename std::vector<Function>::iterator iterator; 
  
  const Space&     space;

  basisset(const Space& space, std::vector<Function> basis) : std::vector<Function>(basis.begin(),basis.end()), 
							      space(space)
  {  }

  basisset(const Space& space, size_t n) : std::vector<Function>(n,space), space(space) {  }
  basisset(const basisset& basis) : std::vector<Function>(basis.begin(),basis.end()), space(basis.space) {  }

  // TODO: Make proper set of operations
  basisset operator *(const scalar_t& s) const;
  friend basisset operator *(const scalar_t& s, const basisset& b)  { return b*s; }
  basisset operator *(const Array2D<scalar_t>& C) const;

  friend basisset operator *(const Array2D<scalar_t>& C, const basisset& phi) {
    assert(C.rows() == phi.size());

    basisset newbasis(phi.space,phi.size());

    for(size_t i=0;i<C.rows();i++){
      Function sum(phi.space);
      for(size_t j=0;j<C.columns();j++)
	sum += phi[j]*C(i,j);

      newbasis[i] = sum;
    }
    
    return newbasis;
  }

  basisset operator +(const basisset& A) const;

  basisset derivative() const;  
  basisset second_derivative() const;
  
  Array2D<scalar_t> laplacian_matrix() const { return derivative().overlap_matrix(); }

  Array2D<scalar_t> overlap_matrix() const;
  Array2D<scalar_t> overlap_matrix(const basisset& phi1) const;

  basisset normalise() const;
  std::vector<scalar_t> norm() const;

  Array2D<scalar_t> calc_error(basisset& phi1) const;

  std::vector<scalar_t> Wksql() const;
  std::vector<scalar_t> hkq(const scalar_t& D, const scalar_t& D2, const Function& V) const;

  std::vector<value_t> get_data_vector() const;

  basisset set_data_vector(const Array2D<value_t>& vect) const;

#if 0
  basisset propagate(const term_t& a0, const term_t& a1, const term_t& a2, const term_t& anl,
		     const scalar_t& direction, const Array2D<scalar_t>& H_nonlin, 
		     const Array2D<scalar_t>& overlap_matrix_inv) const;
#endif

  // TODO: Rearrange parameters
  basisset propagate(const scalar_t& direction, const scalar_t& a1, const scalar_t& a2, const scalar_t& g, const Function& a0, const Array2D<scalar_t>& H_nonlin, const Array2D<scalar_t>& overlap_matrix_inv) const;

  Array2D<value_t> Wsl() const;
  
  basisset nonlin(const Array2D<scalar_t>& H_nonlin) const;
  basisset orthonormalise_advanced(const Array2D<scalar_t>& overlap_matrix_inv, const basisset& phi) const;

};

#endif
