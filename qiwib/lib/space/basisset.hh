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
  
  //  typedef ublas::hermitian_matrix<scalar_t,ublas::upper> hermitian_matrix;

  const Space&     space;
  //  hermitian_matrix overlap;

  basisset(const Space& space, std::vector<Function> basis) : std::vector<Function>(basis.begin(),basis.end()), 
							      space(space)
  {  }

  basisset(const Space& space, size_t n) : std::vector<Function>(n,space), space(space) {  }

  basisset(const basisset& basis) : std::vector<Function>(basis.begin(),basis.end()), space(basis.space) {  }

  // TODO: Make proper set of operations
  basisset operator *(const scalar_t& s){ 
    const basisset& phi(*this);
    basisset newbasis(space,phi.size());

    for(size_t i=0;i<phi.size();i++) 
      newbasis[i] = phi[i]*s;

    return newbasis;
  }

  basisset operator *(const Array2D<scalar_t>& C) const
  {
    assert(C.rows() == this->size());

    const basisset& phi(*this);
    basisset newbasis(space,this->size());

    for(size_t i=0;i<C.rows();i++){
      Function sum(space);
      for(size_t j=0;j<C.columns();j++)
	sum += phi[j]*C(i,j);

      newbasis[i] = sum;
    }
    
    return newbasis;
  }

  basisset derivative(bool periodic = true) const {
    const basisset& phi(*this);

    basisset derivative_basis(space, phi.size());
    
    for(size_t i=0;i<phi.size();i++)
      derivative_basis[i] = space.derivative(phi[i], periodic);

    return derivative_basis;
  }
  
  basisset second_derivative(bool periodic = true) const {
    const basisset& phi(*this);

    basisset derivative_basis(space, phi.size());
    
    for(size_t i=0;i<phi.size();i++)
      derivative_basis[i] = space.second_derivative(phi[i], periodic);

    return derivative_basis;
  }
  
  Array2D<scalar_t> laplacian_matrix() const {
    return derivative().overlap_matrix();
  }

  Array2D<scalar_t> overlap_matrix() const
  {
    Array2D<scalar_t> overlap(this->size(),this->size());
    unsigned int i = 0, j =0;

    for(const_iterator f(this->begin()); f!=this->end();f++,i++){
      j = i;
      for(const_iterator g(f); g != this->end(); g++,j++){
	const scalar_t& ol(space.inner(*f,*g));
	overlap(i,j) = ol;
	overlap(j,i) = Space::conj(ol);
      }
    }
    return overlap;
  }

  Array2D<scalar_t> overlap_matrix(const basisset& phi1) const
  {
    Array2D<scalar_t> overlap(this->size(),this->size());
    unsigned int i = 0, j =0;

    for(const_iterator f(this->begin()); f!=this->end();f++,i++){
      for(const_iterator g(phi1.begin()); g != phi1.end(); g++,j++){
	const scalar_t& ol(space.inner(*f,*g));
	overlap(i,j) = ol;
      }
    }
    return overlap;
  }

  basisset normalise() const
  {
    basisset phi_new(space,this->size());

    for(unsigned int i=0; i<this->size();i++){
      const scalar_t& ol(space.inner((*this)[i],(*this)[i]));
      phi_new[i] = (*this)[i]*pow(ol,-0.5);
    }
    return phi_new;
  }

  std::vector<scalar_t> norm() const
  {
    std::vector<scalar_t> norm(this->size());

    for(unsigned int i=0; i<this->size();i++){
      const scalar_t& ol(space.inner((*this)[i],(*this)[i]));
      norm[i] = pow(ol,0.5);
    }
    return norm;
  }

  Array2D<scalar_t> calc_error(basisset& phi1) const
  {
    Array2D<scalar_t> overlap(this->size(),this->size());
    basisset phi(space,this->size());
    const basisset& phi0(*this);
    std::vector<scalar_t> norm0(this->size());
    std::vector<scalar_t> norm1(this->size());

    norm0 = phi0.norm();
    norm1 = phi1.norm();

    for (unsigned int i=0;i<this->size();i++)
      phi[i] = phi0[i]/norm0[i] - phi1[i]/norm1[i];

    return phi.overlap_matrix();
  }

  std::vector<scalar_t> Wksql() const
  {
    std::vector<scalar_t> Wksql(this->size());
    unsigned int i=0;

    for(unsigned int k=0; k<this->size();k++)
      for(unsigned int s=0; s<this->size();s++)
	for(unsigned int q=0; q<this->size();q++)
	  for(unsigned int l=0; l<this->size();l++){
	    const scalar_t& ol(space.inner((*this)[k],(*this)[s],(*this)[q],(*this)[l]));
	    Wksql[i] = ol;
	    i++;
	  }
    return Wksql;
  }

  std::vector<scalar_t> hkq(const scalar_t& D, scalar_t& D2, std::vector<scalar_t>& V, bool periodic = true) const
  {
    std::vector<scalar_t> hkq(this->size());
    unsigned int i=0;
    const basisset& d1phi(this->derivative(periodic));
    const basisset& d2phi(this->second_derivative(periodic));
    
    for(unsigned int k=0; k<this->size();k++)
	for(unsigned int q=0; q<this->size();q++){
	    const scalar_t& ol1(space.inner(d1phi[k],(*this)[q]));
	    const scalar_t& ol2(space.inner(d2phi[k],(*this)[q]));
	    const scalar_t& olV(space.inner((*this)[k],V,(*this)[q]));
	    hkq[i] = D*ol1+D2*ol2+olV;
	    i++;
	}
    return hkq;
  }
  
};

#endif
