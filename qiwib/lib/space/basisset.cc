#include "basisset.hh"

#define basisset_member(T) template <typename Function> T basisset<Function>::
#define basisset_t template <typename Function> basisset<Function> basisset<Function>::
#define scalar_t typename basisset<Function>::scalar_t
#define value_t typename basisset<Function>::value_t

basisset_t operator *(const scalar_t& s) const { 
  const basisset& phi(*this);
  basisset newbasis(space,phi.size());
  
  for(size_t i=0;i<phi.size();i++) 
    newbasis[i] = phi[i]*s;
  
  return newbasis;
}

basisset_t operator *(const Array2D<scalar_t>& C) const
{
  assert(C.rows() == this->size());

  const basisset& phi(*this);
  basisset newbasis(space,phi.size());

  for(size_t i=0;i<C.rows();i++){
    Function sum(space);
    for(size_t j=0;j<C.columns();j++)
      sum += phi[j]*C(j,i);

    newbasis[i] = sum;
  }
    
  return newbasis;
}

basisset_t operator +(const basisset& A) const
{
  basisset newbasis(space,this->size());

  for(unsigned int i=0;i<this->size();i++) newbasis[i] = (*this)[i]*A[i];  
    
  return newbasis;
}
  
basisset_t derivative() const {
  const basisset& phi(*this);

  basisset derivative_basis(space, phi.size());
    
  for(size_t i=0;i<phi.size();i++)
    derivative_basis[i] = space.derivative(phi[i]);

  return derivative_basis;
}
  
basisset_t second_derivative() const {
  const basisset& phi(*this);

  basisset derivative_basis(space, phi.size());
    
  for(size_t i=0;i<phi.size();i++)
    derivative_basis[i] = space.second_derivative(phi[i]);

  return derivative_basis;
}
  
basisset_member(Array2D<scalar_t>) overlap_matrix() const
{
  Array2D<scalar_t> overlap(this->size(),this->size());
  unsigned int i = 0, j =0;

  for(const_iterator f(this->begin()); f!=this->end();f++,i++){
    j = i;
    for(const_iterator g(f); g != this->end(); g++,j++){
      const scalar_t& ol(space.inner(*f,*g));
      overlap(i,j) = ol;
      overlap(j,i) = FieldTraits<scalar_t>::conj(ol);
    }
  }
  return overlap;
}

basisset_member(Array2D<scalar_t>) overlap_matrix(const basisset& phi1) const
{
  Array2D<scalar_t> overlap(this->size(),phi1.size());
  unsigned int i = 0, j =0;

  for(const_iterator f(this->begin()); f!=this->end();f++,i++){
    j=0;
    for(const_iterator g(phi1.begin()); g != phi1.end(); g++,j++){
      const scalar_t& ol(space.inner(*f,*g));
      overlap(i,j) = ol;
    }
  }
  return overlap;
}

basisset_t normalise() const
{
  basisset phi_new(space,this->size());

  for(unsigned int i=0; i<this->size();i++){
    const scalar_t& ol(space.inner((*this)[i],(*this)[i]));
    phi_new[i] = (*this)[i]*pow(ol,-0.5);
  }
  return phi_new;
}

basisset_member(std::vector<scalar_t>) norm() const
{
  std::vector<scalar_t> norm(this->size());

  for(unsigned int i=0; i<this->size();i++){
    const scalar_t& ol(space.inner((*this)[i],(*this)[i]));
    norm[i] = pow(ol,0.5);
  }
  return norm;
}

basisset_member(Array2D<scalar_t>) calc_error(basisset& phi1) const
{
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

basisset_member(std::vector<scalar_t>) Wksql() const
{
  const basisset& phi(*this);
  std::vector<scalar_t> Wksql(phi.size()*phi.size()*phi.size()*phi.size());
  unsigned int i=0;

  for(unsigned int k=0; k<phi.size();k++)
    for(unsigned int s=0; s<phi.size();s++)
      for(unsigned int q=0; q<phi.size();q++)
	for(unsigned int l=0; l<phi.size();l++){
	  const scalar_t& ol(space.inner(phi[k],phi[s],phi[q],phi[l]));
	  Wksql[i] = ol;
	  i++;
	}
  return Wksql;
}

basisset_member(std::vector<scalar_t>) hkq(const scalar_t& D, const scalar_t& D2, const Function& V) const
{
  unsigned int i=0;
  const basisset& phi(*this);
  const basisset& d1phi(phi.derivative());
  const basisset& d2phi(phi.second_derivative());
  std::vector<scalar_t> hkq(phi.size()*phi.size());
    
  for(unsigned int k=0; k<phi.size();k++)
    for(unsigned int q=0; q<phi.size();q++){
      const scalar_t& ol1(space.inner(phi[k],d1phi[q]));
      const scalar_t& ol2(space.inner(phi[k],d2phi[q]));
      const scalar_t& olV(space.inner(phi[k],V,phi[q]));
      hkq[i] = D*ol1+D2*ol2+olV;
      i++;
    }
  return hkq;
}

basisset_member(std::vector<value_t>) get_data_vector() const
{
  std::vector<value_t> vect(this->size()*(*this)[0].size());
  unsigned int i=0;

  for(unsigned int k=0; k<this->size();k++)
    for(unsigned int s=0; s<space.Nx;s++){
      vect[i] = (*this)[k][s];
      i++;
    }
  return vect;
}

basisset_t set_data_vector(const Array2D<value_t>& vect) const
{
  basisset phi(space,this->size());
  unsigned int i=0;

  for(unsigned int k=0; k<this->size();k++)
    for(unsigned int s=0; s<space.Nx;s++){
      phi[k][s] = vect[i];
      i++;
    }
  return phi;
}

basisset_t propagate(const term_t& a0, const term_t<scalar_t>& a1, const term_t<scalar_t>& a2, const term_t<scalar_t>& anl, const scalar_t<scalar_t>& direction, const Array2D<scalar_t>& H_nonlin, const Array2D<scalar_t>& overlap_matrix_inv) const
{
  const basisset& phi(*this);

  basisset F(space,phi.size());
  const basisset& d1phi(phi.derivative());
  const basisset& d2phi(phi.second_derivative());
  const basisset& Hnl(phi.nonlin(H_nonlin));
    
  for (unsigned int i=0;i<phi.size(); i++){
    F[i]  = a0*phi[i] + a1*d1phi[i] + a2*d2phi[i] + anl*Hnl[i];
    F[i] *= direction;
  }
  return F.orthonormalise_advanced(overlap_matrix_inv,phi);
}  


basisset_t propagate(const Function& a0, const scalar_t& a1, const scalar_t& a2, const scalar_t& anl, const scalar_t& direction, const Array2D<scalar_t>& H_nonlin, const Array2D<scalar_t>& overlap_matrix_inv) const
{
  const basisset& phi(*this);

  basisset F(space,phi.size());
  const basisset& d1phi(phi.derivative());
  const basisset& d2phi(phi.second_derivative());
  const basisset& Hnl(phi.nonlin(H_nonlin));
    
  for (unsigned int i=0;i<phi.size(); i++){
    F[i] = phi[i]*a0 + d1phi[i]*a1 + d2phi[i]*a2 + Hnl[i]*anl;
    F[i] *= direction;
  }
  return F.orthonormalise_advanced(overlap_matrix_inv,phi);
}  


  
basisset_member(Array2D<value_t>) Wsl() const
{
  const basisset& phi(*this);
  Array2D<value_t> Wsl(space.Nx,phi.size()*phi.size());
  unsigned int M(phi.size());    
		
  for(unsigned int s=0; s<M; s++)
    for(unsigned int l=0; l<M; l++)
      for(unsigned int x=0; x<space.Nx; x++)
	Wsl(x,l+s*M) = space.conj(phi[s][x]) * phi[l][x];
    
  return Wsl;
}  
  
basisset_t nonlin(const Array2D<scalar_t>& H_nonlin) const
{
  const basisset& phi(*this);
  Array2D<value_t> Wsl(this->Wsl());
  basisset H_nl(space,phi.size());
  unsigned int M(phi.size());

  // Default constructor should initialize to "0". 
  // for(unsigned int j=0; j<M; j++)
  // 	for(unsigned int x=0; x<space.Nx; x++)
  // 	  H_nl[j][x] = 0.0;
	  
  for(unsigned int j=0; j<M; j++)
    for(unsigned int q=0; q<M; q++)
      for(unsigned int x=0; x<space.Nx; x++){
	value_t sum(0);
	for(unsigned int s=0; s<M; s++)
	  for(unsigned int l=0; l<M; l++){
	    sum += Wsl(x,l+s*M) * H_nonlin[l+q*M+s*M*M+j*M*M*M];
	  }
	H_nl[j][x] += phi[q][x] * sum;
      }
    
  return H_nl;
}    

basisset_t orthonormalise_advanced(const Array2D<scalar_t>& overlap_matrix_inv, const basisset& phi) const
{
  basisset F(*this);
  unsigned int M = this->size();
  scalar_t ovrlp = 0;
    
  for(unsigned int j=0; j<M; j++)
    for(unsigned int q=0; q<M; q++)
      {     	
	ovrlp = space.inner(phi[q], F[j]);
	for(unsigned int k=0; k<M; k++)
	  F[j] -= phi[k] * ovrlp * overlap_matrix_inv(k,q);
      }
  return F;
}
