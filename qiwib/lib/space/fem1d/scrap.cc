  ublas::vector<double> linear_solve(const ublas::banded_matrix<double>& A, 
			      const ublas::vector<double>& b) const
  {
    const unsigned int n = A.size1();

    const char UPLO('U');
    const int N(n-2), KD(1),NRHS(1), LDAB(2);
    int INFO(0);
    double cholesky[2*(n-1)];
    double rhs[n-2];
						
    for(size_t i=1;i<n-1;i++){
      cholesky[2*i-1] = A(i,i);      
      if(i>=2) cholesky[2*i-2] = A(i,i-1);
    }
    std::copy(b.begin(),b.end(),rhs);

    dpbtrf(&UPLO,&N,&KD,cholesky,&LDAB,&INFO);

    if(INFO){
      fprintf(stderr,"LAPACK dbptrf (symmetric banded Cholesky factorisation) failed with INFO = %d\n",INFO);
      abort();
    }
    
    dpbtrs(&UPLO,&N,&KD,&NRHS,cholesky,&LDAB,rhs,&N,&INFO);

    if(INFO){
      fprintf(stderr,"LAPACK dbptrs (symmetric banded linear solver) failed with INFO = %d\n",INFO);
      abort();
    }

    ublas::vector<double> result(n);
    result[0]   = 0;		// Homogeneous Dirichlet
    result[n-1] = 0;
    copy(rhs,rhs+n-2,result.begin()+1);

    return result;
  };
