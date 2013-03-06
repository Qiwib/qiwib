/*
 * Copyright (C) 2012 by Thomas Ernst
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
*/


     #include <octave/oct.h>
     
     DEFUN_DLD (calc_rho_C, args, , "No help available")
     {
       mlock(); 
       int nargin = args.length ();
       octave_value_list retval;
       if (nargin != 5)
         print_usage ();
       else
         {
	   Matrix basis = args(0).matrix_value ();
	   Matrix basis_diff = args(1).matrix_value ();
	   ComplexColumnVector C = args(2).complex_column_vector_value ();
	   int M = args(3).int_value ();
	   int nmax_h = args(4).int_value ();
	   
	   Complex temp = 0;
	   Complex temp2 = 0;
	   Complex CC;
	   int nCr = 0;
	   int nCl = 0;
	   int nindex[M];
	   int k = 0;
	   int s = 0;
	   int q = 0;
	   int l = 0;
	   double basisnCl[M];
	   ComplexColumnVector rho_kq(M*M);
	   ComplexColumnVector rho_ksql(M*M*M*M);
	   //const SparseMatrix basis(basis1);
   	    
           if (! error_state)
           	{ 	
        		for(int i=0; i<M*M; i++) rho_kq(i) = 0;
				for(int i=0; i<M*M*M*M; i++) rho_ksql(i) = 0;
				
				
				for(int n=0; n<nmax_h; n++)
				{
					// nCr and nCl are the indices of the bra and ket of the n-th non zero element H_ijksql
					nCr = basis_diff(2,n); // this is correct as it is supposed to be the other way round
					nCl = basis_diff(1,n);
					for(int i=0; i<M; i++)
					{
						// basisnCl is the state in Fock representation pointed by nCl
						basisnCl[i] = basis(nCl,i);
						
						// nindex stores the index of the orbitals ordered according to the occupation number difference
						// between the ket and the bra (the order is from the largest difference to the smallest)
						for(int i=0; i<M; i++) nindex[i] = int(basis_diff(i+3,n));
					}
					temp = 0;
					temp2 = 0;
	
					// multiplication of the expansion coefficient for the left and right side
					CC = conj(C(nCl)) * C(nCr);

							// and here we multiply CC by the result of < basisnCl | a^{\dagger}_{k} a^{\dagger}_{s} a_{q} a_{l} | basisnCr >
							switch(int(basis_diff(0,n)))
							{
															
								case 0:
									for(k=0; k<M; k++)
									{
										// diagonal terms
										rho_kq(k*M+k) = rho_kq(k*M+k) + CC * basisnCl[k];
										rho_ksql(k*M*M*M+k*M*M+k*M+k) = rho_ksql(k*M*M*M+k*M*M+k*M+k) + CC * basisnCl[k]*(basisnCl[k]-1);
										
										// terms: a^{\dagger}_{k} a^{\dagger}_{s} a_{s} a_{k}
										// we split the loop over s in two because we don't want to include the case s = k which we already considered
										for(s=0; s<k; s++)
										{	
											temp = CC * basisnCl[k] * basisnCl[s];
											rho_ksql(k*M*M*M+s*M*M+k*M+s) = rho_ksql(k*M*M*M+s*M*M+k*M+s) + temp;
											// symmetry of the density matrix
											rho_ksql(k*M*M*M+s*M*M+s*M+k) = rho_ksql(k*M*M*M+s*M*M+k*M+s);
										}
										for(s=k+1; s<M; s++)
										{	
											temp = CC * basisnCl[k] * basisnCl[s];
											rho_ksql(k*M*M*M+s*M*M+k*M+s) = rho_ksql(k*M*M*M+s*M*M+k*M+s) + temp;
											// symmetry of the density matrix
											rho_ksql(k*M*M*M+s*M*M+s*M+k) = rho_ksql(k*M*M*M+s*M*M+k*M+s);
										}										
									}
									break;

								case 21:
									k = nindex[M-1]; s = nindex[0]; l = nindex[0]; q = nindex[0];
									
									// this rely on the fact that q is going to span on every index 0:M-1 different from k
									// a^{\dagger}_{k}a_{q} acting on the left
									rho_kq(k*M+q) = rho_kq(k*M+q) + CC * sqrt( basisnCl[k]*(basisnCl[q]+1) );
									
									// a^{\dagger}_{k} a^{\dagger}_{k} a_{l} a_{k}	
									temp = CC * (basisnCl[k]-1) * sqrt( basisnCl[k]*(basisnCl[l]+1) );
									rho_ksql(k*M*M*M+k*M*M+k*M+l) = rho_ksql(k*M*M*M+k*M*M+k*M+l) + temp;
									rho_ksql(k*M*M*M+k*M*M+l*M+k) = rho_ksql(k*M*M*M+k*M*M+k*M+l);
									
									// a^{\dagger}_{s} a^{\dagger}_{k} a_{s} a_{s}
									temp = CC * basisnCl[s] * sqrt( basisnCl[k]*(basisnCl[s]+1) );
									rho_ksql(k*M*M*M+s*M*M+s*M+s) = rho_ksql(k*M*M*M+s*M*M+s*M+s) + temp;
									rho_ksql(s*M*M*M+k*M*M+s*M+s) = rho_ksql(k*M*M*M+s*M*M+s*M+s);
									
									for(int i=1; i<M-1; i++)
									{
										s = nindex[i];
										// a^{\dagger}_{k} a^{\dagger}_{s} a_{s} a_{l}
										temp = CC * basisnCl[s] * sqrt( basisnCl[k]*(basisnCl[l]+1) );
										rho_ksql(k*M*M*M+s*M*M+s*M+l) = rho_ksql(k*M*M*M+s*M*M+s*M+l) + temp;
										rho_ksql(k*M*M*M+s*M*M+l*M+s) = rho_ksql(k*M*M*M+s*M*M+s*M+l);
										rho_ksql(s*M*M*M+k*M*M+s*M+l) = rho_ksql(k*M*M*M+s*M*M+s*M+l);
										rho_ksql(s*M*M*M+k*M*M+l*M+s) = rho_ksql(k*M*M*M+s*M*M+s*M+l);
									}
									break;
									
								case 22:			
									k = nindex[M-1]; q = nindex[0];
									// a^{\dagger}_{k} a^{\dagger}_{k} a_{q} a_{q}
									temp = CC * sqrt( (basisnCl[k]-1)*basisnCl[k]*(basisnCl[q]+1)*(basisnCl[q]+2) );
									rho_ksql(k*M*M*M+k*M*M+q*M+q) = rho_ksql(k*M*M*M+k*M*M+q*M+q) + temp;
									break;
													
								case 31:
									k = nindex[M-1]; q = nindex[0]; l = nindex[1];
									// a^{\dagger}_{k} a^{\dagger}_{k} a_{q} a_{l}
									temp = CC * sqrt( (basisnCl[k]-1)*basisnCl[k]*(basisnCl[q]+1)*(basisnCl[l]+1) );
									rho_ksql(k*M*M*M+k*M*M+q*M+l) = rho_ksql(k*M*M*M+k*M*M+q*M+l) + temp;
									rho_ksql(k*M*M*M+k*M*M+l*M+q) = rho_ksql(k*M*M*M+k*M*M+q*M+l);
									break;
																	
								case 32:
									k = nindex[M-2]; s = nindex[M-1]; q = nindex[0];
									// a^{\dagger}_{k} a^{\dagger}_{s} a_{q} a_{q}
									temp = CC * sqrt( basisnCl[k]*basisnCl[s]*(basisnCl[q]+1)*(basisnCl[q]+2) );
									rho_ksql(k*M*M*M+s*M*M+q*M+q) = rho_ksql(k*M*M*M+s*M*M+q*M+q) + temp;
									rho_ksql(s*M*M*M+k*M*M+q*M+q) = rho_ksql(k*M*M*M+s*M*M+q*M+q);
									break;
																					
								case 41:
									// a^{\dagger}_{k} a^{\dagger}_{s} a_{q} a_{l}
									k = nindex[M-2]; s = nindex[M-1]; q = nindex[0]; l = nindex[1];
									temp = CC * sqrt( basisnCl[k]*basisnCl[s]*(basisnCl[q]+1)*(basisnCl[l]+1) );
									rho_ksql(k*M*M*M+s*M*M+q*M+l) = rho_ksql(k*M*M*M+s*M*M+q*M+l) + temp;
									rho_ksql(k*M*M*M+s*M*M+l*M+q) = rho_ksql(k*M*M*M+s*M*M+q*M+l);
									rho_ksql(s*M*M*M+k*M*M+q*M+l) = rho_ksql(k*M*M*M+s*M*M+q*M+l);
									rho_ksql(s*M*M*M+k*M*M+l*M+q) = rho_ksql(k*M*M*M+s*M*M+q*M+l);
									break;
																
							}
				}
						
							
           	 retval(0) = rho_kq;
           	 retval(1) = rho_ksql;
           	}
         }	
         
      return retval;
     }

