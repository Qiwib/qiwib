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

// This initialise the matrix dH_mn used to speed the computation of H_ijksql
// See Thomas Ernst thesis (section B.3.5)

#include <octave/oct.h>
	 
     DEFUN_DLD (calc_hilbert_diff_C, args, , "No help available")
     {
       mlock(); 
       int nargin = args.length ();
       octave_value_list retval;
       if (nargin != 3)
         print_usage ();
       else
         {
		Matrix basis = args(0).matrix_value (); 		
	   int M = args(1).int_value ();
	   int nmax = args(2).int_value ();
   	   
	   int temp = 0;
	   int difflnz = 0;
	   int diff[M];
	   int nindex[M];
	   int diffmax = 0;
	   int difftemp = 0;
	   int hilbert_nz = 0;
	   int hilbert_nz_max = nmax;
	   int basisnr[M];
	   Matrix hilbert_diff(M+3,nmax);
   	    
           if (! error_state)
           	{ 	
				
				// nr is an index (from 0 to nmax) labelling the ket
				for(int nr=0; nr<nmax; nr++)//for(int nr=0; nr<nmax; nr++)
				{
				
					// Copy over the basis vector relative to nr
					for(int i=0; i<M; i++) basisnr[i] = basis(nr,i);
					
					// nl is an index (from 0 to nmax) labelling the bra
					for(int nl=0; nl<nmax; nl++)
					{	

						difflnz = 0;
						diffmax = 0;
						for(int i=0; i<M; i++)
						{
							// Compute the difference in the occupation number of every orbital between bra and ket
							diff[i] = basis(nl,i) - basisnr[i];							
							if(diff[i]!=0) difflnz++; // difflnz is a measure of how many orbitals have different occupation numbers
							if(abs(diff[i])>diffmax) diffmax = abs(diff[i]); // diffmax is a measure of the maximum difference in occupation number
						}
						
						// The only non-zero elements in H_ijksql are the ones where at most four orbitals have different 
						// occupation number (difflnz<5) and at most two particles are shifted from one orbital (diffmax<3)
						if(difflnz<5 && diffmax<3)
						{

							for(int i=0; i<M; i++) nindex[i] = i;					
							
							// Here we sort the vector of the difference of the occupation number from the biggest to the smallest
							// we can keep track of the ordering index with the array nindex
							for(int i=0; i<M; i++)
								for(int j=i; j<M; j++)
								{							
									if(diff[j]>diff[i])
									{
										difftemp = diff[i];
										diff[i] = diff[j];
										diff[j] = difftemp;
										difftemp = nindex[i];
										nindex[i] = nindex[j];
										nindex[j] = difftemp;
									}
								}
							
							// temp is just a convenience index useful to tag the case we are based on the
							// values of difflnz and the highest difference in occupation number (diff[0])
							// Much like hashing 
							temp = 10*difflnz+diff[0];

							switch(temp)
							{
								// the structure of hilbert_diff is:
								// hilbert_diff(0,:) is the "temp" tag we compute before to distinguish the cases
								// hilbert_diff(1,:) is the index of the ket
								// hilbert_diff(2,:) is the index of the bra
								// hilbert_diff(3:M,:) is the list of orbitals with different occupation number from the largest difference to the smallest
								// the second dimension of hilbert_diff simply number the case of non-zero elements of H_ijksql
								// sequentially
							
								// this is the case where bra and ket are equal
								case 0:
									hilbert_diff(0,hilbert_nz) = 0;
									hilbert_diff(1,hilbert_nz) = nr;
									hilbert_diff(2,hilbert_nz) = nl;
									for(int i=3; i<M+3; i++)
									{
										hilbert_diff(i,hilbert_nz) = nindex[i-3];
									}
									hilbert_nz++;	
									break;

								// this is the case where one particle is moved from one orbital to the other
								case 21:
									hilbert_diff(0,hilbert_nz) = 21;
									hilbert_diff(1,hilbert_nz) = nr;
									hilbert_diff(2,hilbert_nz) = nl;
									for(int i=3; i<M+3; i++)
									{
										hilbert_diff(i,hilbert_nz) = nindex[i-3];
									}
									hilbert_nz++;								
									break;
									
								// this is the case where two particles are moved from one orbital to the other
								case 22:			
									hilbert_diff(0,hilbert_nz) = 22;
									hilbert_diff(1,hilbert_nz) = nr;
									hilbert_diff(2,hilbert_nz) = nl;
									for(int i=3; i<M+3; i++)
									{
										hilbert_diff(i,hilbert_nz) = nindex[i-3];
									}
									hilbert_nz++;
									break;
													
								// This term doesn't make any sense
								// case 31:
// 									hilbert_diff(0,hilbert_nz) = 31;
// 									hilbert_diff(1,hilbert_nz) = nr;
// 									hilbert_diff(2,hilbert_nz) = nl;
// 									for(int i=3; i<M+3; i++)
// 									{
// 										hilbert_diff(i,hilbert_nz) = nindex[i-3];
// 									}
// 									hilbert_nz++;
// 									break;
																
								// this is the case where two particles are moved from one orbital
								// into two different orbitals (parametric down conversion)
								case 32:
									hilbert_diff(0,hilbert_nz) = 32;
									hilbert_diff(1,hilbert_nz) = nr;
									hilbert_diff(2,hilbert_nz) = nl;
									for(int i=3; i<M+3; i++)
									{
										hilbert_diff(i,hilbert_nz) = nindex[i-3];
									}
									hilbert_nz++;
									break;
									
								// This is the case where we transfer one particle each from two
								// different orbitals and put the in another two different orbitals													
								case 41:
									hilbert_diff(0,hilbert_nz) = 41;
									hilbert_diff(1,hilbert_nz) = nr;
									hilbert_diff(2,hilbert_nz) = nl;
									for(int i=3; i<M+3; i++)
									{
										hilbert_diff(i,hilbert_nz) = nindex[i-3];
									}
									hilbert_nz++;		
									break;
																
							}
							
							// if we hit the maximum size of the hilbert_diff matrix just allocate some more
							if (hilbert_nz==hilbert_nz_max)
								{
									hilbert_diff.resize(M+3,2*hilbert_nz);
									hilbert_nz_max = 2*hilbert_nz;
								}				

						}
					
					}
				}
			
			 hilbert_diff.resize(M+3,hilbert_nz);
           	 retval(0) = hilbert_diff;
           	 retval(1) = hilbert_nz;
           	}
         }	  

      return retval;
     }
