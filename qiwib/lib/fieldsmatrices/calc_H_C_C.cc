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
     	 
     DEFUN_DLD (calc_H_C_C, args, , "No help available")
     {
       mlock(); 
       int nargin = args.length ();
       octave_value_list retval;
       if (nargin != 7)
         print_usage ();
       else
         {
	   Matrix basis = args(0).matrix_value ();
	   Matrix basis_diff = args(1).matrix_value ();
	   ComplexColumnVector h = args(2).complex_column_vector_value ();
	   ComplexColumnVector w = args(3).complex_column_vector_value ();
	   int M = args(4).int_value ();
	   int nmax_h = args(5).int_value ();
	   double eps = args(6).double_value ();
   	   
	   Complex temp = 0;
	   Complex temp2 = 0;
	   int nindex[M];
	   int k = 0;
	   int s = 0;
	   int q = 0;
	   int l = 0;
	   int nr = 0;
	   int nl = 0;
	   double basisnr[M];
	   ComplexMatrix H_C(3,nmax_h);
   	    
           if (! error_state)
           	{ 	
				
				for(int n=0; n<nmax_h; n++)
				{

						nr = basis_diff(1,n);
						nl = basis_diff(2,n);
						for(int i=0; i<M; i++)
						{
							nindex[i] = int(basis_diff(i+3,n));
							basisnr[i] = basis(nr,i);
						}
						temp = 0;
						temp2 = 0;
						
							switch(int(basis_diff(0,n)))
							{
															
								case 0:
									for(int j=0; j<M; j++)
									{
										temp2 = h(j*M+j) + 0.5 * (basisnr[j]-1.) * w(j*M*M*M+j*M*M+j*M+j);
										for(int i=0; i<j; i++)
											temp2 += 0.5 * basisnr[i] * 2.0 * w(j*M*M*M+i*M*M+j*M+i);
										for(int i=j+1; i<M; i++)
											temp2 += 0.5 * basisnr[i] * 2.0 * w(j*M*M*M+i*M*M+j*M+i);
										temp = temp + basisnr[j] * temp2;											
									}
									break;

								case 21:
									k = nindex[0]; q = nindex[M-1];
									temp = h(k*M+q) + basisnr[k] * w(k*M*M*M+k*M*M+k*M+q);
									temp += (basisnr[q]-1.) * w(k*M*M*M+q*M*M+q*M+q);
									for(int i=1; i<M-1; i++)
									{
										s = nindex[i];
										temp = temp + basisnr[s] * 2.0 * w(k*M*M*M+s*M*M+q*M+s);
									}
									temp = sqrt( (basisnr[k]+1.)*basisnr[q] ) * temp;
									break;
									
								case 22:			
									k = nindex[0]; q = nindex[M-1];
									temp = (basisnr[k]+1.) * (basisnr[k]+2.) * (basisnr[q]-1.) * basisnr[q];
									temp = 0.5 * sqrt(temp) * w(k*M*M*M+k*M*M+q*M+q);
									break;
													
								case 31:
									k = nindex[0]; s = nindex[1]; q = nindex[M-1];
									temp = (basisnr[k]+1.) * (basisnr[s]+1.) * (basisnr[q]-1.) * basisnr[q];
									temp = sqrt(temp) * w(k*M*M*M+s*M*M+q*M+q);
									break;
																	
								case 32:
									k = nindex[0]; q = nindex[M-2]; l = nindex[M-1];
									temp = (basisnr[k]+1.) * (basisnr[k]+2.) * basisnr[q] * basisnr[l];
									temp = sqrt(temp) * w(k*M*M*M+k*M*M+q*M+l);
									break;
																					
								case 41:
									k = nindex[0]; s = nindex[1]; q = nindex[M-2]; l = nindex[M-1];
									temp = (basisnr[k]+1.) * (basisnr[s]+1.) * basisnr[q] * basisnr[l];
									temp = sqrt(temp) * 2.0 * w(k*M*M*M+s*M*M+q*M+l);		
									break;
																
							}
							H_C(0,n) = nl+1;
							H_C(1,n) = nr+1;
							if (abs(temp)>eps) {H_C(2,n) = temp;}					
							else {H_C(2,n) = 0;}
				}
			
           	 retval(0) = H_C;
           	}
         }

      return retval;
     }

