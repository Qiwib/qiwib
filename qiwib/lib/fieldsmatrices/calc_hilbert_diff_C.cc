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
				
				for(int nr=0; nr<nmax; nr++)//for(int nr=0; nr<nmax; nr++)
				{
					for(int i=0; i<M; i++) basisnr[i] = basis(nr,i);
					for(int nl=0; nl<nmax; nl++)
					{	

						difflnz = 0;
						diffmax = 0;
						for(int i=0; i<M; i++)
						{
							diff[i] = basis(nl,i) - basisnr[i];							
							if(diff[i]!=0) difflnz++;
							if(abs(diff[i])>diffmax) diffmax = abs(diff[i]);
						}
						
						if(difflnz<5 && diffmax<3)
						{

							for(int i=0; i<M; i++) nindex[i] = i;					
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
							
							temp = 10*difflnz+diff[0];

							switch(temp)
							{
								
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
													
								case 31:
									hilbert_diff(0,hilbert_nz) = 31;
									hilbert_diff(1,hilbert_nz) = nr;
									hilbert_diff(2,hilbert_nz) = nl;
									for(int i=3; i<M+3; i++)
									{
										hilbert_diff(i,hilbert_nz) = nindex[i-3];
									}
									hilbert_nz++;
									break;
																	
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
