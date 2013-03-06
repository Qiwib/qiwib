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

// This creates the basis of the Hilbert space
// This can be made more efficient see: European Journal of Physics, 31(3), 591–602.
// You can create the basis of the Hilbert space already in the right order saving the ordering call later

     #include <octave/oct.h>
	 
     DEFUN_DLD (calc_hilbert_C, args, , "No help available")
     {
       mlock(); 
       int nargin = args.length ();
       octave_value_list retval;
       if (nargin != 3)
         print_usage ();
       else
         {
		   int N = args(0).int_value (); 
		   int M = args(1).int_value (); 
		   int nmax = args(2).int_value ();
		   
		   Matrix basis(nmax,M);
		   int a[N];
		   int b[M][M];
		   int row = 0;
		   int basis_len = 0;
		   
           if (! error_state)
           	{ 
					
	           	for(int i=0; i<nmax; i++)
						for(int j=0; j<M; j++) basis(i,j) = 0;
	           	for(int i=1; i<N; i++) a[i] = 1;
	           	a[0] = M;
	           	
	           	while(a[N-1]<M)
	           	{
					for(int i=0; i<a[0]; i++)
						for(int j=0; j<M; j++) b[i][j] = 0;
					
					for(int i=0; i<a[0]; i++) b[i][i+(M-a[0])] = 1;
					
					for(int s=1; s<N; s++)
					{
						row = a[s]-1;
						for(int i=0; i<a[0]; i++) b[i][row] += 1;
					}
					for(int i=0; i<a[0]; i++)
						for(int j=0; j<M; j++) basis(i+basis_len,j) = b[i][j];
						
					basis_len += a[0];
					
					a[1]++;
					
					for(int r=1; r<N-1; r++)
						if(a[r]==M+1)
						{
							a[r+1] += 1;
							for(int i=1; i<r+1; i++) a[i] = a[r+1];
						}
					a[0] = M - a[1] + 1;		
				}
				basis(basis_len,M-1) = N;
	           	
	           	 
	           	retval(0) = basis;
           	}
         }		   

      return retval;
     }
