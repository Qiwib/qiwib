/*
 * Copyright (C) 2011 by Thomas Ernst
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
	 
     DEFUN_DLD (calc_fields_C, args, , "No help available")
     {
       int nargin = args.length ();
       octave_value_list retval;
       if (nargin != 3)
         print_usage ();
       else
         {
		   int Ng = args(0).int_value (); 
		   int M = args(1).int_value (); 
		   ComplexMatrix phi = args(2).complex_matrix_value ();
		   
		   ComplexColumnVector w(M*M*M*M);
		   Complex sum;
		   
			if (! error_state)
           	{ 
				
				for(int k=0; k<M; k++)	
					for(int s=0; s<M; s++)
						for(int q=0; q<M; q++)
							for(int l=0; l<M; l++)
							{
								sum = 0;
								for(int n=0; n<Ng; n++) sum += conj(phi(n,k)*phi(n,s)) * phi(n,q)*phi(n,l);
								w(l+q*M+s*M*M+k*M*M*M) = sum;
							}
	           	 
	           	retval(0) = w;
           	}
         }
      return retval;
     }
