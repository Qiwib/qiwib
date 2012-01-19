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
	 
     DEFUN_DLD (Calc_H_phi_nl_C, args, , "No help available")
     {
       mlock(); 
       int nargin = args.length ();
       octave_value_list retval;
       if (nargin != 3)
         print_usage ();
       else
         {
		   int M = args(0).int_value ();
		   ComplexColumnVector rho_ksql = args(1).complex_column_vector_value ();
		   ComplexColumnVector rho_inv = args(2).complex_column_vector_value ();
		   
		   ComplexColumnVector H_int(M*M*M*M);
		   Complex temp = 0;
		   
			if (! error_state)
           	{ 
				
				for(int j=0; j<M; j++)	
					for(int s=0; s<M; s++)
						for(int q=0; q<M; q++)
							for(int l=0; l<M; l++)
							{
								temp = 0;
								for(int k=0; k<M; k++)
									temp = temp + rho_inv(k+j*M) * rho_ksql(l+q*M+s*M*M+k*M*M*M);
								H_int(l+q*M+s*M*M+j*M*M*M) = temp;
							}
	           	 
	           	retval(0) = H_int;
           	}
         }
      return retval;
     }
