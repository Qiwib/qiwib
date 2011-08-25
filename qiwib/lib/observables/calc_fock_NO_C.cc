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
	 
     DEFUN_DLD (calc_fock_NO_C, args, , "No help available")
     {
       mlock(); 
       int nargin = args.length ();
       octave_value_list retval;
       if (nargin != 7)
         print_usage ();
       else
         {
		   int N = args(0).int_value (); 
		   ComplexColumnVector C = args(1).complex_column_vector_value ();
		   Complex t_LL = args(2).complex_value();
		   Complex t_LR = args(3).complex_value();
		   Complex t_RL = args(4).complex_value();
		   Complex t_RR = args(5).complex_value();
		   Matrix bin_c = args(6).matrix_value ();
			
		   ComplexColumnVector C_NO(N+1);
		   Complex temp_C = 0;

		   
           if (! error_state)
           	{ 

				for (int mu=0; mu<=N; mu++)
				{
					temp_C = 0;
					for (int k=0; k<=mu; k++)
						for (int n=k; n<=N-mu+k; n++) temp_C += sqrt(bin_c(N,n)/bin_c(N,mu))*bin_c(n,k)*bin_c(N-n,mu-k) * pow(t_LL,n-k) * pow(t_RL,k) * pow(t_LR,N-n-mu+k) * pow(t_RR,mu-k) * C(n);

					C_NO(mu) = temp_C;
				}
           	 
	           	retval(0) = C_NO;
           	}
         }		   

      return retval;
     }

