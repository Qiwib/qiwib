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
     
     DEFUN_DLD (Calc_F_phi_C, args, , "No help available")
     {
       mlock(); 
       int nargin = args.length ();
       octave_value_list retval;
       if (nargin != 8)
         print_usage ();
       else
         {
	   ComplexColumnVector psi = args(0).complex_column_vector_value ();
	   ComplexColumnVector H_sp = args(1).complex_column_vector_value ();
	   ComplexColumnVector H_nonlin = args(2).complex_column_vector_value ();
	   Complex H_d = args(3).complex_value ();
	   int Ng = args(4).int_value ();
	   int M = args(5).int_value ();
	   double g = args(6).double_value ();
	   double dx = args(7).double_value ();
	
	   ComplexMatrix Ovrlp_phi_inv(M,M);
   	   ComplexColumnVector F(M*Ng);
   	   Complex H_lin_diag[Ng];
   	   Complex ovrlp = 0;
   	   Complex temp = 0;
	   double temp_rH = 0;
	   double temp_iH = 0;
	   double temp_rW = 0;
	   double temp_iW = 0;	   
	   double temp_r = 0;
	   double temp_i = 0;
   	   Complex H_nlC[Ng*M];
   	   Complex phi[Ng*M];
   	   Complex phi_cc[Ng*M];
   	   double H_nonlinC[2*M*M*M*M];
	   double W_sl[2*Ng*M*M];
   	   Complex FC[Ng*M];
   	   Complex Ovrlp_phiC[M][M];
   	   Complex Ovrlp_phi_invC[M][M];
   	   Complex H_linC = 0;
   	   Complex H_linCr = 0;
   	   Complex H_linCl = 0;
	   Complex H_linCr_periodic = 0;
   	   Complex H_linCl_periodic = 0;  	   
    	   Complex H_linCrr = 0;
   	   Complex H_linCll = 0;
	   Complex H_linCrr_periodic = 0;
   	   Complex H_linCll_periodic = 0;
	   
           if (! error_state)
           	{
			
		//for(int n=0; n<M*M*M*M; n++) H_nonlinC[n]=H_nonlin(n);	   
		for(int j=0; j<M; j++)
		        for(int q=0; q<M; q++)
				for(int s=0; s<M; s++)
					for(int l=0; l<M; l++)
					{
					  H_nonlinC[0+l+s*M+2*q*M*M+2*j*M*M*M]=real(H_nonlin(l+q*M+s*M*M+j*M*M*M));
					  H_nonlinC[M*M+l+s*M+2*q*M*M+2*j*M*M*M]=imag(H_nonlin(l+q*M+s*M*M+j*M*M*M));
					}
		
					  
           	for(int n=0; n<Ng*M; n++)
           	{
				H_nlC[n]=0;				
				phi[n]=psi(n);
				phi_cc[n]=conj(psi(n));
				FC[n]=0;
		}
		for(int n=0; n<Ng; n++) H_lin_diag[n] = H_sp(n);
		H_linCr = H_sp(Ng);
		H_linCl = H_sp(Ng+1);
		H_linCr_periodic = H_sp(Ng+2);
		H_linCl_periodic = H_sp(Ng+3);
		H_linCrr = H_sp(Ng+4);
		H_linCll = H_sp(Ng+5);
		H_linCrr_periodic = H_sp(Ng+6);
		H_linCll_periodic = H_sp(Ng+7);		
		for(int n=0; n<Ng; n++)
			for(int s=0; s<M; s++)
				for(int l=0; l<M; l++)
				{
				  W_sl[0+l+s*M+2*n*M*M] = real(phi_cc[n+s*Ng] * phi[n+l*Ng]);
				  W_sl[M*M+l+s*M+2*n*M*M] = imag(phi_cc[n+s*Ng] * phi[n+l*Ng]);
				}
		
	        for(int j=0; j<M; ++j)
		        for(int q=0; q<M; ++q)
			        for(int n=0; n<Ng; ++n)
			        {
					temp_r = 0;
					temp_i = 0;
					for(int ls=0; ls<M*M; ++ls)
					{
					  temp_rH = H_nonlinC[ls+2*q*M*M+2*j*M*M*M];
					  temp_iH = H_nonlinC[M*M+ls+2*q*M*M+2*j*M*M*M];
					  temp_rW = W_sl[ls+2*n*M*M];
					  temp_iW = W_sl[M*M+ls+2*n*M*M];
					  temp_r += temp_rH * temp_rW - temp_iH * temp_iW;
					  temp_i += temp_rH * temp_iW + temp_iH * temp_rW;
					}
					H_nlC[n+j*Ng] += Complex(temp_r,temp_i) * phi[n+q*Ng];
				}	
					
		      	 	
           	 for(int j=0; j<M; j++)
		{
			
			H_linC = H_lin_diag[0];
			FC[0+j*Ng] = H_linC * phi[0+j*Ng];
			FC[0+j*Ng] += H_linCr * phi[1+j*Ng] + H_linCl_periodic * phi[Ng-1+j*Ng];
			FC[0+j*Ng] += H_linCrr * phi[2+j*Ng] + H_linCll_periodic * phi[Ng-2+j*Ng];			
			FC[0+j*Ng] = H_d * ( FC[0+j*Ng] + g * H_nlC[0+j*Ng] );
		
			H_linC = H_lin_diag[1];
			FC[1+j*Ng] = H_linC * phi[1+j*Ng];
			FC[1+j*Ng] += H_linCr * phi[2+j*Ng] + H_linCl * phi[0+j*Ng];
			FC[1+j*Ng] += H_linCrr * phi[3+j*Ng] + H_linCll_periodic * phi[Ng-1+j*Ng];			
			FC[1+j*Ng] = H_d * ( FC[1+j*Ng] + g * H_nlC[1+j*Ng] );			
			
			H_linC = H_lin_diag[Ng-1];
			FC[Ng-1+j*Ng] = H_linC * phi[Ng-1+j*Ng];
			FC[Ng-1+j*Ng] += H_linCr_periodic * phi[0+j*Ng] + H_linCl * phi[Ng-2+j*Ng];
			FC[Ng-1+j*Ng] += H_linCrr_periodic * phi[1+j*Ng] + H_linCll * phi[Ng-3+j*Ng];
			FC[Ng-1+j*Ng] = H_d * ( FC[Ng-1+j*Ng] + g * H_nlC[Ng-1+j*Ng] );
			
			H_linC = H_lin_diag[Ng-2];
			FC[Ng-2+j*Ng] = H_linC * phi[Ng-2+j*Ng];
			FC[Ng-2+j*Ng] += H_linCr * phi[Ng-1+j*Ng] + H_linCl * phi[Ng-3+j*Ng];
			FC[Ng-2+j*Ng] += H_linCrr_periodic * phi[0+j*Ng] + H_linCll * phi[Ng-4+j*Ng];
			FC[Ng-2+j*Ng] = H_d * ( FC[Ng-2+j*Ng] + g * H_nlC[Ng-2+j*Ng] );
			
			for(int n=2; n<Ng-2; n++)
				{
				H_linC = H_lin_diag[n];
				FC[n+j*Ng] = H_linC * phi[n+j*Ng];
				FC[n+j*Ng] += H_linCr * phi[1+n+j*Ng] + H_linCl * phi[-1+n+j*Ng];
				FC[n+j*Ng] += H_linCrr * phi[2+n+j*Ng] + H_linCll * phi[-2+n+j*Ng];
				FC[n+j*Ng] = H_d * ( FC[n+j*Ng] + g * H_nlC[n+j*Ng] );
			}
           	 		
           	 }
           	 		
			
		for(int k=0; k<M; k++)
			for(int q=0; q<M; q++)
				{
				Ovrlp_phiC[k][q] = 0;
				for(int n=0; n<Ng; n++) Ovrlp_phiC[k][q] += phi_cc[n+k*Ng] * phi[n+q*Ng] * dx;
				}
	           	   	
	        for(int k=0; k<M; k++)
			   for(int q=0; q<M; q++) Ovrlp_phi_inv(k,q) = Ovrlp_phiC[k][q];
	        Ovrlp_phi_inv = Ovrlp_phi_inv.inverse();	        
	        for(int k=0; k<M; k++)
			   for(int q=0; q<M; q++) Ovrlp_phi_invC[k][q] = Ovrlp_phi_inv(k,q);
			           	
			           	          	 	
		for(int j=0; j<M; j++)
			for(int q=0; q<M; q++)
			{     	
				ovrlp = 0;
				for(int n=0; n<Ng; n++) ovrlp += phi_cc[n+q*Ng] * FC[n+j*Ng];
				ovrlp = ovrlp*dx;
				for(int k=0; k<M; k++)
					for(int n=0; n<Ng; n++) FC[n+j*Ng] -= phi[n+k*Ng] * ovrlp * Ovrlp_phi_invC[k][q];
	           	 }
        	 
	        for(int n=0; n<Ng*M; n++) F(n)=FC[n];
	           	 
	           	 
	          //for(int j=0; j<M; j++)
			   //for(int q=0; q<M; q++)
				//{     	
				//ovrlp = 0;
				//for(int n=0; n<Ng; n++) ovrlp = ovrlp + conj(psi(n+q*Ng)) * F(n+j*Ng);
				//ovrlp = ovrlp;
				//temp = 0;
				//for(int n=0; n<Ng; n++) temp = temp + conj(psi(n+q*Ng)) * psi(n+q*Ng);
				//ovrlp = ovrlp/temp;
	           	//for(int n=0; n<Ng; n++) F(n+j*Ng) = F(n+j*Ng) - psi(n+q*Ng) * ovrlp;
	           	 //}
	           	 
	       	  //for(int j=0; j<M; j++)
			   //for(int q=0; q<M; q++)
				//{     	
				//ovrlp = 0;
				//for(int n=0; n<Ng; n++) ovrlp = ovrlp + conj(psi(n+q*Ng)) * F(n+j*Ng);
				//ovrlp = ovrlp;
				//temp = 0;
				//for(int n=0; n<Ng; n++) temp = temp + conj(psi(n+q*Ng)) * psi(n+q*Ng);
				//ovrlp = ovrlp/temp;
	           	//for(int n=0; n<Ng; n++) F(n+j*Ng) = F(n+j*Ng) - psi(n+q*Ng) * ovrlp;
	           	 //}    	 
	           	 
           	 	      	 	
           	 	
           	 retval(0) = F;
           	}
         }
         
      return retval;
     }

