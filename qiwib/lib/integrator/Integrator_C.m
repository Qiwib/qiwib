## Copyright (C) 2012 by Thomas Ernst
## 
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
## THE SOFTWARE.


function Cx = Integrator_C(dt,t,C)
mlock(); global pa
	
	switch pa.ode_C
	
		case {'RK23'}

			[tarray,v] = ode23( @Calc_F_C,[t t+dt],C, pa.ode_C_opts);
			Cx = v(end,:).';	
			
		case {'RK45'}

			[tarray,v] = ode45( @Calc_F_C,[t t+dt],C, pa.ode_C_opts);
			Cx = v(end,:).';
	
		case {'RK78'}

			[tarray,v] = ode78( @Calc_F_C,[t t+dt],C, pa.ode_C_opts);
			Cx = v(end,:).';
			
		case {'adams' 'bdf'}
			
			lsode_options("absolute tolerance", pa.ode_C_opts(1));
			lsode_options("relative tolerance", pa.ode_C_opts(2));
			lsode_options("maximum order", pa.ode_C_opts(3));
			lsode_options("integration method", pa.ode_C);
			v = lsode( @Calc_F_C,[real(C);imag(C)],[t t+dt]);
			Cx = (v(end,1:pa.nmax)+i*v(end,pa.nmax+1:2*pa.nmax)).';
			
		case {'arnoldi'}

			Cx = Arnoldi(pa.H_C_direction*pa.H_C,C,pa.ode_C_opts(1),pa.ode_C_opts(2),dt);	

		case {'lanczos'}

			Cx = Lanczos(pa.H_C,C,pa.ode_C_opts(1),pa.ode_C_opts(2),dt);
			
		otherwise
			disp('integrator for C does not exist!!!');
			
	end
		
