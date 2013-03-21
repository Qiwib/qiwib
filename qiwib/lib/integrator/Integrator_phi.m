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


function phix = Integrator_phi(dt,t,phi)
mlock(); global pa


	switch pa.ode_phi
			
		case {'RK23'}
		
			[tarray,v] = ode23( @Calc_F_phi,[t t+dt],phi(:),pa.ode_phi_opts);
			phix = reshape(v(end,:).',pa.Ng,pa.M);

		case {'RK45'}
		
			[tarray,v] = ode45( @Calc_F_phi,[t t+dt],phi(:),pa.ode_phi_opts);
			phix = reshape(v(end,:).',pa.Ng,pa.M);		
			
		case {'RK78'}
			
			psi = phi.get_data_vector();
			[tarray,v] = ode78( @Calc_F_phi,[t t+dt],psi(:),pa.ode_phi_opts);
			psi = v(end,:).';
			phix = phi.set_data_vector(pa.scalar_type(psi));

		case {'adams' 'bdf'}
			lsode_options("absolute tolerance", pa.ode_phi_opts(1));
			lsode_options("relative tolerance", pa.ode_phi_opts(2));
			lsode_options("maximum order", pa.ode_phi_opts(3));
			lsode_options("integration method", pa.ode_phi);
			v = lsode( @Calc_F_phi,[real(phi(:));imag(phi(:))],[t t+dt]);
			phix = reshape( (v(end,1:pa.M*pa.Ng)+I*v(end,pa.M*pa.Ng+1:2*pa.M*pa.Ng)).',pa.Ng,pa.M);
										
		otherwise
			disp('integrator for phi does not exist!!!');
			
			
	end
	
