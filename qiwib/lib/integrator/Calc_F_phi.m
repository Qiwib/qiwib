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

	
function F = Calc_F_phi(one,two)
mlock(); global pa space realgrid realfunction realbasis complexgrid complexfunction complexbasis phiCpp VCpp

	if strcmp(pa.ode_phi, 'adams') || strcmp(pa.ode_phi, 'bdf') || strcmp(pa.ode_phi, 'stiff') || strcmp(pa.ode_phi, 'non-stiff')
		Fb = F_function_phi(two,one(1:pa.M*pa.Ng)+I*one(pa.M*pa.Ng+1:2*pa.M*pa.Ng));
		F = [real(Fb);imag(Fb)];	
	else
		F = F_function_phi(one,two);
	end

function F_tmp = F_function_phi(t,psi)
mlock(); global pa space realgrid realfunction realbasis complexgrid complexfunction complexbasis phiCpp VCpp

	%F_tmp = Calc_F_phi_C(psi,pa.H_phi_lin,pa.H_phi_nl,pa.H_phi_direction,pa.Ng,pa.M,pa.g,pa.dx);
%	phiCpp = phiCpp.set_data_vector(psi);

%	FCpp = phiCpp.propagate(pa.H_phi_direction, pa.H_Diff, pa.H_Diff2, pa.g, VCpp, pa.H_phi_nl, eye(pa.M));

	phiCpp_temp = phiCpp.set_data_vector(complex(psi));
	FCpp = phiCpp_temp.propagate(complex(pa.H_phi_direction), complex(pa.H_Diff), complex(pa.H_Diff2), complex(pa.g), VCpp, complex(pa.H_phi_nl), complex(inv(phiCpp_temp.overlap_matrix())) );
	F_tmp = FCpp.get_data_vector();