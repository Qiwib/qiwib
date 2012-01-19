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


function Calc_H_phi_lin()
mlock(); global pa

%% fourth order derivatives!!!

	al = pa.H_Diff/12/pa.dx;
	al2 = pa.H_Diff2/12/pa.dx/pa.dx;
	
	pa.H_phi_lin = zeros(pa.Ng+8,1);
	
	
	pa.H_phi_lin(1:pa.Ng) = -30.0*al2*ones(pa.Ng,1) + pa.V.*ones(pa.Ng,1);
	pa.H_phi_lin(pa.Ng+1) = 16*al2+8*al;
	pa.H_phi_lin(pa.Ng+2) = 16*al2-8*al;
	pa.H_phi_lin(pa.Ng+5) = -al2-al;
	pa.H_phi_lin(pa.Ng+6) = -al2+al;

		
	if strcmp(pa.boundary,'periodic')
		pa.H_phi_lin(pa.Ng+3) = 16*al2+8*al;
		pa.H_phi_lin(pa.Ng+4) = 16*al2-8*al;
		pa.H_phi_lin(pa.Ng+7) = -al2-al;
		pa.H_phi_lin(pa.Ng+8) = -al2+al;
	end

