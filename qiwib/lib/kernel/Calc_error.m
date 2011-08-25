## Copyright (C) 2011 by Thomas Ernst
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

	
function delta = Calc_error(phi1,phi,C1,C)
mlock(); global pa

	for j=1:pa.M, phi_err(:,j) = (phi1(:,j) - phi(:,j)) / sqrt( sum( abs(phi(:,j)).^2 ) ); end
	for j=1:pa.M, for l=1:pa.M, H_O(j,l) = phi_err(:,j)'*phi_err(:,l); end, end
	d_phi = abs( trace(H_O * reshape(pa.rho_kq,pa.M,pa.M).') );

	d_C = (C1 - C) / sqrt( sum(abs(C).^2) );
	d_C = sum( abs(d_C).^2 );
	
	delta = 1/4*d_C + d_phi;