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

function calc_fields()
mlock(); global pa

	al = pa.H_Diff/12/pa.dx;
	al2 = pa.H_Diff2/12/pa.dx/pa.dx;
	
	%% create h_kq
	if strcmp(pa.boundary,'periodic')

		psid = [pa.phi(pa.Ng,:);pa.phi(1:pa.Ng-1,:)];
		psiu = [pa.phi(2:pa.Ng,:);pa.phi(1,:)];
		psidd = [ pa.phi(pa.Ng-1:pa.Ng,:) ; pa.phi(1:pa.Ng-2,:) ];
		psiuu = [ pa.phi(3:pa.Ng,:) ; pa.phi(1:2,:) ];		

	elseif strcmp(pa.boundary,'box')
				
		psid = [zeros(1,pa.M);pa.phi(1:pa.Ng-1,:)];
		psiu = [pa.phi(2:pa.Ng,:);zeros(1,pa.M)];
		psidd = [ zeros(2,pa.M) ; pa.phi(1:pa.Ng-2,:) ];
		psiuu = [ pa.phi(3:pa.Ng,:) ; zeros(2,pa.M) ];
			
	end
	for k=1:pa.M, for q=1:pa.M
			
		psi20 = pa.phi(:,k)' * pa.phi(:,q);
		psi2d = pa.phi(:,k)'*psid(:,q);
		psi2u = pa.phi(:,k)'*psiu(:,q);
		psi2dd = pa.phi(:,k)'*psidd(:,q);
		psi2uu = pa.phi(:,k)'*psiuu(:,q);
		pa.h_kq(q+(k-1)*pa.M) =  -30.0 * al2 * psi20 + (16*al2+8*al) * psi2u + (16*al2-8*al) * psi2d + (-al2-al) * psi2uu + (-al2+al) * psi2dd;
		pa.h_kq(q+(k-1)*pa.M) = pa.h_kq(q+(k-1)*pa.M) + pa.phi(:,k)'*(pa.V .* pa.phi(:,q));
		
	end, end
	
	pa.h_kq = pa.h_kq * pa.dx;

	%% create rho_ksql
	pa.w_ksql = calc_fields_C(pa.Ng,pa.M,pa.phi);
	pa.w_ksql = pa.g * pa.dx * pa.w_ksql;
