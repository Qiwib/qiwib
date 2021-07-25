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

	
function psix = Arnoldi(H,psi0,n,epsil=1e-6,dt)
mlock(); global pa

	n0 = n;
	psix = zeros(length(psi0),1);
	h = zeros(n,n);
	q=zeros(length(psi0),n);
	q(:,1) = psi0/sqrt(abs(psi0'*psi0));
	for k=2:n
		q(:,k) = H*q(:,k-1);
		for j=1:k-1
			h(j,k-1) = q(:,j)'*q(:,k);
			q(:,k) = q(:,k) - h(j,k-1)*q(:,j);
		end
		h(k,k-1) = sqrt( abs(q(:,k)'*q(:,k)) );
		q(:,k) = q(:,k) / h(k,k-1);
		if(k>2)
			epsilon_error = prod(diag(h(1:k-1,1:k-1),1))/factorial((k-1))*dt^(k-1);
			if( abs(epsilon_error) < epsil)
				n = k;
	            break
			end
		end
	end
	%Manthe
	h(n,n) = h(n-1,n-1);
	h(n-1,n) = h(n,n-1);
	
	h = h(1:n,1:n);
	
	if (n0==n) fprintf('Arnoldi: Krylov space too small. Make n larger than %i. The error was %f\n',n,epsilon_error); end
	
	[psi,En] = eig(h); En = diag(En);
	
	psiinv = inv(psi);
	for k=1:n
		a = psi(k,:) * ( exp(En*dt) .* psiinv(:,1) );
		psix = psix + a * q(:,k);
	end
	
