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

	
function psix = Lanczos(H,psi0,n,epsil=1e-6,dt)
mlock(); global pa

	n0 = n;
	psix = zeros(length(psi0),1);
	beta = zeros(n+1,1);
	alpha = zeros(n,1);
	q=zeros(length(psi0),n);
	beta(1) = norm(psi0);
	w = psi0;
	for k=1:n
		q(:,k) = w / beta(k);
		w = H*q(:,k);
		alpha(k) = dot(q(:,k),w);
		if k>1, w = w - beta(k) * q(:,k-1); end
		w = w - alpha(k)*q(:,k);
		beta(k+1) = norm(w);

		if(k>2)
			epsilon_error = prod( beta(1:k) )/factorial((k))*dt^k;
			if( abs(epsilon_error) < epsil), n = k; break; end
		end
	end	

	if (n0==n) fprintf('Lanczos: Krylov space too small. Make n larger than %i. The error was %f\n',n,epsilon_error); end

	h_kr = sparse(diag(alpha(1:n)) + diag(beta(2:n),1) + diag(beta(2:n),-1));
	[psi,En] = eig(h_kr); En = diag(En);

	psiinv = inv(psi);
	for k=1:n
		a = psi(k,:) * (exp(pa.H_C_direction*En*dt) .* psiinv(:,1));
		psix = psix + a * q(:,k);
	end
	
