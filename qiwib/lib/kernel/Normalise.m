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


function Normalise()
mlock(); global pa


	pa.C = pa.C / sqrt(abs(pa.C'*pa.C));
	
	for n1 = 1:pa.M
		for n2 = 1:n1-1
			xx = ( pa.phi(:,n2)'*pa.phi(:,n1) )/sqrt(abs(pa.phi(:,n1)'*pa.phi(:,n1)))/sqrt(abs(pa.phi(:,n2)'*pa.phi(:,n2)));
			pa.phi(:,n1) = pa.phi(:,n1) - xx * pa.phi(:,n2);
		end
		pa.phi(:,n1) = pa.phi(:,n1) /  sqrt( abs(pa.phi(:,n1)'*pa.phi(:,n1)*pa.dx) );
	end
