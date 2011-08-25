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


function calc_H_C()
mlock(); global pa

	nthreads = abs(pa.nthreads);

%% first attempt to parrallelise the code to nthread threads
	if nthreads < 1.5
		
		%% split the matrix into smaller portions to save memory
		n_part = 5E6;
		nparts = ceil(1.0*pa.basis_diff_nz/n_part);
		for n=1:nparts-1, inp{n} = [1+(n-1)*n_part,n*n_part]; end
		inp{nparts} = [1+(nparts-1)*n_part,pa.basis_diff_nz];
		
		%% calculate H_C
		H_C = cellfun(@par_calc_H_C,inp);
		
		pa.H_C = H_C{1};
		for i=2:nparts
			pa.H_C = pa.H_C + H_C{i};
		end	

	else
		
		%% split the matrix into smaller portions to save memory and to parallelise the code
		n_part = 5E6;
		nparts = ceil(1.0*pa.basis_diff_nz/n_part);
		for n=1:nparts-1, inp{n} = [1+(n-1)*n_part,n*n_part]; end
		inp{nparts} = [1+(nparts-1)*n_part,pa.basis_diff_nz];
		
		if nthreads>nparts, nthreads=nparts; end
		
		%% calculate H_C
		H_C = parcellfun(nthreads,@par_calc_H_C,inp,"VerboseLevel",0);
		
		pa.H_C = H_C{1};
		for i=2:nparts
			pa.H_C = pa.H_C + H_C{i};
		end	
	
	end

