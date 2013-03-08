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


function inp_error = hilbert()
mlock(); global pa basis_diff

	cputime=time();
	inp_error = 0;
	
	if pa.hilbert==0 || strcmp(pa.hilbert(1:6),"create")
	
		%% create hilbert space	
		printf("Create Hilbert space:  basis...\r"); fflush(stdout);
		pa.basis = calc_hilbert_C(pa.N,pa.M,bincoeff(pa.N+pa.M-1,pa.N)); if pa.N == 1, pa.basis = eye(pa.M); end
		
		%% The new algorithm for calc_hilbert_c already produces an ordered basis so we can skip this step
		%% if pa.N~=1 && pa.M~=1, pa.basis=sortrows(pa.basis,[pa.M:-1:1]); end
		
		printf("Create Hilbert space:  basis_diff...\r"); fflush(stdout);
		[basis_diff,pa.basis_diff_nz] = calc_hilbert_diff_C(pa.basis,pa.M,bincoeff(pa.N+pa.M-1,pa.N));
		printf("Create Hilbert space:  done in %ds.                                      \n",time()-cputime);
		
	elseif strcmp(pa.hilbert(1:4),"save")
		
		%% create and save hilbert space	
		if !strcmp(pa.hilbert(5),'/'), pa.hilbert = [pa.hilbert(1:4),pa.dir_current,'/',pa.hilbert(5:end)]; end
		printf("Create and save Hilbert space:  basis...\r"); fflush(stdout);
		pa.basis = calc_hilbert_C(pa.N,pa.M,bincoeff(pa.N+pa.M-1,pa.N)); if pa.N == 1, pa.basis = eye(pa.M); end
		
		%% if pa.N~=1 && pa.M~=1, pa.basis=sortrows(pa.basis,[pa.M:-1:1]); end
		
		printf("Create and save Hilbert space:  basis_diff...\r"); fflush(stdout);
		[basis_diff,pa.basis_diff_nz] = calc_hilbert_diff_C(pa.basis,pa.M,bincoeff(pa.N+pa.M-1,pa.N));
		printf("Create and save Hilbert space:  saving...             \r"); fflush(stdout);
		basis_diff_nz = pa.basis_diff_nz; pa.basis = basis;
		save('-z',[pa.hilbert(5:end)],'basis','basis_diff','basis_diff_nz');
		printf("Create and save Hilbert space:  done in %ds.                                      \n",time()-cputime);

	elseif strcmp(pa.hilbert(1:4),"load")
		
		%% load hilbert space
		if !strcmp(pa.hilbert(5),'/'), pa.hilbert = [pa.hilbert(1:4),pa.dir_current,'/',pa.hilbert(5:end)]; end
		printf("Load Hilbert space:  loading...\r"); fflush(stdout);
		load([pa.hilbert(5:end)],'basis');
		if length(basis(:,1)) != bincoeff(pa.N+pa.M-1,pa.N) || length(basis(1,:)) != pa.M
			printf("\n\nERROR: The dimensions of the loaded hilbert matrix do not match pa.N and pa.M from the input file!\n\nQitting Program!\n"); inp_error = 1; return;
		end
		load([pa.hilbert(5:end)],'basis_diff','basis_diff_nz');
		pa.basis = basis; pa.basis_diff_nz = basis_diff_nz;	
		printf("Load Hilbert space:  done in %ds.                                      \n",time()-cputime);
		
	end
	
	
	
	
	%a=[M,ones(1,N-1)]; basis=[];
	%while (a(end)<M)
		%b=[zeros(a(1),M-a(1)),eye(a(1))];
		%for s=2:N
			%row = a(s); len = a(1);
			%b=b+[zeros(len,row-1),ones(len,1),zeros(len,M-row)];
		%end
		%basis=[basis;b];
		
		%a(2)=a(2)+1;
		
		%for r=2:N-1
			%if a(r)==M+1
				%a(r+1)=a(r+1)+1;
				%a(2:r)=a(r+1);
			%end
		%end	
		%a(1)=M-a(2)+1;
		%printf("Create Hilbert space:  %d / %d                  \r", a(end),M); fflush(stdout);
	%end
	%basis=[basis;[zeros(1,M-1),N]];

