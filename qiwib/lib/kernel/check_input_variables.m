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


function inp_error = check_input_variables()
mlock(); global pa basis_diff space grid gridfunction gridbasis

inp_error = 0;
inp_warning = 0;

if pa.relaxation >= 0
	if isempty(pa.ode_C), pa.ode_C='eigs'; end
	if isempty(pa.ode_C_opts), pa.ode_C_opts=20; end
else
 	if isempty(pa.ode_C), pa.ode_C='RK78'; end
 	if isempty(pa.ode_C_opts), odeset("AbsTol",1e-10,"InitialStep",1e-4,"MaxStep",1e-3,"RelTol",1e-10); end
end
if !strcmp(pa.improved_rlx,'imaginary') &&  pa.relaxation>=0
	if length(pa.ode_C_opts)<2, pa.ode_C_opts(2) = 2*pa.ode_C_opts(1); end
	if length(pa.ode_C_opts)<3, pa.ode_C_opts(3) = eps; end
end


%% general
	if isempty(pa.relaxation), printf("ERROR: Please specify pa.relaxation in input file!\n" ); inp_error=1; end
	if isempty(pa.N), printf("ERROR: Please specify pa.N in input file!\n" ); inp_error=1; end 
	if isempty(pa.M), printf("ERROR: Please specify pa.M in input file!\n" ); inp_error=1; end
	if isempty(pa.Ng), printf("ERROR: Please specify pa.Ng in input file!\n" ); inp_error=1; end
	if pa.Ng<2, printf("ERROR: Number of grid points pa.Ng too small!\n" ); inp_error=1; end
	if isempty(pa.L), printf("ERROR: Please specify pa.L in input file!\n" ); inp_error=1; end
	if isempty(pa.xpos0), printf("ERROR: Please specify pa.xpos0 in input file!\n" ); inp_error=1; end
	if length(pa.save_options)!=5, printf("ERROR: Please correct pa.save_options!\n" ); inp_error=1; end
	if isempty(pa.ode_phi)==0 && (strcmp(pa.ode_phi,'lanczos') || strcmp(pa.ode_phi,'arnoldi') || strcmp(pa.ode_phi,'eigs'))
		printf("ERROR: Use of wrong integrator for pa.ode_phi!\n" ); inp_error=1; end	
	if (strcmp(pa.ode_phi,'adams') || strcmp(pa.ode_phi,'bdf')) && isempty(pa.ode_phi_opts), printf("ERROR: Please specify pa.ode_phi_opts in input file!\n" ); inp_error=1; end
	if (strcmp(pa.ode_C,'adams') || strcmp(pa.ode_C,'lanczos') || strcmp(pa.ode_C,'arnoldi') || strcmp(pa.ode_C,'bdf')) && isempty(pa.ode_C_opts), printf("ERROR: Please specify pa.ode_C_opts in input file!\n" ); inp_error=1; end
	if exist('hamiltonian_t')==0 && pa.H_update_step>0, printf("ERROR: No hamiltonian_t specified in input file, but pa.H_update_step>0!\n" ); inp_error=1; end
%% relaxation	
	if pa.relaxation >= 0 && isempty(pa.improved_rlx), printf("ERROR: Please specify pa.improved_rlx in input file!\n" ); inp_error=1; end
	if pa.relaxation >= 0 && (strcmp(pa.improved_rlx,'locked') || strcmp(pa.improved_rlx,'follow') ...
	|| strcmp(pa.improved_rlx,'fixed')) && !strcmp(pa.ode_C,'eigs'),
		printf("ERROR: For improved relaxation please set pa.ode_C='eigs' in the input file!\n" ); inp_error=1;
	end
	if pa.relaxation >= 0 && strcmp(pa.improved_rlx,'imaginary') && strcmp(pa.ode_C,'eigs')
		printf("ERROR: Use of wrong integrator for pa.ode_C!\n" ); inp_error=1; end
	if !strcmp(pa.hilbert(1:4),'save') && !strcmp(pa.hilbert(1:4),'load') && !strcmp(pa.hilbert(1:6),'create') && pa.hilbert(1)!=0, printf("ERROR: Please specify correct pa.hilbert in input file!\n" ); inp_error=1; end
	if (strcmp(pa.hilbert(1:4),'save') || strcmp(pa.hilbert(1:4),'load')) && isempty(pa.hilbert(5:end)), printf("ERROR: Please specify correct pa.hilbert in input file!\n" ); inp_error=1; end
	if !isempty(pa.ode_C) && strcmp(pa.ode_C,'eigs') && bincoeff(pa.N+pa.M-1,pa.N)>10 && bincoeff(pa.N+pa.M-1,pa.N)-2<pa.ode_C_opts(1), pa.ode_C_opts(1)=bincoeff(pa.N+pa.M-1,pa.N)-2; printf("WARNING: Number of calculated eigenvalues in pa.ode_C_opts is too big! New value: %i\n",bincoeff(pa.N+pa.M-1,pa.N)-2); inp_warning=1; end
	if !isempty(pa.ode_C) && strcmp(pa.ode_C,'eigs') && bincoeff(pa.N+pa.M-1,pa.N)>10 && bincoeff(pa.N+pa.M-1,pa.N)-1<pa.ode_C_opts(2), pa.ode_C_opts(2)=bincoeff(pa.N+pa.M-1,pa.N)-1; printf("WARNING: Number of calculated lanczos vectors in pa.ode_C_opts is too big! New value: %i\n",bincoeff(pa.N+pa.M-1,pa.N)-1); inp_warning=1; end
	if !isempty(pa.ode_C) && strcmp(pa.ode_C,'eigs') && bincoeff(pa.N+pa.M-1,pa.N)>10 && pa.ode_C_opts(1)>=pa.ode_C_opts(2), pa.ode_C_opts(1)=pa.ode_C_opts(2)-1; printf("WARNING: Number of calculated eigenvalues vectors in pa.ode_C_opts is larger than number of lanzcos vectors! New value: %i\n",pa.ode_C_opts(1)); inp_warning=1; end
	if !isempty(pa.ode_C) && strcmp(pa.ode_C,'eigs') && bincoeff(pa.N+pa.M-1,pa.N)<=10, printf("WARNING: Hilbert space very small, ignoring pa.ode_C and use the full octave eig function for improved relaxation!\n"); inp_warning=1; end

%% hamiltonian
	if exist('hamiltonian_t')==0 && pa.H_update_step<0, printf("WARNING: No hamiltonian_t specified in input file, using default values!\n" ); inp_warning=1; end

%% initial pa.phi

	if isempty(pa.phi), pa.phi = ones(pa.Ng,pa.M); end
	if length(pa.phi(1,:)) == 1
	    if length(pa.phi(:,1)) == 1, pa.phi = pa.phi * ones(pa.Ng,pa.M);
	    elseif length(pa.phi(:,1)) == pa.Ng && pa.load_phi_C != 1, pa.phi = repmat(pa.phi,1,pa.M);
	    elseif length(pa.phi(:,1)) == pa.Ng && pa.load_phi_C == 1, pa.phi = pa.phi;
	    else printf("ERROR: Definition of pa.phi in input file has wrong dimensions!\n" ); inp_error = 1; return;
	    end
	elseif length(pa.phi(:,1)) == 1
	    if length(pa.phi(1,:)) == pa.Ng && pa.load_phi_C != 1, pa.phi = repmat(pa.phi.',1,pa.M);
	    elseif length(pa.phi(1,:)) == pa.Ng && pa.load_phi_C == 1, pa.phi = pa.phi;
	    else printf("ERROR: Definition of pa.phi in input file has wrong dimensions!\n" ); inp_error = 1; return;
	    end
	elseif length(pa.phi(1,:)) == pa.M
	    if length(pa.phi(:,1)) != pa.Ng, printf("ERROR: Definition of pa.phi in input file has wrong dimensions!\n" ); inp_error = 1; return; end
	elseif length(pa.phi(:,1)) == pa.M
	    if length(pa.phi(1,:)) == pa.Ng, pa.phi = pa.phi.';
	    else printf("ERROR: Definition of pa.phi in input file has wrong dimensions!\n" ); inp_error = 1; return;
	    end
	elseif length(pa.phi(1,:)) < pa.M && length(pa.phi(1,:)) ~= pa.Ng && pa.load_phi_C == 1
	    if length(pa.phi(:,1)) == pa.Ng, pa.phi = pa.phi;
	    else printf("ERROR: Definition of pa.phi in input file has wrong dimensions!\n" ); inp_error = 1; return;
	    end
	elseif length(pa.phi(:,1)) < pa.M && length(pa.phi(:,1)) ~= pa.Ng && pa.load_phi_C == 1
	    if length(pa.phi(1,:)) == pa.Ng, pa.phi = pa.phi.';
	    else printf("ERROR: Definition of pa.phi in input file has wrong dimensions!\n" ); inp_error = 1; return;
	    end
	elseif 
	   printf("ERROR: Definition of pa.phi in input file has wrong dimensions!\n" ); inp_error = 1; return;
	end
	
%% print error message
	if inp_error==1
		printf("\nThere are errors in the input file. Quitting Program!");
	end
	if inp_error==1 || inp_warning==1, printf("\n"); end
