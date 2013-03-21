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

function inp_error = initialize_phi_C();
mlock(); global pa basis_diff space complexgrid complexfunction complexbasis
inp_error = 0;


%%%
%% TODO: grid gridfunction gridbasis instead of (real|complex|realspin2|complexspin2|...)function, etc.

        grid         = complexgrid;
	gridfunction = complexfunction;
	gridbasis    = complexbasis;
	disp("entering initialie_phi_C\n");
%%%%%%%%%%%%%%%%% variables I need %%%%%%%%%%%%%%%%%%%%%%%
	
	printf("Creating initial environment:\r"); fflush(stdout);
	
	%%Initialize all variables
	pa.nmax=bincoeff(pa.N+pa.M-1,pa.N);
	pa.rho_kq = zeros(pa.M*pa.M,1);
	pa.h_kq = zeros(pa.M*pa.M,1);
	if pa.C==0, pa.C = [1;zeros(pa.nmax-1,1)]; end
	pa.H_C = sparse(pa.nmax,pa.nmax);
	
	pa.w_ksql = zeros(pa.M*pa.M*pa.M*pa.M,1);
	pa.rho_ksql = zeros(pa.M*pa.M*pa.M*pa.M,1);
	pa.H_C = sparse(pa.nmax);

	printf("Creating grid: %d,%d,%d\n",pa.xpos0,pa.xpos0+pa.L,pa.Ng);
	space = grid(pa.xpos0,pa.xpos0+pa.L,pa.Ng);
	printf("Setting boundary condition: %d\n",pa.boundary);
	space.set_boundary(pa.boundary);
	printf("Getting x-coordinates\n");
	pa.xpos = space.get_xs(); pa.dx = space.dx;
	printf("Creating basis with %d orbitals\n",pa.M);
	pa.phiCpp = complexbasis(space,pa.M); 
	printf("Done.\n");

	%% grid = realgrid|complexgrid|complexspin2grid
	%% gridfunction = realgridfunction|complexgridfunction|complexspin2gridfunction
	%% gridbasis = (realbasis|complexbasis|complexspin2basis|complexspin3basis)
	
	if isempty(pa.nl), pa.nl=pa.xpos0+pa.L/2; end
	if isempty(pa.nr), pa.nr=pa.xpos0+pa.L/2; end
	[tsts,tsts_i] = min(abs([pa.xpos0:pa.L/(pa.Ng-1):pa.xpos0+pa.L]-pa.nl)); pa.nl = tsts_i;
	[tsts,tsts_i] = min(abs([pa.xpos0:pa.L/(pa.Ng-1):pa.xpos0+pa.L]-pa.nr)); pa.nr = tsts_i;
	
	if exist('hamiltonian_t')!=0, hamiltonian_t(pa.time); end
	if isempty(pa.g), printf("\n\nERROR: Please specify pa.g in input file or pa.g0 as command line parameter!\n" ); inp_error=1; return; end
	if (length(pa.V(:,1)) == 1 || length(pa.V(:,1)) == pa.Ng ) && length(pa.V(1,:)) == 1, pa.V = pa.V .* ones(pa.Ng,1);
	elseif length(pa.V(:,1)) == 1 && length(pa.V(1,:)) == pa.Ng, pa.V = pa.V.';
	else printf("\n\nERROR: wrong dimensions for pa.V!\n"); inp_error = 1; return;
	end
	Calc_H_phi_lin();
	printf("Initializing exernal potential\n");
%%        pa.a0Cpp = gridfunction(pa.V.'); % TODO: Possibility for different spin-components in V
        pa.a0Cpp = gridfunction(pa.V.'); 

	
%%%%%%%%%%%%%%%%% create initial state %%%%%%%%%%%%%%%%%%%%%%%

	printf("Creating initial state. relaxation: %g, load_phi_C: %d\n",pa.relaxation,pa.load_phi_C);
	if pa.relaxation >= 0 && pa.load_phi_C == 0
	        printf("Create initial phi\n");
		%%Create initial phi
		al=pa.H_Diff2/12/pa.dx/pa.dx;
		T = sparse( diag(-30.0*al*ones(pa.Ng,1) + pa.V_build.*ones(pa.Ng,1)) + diag( al*16*ones(pa.Ng-1,1),1) + diag( al*16*ones(pa.Ng-1,1),-1) - diag( al*ones(pa.Ng-2,1),2) - diag( al*ones(pa.Ng-2,1),-2) );

		if strcmp(pa.boundary,'periodic')
			T(1,pa.Ng) = 16*al;
			T(pa.Ng,1) = 16*al;
			T(1,pa.Ng-1) = -al;
			T(pa.Ng,2) = -al;
			T(2,pa.Ng) = -al;
			T(pa.Ng-1,1) = -al;			
		end
		[v,lambda]=eigs(T,pa.M+2,'sa'); [B_dummy,n_eigen] = sort( real(diag(lambda)) );

		pa.phi = v(:,n_eigen(1:pa.M));
		for n=1:pa.M, pa.phi(:,n) = pa.phi(:,n)/sqrt( abs(pa.phi(:,n))'*abs(pa.phi(:,n))*pa.dx); end
		
		%%Create initial C
		printf("Initializing data in orbital basis set\n");
		pa.phiCpp = pa.phiCpp.set_data_vector(complex(pa.phi(:)));
		calc_fields(); calc_H_C();

		
		n=pa.nmax-4; n=max(n,pa.relaxation+1);
		if !strcmp(pa.improved_rlx,'imaginary')
			if length(pa.ode_C_opts)<2, eigs_opts.p = 2*pa.ode_C_opts(1); else, eigs_opts.p=pa.ode_C_opts(2); end
			if length(pa.ode_C_opts)<3, eigs_opts.tol = eps; else, eigs_opts.tol = pa.ode_C_opts(3); end	
			n=min(pa.nmax-4,round(pa.ode_C_opts(1))); n=max(n,pa.relaxation+1);
		end
				
		if isreal(pa.H_C)
			pa.H_C = (pa.H_C+pa.H_C')/2;
			if pa.nmax > 10, [c,dummy,eigs_flag]=eigs(pa.H_C,n,'sa',eigs_opts);
			else [c,dummy]=eig(full(pa.H_C)); end
		else
			if pa.nmax > 10, [c,dummy,eigs_flag]=eigs(pa.H_C,n,'sr',eigs_opts);
			else [c,dummy]=eig(full(pa.H_C)); end		
		end
		[B_dummy,n_eigen] = sort( real(diag(dummy)) );
		pa.C = c(:,n_eigen(pa.relaxation+1));
		pa.C = pa.C / sqrt(sum(abs(pa.C).^2));
		calc_rho();

	elseif pa.relaxation >= 0 && pa.load_phi_C != 1 && pa.load_phi_C != 0
	
		load(pa.load_phi_C);
		[phi,C,inp_error] = AddModes(phi,C); if inp_error>0, return; end
		pa.phi = phi;
		pa.C = C;
		
	elseif pa.relaxation==-1 && pa.load_phi_C != 1 && pa.load_phi_C != 0
	
		load(pa.load_phi_C);
		[phi,C,inp_error] = AddModes(phi,C); if inp_error>0, return; end
		pa.phi = phi .* pa.phi;
		pa.C = C;
	
	elseif pa.relaxation==-2 && pa.load_phi_C != 1 && pa.load_phi_C != 0
	
		load([pa.save_dir_out,'phiC_restart.gz']);
		[phi,C,inp_error] = AddModes(phi,C); if inp_error>0, return; end
		pa.phi = phi;
		pa.t_initial = time;
		pa.C = C;	
		
	elseif pa.load_phi_C == 1
	
		[phi,C,inp_error] = AddModes(pa.phi,pa.C); if inp_error>0, return; end
		pa.phi = phi;
		pa.C = C;	
		
	end
	dbstop("initialize_phi_C.m",156);

	pa.phiCpp = pa.phiCpp.set_data_vector(complex(pa.phi));

	Normalise();
	
	if pa.load_phi_C != 0, calc_fields(); calc_rho(); end

	if pa.relaxation<0
		pa.H_phi_direction = -I; pa.H_C_direction = -I;   %% -i for propagation
	else
		pa.H_phi_direction = -1; pa.H_C_direction = -1;   %% -1 for relaxation
	end
	pa.time = pa.t_initial;
	printf("Creating initial environment: done     \n"); fflush(stdout);
	fflush(stdout);

function [phi,C,inp_error] = AddModes(phi,C)
mlock(); global pa
	
	inp_error = 0;
	cputime=time();	
		
	np = length(pa.phi(1,:)); if pa.load_phi_C == 1, np = pa.M; end
	npl = length(phi(1,:));
	
	%% If the initial files have more single-particle functions than the current simulation then stop
	if np < npl, printf("\nError: Too few single-particle functions for the initial restart file!\n" ); inp_error = 1; return; end
	
	%% Adds additional single-particle functions if the initial set is smaller than the current one
	%% Also reads the initial C and places its values correctly into the larger current C	
	
	if np > npl
		printf("Adding single-particle functions to the initial set: calculating...\r"); fflush(stdout);
		
		nr_eigen = (npl+1)*(np-npl);
		al=pa.H_Diff2/12/pa.dx/pa.dx;
		T = sparse( diag(-30.0*al*ones(pa.Ng,1) + pa.V_build.*ones(pa.Ng,1)) + diag( al*16*ones(pa.Ng-1,1),1) + diag( al*16*ones(pa.Ng-1,1),-1) - diag( al*ones(pa.Ng-2,1),2) - diag( al*ones(pa.Ng-2,1),-2) );
		if strcmp(pa.boundary,'periodic')
			T(1,pa.Ng) = 16*al;
			T(pa.Ng,1) = 16*al;
			T(1,pa.Ng-1) = -al;
			T(pa.Ng,2) = -al;
			T(2,pa.Ng) = -al;
			T(pa.Ng-1,1) = -al;			
		end
		[v,lambda]=eigs(T,nr_eigen,'sa'); [b_dummy,n_eigen] = sort( real(diag(lambda)) );
		phi_tmp = v(:,n_eigen);
		for n=1:nr_eigen
			phi_tmp(:,n) = phi_tmp(:,n)/sqrt( abs(phi_tmp(:,n))'*abs(phi_tmp(:,n))*pa.dx);
			for k=1:npl
				phi_tmp(:,n) = phi_tmp(:,n) - phi(:,k) * (phi(:,k)'*phi_tmp(:,n)*pa.dx)/sqrt( abs(phi(:,k))'*abs(phi(:,k))*pa.dx)/sqrt( abs(phi_tmp(:,n))'*abs(phi_tmp(:,n))*pa.dx);
			end
			phi_tmp(:,n) = phi_tmp(:,n)/sqrt( abs(phi_tmp(:,n))'*abs(phi_tmp(:,n))*pa.dx);
		end
		for n=1:np-npl
			phi(:,n+npl) = zeros(pa.Ng,1);
			for k=1:npl+1
				phi(:,n+npl) = phi(:,n+npl) + phi_tmp(:,k+(n-1)*(npl+1));
			end
		end
		basis_tmp = calc_hilbert_C(pa.N,npl,bincoeff(pa.N+npl-1,pa.N)); if pa.N == 1, basis_tmp = eye(npl); end
		for n=1:length(basis_tmp(:,1))
			n_same(n) = find(sum( abs(pa.basis-[ones(length(pa.basis(:,1)),1)*basis_tmp(n,:),zeros(length(pa.basis(:,1)),np-npl)])' ) == 0 );
		end
		C_tmp = C; C = zeros(length(pa.basis(:,1)),1); C(n_same) = C_tmp;
		
		printf('Adding single-particle functions to the initial set: done in %ds.                          \n',time()-cputime); fflush(stdout);
	end
	
