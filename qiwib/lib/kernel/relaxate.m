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
	
function relaxate()
mlock(); global pa basis_diff space grid gridfunction gridbasis
	printf("\nStart Relaxation for state %d: Initialise...\r",pa.relaxation); fflush(stdout);			
%%initial energy and other variables and display it%%

	printf("Start Relaxation for state %d:                     \n\n",pa.relaxation); fflush(stdout);		
	E = calc_E(); E0=E; dE=0;
	
	cpu_time00 = time(); cpu_time0 = cpu_time00; cpu_time1 = 0; steps_count = 0; 
	C00=pa.C; dE = 0; save_time = pa.save_step+pa.time; H_update_time = pa.H_update_step+pa.time;
	dE_per_second = 2*pa.dE_limit; dE_per_second_old = 2*pa.dE_limit; time_step_old = 0;
	C_corr = 1; phi_corr = 1;
	cpu_time_phi_tot = 0; cpu_time_C_tot = 0; cpu_time_phi_prop_tot = 0; cpu_time_C_prop_tot = 0;
	cpu_time_C = 0; cpu_time_phi = 0; cpu_time_C_prop = 0; cpu_time_phi_prop = 0;
	if strcmp(pa.ode_C,'eigs')
		if length(pa.ode_C_opts)<2, eigs_opts.p = 2*pa.ode_C_opts(1); else, eigs_opts.p=pa.ode_C_opts(2); end
		if length(pa.ode_C_opts)<3, eigs_opts.tol = eps; else, eigs_opts.tol = pa.ode_C_opts(3); end	
	end

%%initial mean fields%%	
	cpu_time_temp = time(); calc_rho(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;
	cpu_time_temp = time(); pa.H_phi_nl = Calc_H_phi_nl(); cpu_time_phi = cpu_time_phi+time()-cpu_time_temp;
	cpu_time_temp = time(); calc_fields(); cpu_time_phi = cpu_time_phi+time()-cpu_time_temp;
	cpu_time_temp = time(); calc_H_C(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;
	save_time_dep(0,save_time,E,dE,C_corr,phi_corr,0,steps_count,cpu_time1,cpu_time0,cpu_time00);

%%the propagation%%

while(pa.time<pa.endtime || pa.endtime<0)

		phi0 = pa.phiCpp;
		C0 = pa.C; E0 = E; h_kq0 = pa.h_kq; w_ksql0 = pa.w_ksql; H_phi_nl0 = pa.H_phi_nl; H_C0 = pa.H_C;
		rho_kq0 = pa.rho_kq; rho_ksql0 = pa.rho_ksql; rho_kq_inv0 = pa.rho_kq_inv;

		%%Calculating C
		cpu_time_temp = time();
		if strcmp(pa.improved_rlx,'follow') && strcmp(pa.ode_C,'eigs')
		
			eigs_opts.v0 = pa.C;
			n=min(pa.nmax-4,round(pa.ode_C_opts(1))); n=max(n,pa.relaxation+1);
			if isreal(pa.H_C)
				pa.H_C = (pa.H_C+pa.H_C')/2;
				if pa.nmax > 10, [c,dummy,eigs_flag]=eigs(pa.H_C,n,'sa',eigs_opts);
				else [c,dummy]=eig(full(pa.H_C)); eigs_flag=0; end
			else
				if pa.nmax > 10, [c,dummy,eigs_flag]=eigs(pa.H_C,n,'sr',eigs_opts);
				else [c,dummy]=eig(full(pa.H_C)); eigs_flag=0; end		
			end
			ovrlp0 = 0; ovrlp1 = ovrlp0;
			for i=1:n
				ovrlp0 = abs(c(:,i)'*C0);
				if ovrlp0 >= ovrlp1
					ovrlp_i = i;
					ovrlp1 = ovrlp0;
				end
			end
			pa.C = c(:,ovrlp_i);
		
		elseif strcmp(pa.improved_rlx,'lock') && strcmp(pa.ode_C,'eigs')
		
			eigs_opts.v0 = pa.C;
			n=min(pa.nmax-4,round(pa.ode_C_opts(1))); n=max(n,pa.relaxation+1);
			if isreal(pa.H_C)
				pa.H_C = (pa.H_C+pa.H_C')/2;
				if pa.nmax > 10, [c,dummy,eigs_flag]=eigs(pa.H_C,n,'sa',eigs_opts);
				else [c,dummy]=eig(full(pa.H_C)); eigs_flag=0; end
			else
				if pa.nmax > 10, [c,dummy,eigs_flag]=eigs(pa.H_C,n,'sr',eigs_opts);
				else [c,dummy]=eig(full(pa.H_C)); eigs_flag=0; end		
			end
			ovrlp0 = 0; ovrlp1 = ovrlp0;
			for i=1:n
				ovrlp0 = abs(c(:,i)'*C00);
				if ovrlp0 >= ovrlp1
					ovrlp_i = i;
					ovrlp1 = ovrlp0;
				end
			end
			pa.C = c(:,ovrlp_i);

		elseif strcmp(pa.improved_rlx,'fixed') && strcmp(pa.ode_C,'eigs')

			eigs_opts.v0 = pa.C;
			n=min(pa.nmax-4,round(pa.ode_C_opts(1))); n=max(n,pa.relaxation+1);
			if isreal(pa.H_C)
				pa.H_C = (pa.H_C+pa.H_C')/2;
				if pa.nmax > 10, [c,dummy,eigs_flag]=eigs(pa.H_C,n,'sa',eigs_opts);
				else [c,dummy]=eig(full(pa.H_C)); eigs_flag=0; end
			else
				if pa.nmax > 10, [c,dummy,eigs_flag]=eigs(pa.H_C,n,'sr',eigs_opts);
				else [c,dummy]=eig(full(pa.H_C)); eigs_flag=0; end		
			end
			[B_dummy,n_eigen] = sort( real(diag(dummy)) );
			pa.C = c(:,n_eigen(pa.relaxation+1));
			
		elseif strcmp(pa.improved_rlx,'imaginary');
					
			pa.C = Integrator_C(pa.dt,pa.time,pa.C);
			eigs_flag = [];
			
		elseif (strcmp(pa.improved_rlx,'fixed')==0 && strcmp(pa.improved_rlx,'lock')...
		 	&& strcmp(pa.improved_rlx,'follow')) || strcmp(pa.ode_C,'eigs')==0
		
			disp('!!could not find improved relaxation integrator for C!!');
			return					
		end		
		cpu_time_C_prop = cpu_time_C_prop+time()-cpu_time_temp;

		%%Calculating phi
		if pa.no_prop_phi == 0
			cpu_time_temp = time(); 
			pa.phiCpp = Integrator_phi(pa.dt,pa.time,phi0);
			cpu_time_phi_prop = cpu_time_phi_prop+time()-cpu_time_temp;

			cpu_time_temp = time(); calc_rho(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;
			cpu_time_temp = time(); pa.H_phi_nl = Calc_H_phi_nl(); cpu_time_phi = cpu_time_phi+time()-cpu_time_temp;

			cpu_time_temp = time(); 
			phi1 = Integrator_phi(pa.dt,pa.time,phi0);
			cpu_time_phi_prop = cpu_time_phi_prop+time()-cpu_time_temp;
		else
			cpu_time_temp = time(); calc_rho(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;
			phi1 = pa.phiCpp;
		end

	%%Calculate observables and output or reset time-step%%	

		delta = Calc_error(phi1,pa.phiCpp,1,1);
		Normalise();

		phi_corr = min(diag( abs( (phi0.overlap_matrix(pa.phiCpp)).^2) ));
		C_corr = abs(pa.C'*C0)^2;
		if pa.qiwib_verbose_output==1, disp([pa.dt,delta]); end
		
		if delta < pa.CMF_error

			%%calc_rho and Calc_H_phi_nl done before to calculate error estimate
			cpu_time_temp = time(); calc_fields(); cpu_time_phi = cpu_time_phi+time()-cpu_time_temp;
			cpu_time_temp = time(); calc_H_C(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;

			E = calc_E(); dE = E-E0;
			pa.time = pa.time + pa.dt; %% for relaxation I actually only propagate dt each time step
			
			steps_count = steps_count+1; cpu_time1 = time(); cpu_time_phi_tot = cpu_time_phi_tot + cpu_time_phi; cpu_time_C_tot = cpu_time_C_tot + cpu_time_C; cpu_time_phi_prop_tot = cpu_time_phi_prop_tot + cpu_time_phi_prop; cpu_time_C_prop_tot = cpu_time_C_prop_tot + cpu_time_C_prop;
			save_time_dep(pa.dt,save_time,E,dE,C_corr,phi_corr,eigs_flag,steps_count,cpu_time1,cpu_time0,cpu_time00,cpu_time_phi,cpu_time_phi_tot,cpu_time_C,cpu_time_C_tot,cpu_time_phi_prop,cpu_time_phi_prop_tot,cpu_time_C_prop,cpu_time_C_prop_tot);
			cpu_time0=cpu_time1;
			
			if (abs(pa.time-save_time)<=1E-12), save_time =  save_time + pa.save_step; end
			if (abs(H_update_time-pa.time)<=1E-12) && pa.H_update_step>0
				H_update_time = H_update_time + pa.H_update_step;
				hamiltonian_t(pa.time);
				Calc_H_phi_lin();
			end
			if pa.H_update_step == 0, hamiltonian_t(pa.time); Calc_H_phi_lin(); end
			
			dE_per_second = dE/pa.dt;
			pa.dt = 1.2*pa.dt; %3/4 * pa.dt * (pa.CMF_error/delta)^(1/4);
			if (time_step_old>pa.dt) pa.dt = time_step_old; end
			time_step_old = pa.dt;
			if pa.dt > pa.max_CMF_step, pa.dt = pa.max_CMF_step; end		
			if pa.H_update_step > 0 && pa.time+pa.dt > H_update_time && H_update_time<=save_time, pa.dt = H_update_time-pa.time;
			elseif pa.time+pa.dt > save_time, pa.dt = save_time-pa.time; end
			time_step_old = time_step_old-pa.dt;
			cpu_time_C = 0; cpu_time_phi = 0; cpu_time_C_prop = 0; cpu_time_phi_prop = 0;

		else
			cpu_time_phi_tot = cpu_time_phi_tot + cpu_time_phi; cpu_time_C_tot = cpu_time_C_tot + cpu_time_C;
			pa.phiCpp = phi0; pa.C = C0; E = E0;
			pa.h_kq = h_kq0; pa.w_ksql = w_ksql0;
			pa.rho_kq = rho_kq0; pa.rho_ksql = rho_ksql0; pa.rho_kq_inv = rho_kq_inv0;
			pa.dt = 0.5*pa.dt; %pa.dt = pa.dt * (pa.CMF_error/delta)^(1/3);
			if pa.dt > pa.max_CMF_step, pa.dt = pa.max_CMF_step; end
			
		end	
	
	if pa.dE_limit > abs(dE_per_second)+abs(dE_per_second_old) && pa.dE_limit > 0
		break;
	end
	dE_per_second_old = dE_per_second;
	
	if exist([pa.save_dir_out,".continue"]) != 2
		break;
	end
	
end
printf("\nRelaxation done.                              \n\n"); fflush(stdout);
system(["rm -f ",pa.save_dir_out,".continue"]);
