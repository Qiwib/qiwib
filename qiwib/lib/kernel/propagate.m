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

	
function propagate()
mlock(); global pa basis_diff space grid gridfunction gridbasis

	printf("\nStart Propagation: Initialise...\r"); fflush(stdout);		
%%initial energy and other variables and display it%%
	
	printf("Start Propagation: Running...                    \n\n"); fflush(stdout);	
	E = calc_E(); E0=E; dE=0;

	cpu_time00 = time(); cpu_time0 = cpu_time00; cpu_time1 = 0; steps_count = 0;
	save_time = pa.save_step+pa.time; H_update_time = pa.H_update_step+pa.time; time_step_old = 0;
	C_corr = 1; phi_corr = 1;
	cpu_time_phi_tot = 0; cpu_time_C_tot = 0; cpu_time_phi_prop_tot = 0; cpu_time_C_prop_tot = 0;
	cpu_time_C = 0; cpu_time_phi = 0; cpu_time_C_prop = 0; cpu_time_phi_prop = 0;


%%initial mean fields%%	
	cpu_time_C = time(); calc_rho(); cpu_time_C = time()-cpu_time_C;
	cpu_time_phi = time(); pa.H_phi_nl = Calc_H_phi_nl(); cpu_time_phi = time()-cpu_time_phi;
	cpu_time_temp = time(); calc_fields(); cpu_time_phi = cpu_time_phi+time()-cpu_time_temp;
	cpu_time_temp = time(); calc_H_C(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;
	
	save_time_dep(0,save_time,E,dE,C_corr,phi_corr,0,steps_count,cpu_time1,cpu_time0,cpu_time00);


%%the propagation%%

while(pa.time<pa.endtime || pa.endtime<0)

		phi0 = pa.phiCpp; C0 = pa.C; E0 = E; h_kq0 = pa.h_kq; w_ksql0 = pa.w_ksql; H_phi_nl0 = pa.H_phi_nl; H_C0 = pa.H_C;
		rho_kq0 = pa.rho_kq; rho_ksql0 = pa.rho_ksql; rho_kq_inv0 = pa.rho_kq_inv;

	%%the next few lines is the propagation

		if pa.no_prop_phi == 0,		
			cpu_time_temp = time(); 
			C12 = Integrator_C(pa.dt/2,pa.time,C0);
			cpu_time_C_prop = cpu_time_C_prop+time()-cpu_time_temp;
		
			cpu_time_temp = time(); 
			phi12b = Integrator_phi(pa.dt/2,pa.time,phi0);
			cpu_time_phi_prop = cpu_time_phi_prop+time()-cpu_time_temp;

			pa.C = C12;
			cpu_time_temp = time(); calc_rho(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;
			cpu_time_temp = time(); pa.H_phi_nl = Calc_H_phi_nl(); cpu_time_phi = cpu_time_phi+time()-cpu_time_temp;
			
			cpu_time_temp = time();
			phi12 = Integrator_phi(pa.dt/2,pa.time,phi0);
			pa.phiCpp = Integrator_phi(pa.dt/2,pa.time+pa.dt/2,phi12);
			cpu_time_phi_prop = cpu_time_phi_prop+time()-cpu_time_temp;

			cpu_time_temp = time(); calc_fields(); cpu_time_phi = cpu_time_phi+time()-cpu_time_temp;
			cpu_time_temp = time(); calc_H_C(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;

			cpu_time_temp = time(); 
			pa.C = Integrator_C(pa.dt/2,pa.time+pa.dt/2,C12);
			pa.H_C_direction = i; C0b = Integrator_C(pa.dt/2,pa.time+pa.dt/2,C12); pa.H_C_direction = -i;
			cpu_time_C_prop = cpu_time_C_prop+time()-cpu_time_temp;
		else
			cpu_time_temp = time(); calc_fields(); cpu_time_phi = cpu_time_phi+time()-cpu_time_temp;
			phi12 = pa.phiCpp; phi12b = phi12;			
			cpu_time_temp = time();
			pa.C = Integrator_C(pa.dt,pa.time,C0);
			pa.H_C_direction = i; C0b = Integrator_C(pa.dt,pa.time+pa.dt,pa.C); pa.H_C_direction = -i;
			cpu_time_C_prop = cpu_time_C_prop+time()-cpu_time_temp;
		end



	%%Calculate observables and output or reset time-step%%
		
		delta = Calc_error(phi12,phi12b,C0,C0b);
		%phi_corr = min(abs( (sum(conj(phi0) .* pa.phi) * pa.dx).^2));
		phi_corr = min(abs( diag(phi0.overlap_matrix(pa.phiCpp)).^2));
		C_corr = abs(pa.C'*C0)^2;
		if pa.qiwib_verbose_output==1, disp([pa.dt,delta]); end		
	
		if delta <= pa.CMF_error

			cpu_time_temp = time(); calc_rho(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;
			cpu_time_temp = time(); pa.H_phi_nl = Calc_H_phi_nl(); cpu_time_phi = cpu_time_phi+time()-cpu_time_temp;
			%% calc_fields has already been calculated before for propagation of C
			cpu_time_temp = time(); calc_H_C(); cpu_time_C = cpu_time_C+time()-cpu_time_temp;
			
			E = calc_E(); dE = E-E0;
			pa.time = pa.time + pa.dt;
			
			steps_count = steps_count+1; cpu_time1 = time(); cpu_time_phi_tot = cpu_time_phi_tot + cpu_time_phi; cpu_time_C_tot = cpu_time_C_tot + cpu_time_C; cpu_time_phi_prop_tot = cpu_time_phi_prop_tot + cpu_time_phi_prop; cpu_time_C_prop_tot = cpu_time_C_prop_tot + cpu_time_C_prop;
			save_time_dep(pa.dt,save_time,E,dE,C_corr,phi_corr,0,steps_count,cpu_time1,cpu_time0,cpu_time00,cpu_time_phi,cpu_time_phi_tot,cpu_time_C,cpu_time_C_tot,cpu_time_phi_prop,cpu_time_phi_prop_tot,cpu_time_C_prop,cpu_time_C_prop_tot);
			cpu_time0=cpu_time1;
			
			if ( abs(pa.time-save_time)<=1E-12 ), save_time =  save_time + pa.save_step; end			
			if (abs(H_update_time-pa.time)<=1E-12) && pa.H_update_step>0
				H_update_time = H_update_time + pa.H_update_step;
				hamiltonian_t(pa.time);
				Calc_H_phi_lin();
			end
			if pa.H_update_step == 0, hamiltonian_t(pa.time); Calc_H_phi_lin(); end		
			
			pa.dt = 0.8 * pa.dt * (pa.CMF_error/delta)^(1/4);
			if (time_step_old>pa.dt) pa.dt = time_step_old; end
			time_step_old = pa.dt;
			if pa.dt > pa.max_CMF_step, pa.dt = pa.max_CMF_step; end
			if pa.H_update_step > 0 && pa.time+pa.dt > H_update_time && H_update_time<=save_time, pa.dt = H_update_time-pa.time;
			elseif pa.time+pa.dt > save_time, pa.dt = save_time-pa.time; end
			time_step_old = time_step_old-pa.dt;
			cpu_time_C = 0; cpu_time_phi = 0; cpu_time_C_prop = 0; cpu_time_phi_prop = 0;
			
		else
		
			pa.phiCpp = phi0; pa.C = C0; E = E0; pa.h_kq = h_kq0; pa.w_ksql = w_ksql0; pa.H_phi_nl = pa.H_phi_nl0; pa.H_C = H_C0;
			pa.rho_kq = rho_kq0; pa.rho_ksql = rho_ksql0; pa.rho_kq_inv = rho_kq_inv0;
			pa.dt = 0.5 * pa.dt * (pa.CMF_error/delta)^(1/4);
			if pa.dt>pa.max_CMF_step, pa.dt=pa.max_CMF_step; end
			
		end

	if exist([pa.save_dir_out,".continue"]) != 2
		break;
	end

end
printf("\nPropagation done.                              \n\n"); fflush(stdout);
system(["rm -f ",pa.save_dir_out,".continue"]);
