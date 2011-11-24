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

	
function save_time_dep(dt,save_time,E,dE,C_corr,phi_corr,eigs_flag,steps_count,cpu_time1,cpu_time0,cpu_time00,cpu_time_phi,cpu_time_phi_tot,cpu_time_C,cpu_time_C_tot,cpu_time_phi_prop,cpu_time_phi_prop_tot,cpu_time_C_prop,cpu_time_C_prop_tot)
mlock(); global pa basis_diff space realgrid realfunction realbasis phiCpp

	calc_g1(); n_p_phi = sum(pa.g1_diag) * pa.dx; n_p_C = abs(pa.C'*pa.C);
	nl = sum(pa.g1_diag(1:pa.nl))*pa.dx; nm = sum(pa.g1_diag(pa.nl+1:pa.nr-1))*pa.dx; nr = sum(pa.g1_diag(pa.nr:end))*pa.dx;
	if pa.qiwib_output && strcmp(pa.ode_C,'eigs') && pa.t_initial==pa.time, disp(real([pa.time,dt,E,dE,0,C_corr,phi_corr,n_p_C,n_p_phi,eigs_flag]));
	elseif pa.qiwib_output && !strcmp(pa.ode_C,'eigs') && pa.t_initial==pa.time, disp(real([pa.time,dt,E,dE,0,C_corr,phi_corr,n_p_C,n_p_phi]));
	elseif pa.qiwib_output && strcmp(pa.ode_C,'eigs') && pa.t_initial!=pa.time, disp(real([pa.time,dt,E,dE,dE/dt,C_corr,phi_corr,n_p_C,n_p_phi,eigs_flag]));
	elseif pa.qiwib_output && !strcmp(pa.ode_C,'eigs') && pa.t_initial!=pa.time, disp(real([pa.time,dt,E,dE,dE/dt,C_corr,phi_corr,n_p_C,n_p_phi]));
	elseif !pa.qiwib_output && pa.t_initial==pa.time, printf ('	Running...\n	(see log file for output)\n'); end
	fflush(stdout);
	
	if pa.show_plot == 1
		plot(pa.xpos,pa.g1_diag,pa.xpos,real(ones(pa.Ng,1).*pa.V*max(pa.g1_diag)/max(eps+abs(pa.V))),pa.xpos,imag(ones(pa.Ng,1).*pa.V*max(pa.g1_diag)/max(eps+abs(pa.V)))); drawnow;
	end

	if pa.relaxation!=-2 || pa.t_initial!=pa.time
    
		for i=1:pa.M, pa.phi(:,i) = (phiCpp(i-1).get_data()).'; end

		save("-z",[pa.save_dir_out,'phiC_restart.gz'],'-struct','pa','N','M','time','phi','xpos','C');
		
		if pa.t_initial==pa.time, save_file = fopen ([pa.save_dir_out,'log'], "w"); fprintf(save_file,"# time (1)\n# dt (2)\n# E (3)\n# dE (4)\n# dE/dt (5)\n# C_corr (6)\n# phi_corr (7)\n# sum(|C|^2) (8)\n# N (9)\n# [eigs converted (10)]\n##\n");
		else save_file = fopen ([pa.save_dir_out,'log'], "a"); end
			if strcmp(pa.ode_C,'eigs') && pa.t_initial==pa.time, fprintf(save_file,"%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e\n",pa.time,dt,E,dE,0,C_corr,phi_corr,n_p_C,n_p_phi,eigs_flag);
			elseif !strcmp(pa.ode_C,'eigs') && pa.t_initial==pa.time, fprintf(save_file,"%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e\n",pa.time,dt,E,dE,0,C_corr,phi_corr,n_p_C,n_p_phi);
			elseif strcmp(pa.ode_C,'eigs') && pa.t_initial!=pa.time, fprintf(save_file,"%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e\n",pa.time,dt,E,dE,dE/dt,C_corr,phi_corr,n_p_C,n_p_phi,eigs_flag);
			elseif !strcmp(pa.ode_C,'eigs') && pa.t_initial!=pa.time, fprintf(save_file,"%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e	%15.14e\n",pa.time,dt,E,dE,dE/dt,C_corr,phi_corr,n_p_C,n_p_phi); end
		fclose(save_file);
		
		if pa.t_initial==pa.time, save_file = fopen ([pa.save_dir_out,'steps'], "w"); fprintf(save_file,"# steps_count (1)\n# pa.dt (2)\n# t (3)\n# CPU times (in seconds) for\n# - last time step (4)\n# - total time so far (5)\n# - last time step (W_ksql, h_kq, H_nonlinear) (6)\n# - total time (Wksql, hkq, Hnonlinear) (7)\n# - last time step (Hc,rho_ksql,rho_kq) (8)\n# - total time (Hc,rho_ksql,rho_kq) (9)\n# - last time step (phi propagation) (10)\n# - total time (phi propagation) (11)\n# - last time step (C propagation) (12)\n# - total time (C propagation) (13)\n# - date and time when time step was finished (14)\n##\n");
		else save_file = fopen ([pa.save_dir_out,'steps'], "a"); end
			if pa.t_initial!=pa.time, fprintf(save_file,"%i	%8.6f	%8.6f	%8.3f	%8.3f	%8.3f	%8.3f	%8.3f	%8.3f	%8.3f	%8.3f	%8.3f	%8.3f	%s\n",steps_count,dt,pa.time,cpu_time1-cpu_time0,cpu_time1-cpu_time00,cpu_time_phi,cpu_time_phi_tot,cpu_time_C,cpu_time_C_tot,cpu_time_phi_prop,cpu_time_phi_prop_tot,cpu_time_C_prop,cpu_time_C_prop_tot,datestr(clock)); end
		fclose(save_file);	
				
		if pa.t_initial==pa.time, save_file = fopen ([pa.save_dir_out,'pop_nat'], "w"); fprintf(save_file,"# time, population of natural orbitals as a percentage\n");
		else save_file = fopen ([pa.save_dir_out,'pop_nat'], "a"); end
			temp = [pa.time;sort(abs(eig(reshape(pa.rho_kq,pa.M,pa.M))/n_p_phi),'descend')*100]; fprintf(save_file,"%8.6f	",temp); fprintf(save_file,"\n");
		fclose(save_file);
		
		if pa.t_initial==pa.time, save_file = fopen ([pa.save_dir_out,'spatial_population'], "w"); fprintf(save_file,"# time, N, N_left, N_middle, N_right\n");
		else save_file = fopen ([pa.save_dir_out,'spatial_population'], "a"); end
			fprintf(save_file,"%15.14e	%15.14e	%15.14e	%15.14e	%15.14e\n",pa.time,n_p_phi,nl,nm,nr);
		fclose(save_file);	
		
		if (abs(pa.time-save_time)<=1E-12) || pa.t_initial==pa.time
			if pa.save_options(1) == 1
				save("-z",[pa.save_dir_out,'time_dep/',num2str(pa.time),'_density.gz'],"-struct","pa","time","g1_diag","xpos","V");
			end
			if pa.save_options(2) == 1
				save("-z",[pa.save_dir_out,'time_dep/',num2str(pa.time),'_g1_0x.gz'],"-struct","pa","time","g1_slice","xpos");
			end
			if pa.save_options(3) == 1
				save("-z",[pa.save_dir_out,'time_dep/',num2str(pa.time),'_phiC.gz'],"-struct","pa","N","M","time","phi","xpos","C");
			end
			if pa.save_options(4) == 1
				[AA,BB] = eig(reshape(pa.rho_kq,pa.M,pa.M));
				[BB,nBB] = sort(diag(BB),'descend'); AA = AA(:,nBB); pa.phiNO_n = diag(BB);
				pa.phiNO = zeros(pa.Ng, pa.M);
				for i=1:pa.M, for j=1:pa.M
					pa.phiNO(:,i) += AA(j,i) * pa.phi(:,j);
				end, end
				save("-z",[pa.save_dir_out,'time_dep/',num2str(pa.time),'_phiNO.gz'],"-struct","pa","time","phiNO","xpos","phiNO_n");
			end
			if pa.save_options(5) == 1
				save("-z",[pa.save_dir_out,'time_dep/',num2str(pa.time),'_rho.gz'],"-struct","pa","rho_kq","rho_ksql");
			end	
		end
	
	end	

