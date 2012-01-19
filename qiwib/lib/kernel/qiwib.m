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
	
function qiwib(input_file, out_arg, ow=[], vv=1, gg=1, print_qiwib_version='unknown',qiwib_output=0,qiwib_verbose_output,dir_current=0,command_args)

%Example output:
%qiwib('inputvariables_relax.m','results/test_relax_15p4mtest/','ow')
 
%global structure, to save memory basis_diff is seperate
mlock(); space -global; global pa basis_diff space realgrid realfunction realbasis

%thats my global structure
pa.basis = 0;
basis_diff = 0;
pa.basis_diff_nz = 0;
pa.hilbert = ['create'];
pa.hilbert_file = 0;
pa.N  = [];
pa.M  = [];
pa.L  = [];
pa.phi  = [];
pa.C  = 0;
pa.Ng  = [];
pa.dx  = 0;
pa.nmax = 0;
pa.V = 0;
pa.V_build = 0;
pa.relaxation = [];
pa.endtime = -1;
pa.dt = 1E-4;
pa.g = 0;
pa.h_kq = 0;
pa.w_ksql = 0;
pa.rho_kq = 0;
pa.rho_ksql = 0;
pa.rho_kq_inv= 0;
pa.CMF_error= 1E-6;
pa.max_CMF_step = 1.0;
pa.xpos = 0;
pa.xpos0 = [];
pa.boundary = 'box';
pa.improved_rlx = 'fixed';
pa.save_dir_out = "";
pa.load_phi_C = 0;
pa.ode_phi = 'RK78';
pa.ode_phi_opts = odeset("AbsTol",1e-10,"InitialStep",1e-4,"MaxStep",1e-3,"RelTol",1e-10);
pa.ode_C = [];
pa.ode_C_opts = [];
pa.H_phi_lin = 0;
pa.H_phi_nl0 = 0;
pa.H_phi_nl1 = 0;
pa.H_phi_nl = 0;
pa.H_phi_direction = 0;
pa.H_C_direction = 0;
pa.Ov_inv = 0;
pa.psi = 0;
pa.H_nl = 0;
pa.nthreads = 0;
pa.nl = [];
pa.nr = [];
pa.dE_limit = -1;
pa.show_plot = 0;
pa.save_step = 1;
pa.H_update_step = -1;
pa.t_initial = 0;
pa.H_Diff = 0;
pa.H_Diff2 = -1/2;
pa.save_options = [1,0,0,0,0];
pa.time = 0;
pa.g1 = 0;
pa.g1_slice = 0;
pa.g1_diag = 0;
pa.C = 0;
pa.qiwib_output = qiwib_output;
pa.dir_current = dir_current;
pa.no_prop_phi = 0;
pa.qiwib_verbose_output = qiwib_verbose_output;
pa.print_qiwib_version = print_qiwib_version;

%% here we have the two variables from the input parameters that can be used in the input-file
pa.V0 = vv;
pa.g0 = gg;

%dir_qiwib = pwd;
%ignore_function_time_stamp ("all");
more off;
%addpath(dir,[dir,"/fieldsmatrices"],[dir,"/integrator"],[dir,"/kernel"],[dir,"/observables"]);

%%nice logo
	printf("\n");
	printf("  ---------\n");
	printf("  | QiwiB |\n");
	printf("  ---------(%s)\n\n",pa.print_qiwib_version);
	printf("     ()\n");
	printf("  --O//)-\n");
	printf("     ()\n\n\n");

%mlock("calc_H_C_C");
	
	%%read inputfile and print some of the variables
	if exist(input_file) != 2
		printf("Input file does not exist! Please check path and name.\n\n");
		return
	end
	source(input_file); if check_input_variables(), return; end;
	
	printf("Additional Information:\n"); 
	printf("  	Length of hilbert space:	%15.14g\n", bincoeff(pa.N+pa.M-1,pa.N) );
	printf("  	Output directory: 		%s\n",out_arg);
	printf("\n"); fflush(stdout);

	if !strcmp(out_arg(end),'/'), out_arg=[out_arg,'/']; end
	pa.save_dir_out=out_arg;
	[s, mess, messid] = mkdir(pa.save_dir_out);
	
	if pa.relaxation!=-2
		if strcmp(mess,'')==0 && strcmp(ow,'ow')==0 && s==1
			disp('!!! Directory exists, use the option overwite (-ow) !!!'); return
		end
		system(['mkdir -p ',pa.save_dir_out]);
		confirm_recursive_rmdir(0); rmdir([pa.save_dir_out,'time_dep'],"s"); confirm_recursive_rmdir(0);
		[s, mess, messid] = mkdir([pa.save_dir_out,'time_dep']);
		copyfile(input_file,[pa.save_dir_out,'input']);
		my_input_file = fopen ([pa.save_dir_out,'input'], "a");
			fprintf(my_input_file,'\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n%%%% QiwiB (%s) using Octave %s was called with the following arguments:\n%%%% qiwib %s\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',pa.print_qiwib_version,version,command_args);
		fclose(my_input_file);

	elseif pa.relaxation==-2
		if strcmp(mess,'')!=0 && s==1
			disp('!!! Output directory does not exist !!!');
			rmdir(pa.save_dir_out);
			return
		end		
	end
	system([["touch ",pa.save_dir_out,".continue"]]);

%%create/save/load basis
	if hilbert() == 1, return; end	

%%initialise and then let's go

	format short e;
	if initialize_phi_C(), return; end;
	if pa.relaxation<0
		propagate();
	else
		relaxate();
	end
	
	more on; format;
