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
mlock(); space -global; 
global pa basis_diff space grid gridfunction gridbasis realgrid realfunction realbasis complexgrid complexfunction complexbasis realspin2grid realspin2function realspin2basis complexspin2grid complexspin2function complexspin2basis


%dir_qiwib = pwd;
%ignore_function_time_stamp ("all");
more off;
%addpath(dir,[dir,"/fieldsmatrices"],[dir,"/integrator"],[dir,"/kernel"],[dir,"/observables"]);
source("./lib/kernel/defaults.input")

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
	disp("Reached initialize_phi_C()\n");
        [grid,gridfunction,gridbasis] = select_grid(pa.scalar_type,pa.Ncomponents);
	if initialize_phi_C(), return; end;
	if pa.relaxation<0
		propagate();
	else
		relaxate();
	end
	
	more on; format;
endfunction


%% This can be made shorter and more general, but not really necessary until
%% we get more grid types.
function [grid,gridfunction,gridbasis] = select_grid(scalar_type,Ncomponents)
global pa realgrid realfunction realbasis complexgrid complexfunction complexbasis realspin2grid realspin2function realspin2basis complexspin2grid complexspin2function complexspin2basis
	printf("Orbital space is %s with %d components\n",ifelse(scalar_type==@real,"real","complex"),Ncomponents)
	if scalar_type == @real
		if Ncomponents == 1
			grid         = realgrid;
			gridfunction = realfunction;
			gridbasis    = realbasis;
		elseif Ncomopnents == 2
			grid         = realspin2grid;
			gridfunction = realspin2function;
			gridbasis    = realspin2basis;
		else
			disp("Ncomponents > 2. Rebuild Octave module for libspace with support for higher number of components.")
			exit;
		end
	elseif scalar_type == @complex
		if Ncomponents == 1
			grid         = complexgrid;
			gridfunction = complexfunction;
			gridbasis    = complexbasis;
		elseif Ncomopnents == 2
			grid         = complexspin2grid;
			gridfunction = complexspin2function;
			gridbasis    = complexspin2basis;
		else
			disp("Ncomponents > 2. Rebuild Octave module for libspace with support for higher number of components.")
			exit;
		end
	else
	    printf("Unsupported scalar type '%s' selected.\n",scalar_type)
	    exit;
	end        
 endfunction		       
		       

