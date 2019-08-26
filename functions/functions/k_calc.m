function [cells,truncation_info]=k_calc(gencalc,cells,cells_to_calc,gases_to_calc,gasNames, disp_log)

% [cells,truncation_info]=K_CALC(gencalc,cells,cells_to_calc,gases_to_calc,disp_log) 
%
% calculates the absorption k-vectors on a cell-by-cell basis using the gas 
% cell data in 'cells' and the calculation information in structure 'gencalc'.
%
% Inputs:
% GENCALC is a structure containing calculation grids and indexing calculation (struct)
% CELLS is a structure containing atmospheric cell information
% CELLS_TO_CALC specifies the order of the cells in the calculation (row vector)
% GASES_TO_CALC numbers the gases to be used in the calculation (row vector)
% DISP_LOG is set to 'off' to supress all non-error output (for www version)
%
% Outputs:
% CELLS is a structure containing the k-vector information in [cm^2 molecule^-1]
% 		  in (wavnumber,cell number,gas) format  
%
% K_CALC calls GENSPECT functions HITLOADR,LSTRENGTH,WIDTHL,WIDTHD
%
% (C) Ben Quine 23-MAY-2000.

% record time to process
cpu_time = cputime;
% y_range records the maxiumium and minimium values of the voigt parameter y that are encoutered [DIAGNOSTIC ONLY]
%y_range=[100 0];
%dx_range=[1e10 0];

% line_calcs record the number of line calculations that are used during the call to k_calc
cells.num_line_calcs=0;
% set 'model' is a lower case string
gencalc.model=lower(setstr(gencalc.model)); 
% setup truncation information string (www version)
truncation_info = '';

% set display log setting default
if exist('disp_log')~=1
   disp_log='on';
end

% Define Cells to include in calculations
if exist('cells_to_calc')
   % check that cells listed exist
   if (max(cells_to_calc) > cells.ncells) | (min(cells_to_calc) < 1)
      error('GENSPECT K_CALC: Invalid Cell Range Specified')
   end
   % define range
	numcells=length(cells_to_calc);
else
   numcells=cells.ncells;
	cells_to_calc = [1:numcells];   
end

% Define Gases to include in calculations
if ~exist('gases_to_calc')
   % define gas range
	gases_to_calc = [1:cells.ngases];
else
   % check for valid range of gases
   if (max(gases_to_calc) > cells.ngases) || (min(gases_to_calc) < 1)
      error('GENSPECT K_CALC: Invalid Gas Range Specified')
   end
end
% check mixing ratio ranges
if min(cells.mixrat(cells_to_calc,gases_to_calc))<0 || max(cells.mixrat(cells_to_calc,gases_to_calc))>1
    error('GENSPECT K_CALC: Mixing ratio out of range [0,1]. Cannot compute absorption coefficient')
end

% square root of ln2
sqrtln2=0.83255461115770;

% Load Grid division data
global grid_division
if isempty(grid_division)
	disp('GENSPECT K_CALC: Loading Grid Division Data')
   % This data is built by function build_gengrid_splits.m
   if gencalc.accuracy<=0.01
      load('grid_division_data_0_01.mat');
      disp('GENSPECT K_CALC: Interpolation accuracy set to better then 0.01%')
   elseif gencalc.accuracy<=0.1
      load('grid_division_data_0_1.mat');
      disp('GENSPECT K_CALC: Interpolation accuracy set to better then 0.1%')
   else
		load('grid_division_data_1_0.mat');
	   disp('GENSPECT K_CALC: Interpolation accuracy set to better then 1%')
	end
end
% determine which line function to use
if strcmp(lower(gencalc.typ),'voigt')
   default_line_function=1;
else
   default_line_function=0;
end
% record the minimium grid interval
dv=gencalc.grid_interval(1);
% initialise calculation grids
n = gencalc.grid_pts(1);
n_gases = length(gases_to_calc);
% initialise k_vector storage arrays
cells.k_vector=zeros(n,numcells,n_gases);
% initalise variable for storing previous aL aD and v parameters to avoid indexing routine
aL_prev=0;
aD_prev=0;
v_prev=0;
y_prev=NaN;
dx_prev=NaN;

% do main loop for gases
for jj = 1:n_gases;
   j = gases_to_calc(jj);
   %Load the hitran line data from file 
   [lins,trunc_info]=hitloadr(cells.param_filename{j},cells.gasname{j},cells.iso(j), ...
      gencalc.start_wavnum-gencalc.wingcutoff,gencalc.stop_wavnum+gencalc.wingcutoff,disp_log, gasNames);
%	% Subbed code to remove wing lines [FOR CODE TESTING ONLY]
%	[lins,trunc_info]=hitloadr(cells.param_filename{j},cells.gasname{j},cells.iso(j),gencalc.start_wavnum,gencalc.stop_wavnum,disp_log);
	% store truncation info (www version)
   truncation_info = [truncation_info trunc_info];
   % do main loop for cells
   ii = 0;
   while ii<numcells
      ii=ii+1;
	   % index of path to calculate
	   cell_num = cells_to_calc(ii);
		%Correct line strength for temperature
		if(gencalc.model=='hitran')
			linestrength=lstrength(lins,cells.tmean(cell_num));  
		elseif(gencalc.model=='elsass')
			linestrength=lins.intens;
		else
			disp(['GENSPECT K_CALC: Incorrect line data model'])
		end
      %Correct lorentz line width
	  aL_lines=widthl(cells.tmean(cell_num),lins.temp_coeff,cells.pmean(cell_num),cells.ppartial(cell_num),lins.air_broad,lins.self_broad);
      %Correct doppler line width
      aD_lines=widthd(cells.tmean(cell_num),lins.wavenum,lins.iso,cells.gas(j));
      for i=1:gencalc.num_grids
         k_vec{i}=zeros(1,gencalc.grid_pts(i));
      end
      for line_i = 1:lins.n_lines
      	% record current line's parameters
         aL=aL_lines(line_i);
         aD=aD_lines(line_i);
         v=lins.wavenum(line_i);
         % Calculate indexes for calculation division 
         y=sqrtln2*(aL/aD);
			% calculate base grid resolution in x-space
         dx=dv*sqrtln2/aD;
         if y~=y_prev
	         % determine which line of data to use
				ny = floor((length(grid_division.y)-1)/(grid_division.y1-grid_division.y0)*(log10(y)+grid_division.y1)+1);
				% ensure that index is within precomputed domain
				ny = max(min(ny,grid_division.n_pts),1);
			end
			if y~=y_prev | dx~=dx_prev
				%% PHASE I: Compute the bounds of the calculation in wavenumber space applying
				%%          limits imposed by the wing cut-off and the required computational range.
				%%				For the remaining ranges compute the gencalc to use.
				% determine the base grid multiplier [conservative rounding down]
				multiplier = max(ceil(log2(dx/grid_division.dx0)+1),1);
				% compute the calculation lower bondary in v space
				max_v = aD/sqrtln2*grid_division.max_x{ny};
				% find and remove grid splits that exceed or equal base grid multiplier
				% limit number of grid divisions to number to be used in calculation
				% limit number of divisions by size of wing computation [symmetric limitations all applied here].
				grid_no=grid_division.grid_no{ny}-multiplier+1;
				% find grid points within the wing cut-off distance [greater than one division]
				ind=find(grid_no>0 & grid_no<gencalc.num_grids & max_v<gencalc.wingcutoff);
				% Generate series of computation boundaries in wavenumber space and a list of grid numbers to use
				% [note that for n computation boundaries there are n-1 grids].
				if length(ind)>0
			   	max_v = [-gencalc.wingcutoff -max_v(ind(end:-1:1)) +max_v(ind) +gencalc.wingcutoff];
					% add additional one for the boundary grid
			   	grid_no = grid_no([ind(end:-1:1)+1 ind(1:end) ind(end)+1]);					% for multiple calculation ranges form +- grids
				else		% only 1 grid to be used
				  	max_v = [-gencalc.wingcutoff +gencalc.wingcutoff];
				   grid_no=1;
				end
			end
         % apply non-symmetric calculation limitations  - user imposed calculation range limits.
			% upper bound
			ind = find((v+max_v)<gencalc.stop_wavnum);
			n_ind = length(ind);
			if n_ind==0
   			% non of the calculation is in required range
			   grid_no_nonsym=[];
			   max_v_nonsym=[];
			elseif n_ind<length(max_v)
			   % some of the calculation is in range
			   grid_no_nonsym=grid_no(ind);
            max_v_nonsym = [max_v(ind) gencalc.stop_wavnum-v];
         else
			   grid_no_nonsym=grid_no;
            max_v_nonsym = max_v;
  			end   
			% lower bound
			ind = find((v+max_v)>gencalc.start_wavnum);
			n_ind = length(ind);
			if n_ind==0
			   % non of the calculation is in required range
   			grid_no_nonsym=[];
			   max_v_nonsym=[];
			elseif n_ind<length(max_v)
			   % some of the calculation is in range
			   grid_no_nonsym=grid_no(ind-1);
				max_v_nonsym = [gencalc.start_wavnum-v max_v(ind)];   
         else
			   grid_no_nonsym=grid_no;
            max_v_nonsym = max_v;
			end   
			% store the grid resolutions 
         grid_int = gencalc.grid_interval(grid_no_nonsym);
         % intervals need to be computed? move to phase II
         wav=[];
         wav{1}=[];
         if length(grid_int)>0		   
				%% Phase II: Calculate the closest [convervative] grid division point corresponding to the
				%%           wavenumber ranges specified
				% calculate conservative grid divisions based on the coarsest neighbouring interval
				v_bounds = max_v_nonsym + (v-gencalc.start_wavnum);
				grid_round = grid_int([1 1:end end]);
				i=0;
				prev_bounds=0;
				while i<length(max_v_nonsym)
				   i=i+1;
				   % Grid rounding code
				  	int_round=max(grid_round(i),grid_round(i+1));
				   if max_v_nonsym(i)<0 
				      v_bounds(i) = floor(v_bounds(i)/int_round)*int_round;
				   else
				      v_bounds(i) = ceil(v_bounds(i)/int_round)*int_round;
				   end
   	      end
				% generate corresponding indices for each interval [required max/min as there is only one point at calc bounds]
				n=length(grid_no_nonsym);
				m=0;					% number of index elements recorded
				p=1;					% pointer to wavenumber vector
				for i=1:n
		 			if grid_no_nonsym(i)==1
   	   				% compute indices [round corrects for small computational errors]
					   	low = round(min(v_bounds(i)/grid_int(i)+1,gencalc.grid_pts(1))); 
					  		upp = round(min(max(v_bounds(i+1)/grid_int(i),1),gencalc.grid_pts(1))); 
					  	else
		  				   % compute indices [round corrects for small computational errors]
		   				low = round(min(v_bounds(i)/grid_int(i)*2+1,gencalc.grid_pts(grid_no_nonsym(i)))); 
		   				upp = round(min(max(v_bounds(i+1)/grid_int(i)*2,1),gencalc.grid_pts(grid_no_nonsym(i)))); 
   	            end
				     	% Record grid for calculation if it contains entries [singleton entries are removed as they cause interpolation errors]
   	            if low<upp			
   	         		m=m+1;
% 		   				% Not required, but usefull for checking computation bounds
%					      lower_v(m) = gencalc.grid{grid_no_nonsym(i)}(lower_ind);
%					      upper_v(m) = gencalc.grid{grid_no_nonsym(i)}(upper_ind);
   				      % compile wavenumber vector and index pointers
                     grid_number(m)=grid_no_nonsym(i);
                     lower_ind(m)=low;
                     upper_ind(m)=upp;
  		               stop=p+upp-low;
     	   				wav_ind(m,1:2)=[p stop];
							wav{m}=gencalc.grid{grid_number(m)}(low:upp);         
							p=stop+1;
    					end
				end
				% record final number of entries
      	   n=m;
      	else					% no lines in index
      	   n=0;
			end   
      	% record values of aL aD and v used to build indexing
      	aL_prev=aL;
      	aD_prev=aD;
      	v_prev=v;
  			% concatinate wavnumber ranges together into one vector for processing
      	wav=[wav{:}];
%%    	y_range=[min(y_range(1),index.y) max(y_range(2),index.y)];
%%    	dx_range=[min(dx_range(1),index.dx) max(dx_range(2),index.dx)];
			% loop through grid calculations [as required] 
			if n>0
				if default_line_function
	   	      k = voigt(wav,v,aL,aD,[]).*linestrength(line_i);
	   	   else
	   	      k = lineshap(wav,v,aL,aD,cells.tmean(cell_num),gencalc.typ).*linestrength(line_i);
            end
            %  record the number of compuations
	      	cells.num_line_calcs=cells.num_line_calcs + length(wav);
	      	% index line computations to appropraite grids
	      	for i=1:n
	  	   	   k_vec{grid_number(i)}(lower_ind(i):upper_ind(i)) = k_vec{grid_number(i)}(lower_ind(i):upper_ind(i)) + k(wav_ind(i,1):wav_ind(i,2));
		  		end
      	end
   	end
   	% Interopolate and add contributions from wing part of calculation
   	cells.k_vector(:,cell_num,j) = k_vec{1}';
   	if gencalc.num_grids>1		% do interpolation
		   for i=2:gencalc.num_grids
		      cells.k_vector(:,cell_num,j) = cells.k_vector(:,cell_num,j) + interp1(gencalc.grid{i},k_vec{i},gencalc.grid{1})';
   	   end
   	end
   end
end
% add a wavelength scale
cells.wavnum = gencalc.grid{1}';
% add a legend
cells.k_legend1 = 'k vector units: [cm^2 molecule^-1].';
%cells.y_range=y_range;
%cells.dx_range=dx_range;
cells.cputime = cputime-cpu_time;
