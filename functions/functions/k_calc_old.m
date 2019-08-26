function [cells,truncation_info]=k_calc_old(gencalc,cells,cells_to_calc,gases_to_calc,disp_log)

% [cells,truncation_info]=K_CALC_OLD(gencalc,cells,cells_to_calc,gases_to_calc,disp_log) 
%
% calculates the absorption k-vectors on a cell-by-cell basis using the gas cell data 
% in 'cells' and the calculation information in structure 'gencalc'.
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
% K_CALC_OLD calls GENSPEC functions HITLOADR,LSTRENGTH,WIDTHL,WIDTHD
%
% THIS FUNCTION IS NOW OBSOLETE - USE K_CALC.M
%
% (C) Ben Quine 29-FEB-2000

% record time to process
cpu_time = cputime;
% line_calcs record the number of line calculations that are used during the call to k_calc
line_calcs=0;
gencalc.model=lower(setstr(gencalc.model)); %model is a lower case string
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
   if (max(gases_to_calc) > cells.ngases) | (min(gases_to_calc) < 1)
      error('GENSPECT K_CALC: Invalid Gas Range Specified')
   end
end

% initialise calculation grid
n = gencalc.calc_points;
n_gases = length(gases_to_calc);
cells.k_vector=zeros(n,numcells,n_gases);
cells.k_wing=zeros(gencalc.wing_points,numcells,n_gases);

% do main loop for gases
for jj = 1:n_gases;
   j = gases_to_calc(jj);
   %Load the hitran line data from file 
%   [lins,trunc_info]=hitloadr(cells.param_filename{j},cells.gasname{j},cells.iso(j),gencalc.lower_calc_limit,gencalc.upper_calc_limit,disp_log);
  [lins,trunc_info]=hitloadr(cells.param_filename{j},cells.gasname{j},cells.iso(j),gencalc.start_wavnum,gencalc.stop_wavnum,disp_log);
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
		loren=widthl(cells.tmean(cell_num),lins.temp_coeff,cells.pmean(cell_num),...
	      cells.ppartial(cell_num),lins.air_broad,lins.self_broad);
		%Correct doppler line width
      doppl=widthd(cells.tmean(cell_num),lins.wavenum,lins.iso,cells.gas(j));
      
      % Start with first line and first calculation interval
      line_i = 0;
      calc_int = 0;      
      interval_stop=0;			% marks upper interval limit for calculation
      % Loop through normal calculation intervals
      for line_i = 1:lins.n_lines
         % check to see if interval number needs [should never exceed interval number as line data is wavenumber truncated]
         while lins.wavenum(line_i)>interval_stop
            calc_int = calc_int + 1;
            % define new upper interval limit
 		      interval_stop = gencalc.interval_stop(calc_int);
            % record calculation range for this interval [zero entry means no calculation required for specified range]
            if lins.wavenum(line_i)<=interval_stop	% set indexes if line is now within interval
	            fine_ind = [gencalc.fine_index_start(calc_int):gencalc.fine_index_stop(calc_int)];
		  	      wing_ind = [gencalc.wing1_start(calc_int):gencalc.wing1_stop(calc_int) gencalc.wing2_start(calc_int):gencalc.wing2_stop(calc_int)];
            end
         end
         % do fine calculation [add it straight to vector]
         if ~isempty(fine_ind)
            line_calcs=line_calcs + length(fine_ind);			% record the number of compuations
	         cells.k_vector(fine_ind,cell_num,j) = cells.k_vector(fine_ind,cell_num,j)' + ...
               lineshap(gencalc.calc_grid(fine_ind),lins.wavenum(line_i), loren(line_i),doppl(line_i),cells.tmean(cell_num),gencalc.typ).*linestrength(line_i);
         end
         % do the wing calculation and add it to the wing vector
         if ~isempty(wing_ind)
            line_calcs=line_calcs + length(wing_ind);			% record the number of compuations
	         cells.k_wing(wing_ind,cell_num,j)=cells.k_wing(wing_ind,cell_num,j)' + lineshap(gencalc.wing_grid(wing_ind), ...
	            lins.wavenum(line_i), loren(line_i),doppl(line_i),cells.tmean(cell_num),gencalc.typ).*linestrength(line_i);
         end
      end
      % Interopolate and add contributions from wing part of calculation
      cells.k_vector(:,cell_num,j) = cells.k_vector(:,cell_num,j) + interp1(gencalc.wing_grid,cells.k_wing(:,cell_num,j),gencalc.calc_grid)';
   end
end
% add a wavelength scale
cells.wavnum = gencalc.calc_grid';
% add a legend
cells.k_legend1 = 'k vector units: [cm^2 molecule^-1]';
cells.num_line_calcs = line_calcs;
cells.cputime = cputime-cpu_time;
