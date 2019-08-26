function cells=exptcell(P,T,length_of_cell,gas,iso,mixrat,database_name)

% cells = EXPTCELL(P,T,length_of_cell,gas,iso,mixrat,database_name) 
%
% Generates a structure 'cells' corresponding to a simple laboratory
% cell with specified temperature, pressure and length containing a 
% specified gas mixture.
%
% P is the pressure of the cell (scalar)(Pa)
% T is the cell temperature (scalar)(K)
% LENGTH_OF_CELL is the length of the experimental cell (scalar)(m) 
% GAS is a cell array of gas names (cell array)(string entries)
% ISO is an array of isotopes numbers (double array)
% MIXRAT is an array of the concentrations of each gas in the cell (1 x ngases)
% DATABASE_NAME (optional) String identifying database to use ('hitran96', 'hitran2000') or structure array
%       of alternate gas parameter input filenames (string), one per gas [Defaults to Hitran2000]. 
%
%
% (C) Ben Quine, 2000.

% Record parameter datafiles to use
% set default database if none specified
if ~exist('database_name','var')
    database_name='hitran2000';
end
[cells.param_filename, cells.database] = hitfilenames(database_name,gas);

cells.thickness=length_of_cell;
cells.tmean=T;
cells.pmean=P;
cells.zmean=cells.thickness/2;
cells.umean=P/(T*1.380662e-23);
cells.zbottom=0;
cells.ztop=cells.thickness;

% list gas in each cell
cells.gasname=gas;
cells.gas = gasnum(gas);
cells.iso = iso;
cells.ncells = 1;
cells.ngases = length(cells.gas);

% error checking for gas/iso to ensure the same number of entries
if cells.ngases~=length(cells.iso)
   error('GENSPECT EXPTCELL: The number of gases and isotopes specified are different.')
end

% For each gas, attribute mixing ratio associated with each cell
if exist('mixrat')==1
   data_type=whos('mixrat');
	if strcmp(data_type.class,'cell')
			for i=1:cells.ngases
		   cells.mixrat(:,i)=mixrat{i}; 
   		end      
   elseif strcmp(data_type.class,'double')
      cells.mixrat=mixrat; 						
   else														
      warning('GENSPECT EXPTCELL: Unrecognised Data Class for mixrat: Mixing ratios set to zero');
	   cells.mixrat=zeros(ncells,cells.ngases);
   end
else			% set mixing ratios to zero if they don't exist
   cells.mixrat=zeros(ncells,cells.ngases);
end

% record partial pressure
if size(cells.mixrat,2)==cells.ngases
	cells.ppartial=cells.pmean*ones(1,cells.ngases).*cells.mixrat;
else
   warning('GENSPECT EXPTCELL: Mixing ratio undefined for multi gas problem. cells.ppartial not computed.')
end
