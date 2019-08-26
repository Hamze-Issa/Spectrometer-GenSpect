function paths=path_limb(cells,cell_index,gas_index,mixrat,tangent_height,ang_size,paths)

% paths = PATH_LIMB(cells,cell_index,gas_index,mixrat,tangent_height,ang_size,paths)
%
% Generates a limb veiwing path through a set of gas cells defined in CELLS. A limb
% viewing geometry may be entirely described in terms of a single parameter TANGENT_HEIGHT
% which is equivilent to the lowest point in the path. The Earth radius is assumed to be
% 6730000 m. The function make no correction for refraction (as yet).
%
% Inputs:
% CELLS is a structure defining the path cells to be used (structure)
% CELL_INDEX lists the gas cells to be included (row vector). Pass as [] to include 
%            all cells. [Note cells are sorted into altitude order by function].
% GAS_INDEX lists the gases to be included (row vector). Pass as [] to include all gases.
% MIXRAT defines the mixing ratio MATRIX elements for the gases and cells corresponding 
%        to cell_index and gas_index entries. Gases by column, cells are by row.
% TANGENT_HEIGHT is the tangent height of the limb view (m).
% ANG_SIZE is the angular size of the source. Defaults to sr-1 [sr]
% PATHS  is an optional parameter. Pass an existing PATHS structure to append this 
%        path segment to it OR pass an index of CELL wavenumber elements to include in 
%        the path.
%
% Outputs:
% The function calls on the structure cells and returns a structure PATHS 
% containing all the requisit information for radiance calculations. It 
% contains the following subvariables:
%
% N_ELEMENTS is the number of segments in the path
% TMEAN is the mean temperature of each gas cell
% PMEAN is the mean pressure of each gas cell
% ZMEAN is the mean altitude of each gas cell
% UMEAN is the number of molecules in each gas cell (m^{-3}) 
% MIXRAT is the mixing ratio of the absorbing gases
%			for multi gas problems it is a matrix of type [gases by column]
% UPARTIAL is the number density of each gas (m-3)
% NTOTAL is the number of molicules by gas (m-2 of path)
% ZBOTTOM and ZTOP are the heights of the top and bottom of the cells. 
% THETA is the rediation incidence angle (Radians from zenith)
% ANG_SIZE is the angular size of the path through the atmosphere [defaults to 1 or (sr-1)]
% LOG_TRANSMISSION is the natural log of the transmission (transmission = exp(log_transmission))
% EMISSIVITY is the emissivity by frequency for each layer
% WAVNUM is the wavenumber grid used in the calculations
% THICKNESS is the layer thickness
% GAS is a structure listing the Hitran numbers of each gas considered
% SEGMENT_TYPE is a structure describing each path segment (eg: 'GAS CELL')
% TYPE_INDEX is an index describing each path segment type
%            0 = SOURCE type segment
%            1 = GAS CELL type segment
%            2 = REFLECTED PATH type segment
% LAST_MODIFIED gives the date when PATHS was last modified
%
% (C) Ben Quine, 2001.

% check if cell_index is defiend
if isempty(cell_index)
   cell_index=[1:cells.ncells]';
else		% force index to be column vector
   cell_index = reshape(cell_index,length(cell_index),1);
end
% check if gas_index is defiend
if isempty(gas_index)
   gas_index=[1:cells.ngases]';
end
% Determine whether existing path has been defined
max_ind = length(cells.wavnum);
if ~exist('paths')
   paths.n_elements = 0;
   wav_index=[1:max_ind];
elseif isempty(paths)
   paths.n_elements = 0;
   wav_index=[1:max_ind];
elseif isa(paths,'struct')			% check that wavnumber ranges are the same
   if length(paths.wavnum)~=length(cells.wavnum)
      error('GENSPECT PATH_LIMB: Wavenumber ranges of PATHS and CELLS are different')
   end
   wav_index=[1:max_ind];
elseif isa(paths,'double')			% used to limit range of path
   if min(size(paths))==1
      wav_index = paths;
      clear paths;
	   paths.n_elements = 0;
   else
      error('GENSPECT PATH_LIMB: PATHS variable passed to function unrecognised')
   end
   if max(wav_index)>length(cells.wavnum) | min(wav_index)<1
      error('GENSPECT PATH_LIMB: PATHS index variable has values not consistant with CELLS')
   end
end
% check mixing ratio definintion
if size(mixrat)~=[length(cell_index),length(gas_index)]
      error('GENSPECT PATH_LIMB: Mixing Ratio not correctly defined. Cannot calculate concentrations')
end

% Find which cell the tangent_height falls within [reindexing to full cells structure]
tangent_cell_sub_ind = find(cells.ztop(cell_index)>tangent_height & cells.zbottom(cell_index)<=tangent_height);
tangent_cell_ind=cell_index(tangent_cell_sub_ind);
% check tangent_height only occurs within one cell and find that cell
if length(tangent_cell_ind)==0
   error('GENSPECT PATH_LIMB: Tangent height specified is not within the cell range supplied')
elseif length(tangent_cell_ind)>1
   error('GENSPECT PATH_LIMB: Tangent height specified occurs within multiple cells. Problem not unique!')
end
% find indexes of all the cells above the tanget height and order them by increasing altitude
cells_above_tangent_sub_ind = find(cells.zbottom(cell_index)>tangent_height);
[nul,ind] = sort(cells.zbottom(cell_index(cells_above_tangent_sub_ind)));
% reindex to section of cells structure [required later for ref to mixrat]
cells_above_tangent_sub_ind = cells_above_tangent_sub_ind(ind);
% reindexing to full cells structure
cells_above_tangent_ind = cell_index(cells_above_tangent_sub_ind);

% compute path length through each one of the atmospheric shells
% Path length is computed by defining the radiation path in simple 2D geometry
% [currently a straight line - no refraction] with X-perpendicular to zenith tangent
% and Y-vertical along zenith tangent.

% Define Earth radius
%R = 6730e3;						
R = 3397.2e3;					% for Mars
% do lowest cell first
tangent_cell_pathl = 2*sqrt((cells.ztop(tangent_cell_ind)+R)^2 - (R+tangent_height)^2);
% do rest of cells
other_cells_pathl = sqrt((cells.ztop(cells_above_tangent_ind)+R).^2 - (R+tangent_height)^2)- ...
   			sqrt((cells.zbottom(cells_above_tangent_ind)+R).^2 - (R+tangent_height)^2);

% number of gases
ngases = length(gas_index);
% record size of initial path
n_elements_initial = paths.n_elements;
% setup index to cells for each path element added
path_ind = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate structure array paths [first path section to but not including lower tangent cell]
start_ind = paths.n_elements+1;

% find out how many elements to add to path
n_elements_2add = length(cells_above_tangent_ind);

if n_elements_2add>0

	stop_ind = paths.n_elements + n_elements_2add;
	path_ind = [path_ind; cells_above_tangent_ind(end:-1:1)];		% record indexes for k-calc reference
	paths.n_elements=stop_ind;
	
	paths.zbottom(start_ind:stop_ind)=cells.zbottom(cells_above_tangent_ind(end:-1:1));
	paths.ztop(start_ind:stop_ind)=cells.ztop(cells_above_tangent_ind(end:-1:1));
	paths.pmean(start_ind:stop_ind)=cells.pmean(cells_above_tangent_ind(end:-1:1));
	paths.tmean(start_ind:stop_ind)=cells.tmean(cells_above_tangent_ind(end:-1:1));
	paths.umean(start_ind:stop_ind)=cells.umean(cells_above_tangent_ind(end:-1:1));
	paths.zmean(start_ind:stop_ind)=cells.zmean(cells_above_tangent_ind(end:-1:1));
	paths.thickness(start_ind:stop_ind) = other_cells_pathl(end:-1:1);

	upartial=cells.umean(cells_above_tangent_ind(end:-1:1))*ones(1,ngases).*mixrat(cells_above_tangent_sub_ind(end:-1:1),:); 	%the number of absorbing gas molecules m-3
	ntotal=other_cells_pathl(end:-1:1)*ones(1,ngases).*upartial/1e4;
    mixrat_down = mixrat(cells_above_tangent_sub_ind(end:-1:1),:);
	% loop to record molicule numbers and mixing ratios by layer
	j=0;
	for i=start_ind:stop_ind
	   j=j+1;
	   paths.mixrat{i}=mixrat_down(j,:);
	   paths.upartial{i}=upartial(j,:);
	   paths.ntotal{i}=ntotal(j,:);
	   paths.gas{i}=cells.gas(gas_index);
	   paths.theta{i}=0;
	   paths.segment_type{i} = 'GAS CELL';
	end
   paths.type_index(start_ind:stop_ind)=ones(1,n_elements_2add);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate structure array paths [mid path section]
mid_ind = paths.n_elements+1;
paths.n_elements=mid_ind;
path_ind = [path_ind; tangent_cell_ind];		% record indexes for k-calc reference
	
paths.zbottom(mid_ind)=cells.zbottom(tangent_cell_ind);
paths.ztop(mid_ind)=cells.ztop(tangent_cell_ind);
paths.pmean(mid_ind)=cells.pmean(tangent_cell_ind);
paths.tmean(mid_ind)=cells.tmean(tangent_cell_ind);
paths.umean(mid_ind)=cells.umean(tangent_cell_ind);
paths.zmean(mid_ind)=cells.zmean(tangent_cell_ind);
paths.thickness(mid_ind) = tangent_cell_pathl;

upartial=cells.umean(tangent_cell_ind)*ones(1,ngases).*mixrat(tangent_cell_sub_ind,:); 	%the number of absorbing gas molecules m-3
ntotal=tangent_cell_pathl*ones(1,ngases).*upartial/1e4;
paths.mixrat{mid_ind}=mixrat(tangent_cell_sub_ind,:);
paths.upartial{mid_ind}=upartial;
paths.ntotal{mid_ind}=ntotal;
paths.gas{mid_ind}=cells.gas(gas_index);
paths.theta{mid_ind}=0;
paths.segment_type{mid_ind} = 'GAS CELL';
paths.type_index(mid_ind)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate structure array paths [first path section to but not including lower tangent cell]
start_ind = paths.n_elements+1;

% find out how many elements to add to path
n_elements_2add = length(cells_above_tangent_ind);

if n_elements_2add>0

	stop_ind = paths.n_elements + n_elements_2add;
	path_ind = [path_ind; cells_above_tangent_ind];		% record indexes for k-calc reference
	paths.n_elements=stop_ind;
	
	paths.zbottom(start_ind:stop_ind)=cells.zbottom(cells_above_tangent_ind);
	paths.ztop(start_ind:stop_ind)=cells.ztop(cells_above_tangent_ind);
	paths.pmean(start_ind:stop_ind)=cells.pmean(cells_above_tangent_ind);
	paths.tmean(start_ind:stop_ind)=cells.tmean(cells_above_tangent_ind);
	paths.umean(start_ind:stop_ind)=cells.umean(cells_above_tangent_ind);
	paths.zmean(start_ind:stop_ind)=cells.zmean(cells_above_tangent_ind);
	paths.thickness(start_ind:stop_ind) = other_cells_pathl;

	upartial=cells.umean(cells_above_tangent_ind)*ones(1,ngases).*mixrat(cells_above_tangent_sub_ind,:); 	%the number of absorbing gas molecules m-3
	ntotal=other_cells_pathl*ones(1,ngases).*upartial/1e4;
	% loop to record molicule numbers and mixing ratios by layer
	j=0;
	for i=start_ind:stop_ind
	   j=j+1;
	   paths.mixrat{i}=mixrat(cells_above_tangent_sub_ind(j),:);
	   paths.upartial{i}=upartial(j,:);
	   paths.ntotal{i}=ntotal(j,:);
	   paths.gas{i}=cells.gas(gas_index);
	   paths.theta{i}=0;
	   paths.segment_type{i} = 'GAS CELL';
	end
   paths.type_index(start_ind:stop_ind)=ones(1,n_elements_2add);
end

% add wavnumber scale
paths.wavnum=cells.wavnum(wav_index);
 
% loop to calculate log transmission
start_ind = n_elements_initial+1;
stop_ind = paths.n_elements;
log_trans=zeros(length(wav_index),stop_ind-start_ind+1);
% extract number densities from structure
ntotal = reshape([paths.ntotal{start_ind:stop_ind}],length(gas_index),stop_ind-start_ind+1)';
for i=1:ngases
   % note 1/1e4 factor to convert units from SI to hitran prefered 
	log_trans = log_trans-ones(length(wav_index),1)*ntotal(:,i)'.*cells.k_vector(wav_index,path_ind,gas_index(i));
end
paths.log_transmission(:,start_ind:stop_ind)=log_trans;
paths.emissivity(:,start_ind:stop_ind) = 1-exp(log_trans);
paths.ang_size(start_ind:stop_ind)=ang_size;

paths.upartial_legend = 'upartial is number density of gas m^{-3}';
paths.ntotal_legend = 'ntotal: number of active gas molecules per square cm of path (cm^{-2})';
% add time to path
paths.last_modified = datestr(now);

