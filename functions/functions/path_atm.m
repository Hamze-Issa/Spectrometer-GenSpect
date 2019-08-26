function paths=path_atm(cells,cell_index,gas_index,mixrat,theta,ang_size,paths,wav_index)

% paths = PATH_ATM(cells,cell_index,gas_index,mixrat,theta,ang_size,paths,wav_index)
%
% Generates a nardir veiwing path through a set of gas cells defined in CELLS. 
% The path is a linear  path through the set of gas cells (the order and cells 
% used is determined by the start and stop altitudes) at an incidence angle 
% THETA (measured from zenith in degrees).
%
% CELLS is a structure defining the path cells to be used (structure)
% CELL_INDEX lists the gas cells to be included (row vector). Pass as [] to include 
%            all cells.
% GAS_INDEX lists the gases to be included (row vector). Pass as [] to include all gases.
% MIXRAT defines the mixing ratio MATRIX elements for the gases and cells corresponding 
%        to cell_index and gas_index entries. Gases by column, cells are by row.
% THETA  is the solar theta angle (scalar)(radians).
% ANG_SIZE is the angular size of the path through the atmosphere [defaults to 1 or (sr-1)]
% PATHS  is an optional parameter. Pass an existing PATHS structure to append this 
%        path segment to it OR empty array [].
% WAV_INDEX is an optional parameter. Pass an index of CELL wavenumber elements to include in 
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
% UMEAN is the number of molecules in each gas cell [m^{-3}] 
% MIXRAT is the mixing ratio of the absorbing gases
%			for multi gas problems it is a matrix of type [gases by column]
% UPARTIAL is the number density of each gas [m^{-3}]
% NTOTAL is the number of molicules by gas [cm^{-2} of path]
% ZBOTTOM and ZTOP are the heights of the top and bottom of the cells. 
% THETA is the rediation incidence angle (Radians from zenith)
% LOG_TRANSMISSION is the natural log of the transmission (transmission = exp(log_transmission))
% EMISSIVITY is the emissivity by frequency for each layer
% ANG_SIZE is the angular size of the path through the atmosphere [defaults to 1 or (sr-1)]
% WAVNUM is the wavenumber grid used in the calculations
% THICKNESS is the layer thickness
% GAS is a structure listing the Hitran numbers of each gas considered
% SEGMENT_TYPE is a structure describing each path segment (eg: 'GAS CELL')
% TYPE_INDEX is an index describing each path segment type
%            0 = SOURCE type segment
%            1 = GAS CELL type segment
%            2 = REFLECTED PATH type segment
%            3 = User defined SOURCE
% LAST_MODIFIED gives the date when PATHS was last modified
%
% (C) Ben Quine, 31-OCT-2003.

% check if cell_index is defiend
if isempty(cell_index)
   cell_index=[1:cells.ncells]';
end
% check if gas_index is defiend
if isempty(gas_index)
   cell_index=[1:cells.ngases]';
end
% recover frequency elements to compute
max_ind = length(cells.wavnum);
if exist('wav_index')
    if isempty(wav_index)
       wav_index=[1:max_ind];
    end
    if max(wav_index)>length(cells.wavnum) | min(wav_index)<1
         error('GENSPECT PATH_ATM: PATHS index variable has values not consistant with CELLS')
    end
else
   wav_index=[1:max_ind];
end
% Determine whether existing path has been defined
if ~exist('paths')
   paths.n_elements = 0;
elseif isempty(paths)
   paths.n_elements = 0;
elseif isa(paths,'struct')			% check that wavnumber ranges are the same
   if length(paths.wavnum)~=length(cells.wavnum(wav_index))
      error('GENSPECT ATMPATH: Wavenumber ranges of PATHS and CELLS are different')
   end
else
      error('GENSPECT PATH_ATM: PATHS variable passed to function unrecognised')
end
% check mixing ratio definintion
if size(mixrat)~=[length(cell_index),length(gas_index)]
      error('GENSPECT PATH_GAS: Mixing Ratio not correctly defined. Cannot calculate concentrations')
end

% find out how many elements to add to path
n_elements_2add = length(cell_index);
ngases = length(gas_index);

start_ind = paths.n_elements+1;
stop_ind = paths.n_elements + n_elements_2add;
paths.n_elements=stop_ind;

%Generate structure array paths
paths.zbottom(start_ind:stop_ind)=cells.zbottom(cell_index);
paths.ztop(start_ind:stop_ind)=cells.ztop(cell_index);
paths.pmean(start_ind:stop_ind)=cells.pmean(cell_index);
paths.tmean(start_ind:stop_ind)=cells.tmean(cell_index);
paths.umean(start_ind:stop_ind)=cells.umean(cell_index);
paths.zmean(start_ind:stop_ind)=cells.zmean(cell_index);
paths.thickness(start_ind:stop_ind) = (paths.ztop(start_ind:stop_ind) - paths.zbottom(start_ind:stop_ind));

upartial=cells.umean(cell_index')*ones(1,ngases).*mixrat; 	%the number of absorbing gas molecules m-3
ntotal=(cells.ztop(cell_index) - cells.zbottom(cell_index))*ones(1,ngases).*upartial.*sec(theta)/1e4;
% loop to record molicule numbers and mixing ratios by layer
j=0;
for i=start_ind:stop_ind
   j=j+1;
   paths.mixrat{i}=mixrat(j,:);
   paths.upartial{i}=upartial(j,:);
   paths.ntotal{i}=ntotal(j,:);
   paths.gas{i}=cells.gas(gas_index);
   paths.theta{i}=theta;
   paths.segment_type{i} = 'GAS CELL';
end
paths.type_index(start_ind:stop_ind)=ones(1,n_elements_2add);
paths.ang_size(start_ind:stop_ind)=ang_size;
paths.wavnum=cells.wavnum(wav_index);

% loop to calculate log transmission
log_trans=zeros(length(wav_index),n_elements_2add);
for i=1:ngases
   % note 1/1e4 factor to convert units from SI to hitran prefered 
	log_trans = log_trans-ones(length(wav_index),1)*ntotal(:,i)'.*cells.k_vector(wav_index,cell_index,gas_index(i));
end
paths.log_transmission(:,start_ind:stop_ind)=log_trans;
paths.emissivity(:,start_ind:stop_ind) = 1-exp(log_trans);

paths.upartial_legend = 'upartial is number density of gas [m^{-3}]';
paths.ntotal_legend = 'ntotal: number of active gas molecules per square cm of path [cm^{-2}]';
% add time to path
paths.last_modified = datestr(now);

