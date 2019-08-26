function paths=path_source(T,paths,ang_size)

% paths = PATH_SOURCE(T,paths,ang_size) defines a blackbody source.
%
% Inputs:
% T is the temperature of the source in Kelvin. If no value is passed, the temperature
%		is taken to be 6000K. (scalar) [K]
% PATHS is either an existing path segment PATHS (struct) or a vector of wavenumbers
%			(row vector) [cm^-1]
% ANG_SIZE is the angular size of the source. Defaults to sr-1 [sr]
%
% Outputs a structure PATHS consistant with PATH_ATM and PATH_REFLECT
%
% (C) Ben Quine, 31-OCT-2003.

% Determine whether existing path has been defined
if ~exist('ang_size')
   ang_size=1;							% set to units sr-1
end
if ~exist('paths')
   error('GENSPECT PATH_SOURCE: No existing path segment or wavenumber range defined')
end
if isa(paths,'struct')
    % path structure already defined. Append new info
elseif isa(paths,'double')
    % wavnumber range passed. Create path
   wavenum = paths;
   clear paths;
   paths.wavnum = wavenum;
   paths.n_elements = 0;
else			% check that wavnumber ranges are the same
      error('GENSPECT PATH_SOURCE: variable PATHS is not recognised as existing path segment or wavenumber range')
end

n_el = paths.n_elements+1;
n_wav=length(paths.wavnum);
paths.n_elements = n_el;

paths.theta{n_el}=0;
paths.segment_type{n_el} = ['SOURCE'];
paths.type_index(n_el)=0;
paths.ang_size(n_el)=ang_size;
paths.thickness(n_el) = 0;
paths.zbottom(n_el)=0;
paths.ztop(n_el)=0;
paths.zmean(n_el)=0;
paths.tmean(n_el) = T;

% set emissivity to 1
paths.emissivity(:,n_el) = ones(n_wav,1);
% set transmission to 1 (log_transmission to 0)
paths.log_transmission(:,n_el) = zeros(n_wav,1);

% add time to path
paths.last_modified = datestr(now);

