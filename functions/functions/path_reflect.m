function paths=path_reflect(r,e,inc_angle,reflect_angle,ang_size,surface_type,T,paths)

% paths = PATH_REFLECT(r,e,inc_angle,reflect_angle,ang_size,surface_type,T,paths)  
%
% Defines a path segment corresponding to a reflecting surface
%
% R is the reflectivity of the earth's surface for 'smooth' only [scalar]
% E is the emissivity of the earth's surface (usually set to 1-R) [scalar]
% INC_ANGLE is the incident radiadion angle [rad where 0 is normal incidence]
% REFLECT_ANGLE is the reflected angle [rad]
% ANG_SIZE is the angular size of the reflection. Defaults to sr-1 [sr]
% SURFACE_TYPE is either 'smooth' or 'lambertian'
% T is the temperature of the surfact (eg: 298 for the earth). Scalar [K]
% PATHS is a structure containing an existing path segment OR a vector 
%       of wavenumbers (row vector) [cm^{-1}]
%
% Appends or creates a PATH structure. See PATH_ATM.M for details
%
% (C) Ben Quine, University of Toronto, 02-OCT-2001.

% Determine whether existing path has been defined
if ~exist('paths')
   error('GENSPECT PATH_REFLECT: No existing path segment or wavenumber range defined')
end
if isa(paths,'struct')
   wav=paths.wavnum;
elseif isa(paths,'double')
   wav = paths;
   clear paths;
   paths.n_elements = 0;
else			% check that wavnumber ranges are the same
      error('GENSPECT PATH_REFLECT: variable PATHS is not recognised as existing path segment or wavenumber range')
end

n_el = paths.n_elements+1;
paths.n_elements = n_el;
paths.theta{n_el}=[inc_angle,reflect_angle];

%Compute surface albedo
if strcmp(surface_type,'smooth'),
   log_trans = log(r/(sec(inc_angle)*sec(reflect_angle)));
elseif strcmp(surface_type,'lambertian'),
   log_trans = log(r/(2*pi*sec(reflect_angle)));
else
   error('GENSPECT PATH_REFLECT: Invalid surface type, must be ''smooth'' or ''lambertian''')
end

paths.segment_type{n_el} = [upper(surface_type) ' REFLECTION'];
paths.type_index(n_el)=2;
paths.tmean(n_el) = T;
paths.thickness(n_el) = 0;
paths.ang_size(n_el) = ang_size;

%Compute surface blackbody emission
n_wav=length(wav);
paths.emissivity(:,n_el) = ones(n_wav,1).*e;
paths.log_transmission(:,n_el) = ones(n_wav,1).*log_trans;

% add time to path
paths.last_modified = datestr(now);
