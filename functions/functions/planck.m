function radiance=planck(wavenum,temp)

% radiance=PLANCK(wavenum,temp) 
%
% Computes the planck or blackbody radiance of a body given its radiating 
% temperature [K] and the wavenumber [cm^{-1}] of radiation.
%
% If WAVENUM and TEMP are column and row vectors, respectively, then PLANCK
% will return the corresponding matrix (n_v .* n_T) of planck function values.
% It does so according to the formula B(v,T)=2h*c^2*v^3/(exp(hcv/kT)-1)
% The constants 2hc^2 and hc/kT have 7 significant digits of precision.
%
% WAVENUM is the wavenumber of radiation (ROW vector) [cm^{-1}]
% TEMP    is the radiating temperature of the blackbody (COLUMN vector) [K]
% RADIANCE is in units [Wm^{-2}sr^{-1}(cm^{-1})^{-1}].
%
% (C) Ben Quine, University of Toronto, 02-OCT-2001.

% reference: J T Houghton, The Physics of Atmospheres, 2nd ed., 1997.

% Define Constants 
c1=1.191062e-8; %2hc^2 [Wm^-2sr^-1(cm^-1)^4]
c2=1.438786;    %hc/kT [K(cm^-1)^-1]

if size(wavenum,2)~=1 
   error('GENSPECT PLANCK: Wavenumber input must be a column vector')
end
if size(temp,1)~=1
   error('GENSPECT PLANCK: Temperature input must be a row vector')
end

% Do matrix calculation
tiny=1e-30;
radiance=((c1.*(wavenum.^3))*ones(size(temp)))./(exp(wavenum*(c2./(tiny+temp)))-1);
