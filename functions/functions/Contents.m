% University of Toronto GENSPECT v1.1
%
% (C) Ben Quine, 15-APR-2000.
%     Email: ben@genspect.com
%
%  Atmosphere Resolution
%
%	ATMREAD  - reads and saves an atmospheric profile file
%	ATMCELL  - resolves an atmospheric profile into cells
%	CONCREAD - reads and saves a concentration profile file
%	PATH_ATM  - generates path segments through a set of gas cells
%	PATH_REFLECT  - calculates a reflection path segment
%	PATH_SOURCE  - calculates a source path segment
%   RADIANCE - calculates radiance and emission for a defined path
%	CURTISGOD - performs the curtis-godson approximation
%	MLWEIGHT - gives the molecular mass (in kg) of a given atmospheric gas
%	QTWEIGHT - table lookup Q values for a given T
%
%  Line Calculations
%
%   EXPTCELL    - Defines a laboratory cell for line-by-line calculation
%	HITLOADR 	- loads and saves hitran data -mat files over a specified range
%   HITFILENAMES- generates datafile names for Hitran format databases
%	GASNUM   	- gives hitran's number code for a gas
%	GASNAME   	- gives gas name according to hitran's number code
%	LSTRENGTH	- corrects line strength for temperature
%	WIDTHL   	- generates lorentz half-width
%	WIDTHD   	- generates doppler half-width
%	LINESHAP 	- calculates a line profile (lorentz,doppler,vanveck,voigt)
%	VOIGT    	- calculates a voigt profile
%	LORENTZ  	- calculates a lorentz profile
%	DOPPLER  	- calculates a doppler profile
%	VANVLECK  	- calculates a vanvleck profile
%	K_CALC   	- calculates wing and centre line absorption coeffs cell-by-cell
%	PLANCK   	- generates a blackbody function for an object
%
%  Calculation Scripts
%
%  BENCHTST - sample calculation script
%  DOWNWARD_PATH_EXAMPLE - downward transmission example script
%  RELFECTED_PATH_EXAMPLE - sample reflected path calculation script
%  LIMB_VIEW_EXAMPLE - Limb view along a tangent height example script
%  PRESS_VS_ALT_EXAMPLE - compares a pressure division calculation with an altitude division
%  
%
%  File I/O and Testing
%
%	HITREAD  - reads and saves hitran data for a gas
%   MLWREAD  - reads and save mass data from Hitran file of format molparam.txt
%   QTREAD   - reads and save Total Inturnal Partition data from Hitran ascii file [see qt.dat]
%
%  REBUILD_DATABASE  is a script that reads all available hitran data files, 
%							molparam.txt and qt.dat
%
%  Additions
%
%   INSTCONV - Function to convolve a spectral calculation with an instrument function.
%   K_H2O_CONTIN - computes an estimate the water vapour continium according to Burch 1976
%	FILTREAD - reads and saves out a filter profile file [untested]
%   
%

%  BUILD_GENGRID_SPLITS is a function that generates precomputed grid spliting data required
%								by SPLIT_GENCALC. Calls GENGRID_DIVIDE and GENGRID_INTERP_ERROR
