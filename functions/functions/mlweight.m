function mweight=mlweight(Mol)

% MLWEIGHT(MOL) returns the molecular weight for all isotopes 
% of molecule number Mol (in the hitran database).
%
% MLWEIGHT calls no other GENSPECT functions
%
% MOL is the gas number as referenced in the HITRAN database
%
% (C) Ben Quine, University of Toronto, 1999

% keep data global to this function
global mlw
% load data is not already there
if isempty(mlw)
   load mlw_data.mat
   disp('GENSPECT MLWEIGHT: Loading mlweight data')
end
% Constants [atmoic mass unit]
amu_to_kg=1.6605654e-27;
% extract required weight [all isotomers]
mweight = mlw.mass(mlw.mol_index(Mol)+[1:mlw.ison(Mol)])'*amu_to_kg;
