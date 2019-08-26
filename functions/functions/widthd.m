function aD=widthd(t,v0,Iso,gasnumber)

% WIDTHD(t,v0,iso,gasnumber) computes the doppler half width at temperature T
% and wavenumber v0 for a molecule of molecular mass M (in kg).  
%
% WIDTHD calls no other genasis functions
%
% T is the temperature
% V0 is the wavenumber
% M is the molecular mass of the gas in kg
%
% (C) Ben Quine, 1999

% Physical constants
kB=1.380662e-23; %J/K
c=2.99792458e8; %m/s
const=sqrt(2*log(2)*kB*t/(c^2));

% Calculate molicule mass for each isotope
mol_mass = mlweight(gasnumber);

% formula
aD=v0./sqrt(mol_mass(Iso)).*const;
