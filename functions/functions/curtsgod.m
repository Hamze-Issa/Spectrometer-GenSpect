function pt=curtsgod(layer)

% CURTSGOD(layer) generates a mean pressure and temperature value 
% given vectors of layer pressure and temperature using the Curtis Godson
% approximation.
%
% CURTSGOD calls no other genasis functions
%
% layer.p is a vector of layer pressure
% layer.t is a vector of layer temperature
% layer.z is a vector of layer height
%
% (C) Ben Quine APR-1999

%Boltzmann's constant
kB=1.380662e-23; %J/K

%mean number density in the layer by the ideal gas law
ulayer=layer.p./(layer.t*kB);

%mean number density in the cell
normalise = 1./(sum(ulayer)-0.5*(ulayer(1)+ulayer(end)));
%mean pressure, temperature in the cell
pt.p=normalise.*(sum(ulayer.*layer.p)-0.5*(ulayer(1)*layer.p(1)+ulayer(end)*layer.p(end)));
pt.t=normalise.*(sum(ulayer.*layer.t)-0.5*(ulayer(1)*layer.t(1)+ulayer(end)*layer.t(end)));
pt.z=normalise.*(sum(ulayer.*layer.z)-0.5*(ulayer(1)*layer.z(1)+ulayer(end)*layer.z(end)));
pt.u=(sum(ulayer)-0.5*(ulayer(1)+ulayer(end)))./length(ulayer);
