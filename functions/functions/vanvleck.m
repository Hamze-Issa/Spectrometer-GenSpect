function g=vanvleck(v,v0,aL,aD,T)

% VANVLECK(v,v0,aL,aD,T) generates a lineshape according to the formula
% given by Van Vleck and Huber.
%
% VANVLECK calls no other GENSPECT functions
%
% V is the frequency range
% V0 is the band frequency
% aL is the lorentz line half width
% aD is the doppler line half width
% T is the temperature

%define constants
h=4.135701e-15; %planck's constant (eVs)
c=2.997924580e+8; %speed of light in a vacuum
kB=1.380662e-23; %Boltzmann constant(JK-1)
v=v+v0; %adjust range

%do calc
	f=h*c./(2*kB.*T);
	gLn=aL./pi./(aL.^2+(-v-v0).^2);
	gLp=aL./pi./(aL.^2+(-v-v0).^2);

	g=(v./v0).*(tanh(f.*v)./tanh(f.*v0)).*(gLp+gLn);