function g=voigt(v,v0,aL,aD,T)

% VOIGT(v,v0,aL,aD,T) calculates a voigt profile and
% defaults to the approximation called VOIGTHUM in this
% toolbox
%
% VOIGT calls functions VOIGTWELLS
%
% V is the frequency range [row vector]
% V0 is wavenumber [scalar]
% aL is lorentz line half width [scalar]
% aD is the doppler line half width [scalar]
% T is a dummy variable [not required (used by other lineshape functions)]
%
% (C) Ben Quine, 2001.

%define x, y, and z in the voigt function
c=0.83255461115770;
x=(v-v0).*c./aD;
y=(aL./aD).*c;
% Note ((log(2)/pi).^0.5)=0.46971863934983
g=(0.46971863934983./aD).*voigtwells(x,y);
%g=(0.46971863934983./aD).*voigthum(x,y);
