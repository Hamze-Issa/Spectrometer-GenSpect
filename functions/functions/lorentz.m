function g=lorentz(v,v0,aL,aD,T)

% LORENTZ(v,v0,aL,aD,T) generates a lorentz lineshape given frequency, 
% frequency range, and half width at half max of the lorentz profile
%
% Inputs:
% V is the frequency range about band freq (cm-1)
% V0 is the band frequency (cm-1)
% aL is lorentz line half width (null variable)
% aD need not be passed
% T is temperature (null variable)
%
% (C) Ben Quine, 2000.

g=aL./pi./(aL.^2+(v-v0).^2);

