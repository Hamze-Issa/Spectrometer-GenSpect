function g=lineshape(v,v0,aL,aD,T,typ)

% g=LINESHAPE(v,v0,aL,aD,T,typ) calculates the power of a given line
%
% V is the frequency range
% V0 is the band frequency
% aL is the lorentz line half width
% aD is the doppler line half width
% T is the temperature [Needed by VANVLECK]
% typ is either 'lorentz','doppler','vanvleck',or 'voigt'
%
% LINESHAPE calls functions VOIGT,LORENTZ,DOPPLER,VANVECK,VOIGTHUI
%
% (C) Ben Quine, 2000.

typ=lower(setstr(typ)); 
eval(['g=' typ '(v,v0,aL,aD,T);']);
