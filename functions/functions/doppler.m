function g=doppler(v,v0,aL,aD,T)

% g = DOPPLER(v,v0,aL,aD,T) generates a doppler profile
%
% v is the frequency range
% v0 is the band frequency
% aD is the Doppler line half width
% aL and T are null variables
%
% (c) Ben Quine, 2000.

g=exp((-log2(2)*(v-v0).^2)./(aD.^2))./(aD.*(pi^0.5));