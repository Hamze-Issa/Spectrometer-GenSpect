function g=voigthui(v,v0,aL,aD)

% VOIGTHUI(v,v0,aL,aD) generates a voigt profile using the complex 
%   error function approximation as provided in A.K. Hui et al,
%	"Rapid Computation of the Voigt and Complex Error Functions",
%	JQSRT 1978.
%	This function uses the p=6 approximation and has error less than
%	10^(-6)
%
% VOIGTHUI calls no other genspec functions
%
% v is the frequency range
% v0 is the band frequency
% aL is the lorentz line half width
% aD is the doppler line half width
%
% Modified for speed 23-09-98 by Ben Quine.

c=0.83255461115770;
zh = aL./aD.*c - i*((v-v0).*c./aD);

%define the coefficients for the rational approximation

a= [122.607931777104326 214.382388695706425 181.92853309218549 93.155580458138441 30.180142196210589 5.9126209773153 0.564189583562615];
b= [122.607931773875350 352.730625110963558 457.334478783897737 348.703917719495792 170.354001821091472 53.992906912940207 10.479857114260399];

g=real(((((((a(7)*zh + a(6) ).*zh + a(5)).*zh + a(4)).*zh + a(3)).*...
   zh + a(2)).*zh+a(1))./(((((((zh+b(7)).*zh+b(6)).*zh+b(5)).*...
   zh+b(4)).*zh+b(3)).*zh+b(2)).*zh+b(1)));
