function aL=widthl(t,tcoeff,p,ppartial,air,self)

% WIDTHL(t,tcoeff,p,ppartial,air,self) computes the lorentz 
% half width of a line.
%
% WIDTHL calls on no other GENSPECT functions
%
% T is temperature
% TCOEFF is the temperature coefficient
% P is the pressure
% PPARTIAL is the partial pressure
% AIR is the air broadened half width
% SELF is the self broadened half width
%
% (C) Ben Quine, 2000.

%Reference temperature and pressure
tref=296;
pref=101325;
%formula
aL=((tref./t).^tcoeff).*((air.*(p-ppartial)./pref)+(self.*ppartial./pref));
