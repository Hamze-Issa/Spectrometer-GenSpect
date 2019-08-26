function strength=lstrength(lins,temp)

% LSTRENGTH(lins,temp) calculates the
% line strength as corrected by the temperature of the absorbing gas.
%
% LSTRENGTH calls functions QTWEIGHT
%  
% LINS is the line database
% TEMP is the temperature
%
% (C) Ben Quine, 2000 

Mol=lins.gas;
% define constants

c2=1.438786;					% units are cm K [Kaye And Laby 14th Ed.] 
Tref=296;						% units are K
% load qt data for specified molicule
QT=qtweight(Mol,temp);
% calculate strength
strength=exp(c2.*lins.E_lower.*(1./Tref-1./temp)).* ...
  lins.intens.*QT.ratio(lins.iso).*(1-exp(lins.wavenum.*(-c2/temp)))./(1-exp(lins.wavenum.*(-c2/Tref)));
