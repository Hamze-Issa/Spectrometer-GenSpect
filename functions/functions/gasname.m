function name=gasname(number)

% GASNAME(gas_number) takes the gas number and outputs the HITRAN gas name (string).
%
% NUMBER is a HITRAN gas number (cell array)
% For multi-gas problems 'gas_number' is of class 'cell' with elements 
% gas_number{1}='3', gas_number{2}='14' etc. Returns a vector of gas numbers.
%
% Current gases supported run from numbers 1-38: (An entry of 0 will return
% a list of all the gas names.)
%
% 'All';'H2O';'CO2';'O3';'N2O';'CO';'CH4';'O2';'NO';'SO2';'NO2';'NH3';'HNO3';
% 'OH';'HF';'HCl';'HBr';'HI';'ClO';'OCS';'H2CO';'HOCl';'N2';'HCN';'CH3Cl';
% 'H2O2';'C2H2';'C2H6';'PH3';'COF2';'SF6';'H2S';'HCOOH';'HO2';'O';'ClONO2';
% 'NO+';'HOBR';'C2H4'
%
% (C) Ben Quine, 2000, Rev 1.01.

if nargin ~= 1
   error('GENSPECT GASNAME: Incorrect number of input arguments.')
end

data_type=whos('number');
if strcmp(data_type.class,'cell')
   no_species = size(number,2);
elseif strcmp(data_type.class,'double')
   no_species = 1;
else
   error('GENSPECT GASNAME: Input to function must be ''double'' or ''cell'' class');
end

Mol_id={'All   ';'H2O   ';'CO2   ';'O3    ';'N2O   ';'CO    ';'CH4   ';'O2    '; ...
     'NO    ';'SO2   ';'NO2   ';'NH3   ';'HNO3  ';'OH    ';'HF    '; ...
     'HCl   ';'HBr   ';'HI    ';'ClO   ';'OCS   ';'H2CO  ';'HOCl  '; ...
     'N2    ';'HCN   ';'CH3Cl ';'H2O2  ';'C2H2  ';'C2H6  ';'PH3   '; ...
     'COF2  ';'SF6   ';'H2S   ';'HCOOH ';'HO2   ';'O     ';'ClONO2'; ...
     'NO+   '};             % 'HOBR  ';'C2H4  '}; % These species commented out until QT data released, BQ, Feb 2001

for i=1:no_species
   if strcmp(data_type.class,'cell')
      if number{i} == 0,
         name=Mol_id';
         break
      else
         id_str=number{i};
      end
   else
     if number==0
         name=Mol_id';
         break
     elseif number<length(Mol_id)
         id_str=number;
     end
   end
   if id_str<length(Mol_id)
       name(i)=deblank(Mol_id(id_str+1));
   else
       name(i)={'Unknown'};
   end
end