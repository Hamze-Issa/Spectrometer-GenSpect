function [number, id_mols]=gasnum(species)

% [number, id_mols]=GASNUM(gas_name_str) takes string gas name and outputs the 
% HITRAN Gas number (integer) followed by the 2-digit id string (eg. '01')  
% to which it corresponds in the hitran database. The return of 'id_mols' is 
% optional.
%
% 'gas_name_str' is a string gas name such as 'hno3'
% for multi-gas problems 'gas_name_str' is of class 'cell' with 
% elements gas{1}='co' gas{2}='hno3' etc. Returns a vector of gas numbers.
%
% Current gases supported are: 
%
% 'All';'H2O';'CO2';'O3';'N2O';'CO';'CH4';'O2';'NO';'SO2';'NO2';'NH3';'HNO3';
% 'OH';'HF';'HCl';'HBr';'HI';'ClO';'OCS';'H2CO';'HOCl';'N2';'HCN';'CH3Cl';
% 'H2O2';'C2H2';'C2H6';'PH3';'COF2';'SF6';'H2S';'HCOOH';'HO2';'O';'ClONO2';
% 'NO+';'HOBR';'C2H4'
%
% (C) Ben Quine 2000.
%

% Ben Quine 12-APR-1999: Fixed problem with 'All' specification and added id_mols
% Ben Quine 14-APR-1999: Added HOBr and C2H4 to list of Hitran Gases - consistant with Hitran 96

if ~exist('species')
   error('GENSPECT GASNUM: No gas species string passed')
end
data_type=whos('species');
if strcmp(data_type.class,'cell')
   no_species = size(species,2);
elseif strcmp(data_type.class,'char')
   no_species = 1;
else
   error('GENSPECT GASNUM: String passed to function was not "char" or "cell" class')
end

Mol_id=[ 'All   ';'H2O   ';'CO2   ';'O3    ';'N2O   ';'CO    ';'CH4   ';'O2    ';...
     'NO    ';'SO2   ';'NO2   ';'NH3   ';'HNO3  ';'OH    ';'HF    '; ...
     'HCl   ';'HBr   ';'HI    ';'ClO   ';'OCS   ';'H2CO  ';'HOCl  '; ...
     'N2    ';'HCN   ';'CH3Cl ';'H2O2  ';'C2H2  ';'C2H6  ';'PH3   '; ...
     'COF2  ';'SF6   ';'H2S   ';'HCOOH ';'HO2   ';'O     ';'ClONO2'; ...
     'NO+   '];             % 'HOBR  ';'C2H4  ']; % These species commented out untill QT data available, BQ, Jan 2001
Mol_id=lower(Mol_id);
for i=1:no_species
   if data_type.class=='cell'						% specify cell element
      id_str=species{i};
   else   
      id_str=species;								% specify normal character string    
   end   
   id_str=lower(id_str);
	id_mols=[];
	all_flag=0;
	if findstr([' ' id_str ' '],[' ' deblank(Mol_id(1,:)) ' '])
		all_flag=1;
	end
	for j=2:size(Mol_id,1)
	    flag=isempty(findstr([' ' id_str ' '],[' ' deblank(Mol_id(j,:)) ' ']));
	    if all_flag | ~flag
		id=num2str(j-1);
		if size(id,2)<2
			id_mols=[id_mols; '0' id];
		else
			id_mols=[id_mols; id(1:2)];
		end
	    end
	end
	if ~isempty(id_mols) & ~all_flag
	   number(i)=str2num(id_mols);
   elseif all_flag
	   number=str2num(id_mols);
   else
	   number(i) = NaN;
	end
end