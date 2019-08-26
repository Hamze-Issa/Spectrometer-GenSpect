function mlw=mlwread(filename,savename)

% qt=qtread(filename,savename)
%
% This function reads Hitran96 ascii file molparam.txt containing mass 
% and abdundance parameters for the gas species used by Hitran
% The data is loaded and then saved in matlab format.
%
% (C) Ben Quine, University of Toronto 1999

if ~exist('filename')
   error('GENSPECT MLWREAD: No filename Specified')
end   
if exist(filename)~=2
   error('GENSPECT MLWREAD: No file found')
end
fid=fopen(filename,'r');
gas_to_find=1;
i=0;
while ~feof(fid)
	% look for molicule
   data=fgetl(fid);
   
   if length(data)>=6
      gas_check1 = gasnum(data(1:6));
      gas_check2 = findstr(data,num2str(gas_to_find));
      if gas_check1==gas_to_find & ~isempty(gas_check2)
         % found new gas
         disp(['Found Gas: ' num2str(gas_to_find)]);
         gas_found=gas_to_find;
         gas_to_find=gas_to_find+1;
         mlw.ison(gas_found)=0;
      else
         num_data = str2num(data);
         if length(num_data)==5
            mlw.ison(gas_found)=mlw.ison(gas_found)+1;
            i=i+1;
            mlw.ison_name(i) = num_data(1);
            mlw.abundance(i)= num_data(2);
            mlw.Q296(i) = num_data(3);
            mlw.gj(i) = num_data(4);
            mlw.mass(i) = num_data(5);
         end
      end
   end
end

% Close file
disp('Read in molparam.txt')
fclose(fid);
% Add an index to allow the lookup of required parameters
mlw.mol_index=[0 cumsum(mlw.ison(1:end-1))];
% add note for use of qt.mol_index
mlw.legend='Find parameters by mlw.???(qt.mol_index(Mol)+Iso)';
% Add data tracablility
qt.file = filename;
qt.date_str = datestr(now);
% save out data
save(savename,'mlw','-mat')
disp(['GENSPEC MLWREAD: Data saved in: ' savename]);