function qt=qtread(filename,savename)

% qt=qtread(filename,savename)
%
% This function reads Hitran96 ascii files containing estimation parameters for the
% Total internal partition function. The data is loaded and then saved in matlab
% format.
%
% (C) Ben Quine, 1999

global qt

if ~exist('filename')
   error('GENSPECT QTREAD: No filename Specified')
end   
if exist(filename)~=2
   error('GENSPECT QTREAD: No file found')
end
fid=fopen(filename,'r');
% look for molicule and species counter
data=fgetl(fid);
while isempty(findstr(upper(data),'NMOL'))
   data=fgetl(fid);
end
% read line
data=str2num(fgetl(fid));
if length(data)~=2
   error('GENSPECT QTREAD: Line should contain NMOLS and NSPECIES values')
else
   qt.nmol=data(1);
   qt.nspecies=data(2);
end
% look for TMIN/TMAX
data=fgetl(fid);
while isempty(findstr(upper(data),'TMIN'))
   data=fgetl(fid);
end
% read line
data=str2num(fgetl(fid));
if length(data)~=2
   error('GENSPECT QTREAD: TMIN/TMAX line should be followed by 2 entries')
else
   qt.tmin=data(1);
   qt.tmax=data(2);
end

% look for number of ranges
data=fgetl(fid);
while isempty(findstr(upper(data),'NRANG')) & ~feof(fid)
   data=fgetl(fid);
end
qt.n_coeff=str2num(fgetl(fid));

% look for temperature limits
data=fgetl(fid);
while isempty(findstr(upper(data),'MAXT')) & ~feof(fid)
   data=fgetl(fid);
end
for i=1:qt.n_coeff
	qt.t_limits(i)=str2num(fgetl(fid));
end
% look for isotope division
data=fgetl(fid);
while isempty(findstr(upper(data),'ISO')) & ~feof(fid)
   data=fgetl(fid);
end
% read isotope numbers
for i=1:qt.nmol
   qt.ison(i)=str2num(fgetl(fid));
end

% look for Q296
data=fgetl(fid);
while isempty(findstr(upper(data),'Q296'))  & ~feof(fid)
   data=fgetl(fid);
end
% read Q296 values
for i=1:qt.nspecies
   qt.Q296(i)=str2num(fgetl(fid));
end

% look for QCOEFF
data=fgetl(fid);
while isempty(findstr(upper(data),'QCOEF')) & ~feof(fid)
   data=fgetl(fid);
end
% read Coeff values
for j=1:qt.n_coeff
	for i=1:qt.nspecies
	   qt.coeff(j,i,:)=str2num(fgetl(fid));
	end
   data=fgetl(fid);
	while isempty(findstr(upper(data),'QCOEF')) & ~feof(fid)
	   data=fgetl(fid);
	end
end
% Close file
disp('Read in Total internal partition temperature coefficents')
fclose(fid);
% Add an index to allow the lookup of required parameters
qt.mol_index=[0 cumsum(qt.ison(1:end-1))];
% add note for use of qt.mol_index
qt.legend='Find parameters by qt.coeff(T,[qt.mol_index(Mol)+Iso],1:4)';
% Add data tracablility
qt.file = filename;
qt.date_str = datestr(now);
% save out data
save(savename,'qt','-mat')
disp(['GENSPECT QTREAD: Data saved in: ' savename]);