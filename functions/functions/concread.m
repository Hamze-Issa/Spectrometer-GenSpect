function mixrat=concread(concfile)

% CONCREAD(concfile) reads and saves out a volume mixing ratio profile.
% The profile is placed in a struct array called MIXRAT and saved in
% genfile.con
%
% CONCREAD calls no other genasis functions
%
% CONCFILE is the (string) name and path of the mixing ratio file
%
% (C) Ben Quine, 2000.

%open the input file,read the title,format,and number of data lines 
fid=fopen(concfile);
titl=fgetl(fid);
numlins=str2num(fgetl(fid));

%read the profiles in line by line into an array
profile=zeros(numlins,2);
i=0;
while 1 
   i=i+1;
   line1=fgetl(fid);
   if ischar(line1)==0,break,end
   pf=str2num(line1);
   profile(i,1:2)=pf;
end

fclose(fid);
if i-1~=numlins
   disp('GENSPECT CONCREAD: Number of atmosphere lines does not match value claimed in file')
end

%store the concentration profile in two separate vectors
z=profile(:,1);
ratio=profile(:,2);

%put the vectors into standard form structure
mixrat=struct('z',z,'conc',ratio);
