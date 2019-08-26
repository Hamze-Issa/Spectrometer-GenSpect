function atmos=atmread(zptfile)

% ATMOSP = ATMREAD(zptfile) reads and saves out an atmospheric profile.
% The profile is placed in a struct array called ATMOSP and saved in
% genfile.atm  
%
% ATMREAD calls functions PROFILE and stores it in structure ATMOSP
%
% ZPTFILE is the (string) name and path of the data file
%
% (C) Ben Quine, 2000.

% open the input file,read the title,format,and number of data lines 
fid=fopen(zptfile);
titl=fgetl(fid);
formatf=fgetl(fid);
numlins=str2num(fgetl(fid));

% read the profiles in line by line into an array
profile=zeros(numlins,3);
i=0;
while 1 
   i=i+1;
   line1=fgetl(fid);
   if ischar(line1)==0,break,end
   pf=str2num(line1);
   pf(2) = pf(2)* 100;
   pf(3) = pf(3) + 273.15;
   profile(i,1:3)=pf;
end
fclose(fid);
if i-1~=numlins
   disp(['GENSPECT ATMREAD: Number of atmosphere lines, ' num2str(i) ' does not match value claimed in file, ' num2str(numlins)])
end

%using format, assign the parts of the array (z,p,t) their own vectors
if(findstr(lower(formatf),'ztp'))
   znum=1;,pnum=3;,tnum=2;
else %default to zpt
   znum=1;,pnum=2;,tnum=3;
end

%store the profile in three separate vectors
z1=profile(:,znum);
p1=profile(:,pnum);
t1=profile(:,tnum);

%Put the vectors 'back' into a standard form zpt structure
atmosp=struct('z',z1,'p',p1,'t',t1);
atmosp.z_legend = 'altitude in metres';
atmosp.p_legend = 'pressure in Pascals';
atmosp.T_legend = 'temperature in Kelvin';
atmos=atmosp;
