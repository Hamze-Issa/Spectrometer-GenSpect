function gencalc=gengrid(start_wavnum,stop_wavnum,calc_resolution,wingcutoff,typ,model,no_grids,accuracy,disp_log)

% gencalc = GENGRID(start_wavnum, stop_wavnum, calc_resolution, wingcutoff,type,model,n_grids,accuracy,disp_log) 
% generates a structure gencalc which contains calculation information. The function
% divides the computational task into intervals and then determines where line calculations can
% may be interpolated on to a series of lower resolution grids.
%
%
% start_wavnum is the lowest wavenumber to be calculated
% stop_wavnum is the highest wavenumber to be calculated
% calc_resolution is the required calculation resolution
% wingcutoff is the number of waven numbers to include for the wing calculation
% typ is the lineshape type: 'lorentz' 'doppler' or 'voigt'
% model is the database type to use: 'hitran' [note that the database version is set in EXPTCELL or ATMCELL]
% n_grids (optional) sets the number of binary grid divisions used to reduce the
%			  calculation order. Range [1:10], defaults to 6 [must be set if disp_log is set].
% accuracy is the percentage calculation accuracy required [either 0.01% or 0.1% or 1% (0.1% default)]. 
% disp_log (optional) sets the display settings 'off','www','on'
%
% (C) Ben Quine, 06-APR-2000.

% Generate computation grids
if ~exist('no_grids')
   no_grids=6;
else
   if no_grids<1 or no_grids>10
      disp('GENSPECT GENGRID: Number of grids outside acceptable range')
      no_grids=max(min(no_grids,10),1);
   end
end
% set display log defaults
if exist('disp_log')==1
   if strcmp(lower(disp_log),'off')
      disp_log=0;
   elseif strcmp(lower(disp_log),'www')
      disp_log=2;
   else
      disp_log=1;
   end
else
   disp_log=1;
end

% Code indicator
gencalc.code = 'GENSPECT: www.genspect.com';
% insert version number
gencalc.version = 'Version: 1.2 (12-DEC-2015)';
if disp_log==1
   gencalc
end
if ~exist('accuracy')
   gencalc.accuracy=0.1;
else
   if accuracy==0.1 | accuracy==0.01 | accuracy==1
      gencalc.accuracy=accuracy;
   else
      disp('GENSPECT GENGRID: Accuracy not set to 0.01%, 0.1% or 1%  - Defaulting to 0.1%')
      gencalc.accuracy=0.1;
   end
end
% Line to clear globally defined grid_division data
clear global grid_division

min_resolution = 0.1;
% Define limits for calculation resolution
if calc_resolution>min_resolution
   disp('GENSPECT GENCALC: Calculation Resolution cannot be greater than 0.1 cm-1. Terminating.')
   return
end
% Check specified wavenumber range


% set maximium number of grids
max_no_grids=10;

% check specified calculation ranges
if start_wavnum>stop_wavnum
   error('GENSPECT GENCALC: START wavenmber must be greater than STOP wavenumber')
end
if no_grids>max_no_grids
   no_grids=max_no_grids;
   disp(['GENSPECT GENCALC: Number of Grids requested greater than inturnal limit - setting to ' num2str(max_no_grids)])
end
% compute number of calculation points and force the right number
n_pts = (stop_wavnum-start_wavnum)/calc_resolution + 1;
if fix(n_pts/(2^no_grids))*(2^no_grids)~=n_pts & no_grids>1
   disp('GENSPECT GENCALC: Number of calculation points must be a multiple of 2^(number of grids).')
   disp('GENSPECT GENCALC: Increasing upper limit of the calculation to enforce this.');
   n_pts=ceil(n_pts/(2^no_grids))*(2^no_grids);
end
% Record new limits
gencalc.start_wavnum=start_wavnum;
gencalc.stop_wavnum=start_wavnum+calc_resolution*(n_pts-1);
gencalc.num_grids=no_grids;

% Generate basic computation grid
gencalc.grid{1}=linspace(gencalc.start_wavnum,gencalc.stop_wavnum,n_pts);
gencalc.grid_interval(1)=calc_resolution;
gencalc.grid_pts(1)=length(gencalc.grid{1});
% Generate the grids that sequentially divide the wavenumber range 
if gencalc.num_grids>1
   for i=2:no_grids
      gencalc.grid_interval(i)=(2^i)*calc_resolution;
      % form basic grid
      gridpts=linspace(gencalc.start_wavnum+calc_resolution*((2^i)-1),gencalc.stop_wavnum,n_pts/(2^i));
%      [gencalc.start_wavnum:gencalc.grid_interval(i):gencalc.stop_wavnum];
      % add additional points at inturnal boundaries to ensure accurate interpolation
      gridpts=[gridpts gridpts-calc_resolution*((2^i)-1)];
      gencalc.grid{i}=sort(gridpts);
      gencalc.grid_pts(i)=length(gencalc.grid{i});
   end   
end
% Record grid creation information
gencalc.creation_date=datestr(now);
gencalc.title='GENSPECT GENCALC: calculation grids for line-by-line computations.';

% Calculate the bounds for the wing calculation [must be integer number of wing grid intervals]
gencalc.wingcutoff = wingcutoff;
gencalc.freq_legend = 'All frequencies are in wavenumbers (cm^{-1})';

gencalc.cpu_time=cputime;
gencalc.typ=typ;    % The lineshape function to be used: could be 'doppler','lorentz', 

								% or some user-created '.m' file.

gencalc.model=model; % The line model: only 'hitran' or 'elsass' (elsasser lines) are
								% available.

