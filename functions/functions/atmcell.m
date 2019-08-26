function cells=atmcell(ncells,nlayers,zstart,zstop,atmosp,gas,iso,mixrat,division,database_name)

% cells=ATMCELL(ncells,nlayrs,zstart,zstop,atmosp,gas,iso,mixrat,division,database_name)
%
% ATMCELL generates a set of atmospheric CELLS that represent the atmosphere specified
% by structure ATMOSP. Cells are divided into either 'z', 'p' or 'other' increments 
% according to the input DIVISION. The function then performs the curtis-godson approximation 
% and finds mean pressure, temperature, and number density values for each cell. 
%
% ATMREAD should be run to generate structure ATMOSP before ATMCELL is run.
%
% Inputs:
% NCELLS is the number of coarse cells (scalar)
% NLAYRS is the number of fine layers per coarse cell -- sublayers (scalar)
% ZSTART is the bottom of the atmosphere being considered (scalar unless DIVISION 'other'
%			is chosen, then column vector) [m]
% ZSTOP  is the top of the atmosphere being considered (scalar unless DIVISION 'other'
%			is chosen, then column vector) [m]
% ATMOSP is the structure containing the zpt profile read by ATMREAD. 
%			'atmosp' contains (column vectors) atmosp.z, atmosp.p, and atmosp.t
% GAS is the gas name involved (cell array of strings)
% ISO is an array of isotopes to be used (column vector)
% MIXRAT is a structure describing the mixing ratio profile for the gas
%		for multigas problems MIXRAT becomes a cell array mixrat{1}, mixrat{2}
% DIVISION is string specifying how the cells will be formed: can be by height, pressure
%		or by a specific height vector. Default will be a division by height. (string)
%		(Must be one of 'z', 'p' or 'other') if 'other' is chosen then ZTOP and ZBOTTOM are vectors.
% DATABASE_NAME (optional) String identifying database to use ('hitran96', 'hitran2000') or structure array
%       of alternate gas parameter input filenames (string), one per gas [Defaults to Hitran2000]. 
%
% CELLS is a structure created with fields zbottom, ztop, pmean, tmean, umean, zmean,
%			gasname, gas, iso, ncells, ngases, mixrat, ppartial
%
% (C) Ben Quine, University of Toronto, 07-SEP-1999.

% Modifications: Ben Quine, 01-FEB-1999, Repared Density.
% 			     L.Buckley, 05-JUL-1999, Option to divide cells by pressure or by 'other' height intervals.
%			     Ben Quine  02-SEP-1999, Checked and fixed up, error checking on zstart and zstop.
%                E. McKernan 25-JAN-2001, 'other' option was broken (confused ztop-zstop-zbottom-zstart).
%                                         Different zpt interpolation options allowed (eg spline)
%                                         added param_filename option
%                BQ, 1-FEB-2001, Checked code mods.

% store line parameter database location for each gas

% Record parameter datafiles to use
% set default database if none specified
if ~exist('database_name','var')
    database_name='hitran2000';
end
[cells.param_filename, cells.database] = hitfilenames(database_name,gas);

% interpolation method. 'spline' works well except for zpt files with large gaps. 
% 'linear' is fairly safe but will not capture exponential pressure changes as well.
imethod = 'linear';

% set default
if exist('division')==0,
   division='z';
end
% ensure division is lower case
division = lower(division);
% check that zstart and zstop occur within defined atmospheric range
if max(atmosp.z)<zstop | min(atmosp.z)>zstart
   error('GENSPECT ATMCELL: Calculation start and stop heights fall outside specified atmosphere.');
end

% Code for an atmospheric specified by equal height intervals
if strcmp(division,'z') == 1,
   disp('GENSPECT ATMCELL: Dividing cells by height increments...')
   cells.type = 'Equal Altitude Division';
	cel=linspace(zstart,zstop,ncells+1)';
	cells.zbottom=cel(1:ncells);
	cells.ztop=cel(2:ncells+1);

	%for each cell, generate homogenous layers and interpolate pressure
	%and temperature on them
	for i=1:ncells
		layer.z=linspace(cells.zbottom(i),cells.ztop(i),nlayers)';
		layer.p=interp1(atmosp.z,atmosp.p,layer.z,imethod);
		layer.t=interp1(atmosp.z,atmosp.t,layer.z,imethod);
		%now perform the curtis godson approximation on them
		ptrho=curtsgod(layer);
      cells.pmean(i)=ptrho.p;
      cells.tmean(i)=ptrho.t;
      cells.umean(i)=ptrho.u;
      cells.zmean(i)=ptrho.z;
   end
   
% code for an atmosphere specified by equal pressure intervals
elseif strcmp(division,'p') == 1,
   disp('GENSPECT ATMCELL: Dividing cells by pressure increments...')
   cells.type = 'Equal Pressure Division';
   pstart = interp1(atmosp.z,atmosp.p,zstart,imethod);
	pstop = interp1(atmosp.z,atmosp.p,zstop,imethod);   

	cel=linspace(pstart,pstop,ncells+1)';
	cells.pbottom=cel(1:ncells);
	cells.ptop=cel(2:ncells+1);
	cells.zbottom=interp1(atmosp.p,atmosp.z,cells.pbottom,imethod);
	cells.ztop=interp1(atmosp.p,atmosp.z,cells.ptop,imethod);

	% For each cell, generate homogenous layers and interpolate height
	% and temperature on them.
	for i=1:ncells
		layer.p=linspace(cells.pbottom(i),cells.ptop(i),nlayers)';
		layer.z=interp1(atmosp.p,atmosp.z,layer.p,imethod);
		layer.t=interp1(atmosp.p,atmosp.t,layer.p,imethod);
		%now perform the curtis godson approximation on them
		ptrho=curtsgod(layer);
      cells.pmean(i)=ptrho.p;
      cells.tmean(i)=ptrho.t;
      cells.umean(i)=ptrho.u;
      cells.zmean(i)=ptrho.z;
	end
% code for user defined division by height   
elseif strcmp(division,'other') == 1,
   disp('GENSPECT ATMCELL: Dividing cells by specified height intervals...')
   cells.type='User Specified Divisions';
   if length(zstop)~=ncells | length(zstart)~=ncells
      error('GENSPECT ATMCELL: ZTOP and ZBOTTOM must be vectors of length NCELLS');
   end
  	cells.zbottom=zstart;
   cells.ztop=zstop;
	
	% for each cell, generate homogenous layers and interpolate height and temperature 
	for i=1:ncells
		layer.z=linspace(cells.zbottom(i),cells.ztop(i),nlayers)';
		layer.p=interp1(atmosp.z,atmosp.p,layer.z,imethod);
		layer.t=interp1(atmosp.z,atmosp.t,layer.z,imethod);
		%now perform the curtis godson approximation on them
		ptrho=curtsgod(layer);
      cells.pmean(i)=ptrho.p;
      cells.tmean(i)=ptrho.t;
      cells.umean(i)=ptrho.u;
      cells.zmean(i)=ptrho.z;
	end
   
else
   error('GENSPECT ATMCELL: invalid method for cell DIVISION.')
end
   
%construct the data structure of cell headers
cells.pmean=cells.pmean';
cells.tmean=cells.tmean';
cells.umean=cells.umean';
cells.zmean=cells.zmean';
% list gas in each cell
cells.creation_date=datestr(now);
cells.gasname=gas;
cells.gas = gasnum(gas);
cells.iso = iso;
cells.ncells = ncells;
cells.ngases = length(cells.gas);

% error checking for gas/iso to ensure the same number of entries
if cells.ngases~=length(cells.iso)
  	error('GENSPECT ATMCELL: The number of gases and isotopes specified are different.')
end

% For each gas, attribute mixing ratio associated with each cell
if exist('mixrat')==1	% interpolate mixing ratios at mean density point in layer
  	% determine if mixrat is a strcuture or a cell array
  	data_type=whos('mixrat');
	if strcmp(data_type.class,'cell')
			for i=1:cells.ngases
	   		cells.mixrat(:,i)=interp1(mixrat{i}.z,mixrat{i}.conc,cells.zmean,imethod); 
  			end      
  	elseif strcmp(data_type.class,'struct')
     	cells.mixrat=interp1(mixrat.z,mixrat.conc,cells.zmean,imethod); 
  	else
     	warning('GENSPECT ATMCELL: Unrecognised Data Class for mixrat: Mixing ratios set to zero');
   	cells.mixrat=zeros(ncells,cells.ngases);
  	end
else						% set mixing ratios to zero if they don't exist
  	cells.mixrat=zeros(ncells,cells.ngases);
end
if size(cells.mixrat,2)==cells.ngases
	% record partial pressure
	cells.ppartial=cells.pmean*ones(1,cells.ngases).*cells.mixrat;
else
  	warning('GENSPECT ATMCELL: Mixing ratio undefined for multi gas problem. cells.ppartial not computed.')
end
