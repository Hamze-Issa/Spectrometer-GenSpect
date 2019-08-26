% BENCHTST is a demonstration script file for the GENSPECT line-by-line radiative transfer
% code which generates lines data for CO and H2O and computes atmospheric emission radiance
% from a 15-20km layer using a two cell approximation.
%
% (C) Ben Quine, University of Toronto 05-MAR-2001
%

gas{1}='co';   			% The string name of the hitran96 gas
gas{2}='h2o';   		% The string name of the hitran96 gas
iso=[1 0]';			    % The isotope 'code' for the gas: '1' will be the most abundant
						% Isotope code '0' corresponds to all isotopes

theta=0;				% the angle between reflected radiation and the zenith
interp_acc = 1;              % specify 1% interpolation accuracy

startWaveLength = 200; % in nanometers
endWaveLength = 1700;  % in nanometeres
startWaveNumber = l2nu(endWaveLength)
endWaveNumber = l2nu(startWaveLength)
gencalc=gengrid(2150,endWaveNumber,0.002,100,'voigt','hitran',6,interp_acc); 		%'grid generation'
%gencalc=gengrid(2150,2162,0.001,20,'lorentz','hitran'); 		%'grid generation'

disp('Read in the atmospheric profile file..')
zpt=atmread('zpt_us_std_jun40.prf');
disp('Read the mixing ratio profile file...')

mixrats{1}=concread('mixrat_co.mxr');
mixrats{2}=concread('mixrat_h2o.mxr');

disp('Generate two cell Calculation...')
cells=atmcell(2,100,15000,20000,zpt,gas,iso,mixrats,'z');  % using default database
%cells=atmcell(2,100,15000,20000,zpt,gas,iso,mixrats,'z','hitran96');
%cells=atmcell(2,100,15000,20000,zpt,gas,iso,mixrats,'z','hitran2000');

disp('Calculate the centre-line and the wing-line k-vectors...')
cells_to_calc = [1 2];
cells=k_calc(gencalc,cells,cells_to_calc);

figure(1)
clf
disp('Generate an atmospheric path for radiance calculations...')
gases_to_use = [1 2];
cell_order=[cells.ncells:-1:1];
ang_size=1;                         % set ang_size to unity for units of W/m2/cm-1/sr

paths=path_atm(cells,cell_order,gases_to_use,cells.mixrat(cell_order,gases_to_use),theta,ang_size);
paths=radiance(paths,1,'by segment');
plot(paths.wavnum,paths.radiance,'b'),grid

%axis([2150 2160 0 0.4e-6])
xlabel('Wavenumber (cm-1)')
ylabel('Radiance')
title(['GENSPECT Calculation for ' gas{1} ': Isotope Comparison']);

gencalc.cpu_time = cputime - gencalc.cpu_time;
disp(['The whole run took ' num2str((gencalc.cpu_time)/60) ' minutes'])
