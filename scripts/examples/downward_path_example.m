% DOWNWARD_PATH_EXAMPLE is a demonstration script file for the GENSPECT 
% line-by-line radiative transfer code which calculates transmission through CO2
% over a downward atmospheric path from 25 to 15 km.
%
% (C) Ben Quine, University of Toronto 16-NOV-2001

gas{1}='co2';   			% The string name of the hitran96 gas
gas{2}='h2o';
iso=[0 1]';					% The isotope 'code' for the gas: '1' will be the most abundant
							% Isotope code '0' corresponds to all isotopes

theta=0;					% the angle between reflected radiation and the zenith

startWaveLength = 600; % in nanometers
endWaveLength = 1800;  % in nanometeres
startWaveNumber = l2nu(endWaveLength)
endWaveNumber = l2nu(startWaveLength)
gencalc=gengrid(startWaveNumber,endWaveNumber,0.001,100,'voigt','hitran',6,0.1); 
%gencalc=gengrid(2100,2200,0.001,100,'voigt','hitran',6,0.1); 		%'grid generation'
%gencalc=gengrid(2150,2162,0.001,20,'lorentz','hitran'); 		%'grid generation'

disp('Read in the atmospheric profile file..')
zpt=atmread('zpt_us_std_jun40.prf');
disp('Read the mixing ratio profile file...')

mixrats{1}=concread('mixrat_co2.mxr');

disp('Generate Atmospheric Cells for Calculation...')
cells=atmcell(4,100,15000,25000,zpt,gas,iso,mixrats,'z');

disp('Calculate the centre-line and the wing-line k-vectors...')
cells_to_calc = [1:4];
cells=k_calc(gencalc,cells,cells_to_calc);

disp('Generate an atmospheric path for radiance calculations...')
gases_to_use = [1];
% Define reverse cell order [radiation pasing from top to bottom of atmosphere]
cell_order=[cells.ncells:-1:1];
ang_size=1;
% Generate path
paths=path_atm(cells,cell_order,gases_to_use,cells.mixrat(cell_order,gases_to_use),theta,ang_size);
% Compute Radiance
paths=radiance(paths,1,'by path');

% Plot Results
figure(1)
clf
plot(paths.wavnum,paths.transmission,'b'),grid
xlabel('Wavenumber (cm-1)')
ylabel('Transmission')
title(['GENSPECT Calculation for ' gas{1}]);

figure;
clear waveLengths;
waveLengths = zeros(1, length(paths.wavnum));
for i = 1 : length(paths.wavnum)
    waveLengths(i) = nu2l(paths.wavnum(i));
end
plot(waveLengths,paths.transmission),grid;
xlabel('WaveLength (nm)');
ylabel('Transmission');
title('Transmission');
gencalc.cpu_time = cputime - gencalc.cpu_time;
disp(['The whole run took ' num2str((gencalc.cpu_time)/60) ' minutes'])

