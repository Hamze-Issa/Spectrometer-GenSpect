% PRESS_VS_ALT_EXAMPLE is a demonstration script file which generates lines data for CO
%
% (C) Ben Quine, University of Toronto 22-APR-1999
% 		Revised 18-MAY-2000
%

gas{1}='co';   			% The string name of the hitran96 gas
iso=[1]';					% The isotope 'code' for the gas: '1' will be the most abundant
								% Isotope code '0' corresponds to all isotopes
theta=0;						% the angle between reflected radiation and the zenith
ang_size=1;                     % Angular FOV size [sr] 

gencalc=gengrid(2150,2155,0.01,100,'voigt','hitran',6,0.1); 		%'grid generation'

% load zpt structure containing altitude, pressure and temperature data
disp('Read in the atmospheric profile file..')
zpt=atmread('zpt_us_std_jun40.prf');
disp('Read the mixing ratio profile file...')

% load a mixing ratio structure containing mixrat{}.z and mixrat{}.conc
mixrats{1}=concread('mixrat_co.mxr');

% Generate a set of cells to represent an atmosphere
n_cells=input('How many layers to compute for? ');
disp('Generate Atmospheric Cells for Calculation...')
cells_z=atmcell(n_cells,50,000,45000,zpt,gas,iso,mixrats,'z');
cells_p=atmcell(n_cells,50,000,45000,zpt,gas,iso,mixrats,'p');

% compute the absorption coefficients for these cells
disp('Calculate the centre-line and the wing-line k-vectors...')
cells_to_calc = 1:n_cells;
cells_z=k_calc(gencalc,cells_z,cells_to_calc);
cells_p=k_calc(gencalc,cells_p,cells_to_calc);

% Calculate transmission along a path using absorption coefficients representing the atmosphere
disp('Generate an atmospheric path for radiance calculations...')
gases_to_use = [1];
cell_order=[n_cells:-1:1];
% Altitude Division
paths_z=path_atm(cells_z,cell_order,gases_to_use,cells_z.mixrat(cell_order,gases_to_use),theta,ang_size);
paths_z=radiance(paths_z,1,'by path');
% Pressure division
paths_p=path_atm(cells_p,cell_order,gases_to_use,cells_p.mixrat(cell_order,gases_to_use),theta,ang_size);
paths_p=radiance(paths_p,1,'by path');

% Plot out the results and difference results
figure(1)
clf
subplot(211)
plot(paths_z.wavnum,paths_z.radiance,'b',paths_p.wavnum,paths_p.radiance,'r'),grid
xlabel('Wavenumber (cm-1)')
ylabel('Radiance')
title(['GENSPECT Calculation for ' gas{1} ': Isotope Comparison']);

subplot(212)
plot(paths_z.wavnum,paths_z.radiance-paths_p.radiance,'b'),grid
xlabel('Wavenumber (cm-1)')
ylabel('Radiance Difference')

worst_case_error_pcnt = (max((paths_z.radiance-paths_p.radiance)./paths_z.radiance))*100

n_mols_z = sum(paths_z.umean.*paths_z.thickness)./1e4
n_mols_p = sum(paths_p.umean.*paths_p.thickness)./1e4

n_mol_diff_percent = (n_mols_p-n_mols_z)/n_mols_z*100

% Compute the time for the example
gencalc.cpu_time = cputime - gencalc.cpu_time;
disp(['The whole run took ' num2str((gencalc.cpu_time)/60) ' minutes'])
