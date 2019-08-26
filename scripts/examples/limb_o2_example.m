% LIMB_VIEW_EXAMPLE is a demonstration script file for the GENSPECT line-by-line radiative transfer
% toolbox which generates lines data for O2. This script calculates transmission along a limb path
% at a defined tangent altitude using path_limb
%
% (C) Ben Quine, University of Toronto 15-NOV-2001

gas{1}='o2';   			        % The string name of the hitran96 gas
iso=[1]';					    % The isotope 'code' for the gas: '1' will be the most abundant
								% Isotope code '0' corresponds to all isotopes

theta=0;						% the angle between reflected radiation and the zenith
ang_size=1;                     % Angular FOV size [sr] 

gencalc=gengrid(12960,13000,0.001,100,'voigt','hitran',6,1); 	%'grid generation'
%gencalc=gengrid(2150,2162,0.001,20,'lorentz','hitran'); 		%'grid generation'

disp('Read in the atmospheric profile file..')
zpt=atmread('zpt_us_std_jun40.prf');
disp('Read the mixing ratio profile file...')

mixrats{1}=concread('mixrat_o2.mxr');

disp('Generate Atmospheric Cells for Calculation...')
cells=atmcell(20,100,40000,90000,zpt,gas,iso,mixrats,'z');

disp('Calculate the centre-line and the wing-line k-vectors...')
cells_to_calc = [1:20]';
cells=k_calc(gencalc,cells,cells_to_calc);

disp('Generate an atmospheric path for radiance calculations...')
gases_to_use = [1];
cell_order=[1:cells.ncells]';

tangent_height = 40000;					% limb view tangent height
paths=path_limb(cells,cell_order,gases_to_use,cells.mixrat(cell_order,gases_to_use),tangent_height,ang_size);
paths=radiance(paths,1,'by path');

figure(1)
clf
plot(paths.wavnum,paths.transmission,'b'),grid

%axis([2150 2160 0 0.4e-6])
xlabel('Wavenumber (cm-1)')
ylabel('Transmission')
title(['GENSPECT limb path transmission calculation for ' gas{1} ' and tangent path at ' num2str(tangent_height/1e3) ' km']);

gencalc.cpu_time = cputime - gencalc.cpu_time;
disp(['The whole run took ' num2str((gencalc.cpu_time)/60) ' minutes'])

%% Suggested Excersices

% (1)       Adjust the tangent height of the radiance calculation to see the effect on tansmission.
% (2)       Increase the temperature of the atmosphere by 10% [not physical!] to see effect.
% (3)       Increase the pressure of the atmosphere by 1% [not physical!] to see effect.
% (4)       Use the axis command to zoom in on a line to see lineshape
% (5)       Add a sunlike source and plot radiance transmitted
% (6)       Include an additional gas
