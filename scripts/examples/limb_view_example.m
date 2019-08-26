% LIMB_VIEW_EXAMPLE is a demonstration script file for the GENSPECT line-by-line radiative transfer
% toolbox which generates lines data for CO. This script calculates transmission along a limb path
% at a defined tangent altitude using path_limb
%
% (C) Ben Quine, University of Toronto 22-APR-1999
% 		Revised 14-APR-2000
%

gas{1}='co';   			% The string name of the hitran96 gas
gas{2}='co';   			% The string name of the hitran96 gas
iso=[1 2]';					% The isotope 'code' for the gas: '1' will be the most abundant
								% Isotope code '0' corresponds to all isotopes

% Note: Matgen 1.1 introduces an additional angle to allow calculation in reflected
%	     path problems. Problems with no reflection should set theta to zero. 
                        
theta=0;						% the angle between reflected radiation and the zenith
ang_size=1;                     % Angular FOV size [sr] 

gencalc=gengrid(2150,2162,0.001,100,'voigt','hitran',6,0.1); 		%'grid generation'
%gencalc=gengrid(2150,2162,0.001,20,'lorentz','hitran'); 		%'grid generation'

disp('Read in the atmospheric profile file..')
zpt=atmread('zpt_us_std_jun40.prf');
disp('Read the mixing ratio profile file...')

mixrats{1}=concread('mixrat_co.mxr');
mixrats{2}=concread('mixrat_co.mxr');

disp('Generate Atmospheric Cells for Calculation...')
cells=atmcell(20,100,25000,50000,zpt,gas,iso,mixrats,'z');

disp('Calculate the centre-line and the wing-line k-vectors...')
cells_to_calc = [1:20]';
cells=k_calc(gencalc,cells,cells_to_calc);

disp('Generate an atmospheric path for radiance calculations...')
gases_to_use = [1 2];
cell_order=[1:cells.ncells]';

tangent_height = 25000;					% limb view tangent height
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
