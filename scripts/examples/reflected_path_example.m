% REFLECTED_PATH_EXAMPLE is an example script showing how to incorporate the RADIANCE function
% into GENSPECT calculations. REFLECTED_PATH computes the radiance from the sun, 
% earth and atmospheric cells along the path specified by ATMPATH. If the path is 
% not reflected, the value of the radiance returned is that which is indident on 
% the earth's surface. 
%
% NOTE: This script does not include a diffusivity factor to correctly weight
%       reflected emission and absorption contributions correctly.
%
% (C) University of Toronto, 1999.

gas{1}='co2';   			% The string name of the hitran96 gas
gas{2}='h2o';
gas{3}='o2';
gas{4}='co';

iso=[0 0 0 0]';					% The isotope 'code' for the gas: '1' will be the most abundant
								% Isotope code '0' corresponds to all isotopes
theta=0;						% a viewing angle of 0 indicates a nadir view (the only tested view)
zenith=0;
ang_size=1;                     % Angular FOV size [sr] 

startWaveLength = 800; % 200 in nanometers (wavenumber = 50000 cm-1)
endWaveLength = 1700;  % 2500 in nanometeres (wavenumber = 4000 cm-1)    
startWaveNumber = l2nu(endWaveLength)
endWaveNumber = l2nu(startWaveLength)
%gencalc=gengrid(2150,2160,0.005,20,'voigt','hitran'); 		%'grid generation'
gencalc=gengrid(startWaveNumber,endWaveNumber,0.005,20,'voigt','hitran', 2); 		%'grid generation'

% gencalc = gengrid(start_wavnum, stop_wavnum, calc_resolution, wingcutoff,type,model)
% Output:
% gencalc.calc_grid is the calculation spectral grid vector
% gencalc.wing_grid is the wing calculation spectral grid vector
% gencalc.core_grid determines the division between center and wing calcultions 
% gencalc.wingcutoff limits the number coarse intervals included in the wing calculation

disp('Read in the atmospheric profile file..')
zpt=atmread('zpt_uae_May_2019.prf');
disp('Read the mixing ratio profile file...')

mixrats{1}=concread('mixrat_co.mxr');
mixrats{2}=concread('mixrat_co2.mxr');
mixrats{3}=concread('mixrat_h2o.mxr');
mixrats{4}=concread('mixrat_o2.mxr');

disp('Generate Atmospheric Cells for Calculation...')
cells=atmcell(6,100,27,16590,zpt,gas,iso,mixrats,'z');

disp('Calculate the centre-line and the wing-line k-vectors...')
cells_to_calc = [1 2 3 4];
gases_to_use = [1 2 3 4];
[cells]=k_calc(gencalc,cells,cells_to_calc,gases_to_use);

disp('Generate an atmospheric path for radiance calculations...')
cell_order=[2 1];
% Start with sun-like source path segment
Tsun=6000;
paths=path_source(Tsun,cells.wavnum);
% Include Atmosphere (top-down)
paths=path_atm(cells,cell_order,gases_to_use,cells.mixrat(cell_order,gases_to_use),theta,ang_size,paths);
% Include Reflection
surface='smooth';		%
r=0.4;					% Define surface parameters here
e=0.6;					%
Tearth=300;
paths=path_reflect(r,e,0,0,ang_size,surface,Tearth,paths);
% Include path back through atmosphere
cell_order=[1 2];
paths=path_atm(cells,cell_order,gases_to_use,cells.mixrat(cell_order,gases_to_use),theta,ang_size,paths);

%Compute radiance
paths=radiance(paths,1,'by path');

gencalc.cpu_time = cputime - gencalc.cpu_time;
disp(['The whole run took ' num2str((gencalc.cpu_time)/60) ' minutes'])
figure;
clear waveLengths;
waveLengths = zeros(1, length(paths.wavnum));
for i = 1 : length(paths.wavnum)
    waveLengths(i) = nu2l(paths.wavnum(i));
end
plot(waveLengths,paths.radiance),grid;
xlabel('WaveLength (nm)');
ylabel('Total Radiance (W sr^{-1}m^{-2})');
title('Radiance calculated along a Sun-Atmosphere-Earth-Atmosphere path');
figure;
plot(paths.wavnum,paths.radiance),grid;
xlabel('Wavenumber (cm^{-1})');
ylabel('Total Radiance (W sr^{-1}m^{-2})');
title('Radiance calculated along a Sun-Atmosphere-Earth-Atmosphere path');
