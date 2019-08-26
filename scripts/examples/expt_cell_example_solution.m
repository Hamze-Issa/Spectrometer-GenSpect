% EXPT_CELL_EXAMPLE is a demonstration script file for the GENSPECT line-by-line radiative transfer
% code which generates lines data and radiance for a Sun-like source viewed through a laboratory cell
% of length 1m containing CO and H2O at RTP (101325 Pa, 298 K).
%
% (C) Ben Quine, University of Toronto 05-MAR-2001
%

gas{1}='co';   			% The string name of the hitran gas
gas{2}='n2o';   		% The string name of the hitran gas
iso=[0 0]';				% The isotope 'code' for the gas: '1' will be the most abundant
						% Isotope code '0' corresponds to all isotopes
mixrat{1} = 0.002;       % Define mixing ratios for CO
mixrat{2} = 0.01;        % Define mixing ratios for N2O

theta=0;				% the angle between reflected radiation and the zenith
interp_acc = 1;         % specify 1% interpolation accuracy
ang_size=1;             % Angular FOV [sr]

%% Compute Absorption Coefficient

% define calculation
gencalc=gengrid(2150,2162,0.001,25,'voigt','hitran',6,interp_acc); 		%'grid generation'
% define gas cell
P=101325;
T=298;
cell_length=1.0;
cells=exptcell(P,T,cell_length,gas,iso,mixrat,'hitran2000');
% Do line-by-line calculation
cells=k_calc(gencalc,cells);

%% Compute path radiance

% Define a source
T_source = 6900;
paths = path_source(T_source,cells.wavnum);
% Add gas cells to path
gases_to_use = [1 2];
cell_order = [1];
paths=path_atm(cells,1,gases_to_use,cells.mixrat(cell_order,gases_to_use),theta,ang_size,paths);
% Compute radiance
paths=radiance(paths,1,'by path');

mixrat{1} = 0.001;       % Define mixing ratios for CO
mixrat{2} = 0.01;        % Define mixing ratios for N2O

theta=0;				% the angle between reflected radiation and the zenith
interp_acc = 1;         % specify 1% interpolation accuracy

%% Compute Absorption Coefficient

% define calculation
gencalc=gengrid(2150,2162,0.001,25,'voigt','hitran',6,interp_acc); 		%'grid generation'
% define gas cell
P=101325;
T=298;
cell_length=1.0;
cells=exptcell(P,T,cell_length,gas,iso,mixrat,'hitran2000');
% Do line-by-line calculation
cells=k_calc(gencalc,cells);


% build correlation radiometer to select CO
gases_to_use = [1];
cell_order = [1];
paths2=path_atm(cells,1,gases_to_use,cells.mixrat(cell_order,gases_to_use),theta,ang_size,paths);
% Compute radiance through 'first cell and then correlelation cell
paths2=radiance(paths2,1,'by path');

% Plot results
figure(1)
plot(paths.wavnum,paths.radiance,'b'),grid
xlabel('Wavenumber [cm-1]')
ylabel('Radiance  [Wm^{-2}sr^{-1}(cm^{-1})^{-1}]')
title(['Radiance through laboratory cell']);

figure(3)
plot(paths.wavnum,paths.radiance-paths2.radiance,'b'),grid
xlabel('Wavenumber [cm-1]')
ylabel('Radiance  Difference [Wm^{-2}sr^{-1}(cm^{-1})^{-1}]')
title(['Radiance Difference through laboratory cell']);


gencalc.cpu_time = cputime - gencalc.cpu_time;
disp(['The whole run took ' num2str((gencalc.cpu_time)/60) ' minutes'])

%% Suggested Excersises

% (1)       Change h2o to n2o recompute radiance
% (2)       Adjust mixing ratios to make both gases significant
%
% Building a simple Correleation Radiometer
%
% This this simple example we are going to build a simple test atmosphere comprising 
% two gases and then see if we can select one of the gases absorption using a filter cell.
% To select the gas we employ the correlation radiometer technique of making two measurements
% one with a correlation gas cell in front of the test atmosphere and one without. We then
% difference the radiance measurements to obtain a correlated measurement.
% 
%
% (3)       Form an additional path element comprising only one of the gases
% (4)       Compute the radiance through first cell and second 'correlation' cell
% (5)       Plot the radiance difference between a path with no correlation cell and one with a correlation cell
% (6)       Switch the correlation gas. Does this simple correlation radiometer work?
% (7)       Now double the amount of CO in the atmosphere while maintaining the amount in the correlation cell.
%           Can you detect the change?
% (8)       The difference signal from each line should be symmetric. Why is it asymmetric?


