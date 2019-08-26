
gas{1}='co2';   			    % The string name of the hitran96 gas
gas{2}='h2o';
gas{3}='o2';
gas{4}='co';
gases_to_use = [1 2 3 4];
iso = [0 0 0 0];		        % Isotope code '0' corresponds to all isotopes

theta=0;						% a viewing angle of 0 indicates a nadir view (the only tested view)
zenith=0;
ang_size=1;                     % Angular FOV size [sr] 

% Include Reflection
surface='smooth';		%
r=0.4;					% Define surface parameters here
e=0.6;					%
Tearth=300;
Tsun=6000;

startWaveLength = 800; % 200 in nanometers (wavenumber = 50000 cm-1)
endWaveLength = 1700;  % 2500 in nanometeres (wavenumber = 4000 cm-1)    
startWaveNumber = l2nu(endWaveLength)
endWaveNumber = l2nu(startWaveLength)
gencalc=gengrid(startWaveNumber,endWaveNumber,0.005,20,'voigt','hitran', 2); 		%'grid generation'

disp('Read in the atmospheric profile file..')
zpt=atmread('zpt_uae_May_2019.prf');
disp('Read the mixing ratio profile file...')

mixRatsMesured{1}=concread('mixrat_co.mxr');
mixRatsMesured{2}=concread('mixrat_co2.mxr');
mixRatsMesured{3}=concread('mixrat_h2o.mxr');
mixRatsMesured{4}=concread('mixrat_o2.mxr');

numberOfCells = 1;
numberOfGases = length(gases_to_use);
numberOfMixRatiosPerGas = 11;
cells_to_calc = [1:numberOfCells];

%% Start Loop
clear mixrats
mixrats = cell(1,numberOfGases);
for i=1:numberOfGases
    mixrats{i}.z = variables.mixRatiosZ;
end
pathsMaster.counter = 0;
pathsMaster.paths = cell(2,1);
pathsMaster = recursive(variables, mixrats, 0, constants, gencalc, pathsMaster);

%% test loop
% counter = 0;
% clear mixrats
% mixrats = cell(1,numberOfGases);
% pathsMaster.counter = 0;
% pathsMaster.paths = cell(2,1);
% pathsMaster = test(variables, mixrats, 0, constants, gencalc, pathsMaster);
