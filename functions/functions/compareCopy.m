load('variables.mat');
load('constants.mat');
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

numberOfCells = 1;
numberOfGases = length(gases_to_use);
numberOfMixRatiosPerGas = 11;
cells_to_calc = [1:numberOfCells];

constants.gasNames.genhit2000_02 = load('genhit2000_02.lin', '-mat');
constants.gasNames.genhit2000_01 = load('genhit2000_01.lin', '-mat');
constants.gasNames.genhit2000_07 = load('genhit2000_07.lin', '-mat');
constants.gasNames.genhit2000_05 = load('genhit2000_05.lin', '-mat');

%% Start Loop
clear mixrats
mixrats = cell(1,numberOfGases);
for i=1:numberOfGases
    mixrats{i}.z = variables.mixRatiosZ;
end
combinations.counter = 0;
combinations.mixrats = cell(1,1);
combinations = recursiveCopy(variables, mixrats, 0, constants, combinations);
pathsMaster.counter = 0;
pathsMaster.paths = cell(2,105);
pathsMaster = looper(combinations, constants, gencalc, pathsMaster);
save('Aug_14');