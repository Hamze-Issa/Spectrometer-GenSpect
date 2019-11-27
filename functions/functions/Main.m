%{
The MeznSat Project
Atmospheric Retrieval and Data Processing
Hamzeh Issa

This is the main script for the atmospheric retrieval process. 
This process requires high capability computer, like the Masdar HPC.

This script loads the variables, gases, isotopes, parameters, wavelength
ranges, and mix ratios of the gases. and then uses Looper function to apply
the atmospheric model code on each mix ratio to generate radiance to
wavelength relationship. You can then find the best combination (with the
lease MSE) using errorCalculation function.
%}

% Load variables and constants. They are structures that contain most of the parameters required
% parameters that need tuning are listed in variables while parameters that
% don't are listed in constants.
load('variables.mat');
load('constants.mat');

% Load gases to be considered in the simulation, for the full list of
% gases, check "gasname.m" under functions.
gas{1}='co2';   			    % The string name of the hitran96 gas
gas{2}='h2o';
gas{3}='o2';
gas{4}='co';
gases_to_use = [1 2 3 4];       % all 4 gases
iso = [0 0 0 0];		        % Isotope code '0' corresponds to all isotopes

% declare the angles
theta=0;						% a viewing angle of 0 indicates a nadir view (the only tested view)
zenith=0;
ang_size=1;                     % Angular FOV size [sr] 

% Include surface parameters and temperatures
surface='smooth';		
r=0.4;					% Define surface parameters here
e=0.6;					% Mostly (1 - r)
Tearth=300;
Tsun=6000;

% Enter the wavelength range. It is then turned into wavenumber and fed to
% the grid generation function. You can enter wavenumbers directly, but I thought
% wavelengths might be more convenient.
startWaveLength = 800; % 200 in nanometers (wavenumber = 50000 cm-1)
endWaveLength = 1700;  % 2500 in nanometeres (wavenumber = 4000 cm-1)    
startWaveNumber = l2nu(endWaveLength)
endWaveNumber = l2nu(startWaveLength)
gencalc=gengrid(startWaveNumber,endWaveNumber,0.005,20,'voigt','hitran', 2); 		% grid generation

% Read atmospheric profile data (temperature T and pressure P with altitude Z)from the given file, this file needs to be
% updated for each simulation to get accurate data. It can be acquired from
% online databases.
disp('Read in the atmospheric profile file..')
zpt=atmread('zpt_uae_May_2019.prf');
disp('Read the mixing ratio profile file...')

% number of cells is basically the number of atmospheric layers seen by the
% model, for simplicity and speed, it is set to one, if results are not
% satisfactory, it can be increased.
numberOfCells = 1;
numberOfGases = length(gases_to_use);

% numberOfMixRatiosPerGas is the number of mix ratios per gas to be
% permutated with each other to get all the possible combinations, in this
% case, the number of combinations = numberOfMixRatiosPerGas ^ number of
% gases = 11 ^ 4 = 14,641? combinations, which we will need to apply the
% atmospheric model to each of them and get the radiance/ wavelength
% relationship for each of them. That's why this simulation requires long
% time and lots of computational power. Of course, adding more combinations
% (by increasing numberOfMixRatiosPerGas) increases accuracy and requires
% more time.
numberOfMixRatiosPerGas = 11;
cells_to_calc = [1:numberOfCells];

% Load the gas data for each gas from the Hitran database and add them to
% the constants structure.
constants.gasNames.genhit2000_02 = load('genhit2000_02.lin', '-mat');
constants.gasNames.genhit2000_01 = load('genhit2000_01.lin', '-mat');
constants.gasNames.genhit2000_07 = load('genhit2000_07.lin', '-mat');
constants.gasNames.genhit2000_05 = load('genhit2000_05.lin', '-mat');

%% Start Loop

% mixrats holds the 11 mix ratios for each 4 gases
clear mixrats
mixrats = cell(1,numberOfGases);

% Just a litte loop to add the altitude for the mix ratios for compatibility
% purposes only.
for i=1:numberOfGases
    mixrats{i}.z = variables.mixRatiosZ;
end

% combinations holds all the 11^4 combinations for the possible gas mix
% ratios. They are generated in a recursive loop using generateCombinations
% function.
combinations.counter = 0;
combinations.mixrats = cell(1,1);
combinations = generateCombinations(variables, mixrats, 0, constants, combinations);

% pathsMaster is the big bad structure that holds all the outputs: all the
% paths, their parameters and the radiance/wavelength relationship for each
% mix ratio combination. It is defined, then the big bad function "Looper"
% loops over all the combinations applying the atmospheric model to each
% one of them, and packing the results into pathsMaster.
pathsMaster.counter = 0;
pathsMaster.paths = cell(2,105);
pathsMaster = looper(combinations, constants, gencalc, pathsMaster);
save('November_24');