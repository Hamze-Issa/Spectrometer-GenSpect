
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
cells_to_calc = [1:numberOfCells];

%% Start Loop
counter = 0;
for i=1:size(variables.gasesToUseArray)
    for m=1:10
       gasesToUse = nonzeros(variables.gasesToUseArray(i,:))';
       gases_to_use = 1:length(gasesToUse);
       iso=zeros(length(gasesToUse), 1);		% The isotope 'code' for the gas: '1' will be the most abundant
       clearvars  gas mixrats                   % Isotope code '0' corresponds to all isotopes
       mixrats = cell(1, length(gasesToUse));
       gas = cell(1, length(gasesToUse));
       for k=1:size(gasesToUse, 2)
            mixrats{k}.z = variables.mixRatiosZ;
            if (gasesToUse(k) == 1)
                gas{k} = 'co';
                mixrats{k}.conc = variables.gasesMR.co(:,m);
            end
            if (gasesToUse(k) == 2)
                 gas{k} = 'co2';
                 mixrats{k}.conc = variables.gasesMR.co2(:,m);
            end
            if (gasesToUse(k) == 3)
                 gas{k} = 'h2o';
                 mixrats{k}.conc = variables.gasesMR.h2o(:,m);
            end
            if (gasesToUse(k) == 4)
                 gas{k} = 'o2';
                 mixrats{k}.conc = variables.gasesMR.o2(:,m);
            end 
       end
        disp('Generate Atmospheric Cells for Calculation...')
        cells=atmcell(numberOfCells,100,27,16590,zpt,gas,iso,mixrats,'z');

        disp('Calculate the centre-line and the wing-line k-vectors...')

        [cells]=k_calc(gencalc,cells,cells_to_calc,gases_to_use);

        disp('Generate an atmospheric path for radiance calculations...')

        % Start with sun-like source path segment
        paths=path_source(Tsun,cells.wavnum);
        cell_order=[numberOfCells:1];
        % Include Atmosphere (top-down)
        paths=path_atm(cells,cell_order,gases_to_use,cells.mixrat(cell_order,gases_to_use),theta,ang_size,paths);

        paths=path_reflect(r,e,0,0,ang_size,surface,Tearth,paths);
        % Include path back through atmosphere
        cell_order=[1:numberOfCells];
        paths=path_atm(cells,cell_order,gases_to_use,cells.mixrat(cell_order,gases_to_use),theta,ang_size,paths);

        %Compute radiance
        paths=radiance(paths,1,'by path');

        % figure;
        clear paths.waveLengths;
        paths.waveLengths = zeros(1, length(paths.wavnum));
        for n = 1 : length(paths.wavnum)
            paths.waveLengths(n) = nu2l(paths.wavnum(n));
        end
        counter = counter + 1;
        pathsMaster{counter} = paths;
    end
end


