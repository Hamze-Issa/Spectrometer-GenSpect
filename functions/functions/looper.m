%{
The MeznSat Project
Atmospheric Retrieval and Data Processing
Hamzeh Issa

Looper: loops over all the gas mix ratio combinations applying the atmospheric model to each
one of them, and packing the results (radiance/wavelength relationship for each
mix ratio combination) into pathsMaster.
%}

function pathsMaster = looper(combinations, constants, gencalc, pathsMaster)
    for mixrat = 1 : combinations.counter
        disp('Generate Atmospheric Cells for Calculation...')
        mixrats = combinations.mixrats{mixrat};
        cells = atmcell(constants.numberOfCells,100,27,16590,constants.zpt,constants.gas,constants.iso,mixrats,'z');

        disp('Calculate the center-line and the wing-line k-vectors...')
        [cells] = k_calc(gencalc,cells,constants.cells_to_calc,constants.gases_to_use, constants.gasNames);

        disp('Generate an atmospheric path for radiance calculations...')
        % Start with sun-like source path segment
        paths = path_source(constants.Tsun,cells.wavnum);

        % Include Atmosphere downward path (top - buttom)
        cell_order = [constants.numberOfCells : 1];
        paths = path_atm(cells,cell_order,constants.gases_to_use,cells.mixrat(cell_order,constants.gases_to_use),constants.theta,constants.ang_size,paths);

        % Include path reflection on the ground
        paths = path_reflect(constants.r,constants.e,0,0,constants.ang_size,constants.surface,constants.Tearth,paths);

        % Include path back through atmosphere upward path (buttom - top)
        cell_order = [1 : constants.numberOfCells];
        paths = path_atm(cells,cell_order,constants.gases_to_use,cells.mixrat(cell_order,constants.gases_to_use),constants.theta,constants.ang_size,paths);

        % Compute radiance and store it in "paths" which holds the
        % output paths parameters for a single gas mix ratio combination
        paths = radiance(paths,1,'by path');

        clear paths.waveLengths;
        paths.waveLengths = zeros(1, length(paths.wavnum));
        for n = 1 : length(paths.wavnum)
            paths.waveLengths(n) = nu2l(paths.wavnum(n));
        end
        pathsMaster.counter = pathsMaster.counter + 1;
        disp(mixrat)
        disp(datetime(now,'ConvertFrom','datenum'))
        % Store the resulting pat radiance with its corresponding mix
        % ratio into pathsMaster
        pathsMaster.paths{1, pathsMaster.counter} = paths.radiance;
        pathsMaster.paths{2, pathsMaster.counter} = mixrat;

        % Due to the huge size of data, and to relief the RAM, the data
        % is stored on disk each 100 mix ratio combinations in a file
        % with a specific name. Each 100 mix ratios are termed as "patch"
        % After that, pathsMaster is cleared and redefined to clear the RAM
        if (mod(mixrat, 100) == 0)
            batch = mixrat / 100;
            fname = sprintf('../../scratch/hissa/results/batch%d', batch);
            save(fname, 'pathsMaster');
            clear pathsMaster;
            pathsMaster.counter = 0;
            pathsMaster.paths = cell(2,105);
        end
    end
end