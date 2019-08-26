function pathsMaster = recursive(variables, mixrats, k, constants, gencalc, pathsMaster)
    k = k + 1;
    for i = 1 : constants.numberOfMixRatiosPerGas
        mixrats{k}.conc = chooseGas(variables, k, i);
        if (k < constants.numberOfGases)
            pathsMaster = recursive(variables, mixrats, k, constants, gencalc, pathsMaster);
        end
        if (k == constants.numberOfGases)
            disp('Generate Atmospheric Cells for Calculation...')
            cells=atmcell(constants.numberOfCells,100,27,16590,constants.zpt,constants.gas,constants.iso,mixrats,'z');

            disp('Calculate the centre-line and the wing-line k-vectors...')

            [cells]=k_calc(gencalc,cells,constants.cells_to_calc,constants.gases_to_use);

            disp('Generate an atmospheric path for radiance calculations...')

            % Start with sun-like source path segment
            paths=path_source(constants.Tsun,cells.wavnum);
            cell_order=[constants.numberOfCells:1];
            % Include Atmosphere (top-down)
            paths=path_atm(cells,cell_order,constants.gases_to_use,cells.mixrat(cell_order,constants.gases_to_use),constants.theta,constants.ang_size,paths);

            paths=path_reflect(constants.r,constants.e,0,0,constants.ang_size,constants.surface,constants.Tearth,paths);
            % Include path back through atmosphere
            cell_order=[1:constants.numberOfCells];
            paths=path_atm(cells,cell_order,constants.gases_to_use,cells.mixrat(cell_order,constants.gases_to_use),constants.theta,constants.ang_size,paths);

            %Compute radiance
            paths=radiance(paths,1,'by path');

            % figure;
            clear paths.waveLengths;
            paths.waveLengths = zeros(1, length(paths.wavnum));
            for n = 1 : length(paths.wavnum)
                paths.waveLengths(n) = nu2l(paths.wavnum(n));
            end
            pathsMaster.counter = pathsMaster.counter + 1;
            pathsMaster.paths{1, pathsMaster.counter} = paths;
            pathsMaster.paths{2, pathsMaster.counter} = mixrats;
        end
    end
end