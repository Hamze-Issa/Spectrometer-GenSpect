function pathsMaster = parallelLooper(combinations, constants, gencalc, pathsMaster)
tic
tempArray = cell(2, 100);
    parfor mixrat = 1 : combinations.counter
        disp('Generate Atmospheric Cells for Calculation...')
        mixrats = combinations.mixrats{mixrat};
        cells=atmcell(constants.numberOfCells,100,27,16590,constants.zpt,constants.gas,constants.iso,mixrats,'z');

        disp('Calculate the centre-line and the wing-line k-vectors...')

        [cells]=k_calc(gencalc,cells,constants.cells_to_calc,constants.gases_to_use, constants.gasNames);

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
        tempArray{1, mixrat} = paths;
%         tempArray{2, mixrat} = mixrats;
%         pathsMaster.paths{1, mixrat} = paths;
%         pathsMaster.paths{2, mixrat} = mixrats;
    end
    toc
    pathsMaster.paths
end