%{
The MeznSat Project
Atmospheric Retrieval and Data Processing
Hamzeh Issa

errorCalculation: takes all the batch files that contain the
combinations with thei corresponding radiance/wavelength relationships, and
computes the MS error between that radiance and the measured radiance taken by
the spectrometer. Than, it chooses the path (with its combination) with the least
error as the best path
%}

numberOfMainBatches = 146; % 11 ^ 4 / 100
% nominatedBestPaths is to hold the best path of each batch
nominatedBestPaths = cell(2, numberOfMainBatches);
load measuredPath;
for i = 1 : numberOfMainBatches
    % fname = sprintf('./results/resultsOld/batch%d', i);
    fname = sprintf('../../scratch/hissa/results/batch%d', i); 
    load (fname);
    disp(fname)
    % Each batch has 100 combinations (paths, results, whatever)
    for path = 1 : 100
        pathsMaster.paths{3, path} = (measuredPath - pathsMaster.paths{1, path}).^2;
        pathsMaster.paths{4, path} = sqrt(mean(pathsMaster.paths{3, path}));
    end
    [minimumError, minimumErrorIndex] = min(abs([pathsMaster.paths{4, :}]))
    nominatedBestPaths{1, i} = pathsMaster.paths{1, minimumErrorIndex};
    nominatedBestPaths{2, i} = pathsMaster.paths{2, minimumErrorIndex};
    nominatedBestPaths{3, i} = minimumError;
    clear pathsMaster;
end
% Get the best path from nominatedBestPaths (best path ever!)
[bestPathError, bestPathIndex] = min([nominatedBestPaths{3, :}]);
bestPath.path = nominatedBestPaths{1, bestPathIndex};
bestPath.mixratIndex = nominatedBestPaths{2, bestPathIndex};
bestPath.error = nominatedBestPaths{3, bestPathIndex};

%%
% This is only to get the mix ratio of the best path (because mixed ratios weren't
% saved in pathsMaster to save memory)
load('variables');
load('constants');
numberOfGases = 4;
mixrats = cell(1,numberOfGases);
for i=1:numberOfGases
    mixrats{i}.z = variables.mixRatiosZ;
end
combinations.counter = 0;
combinations.mixrats = cell(1,1);
combinations = generateCombinations(variables, mixrats, 0, constants, combinations);

bestPath.resultMixRatio = combinations.mixrats(bestPath.mixratIndex);
save('bestPath_November_24'); 