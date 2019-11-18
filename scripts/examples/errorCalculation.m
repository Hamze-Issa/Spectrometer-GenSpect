numberOfMainIterations = 146;
nominatedBestPaths = cell(2, numberOfMainIterations);
load measuredPath;
for i = 1 : numberOfMainIterations
    % fname = sprintf('./results/resultsOld/iteration%d', i);
    fname = sprintf('../../scratch/hissa/results/iteration%d', i); 
    load (fname);
    disp(fname)
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
[bestPathError, bestPathIndex] = min([nominatedBestPaths{3, :}]);
bestPath.path = nominatedBestPaths{1, bestPathIndex};
bestPath.mixratIndex = nominatedBestPaths{2, bestPathIndex};
bestPath.error = nominatedBestPaths{3, bestPathIndex};

%%
load('variables');
load('constants');
numberOfGases = 4;
mixrats = cell(1,numberOfGases);
for i=1:numberOfGases
    mixrats{i}.z = variables.mixRatiosZ;
end
combinations.counter = 0;
combinations.mixrats = cell(1,1);
combinations = recursiveCopy(variables, mixrats, 0, constants, combinations);

bestPath.resultMixRatio = combinations.mixrats(bestPath.mixratIndex);
save('bestPath_Aug_14'); 