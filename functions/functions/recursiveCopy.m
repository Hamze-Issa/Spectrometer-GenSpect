function combinations = recursiveCopy(variables, mixrats, k, constants, combinations)
    k = k + 1;
    for i = 1 : constants.numberOfMixRatiosPerGas
        mixrats{k}.conc = chooseGas(variables, k, i);
        if (k < constants.numberOfGases)
            combinations = recursiveCopy(variables, mixrats, k, constants, combinations);
        end
        if (k == constants.numberOfGases)
            combinations.counter = combinations.counter + 1;
            combinations.mixrats{1, combinations.counter} = mixrats;
        end
    end
end