%{
The MeznSat Project
Atmospheric Retrieval and Data Processing
Hamzeh Issa

generateCombinations: generates all mix ratio combinations for all gases in
a recursive manner
%}

function combinations = generateCombinations(variables, mixrats, k, constants, combinations)
    k = k + 1;
    for i = 1 : constants.numberOfMixRatiosPerGas
        mixrats{k}.conc = chooseGas(variables, k, i);
        if (k < constants.numberOfGases)
            combinations = generateCombinations(variables, mixrats, k, constants, combinations);
        end
        if (k == constants.numberOfGases)
            combinations.counter = combinations.counter + 1;
            combinations.mixrats{1, combinations.counter} = mixrats;
        end
    end
end