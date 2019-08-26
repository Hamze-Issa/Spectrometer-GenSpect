function pathsMaster = test(variables, mixrats, k, constants, gencalc, pathsMaster)
    k = k + 1;
    for i = 1 : 2
        mixrats{k}.conc = chooseGas(variables, k, i);
        if (k < 2)
            pathsMaster = test(variables, mixrats, k, constants, gencalc, pathsMaster);
        end
        if (k == 2)
            
            pathsMaster.counter = pathsMaster.counter + 1;
            pathsMaster.paths{1, pathsMaster.counter} = 1;
            pathsMaster.paths{2, pathsMaster.counter} = mixrats;
        end
    end
end