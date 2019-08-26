function conc = chooseGas(variables, gas, m)
    if (gas == 1)
        conc = variables.gasesMR.co(:,m);
    end
    if (gas == 2)
         conc = variables.gasesMR.co2(:,m);
    end
    if (gas == 3)
         conc = variables.gasesMR.h2o(:,m);
    end
    if (gas == 4)
         conc = variables.gasesMR.o2(:,m);
    end 
end