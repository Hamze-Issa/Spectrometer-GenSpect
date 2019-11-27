%{
The MeznSat Project
Atmospheric Retrieval and Data Processing
Hamzeh Issa

chooseGas: returns the gas concentration for the input gas ID
Although unnecessary, this method is used for compatability
%}

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