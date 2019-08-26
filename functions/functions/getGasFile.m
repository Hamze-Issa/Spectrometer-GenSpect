function gasFile = getGasFile(param_filename, gasNames)
    if (strcmp(param_filename, 'genhit2000_02.lin'))
        gasFile = gasNames.genhit2000_02.lins;
    end
    if (strcmp(param_filename, 'genhit2000_01.lin'))
         gasFile = gasNames.genhit2000_01.lins;
    end
    if (strcmp(param_filename, 'genhit2000_07.lin'))
         gasFile = gasNames.genhit2000_07.lins;
    end
    if (strcmp(param_filename, 'genhit2000_05.lin'))
         gasFile = gasNames.genhit2000_05.lins;
    end
end
