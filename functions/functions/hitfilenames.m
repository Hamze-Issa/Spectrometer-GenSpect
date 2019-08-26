function [param_filename, database] = hitfilenames(database_name,gas);

% [param_filename, database] = hitfilenames(database_name,gas);
%
% GENSPECT function to generate datafile names for Hitran format databases.
% 
% GAS is a cell string array of gasnames
% DATABASE_NAME (optional) String identifying database to use ('hitran96', 'hitran2000') or structure array
%       of alternate gas parameter input filenames (string), one per gas [Defaults to Hitran2000]. 
% PARAM_FILENAME is a structure array of database file names corresponding to the chosen option.
% DATABASE is a string describing the database choise ('Hitran96', 'Hitran2000' or 'User Defined')
%
% (C) Ben Quine, 2001

% identify gases and check validity
gasn=gasnum(gas);
if length(gasn)<1 
   error('GENSPECT: No gas names specified');
end
% set default database if none specified
if ~exist('database_name','var')
    database_name='hitran2000';
end
% check to see if userdefined database files specified
if iscell(database_name)
    % check to see if database files are specified and exist, else default
    n_datafiles=length(database_name);
    for whichgas = 1:length(gasn)
        [gasn, id_mol]=gasnum(gas{whichgas});
        if whichgas>n_datafiles
            warning(['GENSPECT: User defined parameter datafile not found for gas ' gas{whichgas} '. Defaulting to Hitran2000']);
            if ~isempty(id_mol)
                param_filename{whichgas} = ['genhit2000' '_' id_mol '.lin'];
            else
                param_filename{whichgas} = ['Non-Hitran Gas'];
            end
        elseif isempty(database_name{whichgas}) | exist(database_name{whichgas})~=2
            warning(['GENSPECT: User defined parameter datafile not found for gas ' gas{whichgas} '. Defaulting to Hitran2000']);
            if ~isempty(id_mol)
                param_filename{whichgas} = ['genhit2000' '_' id_mol '.lin'];
            else
                param_filename{whichgas} = ['Non-Hitran Gas'];
            end
        else 
            param_filename{whichgas} = database_name{whichgas};
        end
    end 
    database='User Defined';
elseif strcmp(lower(database_name),'hitran2000')
    for whichgas = 1:length(gasn)
        [gasn, id_mol]=gasnum(gas{whichgas});
        if ~isempty(id_mol)
            param_filename{whichgas} = ['genhit2000' '_' id_mol '.lin'];
        else
            param_filename{whichgas} = ['Non-Hitran Gas'];
        end
    end 
    database=database_name;
elseif strcmp(lower(database_name),'hitran96')
    for whichgas = 1:length(gasn)
        [gasn, id_mol]=gasnum(gas{whichgas});
        if ~isempty(id_mol)
            param_filename{whichgas} = ['genhit96' '_' id_mol '.lin'];
        else
            param_filename{whichgas} = ['Non-Hitran Gas'];
        end
    end 
    database=database_name;
else
    error('GENSPECT: Invalid parameter database filename specified')
end
