function [fstem_imaging, fname_design, dirname_out, dirname_imaging, datatype, extraArgs] = ...
            FEMA_parseInputs(fname_json)

    %% TO DOs %%
        % - contrasts

    if ~nargin == 1 || ~isfile(fname_json) % are we getting the json file path or the contents? 
        error('Invalid JSON file path.');
    else 
        % --- Case 1: Single JSON file input ---
        fprintf('Loading parameters from JSON file: %s\n', fname_json);
        
        % Load and decode JSON
        configFile = jsondecode(fileread(fname_json));
        
        % Extract required parameters
        fstem_imaging = configFile.params.dependent.name;
        dirname_imaging = configFile.params.dependent.dir_data;
        datatype = strrep(configFile.params.dependent.type_data, 'wise', ''); 
        fname_design = '~/';  % need to write out a config file for makeDesign??
        dirname_out = pwd();  % Use current directory as default
 
        extraArgs = {};
        % transformY 
        if isfield(configFile.params.dependent, 'transform') 
            extraArgs{end+1} = 'transformY';
            extraArgs{end+1} = configFile.params.dependent.transform;
        end
        % outputType 
        extraArgs{end+1} = 'outputType';
        extraArgs{end+1} = setOutputType(datatype);
        % ivnames 
        ivnames = get_invnames(configFile);
        extraArgs{end+1} = 'ivnames';
        extraArgs{end+1} = ivnames;
        % random effects 
        extraArgs{end+1} = 'RandomEffects';
        extraArgs{end+1} = configFile.params.random;
        % GRM          ----
        % preg_file        | - do we always want to set this from hard coded path? 
        % address file ----
        % nperms
        extraArgs{end+1} = 'nperms';
        extraArgs{end+1} = configFile.params.advanced.n_perm;
        % permutation type
        extraArgs{end+1} = 'permtype';
        extraArgs{end+1} = configFile.params.advanced.type_perm;
        % covariance 
        extraArgs{end+1} = 'CovType';
        extraArgs{end+1} = configFile.params.advanced.type_cov;
        % ffx estimation type 
        extraArgs{end+1} = 'FixedEstType';
        extraArgs{end+1} = upper(configFile.params.advanced.type_fixed_est);
        % precision
        extraArgs{end+1} = 'SingleOrDouble';
        extraArgs{end+1} = configFile.params.advanced.precision;
        % binning
        extraArgs{end+1} = 'nbins';
        extraArgs{end+1} = configFile.params.advanced.n_bins;
    end
end

function ivnames = get_invnames(configFile)
    nFFX = length(configFile.params.fixed.vars);
    ivnames = cell(1, nFFX);
    for i = 1:nFFX
        if configFile.params.fixed.vars{i}.of_interest
            ivnames{i} = configFile.params.fixed.vars{i}.name;
        end         
    end
    ivnames = ivnames(~cellfun('isempty', ivnames));
end 

function outputType = setOutputType(datatype)
    switch datatype
        case 'voxel'
            outputType = 'nifti';
        case 'vertex'
            outputType = 'gifti';
        case 'corrmat'
            outputType = 'corrmat';
        case 'roi'
            outputType = 'nifti';
        case 'external'
            outputType = 'tables';
    end
end