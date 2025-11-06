function [fstem_imaging, config_design, dirname_out, dirname_imaging, datatype, dataFile, extraArgs] = FEMA_parseInputs(varargin)

    %% TO DOs %%
    % - contrasts
    % - add gifti to roi output depending on table?
    % inputs must chart or string
    
    idx_config = find(strcmp(varargin, {'config'}));
    idx_data = find(strcmp(varargin, {'data'}));
    idx_out = find(strcmp(varargin, {'output'}));
    
    fname_json = varargin{idx_config+1}; 
    dataFile = varargin{idx_data+1};
    dirname_out = varargin{idx_out+1};
    
    fprintf('Loading parameters from JSON file: %s\n', fname_json);
    
    % Load and decode JSON
    configFile = jsondecode(fileread(fname_json));
    
    % Extract required parameters
    fstem_imaging = configFile.params.dependent.name;
    dirname_imaging = configFile.params.dependent.dir_data;
    datatype = strrep(configFile.params.dependent.type_data, 'wise', ''); 
    config_design = fname_json;  % need to write out a config file for makeDesign??
    
    extraArgs = {};
    % outPrefix
    extraArgs{end+1} = 'outPrefix';
    extraArgs{end+1} = configFile.id;
    % transformY 
    if isfield(configFile.params.dependent, 'transform') 
        extraArgs{end+1} = 'transformY';
        extraArgs{end+1} = configFile.params.dependent.transform;
    end
    % outputType 
    extraArgs{end+1} = 'outputType';
    extraArgs{end+1} = setOutputType(datatype);
    % ivnames - now handled by FEMA_parse_JSON and FEMA_makeDesign
    %ivnames = get_ivnames(configFile);
    %extraArgs{end+1} = 'ivnames';
    %extraArgs{end+1} = ivnames;
    % random effects 
    extraArgs{end+1} = 'RandomEffects';
    extraArgs{end+1} = configFile.params.random;
    % GRM (preg_file and address_file tbd)
    extraArgs{end+1} = 'GRM_file';
    extraArgs{end+1} = configFile.params.dir_grm;
    
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
    extraArgs{end+1} = 'precision';
    extraArgs{end+1} = configFile.params.advanced.precision;
    % binning
    extraArgs{end+1} = 'nbins';
    extraArgs{end+1} = configFile.params.advanced.n_bins;
end

function ivnames = get_ivnames(configFile)
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