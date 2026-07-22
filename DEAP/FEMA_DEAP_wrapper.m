function FEMA_DEAP_wrapper(fstem_imaging, fname_design, dirname_out, dirname_imaging, ...
                           datatype, dirname_cache, dirname_jobs, designid, nfrac, varargin)

if isdeployed
    nfrac = str2double(nfrac);
    if isnan(nfrac)
        error('One or more input arguments could not be converted to a number.');
    end
end

if nargin < 9
    logging('Usage: FEMA_DEAP_wrapper(fstem_imaging,fname_design,dirname_out,dirname_imaging,datatype,dirname_jobs,designid,nfrac)');
    error('Incorrect number of input arguments')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Look into how include the following code as "macro"

inputs = inputParser;

% used variables
addParamValue(inputs,'RandomEffects',"{'F' 'S' 'E'}"); % Default to Family, Subject, and eps
addParamValue(inputs,'synth',0); % AMD - put back synth option
addParamValue(inputs,'colsinterest',"[]");  % if empyt we'll get all columns.
addParamValue(inputs,'contrasts',[]);
addParamValue(inputs,'output','nifti'); % but we only support this format...
addParamValue(inputs,'ivnames','');

parse(inputs,varargin{:})

contrasts = str2num_amd(inputs.Results.contrasts);
if ~isfinite(contrasts)
    fname_contrasts = inputs.Results.contrasts;
    logging('Reading contrast matrix from %s',fname_contrasts);
    contrasts = readtable(fname_contrasts);
end

outputFormat = inputs.Results.output;
if ~(contains(outputFormat,'mat') || contains(outputFormat,'nifti'))
    error('Incorrect output format (%s)',outputFormat)
end

if ~isempty(inputs.Results.ivnames)
    ivnames = split(inputs.Results.ivnames,',');
else
    ivnames = {};
end

RandomEffects = rowvec(split(strrep(strrep(inputs.Results.RandomEffects,'{',''),'}','')));
synth = str2num_amd(inputs.Results.synth);
colsinterest = str2num(inputs.Results.colsinterest);

logging('FEMA_DEAP_wrapper RandomEffects: %s, synth %d, colsinterest %s', strjoin(string(RandomEffects), ', '), synth, strjoin(string(colsinterest), ', '));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember(lower(datatype),{'voxel', 'vertex', 'corrmat', 'roi', 'external'})
    error('Input error: invalid datatype')
end

logging('***Start');

% Print build/compile information
% print_compile_stats();

starttime = now();

% Set random number generator so different every time -- should allow for option to control this
rng shuffle 

dirname_results = sprintf('%s/results/%s_%s_%s', dirname_jobs, designid, datatype, fstem_imaging);

if exist(dirname_results,'dir')
    warning([dirname_results, ' exists; overwriting']);
    % cmd = sprintf('rm -rf %s',dirname_results);
    % [s r e] = jsystem(cmd);
else
    mkdir(dirname_results);
end

% Create analysis job
[colnames_model, X, ivec_mask, mask] = FEMA_DEAP_wrapper_submit(fstem_imaging, fname_design, ...
                                                                dirname_imaging, datatype,   ...
                                                                designid, dirname_jobs);
% dispatch work to workers now that job files exists

% Load design
% fname_volinfo = sprintf('%s/vol_info.mat', dirname_imaging);
% mask          = getfield(load(fname_volinfo,'vol_mask_sub'), 'vol_mask_sub');
% ivec_mask     = find(mask>0.5); % Should save and read ivec_mask instead

% Initialization
zmat     = zeros(size(X,2), length(ivec_mask));
logpmat  = zeros(size(X,2), length(ivec_mask));
beta_hat = zeros(size(X,2), length(ivec_mask));
beta_se  = zeros(size(X,2), length(ivec_mask));
sig2mat  = zeros(length(RandomEffects), length(ivec_mask));
sig2tvec = zeros(1, length(ivec_mask));

% function to process results as each worker completes
    function process_results(fraci)
        MAX_RETRIES = 10;
        fname_results = sprintf('%s/job_%02d_%02d.mat', dirname_results, fraci, nfrac);
        retry = 0;
        tmp = struct;
        while retry < MAX_RETRIES
            try
                tmp = load(fname_results, 'zmat', 'logpmat', 'beta_hat', 'beta_se', 'sig2mat', 'sig2tvec');
                retry = MAX_RETRIES + 1;
                fname_out = sprintf('%s/%s_%s_%02d_%02d.mat', dirname_cache, datatype, ...
                                    fstem_imaging, fraci, nfrac);

                % Should perhaps check to touch file instead, to make sure that worker is done writing?
                cache = load(fname_out,'ivec_frac'); 

                zmat(:,cache.ivec_frac)     = tmp.zmat;
                logpmat(:,cache.ivec_frac)  = tmp.logpmat;
                beta_hat(:,cache.ivec_frac) = tmp.beta_hat;
                beta_se(:,cache.ivec_frac)  = tmp.beta_se;
                sig2mat(:,cache.ivec_frac)  = tmp.sig2mat;
                sig2tvec(:,cache.ivec_frac) = tmp.sig2tvec;
            catch
                retry = retry + 1;
                pause(5);
                logging('Error: Could not load %s, retry %d', fname_results, retry);
            end
        end

        if retry == MAX_RETRIES
            logging('Error: Could not load %s, giving up', fname_results);
        end
    end

% Create client to to dispatch work to workers
client = FEMA_DEAP_client(designid, fstem_imaging, @process_results, nfrac);

% block until all workers are done, will call the process_results function as each worker completes
client.wait_for_completion();

% Save output
% Calling FEMA_save
outputType = {datatype, 'mat'};
info = [];
unstructParams = [];
FEMA_save(outputType, dirname_out, 'outPrefix', designid, ...
          'beta_hat', beta_hat, 'beta_se', beta_se, 'zmat', zmat, 'logpmat', logpmat, ...
          'sig2tvec', sig2tvec, 'sig2mat', sig2mat, 'info', info, ...
          'colnames_model', colnames_model, 'mask', mask, 'unstructParams', unstructParams, ...
          'RandomEffects', RandomEffects, ...
          'colsinterest', colsinterest, ...
          'ymat_names', fstem_imaging);

% save_params = struct('fstem_imaging',fstem_imaging,'datatype',datatype,'outdir',dirname_out,'synth',synth);
% base_variables_to_save = {'X','iid','eid','colnames_model','contrasts','datatype','inputs','zmat','logpmat','beta_hat','beta_se','sig2mat','sig2tvec','save_params','mask'};
% 
% if ~exist(dirname_out,'dir')
%     mkdir(dirname_out);
% end
% 
% fpath_out = sprintf('%s/FEMA_wrapper_output_%s_%s.mat',dirname_out,datatype,fstem_imaging);
% 
% %write column names to json for DEAP
% fname_col = sprintf('%s/FEMA_results_colnames.json',dirname_out);
% out = struct('colnames_model',{colnames_model},'RandomEffects',{RandomEffects});
% jsonStr = jsonencode(out);
% fid = fopen(fname_col,'w');
% fprintf(fid,'%s\n',jsonStr);
% fclose(fid);

% =========================================================================
% Write VOXEL results (mat, nifti, or deap)
% =========================================================================
% if strcmpi(datatype,'voxel')
%     vol_z = zeros([size(mask) size(zmat,1)]);
%     vol_logp = zeros([size(mask) size(zmat,1)]);
%     vol_beta_hat = zeros([size(mask) size(zmat,1)]);
%     vol_beta_se = zeros([size(mask) size(zmat,1)]);
%     for j = 1:size(zmat,1)
%         vol_z(:,:,:,j) = single(fullvol(zmat(j,:),mask));
%         vol_logp(:,:,:,j) = single(fullvol(logpmat(j,:),mask));
%         vol_beta_hat(:,:,:,j) = single(fullvol(beta_hat(j,:),mask));
%         vol_beta_se(:,:,:,j) = single(fullvol(beta_se(j,:),mask));
%     end
% 
%     vol_sig2t = zeros([size(mask) 1]);
%     vol_sig2t(ivec_mask) = single(sig2tvec);
%     vol_sig2 = zeros([size(mask) size(sig2mat,1)]);
%     for j = 1:size(sig2mat,1)
%         vol_sig2(:,:,:,j) = single(fullvol(sig2mat(j,:),mask));
%     end
% 
%     if contains(lower(outputFormat), 'nifti')
%         results = struct('beta_hat',vol_beta_hat,'beta_se',vol_beta_se,'zmat',vol_z,'logpmat',vol_logp,'sig2tvec',vol_sig2t,'sig2mat',vol_sig2);
%         cols = colnames_model;
%         logging("length of colsinterest %d", length(colsinterest))
%         if length(colsinterest) > 0
%             cols = colnames_model(colsinterest);
%         end
%         writeNIFTI(results, dirname_out, fstem_imaging, ivnames, cols);
%     end
% else
%     fprintf(1,'Error: datatype=%s not implemented\n',datatype);
% end

% Need to add code to gather outputs
% FEMA_wrapper_script_bottom;

endtime = now();
logging('***Done*** (%0.2f seconds)',(endtime-starttime)*3600*24);
% ToDos
%   Make sure colsinterest is passed from DEAP
end
