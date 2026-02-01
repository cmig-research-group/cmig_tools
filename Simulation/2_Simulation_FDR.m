%% Package loading
addpath(genpath('F:/Research/Project/FEMA_binary/cmig_tools-main/cmig_tools_utils'));
addpath(genpath('F:/Research/Project/FEMA_binary/cmig_tools-main/FEMA'));
addpath(genpath('F:/Research/Project/FEMA_binary/Code'));

%% Basic settings
nXvars           = 100;
nObservations    = 10000;
nFamMembers      = 1; % 1 individuals per family
% nRepObservations = 5; % 1 observations per individual
RandomEffects    = {'S','E'}; 
eid              = ones(nObservations, 1);
agevec           = zeros(nObservations, 1);
y_RFX            = nan(nObservations, 1);
niter            = 1;
nRepeats         = 1000;

%% File Save path
saveDir = "F:/Research/Project/FEMA_binary/Code/Simulation_2/DATA";
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end


%% Parameter Settings
nRepObservations = [5;10;20];
var_subject = [0.1;0.4;0.7];

% Make a grid
[all_members, all_subject] = ndgrid(nRepObservations, var_subject);

% Overall grid
gridVals = [all_members(:), all_subject(:)];

%% Data generation - Using the different X in different simulation settings

rng(20260128, 'twister');

allSeeds = randi(999999, size(gridVals,1), nRepeats);

if isempty(gcp('nocreate'))
    parpool('local');
end

% Generate binary Y with across different parameters
for nrow = 1:size(gridVals,1)

    % get parameter settings
    nRepObservations = gridVals(nrow, 1);
    nSubject         = nObservations / nRepObservations;
    var_subject = gridVals(nrow, 2);
    var_subject_int = round(var_subject * 10);

    seeds = allSeeds(nrow, :);

    savePath = fullfile(saveDir, sprintf('N%d_M%d_V%d', nObservations, nRepObservations, var_subject_int));
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end

    parfor r = 1:nRepeats

        rng(seeds(r), 'twister');

        % Fixed effects
        mu    = zeros(1, nXvars);
        sigma = eye(nXvars);
        X     = [mvnrnd(mu, sigma, nObservations)];
       
        beta = zeros(nXvars, 1);

        % Random effects
        iid_int = sort(randsample(repmat(1:nSubject, 1, nRepObservations), nObservations));
        fid_int = ceil(iid_int / nFamMembers); 
        
        % Initialize
        iid = cell(size(iid_int)); 
        fid = cell(size(fid_int)); 
        
        % Generate individual and family IDs
        for i=1:length(iid_int)
            iid{i}=sprintf('I%i', iid_int(i)); 
        end
        for i=1:length(fid_int)
            fid{i}=sprintf('F%i', fid_int(i)); 
        end
        
        % Parse family structure
        clusterinfo = FEMA_parse_family(iid, eid, fid, agevec, [], ...
                                       'RandomEffects', RandomEffects);
    
        % Generate random effects
        y_RFX = nan(nObservations, 1);
        nfam = length(clusterinfo);
        for fi = 1:nfam
            % if mod(fi,100)==0
            %     logging('fi=%d/%d',fi,nfam);
            % end
            tmp = sqrt(var_subject) .*                        ...
                mvnrnd(zeros(length(clusterinfo{fi}.jvec_fam), 1), ...
                double(getfield(clusterinfo{fi},'V_S')),1)';
            y_RFX(clusterinfo{fi}.jvec_fam,:) = tmp;
        end
    
        y = y_RFX;
        y_transformed = 1 ./ (1 + exp(-y));  % Logistic transformation

        % y_binary is sampled from a Bernoulli(p), where p = y_transformed.
        % This is equivalent to drawing from a uniform distribution and comparing it to p  
        y_binary = rand(nObservations, 1) < y_transformed;
        test_dat = struct('y',      {y},                'y_binary', {y_binary},...                     
                          'X',      {X},                'fid',      {fid},     ...                          
                          'iid',    {iid},              'eid',      {eid},     ...
                          'Var_IID',  {var_subject});

        % files Saved
        saveName = fullfile(savePath,sprintf('dat_%d.mat',r));
        parsave(saveName, test_dat);
    end

    disp(['Nrow: ', num2str(nrow)]);

end


%% Iterate through the folder

if isempty(gcp('nocreate'))
    parpool('local');
end

SaveResultDir = "F:/Research/Project/FEMA_binary/Code/Simulation_2";

N = 10000;
for M = [5, 10, 20]
    for V = [1,4,7]

        FEMAb_results = zeros(nRepeats, nXvars+1);
        
        tic

        parfor  i = 1:nRepeats
 
            DatDir = fullfile(saveDir, sprintf("N%d_M%d_V%d",[N,M,V]));
            fileName = fullfile(DatDir, sprintf('dat_%d.mat',i));
            loadDat  = load(fileName);
        
            y_binary = loadDat.data.y_binary;
            y        = loadDat.data.y;
            X        = loadDat.data.X;
            
            agevec   = zeros(nObservations, 1);
        
            fid      = loadDat.data.fid;
            iid      = loadDat.data.iid;
            eid      = loadDat.data.eid;
        
            y_transformed = 1 ./ (1 + exp(-y)); 
            W = y_transformed .* (1-y_transformed); 
        
            [beta_hat,      beta_se,        zmat,        logpmat,              ...
             sig2tvec,      sig2mat,        ~,           ~,                    ...
             beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,        ...
             sig2mat_perm,             ~,              ~,           ~] =                               ...
                    FEMA_fit_binary(X, iid, eid, fid, agevec, y_binary, niter, ones(1,nXvars),    ...
                                  [], 'RandomEffects', {'S','E'}, 'returnReusable', true,           ...
                            'RandomEstType','MoM');

            % plot_qq_gwas(logpmat);
 
            FEMAb_results(i, :) = [i, logpmat'];
        
            % display(i);
        
        end
        
        toc

        output_file = fullfile(SaveResultDir, sprintf("N%d_M%d_V%d.txt", [N, M, V]));
        writematrix(FEMAb_results, output_file, 'Delimiter','\t');

    end
end

%% Q-Q plot
FEMAb_results = readmatrix(fullfile(SaveResultDir,"N10000_M20_V7.txt"));
hasInf = any(isinf(FEMAb_results(1:nRepeats, 2:end)), 'all');
if hasInf
    hasInfRows = any(isinf(FEMAb_results), 2);
    infRowIndices = find(hasInfRows);
    fprintf('包含Inf的行号：\n');
    disp(infRowIndices');
else
    logp_vector = reshape(FEMAb_results(1:nRepeats, 2:end), [], 1);
    plot_qq_gwas(logp_vector(:,1));

end
