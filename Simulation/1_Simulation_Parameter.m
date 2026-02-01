%% Package loading
addpath(genpath('F:/Research/Project/FEMA_binary/cmig_tools-main/cmig_tools_utils'));
addpath(genpath('F:/Research/Project/FEMA_binary/cmig_tools-main/FEMA'));
addpath(genpath('F:/Research/Project/FEMA_binary/Code'));

%% Basic settings
nXvars           = 5;
nObservations    = 10000;
nFamMembers      = 1; % up to 5 individuals per family
% nRepObservations = 5; % up to 5 observations per individual
RandomEffects    = {'S','E'}; 
eid              = ones(nObservations, 1);
agevec           = zeros(nObservations, 1);
y_RFX            = nan(nObservations, 1);
niter            = 1;
nsamples         = 10;
nRepeats         = 100;

%% File Save path
saveDir = "F:/Research/Project/FEMA_binary/Code/Simulation_1/DATA";
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end



%% Parameter Settings
nRepObservations = [5;10;20];
var_subject = [0.1;0.4;0.7];

% Make a grid
[all_members, all_subject] = ndgrid(nRepObservations, var_subject);

% Overall grid
gridVals = [all_members(:), all_subject(:), 1 - all_subject(:)];

% range for fixed effects
betaLow     = -0.2;
betaHigh    = 0.2;

%% Data generation
rng(20260131, 'twister'); %0128

allSeeds = randi(999999, size(gridVals,1), nRepeats);

if isempty(gcp('nocreate'))
    parpool('local');
end

for nrow = 1:size(gridVals,1)

    % get parameter settings
    nRepObservations = gridVals(nrow, 1);
    nSubject         = nObservations / nRepObservations;
    var_FFX          = gridVals(nrow, 3);
    var_subject      = gridVals(nrow, 2);
    var_subject_int = round(var_subject * 10);

    seeds = allSeeds(nrow, :);

    savePath = fullfile(saveDir, sprintf('N%d_M%d_V%d', nObservations, nRepObservations, var_subject_int));
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end

    parfor r = 1:nRepeats

        rng(seeds(r), 'twister');

        % Fixed effects with intercept
        mu    = zeros(1, nXvars-1);
        sigma = eye(nXvars-1);
        X     = [ones(nObservations, 1), mvnrnd(mu, sigma, nObservations)];

        % sample beta ~ U[-0.2,0.2]
        beta   = betaLow + (betaHigh-betaLow).*rand(nXvars,1);

        % scale X to ensure that var(y_FFX) = VarFFX instead of changing beta
        y_FFX     = X * beta;
        var_yFFX_current = var(y_FFX);
        scaling_factor = sqrt(var_FFX / var_yFFX_current);
        X_scale = [ones(nObservations, 1), X(:,2:end) * scaling_factor];     
        y_FFX     = X_scale * beta;

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
    
    
        y = y_FFX + y_RFX;
        y_transformed = 1 ./ (1 + exp(-y));  % Logistic transformation

        for i=1:nsamples

            % y_binary is sampled from a Bernoulli(p), where p = y_transformed.
            % This is equivalent to drawing from a uniform distribution and comparing it to p
            y_binary = rand(nObservations, 1) < y_transformed;
            X_tmp = array2table(X_scale);
            test_dat = struct('y',      {y},                'y_binary', {y_binary},...                     
                              'X',      {X_scale},          'fid',      {fid},     ...                          
                              'iid',    {iid},              'eid',      {eid},     ...
                              'beta',   {beta},             'Var_FFX',{var_FFX},    ...
                              'Var_IID',{var_subject});
    
            % files Saved
            saveName = fullfile(savePath,sprintf('dat_%d-%d.mat',[r,i]));
            parsave(saveName, test_dat);

        end

    end
    
    disp(['Nrow: ', num2str(nrow)]);

end



% Iterate through the folder

if isempty(gcp('nocreate'))
    parpool('local');
end

SaveResultDir = "F:/Research/Project/FEMA_binary/Code/Simulation_1/FEMA_binary";


tic

N = 10000;
for M = [5, 10, 20]
    for V = [1,4,7]

        FEMAb_results = zeros(nRepeats, nXvars*2+1);

        variableNames = {'TRUE_BETA1','TRUE_BETA2','TRUE_BETA3','TRUE_BETA4',   ...
                         'TRUE_BETA5','Est_BETA1','Est_BETA2','Est_BETA3',      ...
                         'Est_BETA4','Est_BETA5','Est_Var_IID'};
        
        

        parfor  r = 1:nRepeats

            beta = zeros(nsamples, nXvars);
            beta_hat_sum = zeros(nsamples, nXvars);
            sig2mat_sum = zeros(nsamples, length(RandomEffects));

            for i = 1:nsamples

                DatDir = fullfile(saveDir, sprintf("N%d_M%d_V%d",[N,M,V]));
                fileName = fullfile(DatDir, sprintf('dat_%d-%d.mat',[r,i]));
                loadDat  = load(fileName);
            
                y_binary = loadDat.data.y_binary;
                X        = loadDat.data.X;

                agevec   = zeros(nObservations, 1);

                fid      = loadDat.data.fid;
                iid      = loadDat.data.iid;
                eid      = loadDat.data.eid;
            
                VarFFX   = loadDat.data.Var_FFX;
                VarIID   = loadDat.data.Var_IID;

                [beta_hat,      beta_se,        zmat,        logpmat,              ...
                 sig2tvec,      sig2mat,        ~,           ~,                    ...
                 beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,        ...
                 sig2mat_perm,             ~,              ~,           ~] =                               ...
                        FEMA_fit_binary(X, iid, eid, fid, agevec, y_binary, niter, ones(1,nXvars),    ...
                                      [], 'RandomEffects', {'S','E'}, 'returnReusable', true,           ...
                                'RandomEstType','MoM');
            
                beta(i,:) = loadDat.data.beta';
                beta_hat_sum(i,:) = beta_hat;
                sig2mat_sum(i,:) = sig2mat;

            end

            FEMAb_results(r, :) = [mean(beta), mean(beta_hat_sum), mean(sig2mat_sum(:,1))];

        end

        FEMAb_results = array2table(FEMAb_results, 'VariableNames',variableNames);
        output_file = fullfile(SaveResultDir, sprintf("N%d_M%d_V%d.txt", [N, M, V]));
        writetable(FEMAb_results, output_file, 'Delimiter','\t');
        
    end
end

toc
