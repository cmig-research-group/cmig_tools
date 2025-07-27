%% Package loading
addpath(genpath('./cmig_tools-main/cmig_tools_utils'));
addpath(genpath('./cmig_tools-main/FEMA'));

%% Basic settings
nXvars           = 5;
nObservations    = 10000;
nSubject         = 3000;
nFamMembers      = 5; % up to 5 individuals per family
nRepObservations = 5; % up to 5 observations per individual
RandomEffects    = {'F','S','E'}; 
eid              = ones(nObservations, 1);
agevec           = zeros(nObservations, 1);
y_RFX            = nan(nObservations, 1);
niter            = 1;
nRepeats         = 10;

%% File Save path
saveDir = "./DATA";
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end



%% Parameter Settings
var_FFX    = 0.1:0.1:0.8;
var_family = 0.1:0.1:0.8;
var_subj   = 0.1:0.1:0.8;

% Make a grid
[all_FFX, all_family, all_subj] = ndgrid(var_FFX, var_family, var_subj);

% Overall grid
gridVals = [all_FFX(:), all_family(:), all_subj(:)];

% Only keep values where total variance is 1
gridVals = gridVals(round(sum(gridVals(:,1:3),2),10) == 1, :);

betaLow     = 0.2;
betaHigh    = 0.4;





%% Data generation - Using the same X structure in different simulation settings?

rng(20250622, 'twister');

% Fixed effects with intercept
mu    = zeros(1, nXvars-1);
sigma = eye(nXvars-1);
X     = [ones(nObservations, 1), mvnrnd(mu, sigma, nObservations)];

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


% Generate binary Y with across different parameters
for nrow = 1:size(gridVals,1)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    allSeeds = randi(999999, nRepeats, 1);

    % get parameter settings
    VarFFX = gridVals(nrow, 1);
    sig2mat_true = gridVals(nrow, 2:3);

    for r = 1:nRepeats

        rng(allSeeds(r), 'twister');
    
        % sample beta ~ U[0.2,0.4]
        beta   = betaLow + (betaHigh-betaLow).*rand(nXvars,1);

        % scale X to ensure that var(y_FFX) = VarFFX instead of changing beta
        y_FFX     = X * beta;
        var_yFFX_current = var(y_FFX);
        scaling_factor = sqrt(VarFFX / var_yFFX_current);
        X_scale = [ones(nObservations, 1), X(:,2:end) * scaling_factor];     
        y_FFX     = X_scale * beta;
    
        % Generate random effects
        nfam = length(clusterinfo);
        for fi = 1:nfam
            % if mod(fi,100)==0
            %     logging('fi=%d/%d',fi,nfam);
            % end
            tmp = 0;
            for ri = 1:2
                tmp = tmp + sqrt(sig2mat_true(ri)) .*                        ...
                            mvnrnd(zeros(length(clusterinfo{fi}.jvec_fam), 1), ...
                            double(getfield(clusterinfo{fi},                   ...
                            sprintf('V_%s',RandomEffects{ri}))),1)'; %#ok<GFLD>
            end
            y_RFX(clusterinfo{fi}.jvec_fam,:) = tmp;
        end
    
    
        y = y_FFX + y_RFX;
        y_transformed = 1 ./ (1 + exp(-y));  % Logistic transformation

        % y_binary is sampled from a Bernoulli(p), where p = y_transformed.
        % This is equivalent to drawing from a uniform distribution and comparing it to p  
        y_binary = rand(nObservations, 1) < y_transformed;
        X_tmp = array2table(X_scale);
        test_dat = struct('y',      {y},                'y_binary', {y_binary},...                     
                          'X',      {X_scale},          'fid',      {fid},     ...                          
                          'iid',    {iid},              'eid',      {eid},     ...
                          'beta',   {beta},             'Var_FFX',{VarFFX},    ...
                          'Var_FID',  {sig2mat_true(1)},'Var_IID',{sig2mat_true(2)}, 'r',{r});

        % Saved as Comp_VarFFX_VarFam_VarSubj_Repetition.txt
        save(fullfile(saveDir, sprintf('Comp_%g-%g-%g_%d.mat', [10*gridVals(nrow, 1:3), r])), 'test_dat');

    end

    disp(['Nrow: ', num2str(nrow)]);

end


%% Iterate through the folder 
matFiles = dir(fullfile(saveDir, '*.mat'));

FEMAb_results = zeros(0, 17); % record the estimation
variableNames = {'Var_FFX','Var_FID','Var_IID','r','Est_BETA1','Est_BETA2', ...
                 'Est_BETA3','Est_BETA4','Est_BETA5','Est_SIG1',        ...
                 'Est_SIG2','EST_SIG_TOT','TRUE_BETA1','TRUE_BETA2',    ...
                 'TRUE_BETA3','TRUE_BETA4','TRUE_BETA5'};

for i = 1:length(matFiles)

    fileName = matFiles(i).name;
    loadDat  = load(fullfile(saveDir, fileName));

    y_binary = loadDat.test_dat.y_binary;
    X        = loadDat.test_dat.X;
    fid      = loadDat.test_dat.fid;
    iid      = loadDat.test_dat.iid;
    eid      = loadDat.test_dat.eid;

    VarFFX   = loadDat.test_dat.Var_FFX;
    VarFID   = loadDat.test_dat.Var_FID;
    VarIID   = loadDat.test_dat.Var_IID;
    r        = loadDat.test_dat.r;

    beta = loadDat.test_dat.beta;

    [beta_hat,      beta_se,        zmat,        logpmat,              ...
     sig2tvec,      sig2mat,        ~,           ~,                    ...
     ~,             ~,              ~,           ~,                    ...
     ~,             ~,              binvec_save, ~,                    ...
     ~,             ~,              ~] =                               ...
            FEMA_fit_binary(X, iid, eid, fid, agevec, y_binary, niter, ones(1,nXvars),    ...
                    20, [], 'RandomEffects', {'F','S','E'}, 'returnReusable', true,           ...
                    'RandomEstType','MoM');

    input_row = [VarFFX, VarFID, VarIID, r, beta_hat', sig2mat(1:2,:)', sig2tvec, beta'];
    
    FEMAb_results       = [FEMAb_results; input_row];

end

%% Save the results

FEMAb_results = array2table(FEMAb_results, 'VariableNames',variableNames);
writetable(FEMAb_results, "./FEMA_binary.txt", 'Delimiter','\t');


%% Where u can check the convergence of results
% Enabling the recording in FEMA_fit_binary simutaneously
% -- Line 343-344, 389-390, 454-455 in FEMA_fit_binary.m
% for index = 1:5
%     figure;
%     plot(1:iter,beta_record(1:iter,index),'LineWidth',1);
%     xlabel('Iterations');
%     ylabel(sprintf('Beta%d',index));
% end
% 
% 
% for index = 1:2
%     figure;
%     plot(1:iter,sigmat_record(1:iter,index),'LineWidth',1);
%     % yline(0.01, 'LineWidth', 1.5);
%     xlabel('Iterations');
%     ylabel(sprintf('Normalized Sigma%d',index));
% end
