%% Package loading
addpath(genpath('F:/Research/Project/FEMA_binary/cmig_tools-main/cmig_tools_utils'));
addpath(genpath('F:/Research/Project/FEMA_binary/cmig_tools-main/FEMA'));
addpath(genpath('F:/Research/Project/FEMA_binary/Code'));

%% Basic settings
nCausal          = 1;
nNonCausal       = 100;
nXvars           = nNonCausal+1;
nObservations    = 10000;
nFamMembers      = 1; % up to 5 individuals per family
nRepObservations = 10; % up to 5 observations per individual
nSubject         = 1000;
RandomEffects    = {'S','E'}; 
eid              = ones(nObservations, 1);
agevec           = zeros(nObservations, 1);
y_RFX            = nan(nObservations, 1);
niter            = 1;
nRepeats         = 1000;


%% Parameter Settings
Prevalance = [0.5; 0.4; 0.3; 0.2; 0.1];

var_subject = 0.4;


%% File Save path
saveDir = "F:/Research/Project/FEMA_binary/Code/Simulation_4/DATA";
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end


%% Data generation - Using the same X structure in different simulation settings?

rng(20260221, 'twister'); % 0208

allSeeds = randi(999999, size(Prevalance,1)+1, nRepeats);

if isempty(gcp('nocreate'))
    parpool('local');
end

for nrow = 1:size(Prevalance,1)

    % get parameter settings
    target_rate = Prevalance(nrow);
    target_rate_int = round(target_rate * 100);

    seeds = allSeeds(nrow, :);

    savePath = fullfile(saveDir, sprintf('P%d', target_rate_int));
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end

    parfor r = 1:nRepeats

        rng(seeds(r), 'twister');

        % Fixed effects with intercept
        mu    = zeros(1, nXvars-1);
        sigma = eye(nXvars-1);
        X     = [ones(nObservations,1),mvnrnd(mu, sigma, nObservations)];

        beta   = zeros(nXvars, 1);

        y_FFX     = X * beta;

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
    

        prop_high_risk = target_rate - 0.05;
    
        X_hidden_risk = rand(nObservations, 1) < prop_high_risk;
        
        beta_hidden = 5;
        calc_prob = @(b0) 1 ./ (1 + exp(-(b0 + beta_hidden * X_hidden_risk + y_RFX)));
        obj_fun = @(b0) mean(calc_prob(b0)) - target_rate;
        intercept_adj = fzero(obj_fun, [-10, 10]);

        % obtain y with optimized intercept
        beta_adj = [intercept_adj;  beta_hidden; beta(2:nXvars,:)];
        X = [ones(nObservations,1),X_hidden_risk, mvnrnd(mu, sigma, nObservations)];
        y_FFX     = X * beta_adj;
        y         = y_FFX + y_RFX;
        y_transformed = 1 ./ (1 + exp(-y));  % Logistic transformation

        y_binary = rand(nObservations, 1) < y_transformed;
        test_dat = struct('y',      {y},                'y_binary', {y_binary},    ...                     
                          'X',      {X},                'fid',      {fid},         ...                          
                          'iid',    {iid},              'eid',      {eid},         ...
                          'intercept',  {intercept_adj},'Var_IID',{var_subject},   ...
                          'Prevalence',{target_rate});


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

SaveResultDir = "F:/Research/Project/FEMA_binary/Code/Simulation_4";


tic

for P = [10,20,30,40,50]

    FEMAb_results = zeros(nRepeats, nXvars+2);

    parfor  i = 1:nRepeats

        DatDir = fullfile(saveDir, sprintf("P%d",P));
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
                FEMA_fit_binary(X(:,2:end), iid, eid, fid, agevec, y_binary, niter, ones(1,nXvars),    ...
                              [], 'RandomEffects', {'S','E'}, 'returnReusable', true,           ...
                              'RandomEstType','MoM');

        % plot_qq_gwas(logpmat);

        FEMAb_results(i, :) = [i, logpmat'];

        disp(['i = ', num2str(i)]);

    end

    % toc

    output_file = fullfile(SaveResultDir, sprintf("P%d.txt", P));
    writematrix(FEMAb_results, output_file, 'Delimiter','\t');

end

toc

%% Q-Q plot
FEMAb_results = readmatrix(fullfile(SaveResultDir,"P50.txt"));
hasInf = any(isinf(FEMAb_results(1:nRepeats, 4:end)), 'all');
if hasInf
    hasInfRows = any(isinf(FEMAb_results), 2);
    infRowIndices = find(hasInfRows);
    fprintf('包含Inf的行号：\n');
    disp(infRowIndices');
else
    logp_vector = reshape(FEMAb_results(1:nRepeats, 4:end), [], 1);
    plot_qq_gwas(logp_vector(:,1));

end
