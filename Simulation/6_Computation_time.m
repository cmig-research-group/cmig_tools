 %% Package loading
addpath(genpath('F:/Research/Project/FEMA_binary/cmig_tools-main/cmig_tools_utils'));
addpath(genpath('F:/Research/Project/FEMA_binary/cmig_tools-main/FEMA'));
addpath(genpath('F:/Research/Project/FEMA_binary/Code'));

%% Basic settings
nFamMembers      = 1; % up to 5 individuals per family
nRepObservations = 10; % up to 5 observations per individual
var_subject      = 0.4;
var_FFX          = 1-var_subject;
RandomEffects    = {'S','E'}; 
niter            = 1;
nRepeats         = 50;

%% File Save path
saveDir = "F:/Research/Project/FEMA_binary/Code/Simulation_6/DATA";
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Parameter Settings
nXvars           = [5;10;20;50;100];
nObservations    = [1000;2000;5000;10000;20000;50000;100000];

% Make a grid
[all_nXvars, all_nObs] = ndgrid(nXvars, nObservations);

% Overall grid
gridVals = [all_nXvars(:), all_nObs(:)];

% range for fixed effects
betaLow     = -0.2;
betaHigh    = 0.2;

%% Data generation
rng(20260412, 'twister');

allSeeds = randi(999999, size(gridVals,1), nRepeats);

if isempty(gcp('nocreate'))
    parpool('local');
end

for nrow = 1:size(gridVals,1)

    % get parameter settings
    nXvars           = gridVals(nrow, 1);
    nObservations    = gridVals(nrow, 2);
    nSubject         = nObservations / nRepObservations;
    eid              = ones(nObservations, 1);
    agevec           = zeros(nObservations, 1);
    y_RFX            = nan(nObservations, 1);

    seeds = allSeeds(nrow, :);

    savePath = fullfile(saveDir, sprintf('N%d_X%d', nObservations, nXvars));
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


        % y_binary is sampled from a Bernoulli(p), where p = y_transformed.
        % This is equivalent to drawing from a uniform distribution and comparing it to p
        y_binary = rand(nObservations, 1) < y_transformed;
        X_tmp = array2table(X_scale,'VariableNames', compose('X%d', 1:nXvars)');
        data_table = table(double(y_binary), 'VariableNames', {'y_binary'});
        data_table = [data_table, X_tmp];
        data_table.fid = fid(:);
        data_table.iid = iid(:);
        data_table.eid = eid(:);
        

        % files Saved
        saveName = fullfile(savePath,sprintf('dat_%d.txt',r));
        writetable(data_table, saveName, 'Delimiter','\t');

    end
    
    disp(['Nrow: ', num2str(nrow)]);

end



%% Iterate through the folder

% if isempty(gcp('nocreate'))
%     parpool('local');
% end

SaveResultDir = "F:/Research/Project/FEMA_binary/Code/Simulation_6/FEMA_binary";

results_table = zeros(size(gridVals,1),3);
row = 1;

for Nx = [5,10,20,50,100]
    for N = [1000,2000,5000,10000,20000,50000,100000]

        mem_usage_GB = zeros(nRepeats, 1);

        fprintf("N%d_Nx%d\n",[N,Nx]);

        tic
       
        for  r = 1:1

            DatDir = fullfile(saveDir, sprintf("N%d_X%d",[N,Nx]));
            fileName = fullfile(DatDir, sprintf('dat_%d.txt',r));
            loadDat  = readtable(fileName, 'Delimiter', '\t');
        
            y_binary = loadDat.y_binary;
            X_names = compose('X%d', 1:Nx)';
            X = table2array(loadDat(:, X_names));

            agevec   = zeros(N, 1);

            fid = loadDat.fid;
            iid = loadDat.iid;
            eid = loadDat.eid;
            

            % java.lang.System.gc();
            % 
            % mem_before = memory;
            % mem_before_MB = mem_before.MemUsedMATLAB / (1024^3);


            [beta_hat,      beta_se,        zmat,        logpmat,              ...
             sig2tvec,      sig2mat,        ~,           ~,                    ...
             beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,        ...
             sig2mat_perm,             ~,              ~,           ~] =                               ...
                    FEMA_fit_binary(X(:,2:end), iid, eid, fid, agevec, y_binary, niter, ones(1,Nx),    ...
                                  [], 'RandomEffects', {'S','E'}, 'returnReusable', true,           ...
                            'RandomEstType','MoM');

            % mem_after = memory;
            % mem_after_MB = mem_after.MemUsedMATLAB / (1024^3);
            % mem_usage_GB(r) = mem_after_MB - mem_before_MB;
            % 
            % java.lang.System.gc();

        end

        total_time = toc;

        results_table(row, :) = [Nx, N, total_time];
        row = row + 1;
        
    end
end
%%
saveName = fullfile(SaveResultDir, 'FEMA_results.txt');
results_table = array2table(results_table, 'VariableNames', {'Nx', 'N', 'TotalTime_seconds'});
writetable(results_table, saveName, 'Delimiter', '\t')
