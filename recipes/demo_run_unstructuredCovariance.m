%% Demo for using unstructured covariance
% This demo creates synthetic data that has time-varying random effects and
% then estimates and compares the parameters using FEMA

%% Settings
% Total number of observations
nObservations = 20000;

% Total number of unique individuals
nIndividuals = 5000;

% Total number of visits
numVisits = 5;

% Maximum number of individuals in a family
maxNumInFamily = 5;

% Minimum number of observations per visit pair
minNumObs = 500;

% Largest condition number
maxCondNum = 1000;

% Proportion of overall phenotypic variance explained by:
% Total phenotypic variance is 1
var_noise  = 0.1;
var_subj   = 0.3;
var_family = 0.2;
var_FFX    = 0.4;
var_RFX    = var_family + var_subj;

% Range for beta coefficients
betaLow  = -0.2;
betaHigh = 0.2;

% Number of X variables
nXvars = 15;

% Number of y variables
nyVars = 4;

% For age, mean and standard deviations (in days)
meanAges = [0, 60, 180, 240, 365];
stdAges  = [0, 10, 12,  14,  15];

% Max itertions for nearestSPD converge
maxIter = 1000;

% Random effects
RandomEffects = {'F', 'S', 'E'};

% Minimum and maximum values for variances and correlations
min_var_subject  =  0.2;
max_var_subject  =  0.8;
min_corr_subject = -0.7;
max_corr_subject =  0.7;

min_var_family   =  0.2;
max_var_family   =  0.8;
min_corr_family  = -0.7;
max_corr_family  =  0.7;

% Probability range for missing visits
probVisit_Low  = 0.1;
probVisit_High = 0.8;

% Set seed
rng(20260306, 'twister');

%% Simulate data
while true
    % Generate covariance matrices: subject
    while true
        % Arbitrary covariance matrix - subject effect
        allVariances    = (min_var_subject + (max_var_subject - min_var_subject) .* rand(numVisits, 1));
        allCorrelations = (min_corr_subject + (max_corr_subject - min_corr_subject) .* rand(((numVisits * (numVisits - 1))/2), 1));
    
        % Diagonal and off-diagonal locations
        locsDiagonal = 1:numVisits+1:numVisits^2;
        locsOffDiag  = find(tril(ones(numVisits), -1)); %#ok<NASGU>

        % Create covariance matrix
        count = 1;
        varCovarMatrix_subj = zeros(numVisits, numVisits);
        varCovarMatrix_subj(locsDiagonal) = allVariances;
        for i1 = 1:numVisits
            for i2 = i1+1:numVisits
                varCovarMatrix_subj(i1, i2) = allCorrelations(count) * sqrt(varCovarMatrix_subj(i1,i1) * varCovarMatrix_subj(i2,i2));
                count = count + 1;
            end
        end

        % Make symmetric and positive semidefinite
        varCovarMatrix_subj = varCovarMatrix_subj.' + triu(varCovarMatrix_subj, 1);
        [varCovarMatrix_subj, converge] = nearestSPD_timeout(varCovarMatrix_subj, maxIter);
    
        if converge
            if cond(varCovarMatrix_subj) < maxCondNum
                break;
            end
        end
    end
    
    % Generate covariance matrices: family
    while true
        % Arbitrary covariance matrix - family effect
        allVariances    = (min_var_family + (max_var_family - min_var_family) .* rand(numVisits, 1));
        allCorrelations = (min_corr_family + (max_corr_family - min_corr_family) .* rand(((numVisits * (numVisits - 1))/2), 1));

        % Diagonal and off-diagonal locations
        locsDiagonal = 1:numVisits+1:numVisits^2;
        locsOffDiag  = find(tril(ones(numVisits), -1));

        % Create covariance matrix
        count = 1;
        varCovarMatrix_family = zeros(numVisits, numVisits);
        varCovarMatrix_family(locsDiagonal) = allVariances;
        for i1 = 1:numVisits
            for i2 = i1+1:numVisits
                varCovarMatrix_family(i1, i2) = allCorrelations(count) * sqrt(varCovarMatrix_family(i1,i1) * varCovarMatrix_family(i2,i2));
                count = count + 1;
            end
        end

        % Make symmetric and positive semidefinite
        varCovarMatrix_family = varCovarMatrix_family.' + triu(varCovarMatrix_family, 1);
        [varCovarMatrix_family, converge] = nearestSPD_timeout(varCovarMatrix_family, maxIter);

        if converge
            if cond(varCovarMatrix_family) < maxCondNum
                break;
            end
        end
    end
    
    % Genereate IDs - let all IDs be generated at this point
    iid_int = sort(repmat(1:nIndividuals, 1, numVisits));
    
    % Generate family IDs
    fid_int = ceil(iid_int / maxNumInFamily);
    
    % Generate family and individual IDs
    iid = cell(size(iid_int));
    fid = cell(size(fid_int));
    for i=1:length(iid_int)
        iid{i}=sprintf('I%i', iid_int(i));
    end
    for i=1:length(fid_int)
        fid{i}=sprintf('F%i', fid_int(i));
    end
    
    % Generate event IDs
    allEvents = (1:numVisits)';
    eid       = repmat(allEvents, nIndividuals, 1);
    
    % Now, annihilate vists
    while true
        work_iid = iid;
        work_eid = eid;
        work_fid = fid;
        toDelete = nIndividuals*numVisits - nObservations;
        cumulSum = 0;
        while cumulSum ~= toDelete
            delProb  = zeros(numVisits, 1);
            for ii = 1:numVisits
                delProb(ii) = probVisit_Low + abs(probVisit_High - probVisit_Low).*rand;
            end
            cumulSum = sum(round((delProb .* nIndividuals)));
        end
        delProb  = sort(delProb);
        toDelete = round(delProb .* nIndividuals);
        
        for ii = 1:numVisits
            tmp_allLocs  = find(work_eid == ii);
            tmp_toDelete = tmp_allLocs(randperm(length(tmp_allLocs), toDelete(ii)));
            work_eid(tmp_toDelete) = [];
            work_iid(tmp_toDelete) = [];
            work_fid(tmp_toDelete) = [];
        end

        % Make sure we have enough samples for every pair of visits
        count = 1;
        allCounts = zeros((numVisits * (numVisits + 1))/2, 1);
        for i1 = 1:numVisits
            for i2 = 1:i1
                allCounts(count,1) = length(intersect(work_iid(work_eid == i1), work_iid(work_eid == i2)));
                count = count + 1;
            end
        end

        if allCounts >= minNumObs
            iid = work_iid;
            eid = work_eid;
            fid = work_fid;
            break;
        end
    end

    % Simulate age for each ID, each event
    allIDs = unique(iid, 'stable');
    agevec = zeros(length(iid), 1);
    count  = 1;
    for i  = 1:length(allIDs)
        whichEvents = eid(strcmpi(iid, allIDs{i}));
        howMany     = length(whichEvents);
        agevec(count:count+howMany-1) = floor(stdAges(whichEvents) .* randn(1, 1) + meanAges(whichEvents));
        count = count + howMany;
    end
    
    % Create clusterinfo
    clusterinfo = FEMA_parse_family(iid, eid, fid, agevec, [], 'RandomEffects', RandomEffects);
    
    % Create spline basis
    knots         = meanAges;
    basisFunction = createBasisFunctions(agevec, knots, 'nsk', [], 'svd');
    
    % Simulate fixed effects
    % How many additional X variables do we need to create beyond intercept
    % and the basis functions?
    remX   = nXvars - (size(basisFunction, 2) + 1);
    otherX = randn(length(iid),remX);

    % Put all X variables together
    X         = [ones(length(iid), 1), basisFunction, otherX];
    colnames  = [{'Intercept'}, strrep(strcat({'Basis_'}, num2str((1:size(basisFunction,2))')), ' ', '')', ...
                                strrep(strcat({'Other_'}, num2str((1:size(otherX,2))')), ' ', '')'];

    eff_FFX   = betaLow + (betaHigh - betaLow) .* rand(nXvars, nyVars);
    y_FFX     = X * eff_FFX;
    locIC     = false(size(X,2),1);
    locIC(1)  = true;

    % Mean and standard deviation of y_FFX
    std_FFX  = std(y_FFX);
    mean_FFX = mean(y_FFX);

    % Standardize FFX_y: both mean centering and standardizing
    y_FFX = (y_FFX - mean_FFX)./ std_FFX;
    y_FFX = sqrt(var_FFX)' .* y_FFX;

    % Update effect sizes
    gTruth_FFX            = eff_FFX;
    gTruth_FFX(~locIC, :) = gTruth_FFX(~locIC,:) ./ std_FFX;
    gTruth_FFX(locIC,  :) = (gTruth_FFX(locIC,:) - mean_FFX) ./ std_FFX;
    gTruth_FFX            = sqrt(var_FFX)' .* gTruth_FFX;

    % Repeated measurements - same effect across subjects
    try
        y_RFX   = zeros(length(iid), nyVars);
        for fi  = 1:length(clusterinfo)
        
            % Family effect
            [~, b]   = ismember(eid(clusterinfo{fi}.jvec_fam), allEvents);
            V_family = varCovarMatrix_family(b, b);
            tmp1     = sqrt(var_family/var_RFX) .* mvnrnd(zeros(length(clusterinfo{fi}.jvec_fam), 1), double(V_family) .* clusterinfo{fi}.V_F, nyVars)';
        
            % Subject effect
            V_subj  = varCovarMatrix_subj(b, b);
            tmp2    = sqrt(var_subj/var_RFX) .* mvnrnd(zeros(length(clusterinfo{fi}.jvec_fam), 1), double(V_subj) .* clusterinfo{fi}.V_S, nyVars)';
        
            % Put them together
            y_RFX(clusterinfo{fi}.jvec_fam,:) = sqrt(var_RFX) .* (tmp1 + tmp2);
        end
        break;
    catch
    end
end

% Noise
noise   = rand(length(iid),  nyVars);
noise   = (noise - mean(noise))./std(noise);
noise   = noise .* sqrt(var_noise)';

% Put phenotype together
ymat = y_RFX + y_FFX + noise;

% Updated ground truth information for variance components
gTruth_varCovar_family  = varCovarMatrix_family .* var_family;
gTruth_varCovar_subject = varCovarMatrix_subj .* var_subj + eye(numVisits) .* var_noise;

%% Estimation
% Settings for FEMA_fit
contrasts      = [];
nbins          = 0;
CovType        = 'unstructured';
returnReusable = false;
GRM            = [];

% Fit the model
[beta_hat,      beta_se,        zmat,        logpmat,              ...
 sig2tvec,      sig2mat,        Hessmat,     logLikvec,            ...
 beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,        ...
 sig2mat_perm,  logLikvec_perm, binvec_save, nvec_bins,            ...
 tvec_bins,     FamilyStruct,   coeffCovar,  unstructParams,       ...
 residuals_GLS, info] = FEMA_fit(X, iid, eid, fid, agevec, ymat,   ...
                                 contrasts, nbins, GRM,            ...
                                 'RandomEffects', RandomEffects,   ...
                                 'CovType', CovType, 'returnResiduals', returnReusable);

%% Compare parameters
% Fixed effects: plot everything
figure('Units', 'centimeters', 'Position', [10 10 18 8], 'Name', 'FixedEffects');
subplot(1, 2, 1);
scatter(gTruth_FFX(:), beta_hat(:), 'o', 'filled');
xlabel('True \beta (all X)');
ylabel('Estimated \beta (all X)');

% Fixed effects: colour coded by variable
subplot(1, 2, 2);
scatter(gTruth_FFX, beta_hat, 'o', 'filled');
xlabel('True \beta');
ylabel('Estimated \beta');
ll = legend(strcat({'y'}, num2str((1:nyVars)')), 'Location', 'eastoutside', 'Orientation', 'vertical', 'Box', 'off');
ll.Position(2) = ll.Position(2) - eps;

% Random effects: family effect
for yy = 1:nyVars
    figure('Units', 'centimeters', 'Position', [10 10 18 6], 'Name', ['FamilyEffect y:', num2str(yy)]);
    subplot(1, 2, 1);
    tmpLim = max(max(abs(gTruth_varCovar_family), [], 'all'), max(abs(sig2mat(:,:,1,yy) .* sig2tvec(yy)), [], 'all'));
    imagesc(gTruth_varCovar_family);
    axis xy image
    colormap('hot');
    clim([-tmpLim, tmpLim]);
    colorbar;
    title('True varCovar: family');
    for i1 = 1:numVisits
        for i2 = 1:i1
            text(i1, i2, num2str(gTruth_varCovar_family(i1, i2), '%.2f'), 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
    xticks(1:numVisits);
    yticks(1:numVisits);
    box off
    
    subplot(1, 2, 2);
    toPlot = squeeze(sig2mat(:,:,1,yy)) .* sig2tvec(yy);
    imagesc(toPlot);
    axis xy image
    colormap('hot');
    clim([-tmpLim, tmpLim]);
    colorbar;
    title('Estimated varCovar: family');
    for i1 = 1:numVisits
        for i2 = 1:i1
            text(i1, i2, num2str(toPlot(i1, i2), '%.2f'), 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
    xticks(1:numVisits);
    yticks(1:numVisits);
    box off
end

% Random effects: subject effect
for yy = 1:nyVars
    figure('Units', 'centimeters', 'Position', [10 10 18 6], 'Name', ['Subject and Error Effect y:', num2str(yy)]);
    subplot(1, 2, 1);
    tmpLim = max(max(abs(gTruth_varCovar_subject), [], 'all'), max(abs(sig2mat(:,:,2,yy) .* sig2tvec(yy)), [], 'all'));
    imagesc(gTruth_varCovar_subject);
    axis xy image
    colormap('hot');
    clim([-tmpLim, tmpLim]);
    colorbar;
    title('True varCovar: subject + error');
    for i1 = 1:numVisits
        for i2 = 1:i1
            text(i1, i2, num2str(gTruth_varCovar_subject(i1, i2), '%.2f'), 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
    xticks(1:numVisits);
    yticks(1:numVisits);
    box off
    
    subplot(1, 2, 2);
    toPlot = squeeze(sig2mat(:,:,2,yy)) .* sig2tvec(yy);
    imagesc(toPlot);
    axis xy image
    colormap('hot');
    clim([-tmpLim, tmpLim]);
    colorbar;
    title('Estimated varCovar: subject + error');
    for i1 = 1:numVisits
        for i2 = 1:i1
            text(i1, i2, num2str(toPlot(i1, i2), '%.2f'), 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
    xticks(1:numVisits);
    yticks(1:numVisits);
    box off
end