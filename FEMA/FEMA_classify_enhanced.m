function [results, models] = FEMA_classify_enhanced(fstem_imaging, fname_design, dirname_out, dirname_imaging, datatype, target, varargin)
%
% Enhanced classification wrapper for FEMA pipeline with multiple algorithms
% and improved cross-validation
%
% USAGE: [results, models] = FEMA_classify_enhanced(fstem_imaging, fname_design, dirname_out, dirname_imaging, datatype, target, varargin)
%
% INPUTS
%   fstem_imaging <char>       : name of vertex/voxel-mapped phenotype (e.g., 'thickness-sm16', 'FA')
%   fname_design <cell>        : cell array with path to design matrix files
%   dirname_out <cell>         : cell array with output directory paths
%   dirname_imaging <char>     : path to imaging data directory
%   datatype <char>            : 'voxel','vertex','external', 'corrmat'
%   target <cell>              : target variable names for classification
%
% Optional input arguments:
%   classifier <char>          : 'svm' (default), 'logistic', 'randomforest', 'ensemble'
%   cv_method <char>           : 'kfold' (default), 'stratified', 'family_holdout'
%   n_folds <num>              : number of CV folds (default 10)
%   n_repetitions <num>        : number of CV repetitions (default 10)
%   feature_selection <char>   : 'none' (default), 'svd', 'lasso', 'mutual_info'
%   n_components <num>         : number of components for feature selection (default 250)
%   hyperparameter_tuning <bool>: whether to tune hyperparameters (default true)
%   age_binning <bool>         : whether to analyze by age bins (default true)
%   age_bin_size <num>         : age bin size in years (default 2.5)
%   output_format <char>       : 'mat' (default), 'json', 'csv'
%   save_models <bool>         : whether to save trained models (default false)
%   parallel <bool>            : whether to use parallel processing (default false)
%   verbose <bool>             : whether to display progress (default true)
%   ranknorm <bool>            : whether to apply rank normalization (default false)
%   varnorm <bool>             : whether to apply variance normalization (default false)
%   ico <num>                  : icosahedron order for vertex-based data (default 5)
%   pihat_file <char>          : path to pihat file (default [])
%   preg_file <char>           : path to preg file (default [])
%   address_file <char>        : path to address file (default [])
%
% OUTPUTS
%   results <struct>           : structure containing classification results
%   models <struct>            : structure containing trained models (if save_models=true)
%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logging('***Start FEMA_classify_enhanced v2.0***');
starttime = now();
rng shuffle % Set random number generator so different every time

% Input validation
if nargin < 6
    error('Usage: FEMA_classify_enhanced(fstem_imaging,fname_design,dirname_out,dirname_imaging,datatype,target,varargin)');
end

% Parse optional arguments
p = inputParser;
addParameter(p, 'classifier', 'svm', @ischar);
addParameter(p, 'cv_method', 'kfold', @ischar);
addParameter(p, 'n_folds', 10, @isnumeric);
addParameter(p, 'n_repetitions', 10, @isnumeric);
addParameter(p, 'feature_selection', 'svd', @ischar);
addParameter(p, 'n_components', 250, @isnumeric);
addParameter(p, 'hyperparameter_tuning', true, @islogical);
addParameter(p, 'age_binning', true, @islogical);
addParameter(p, 'age_bin_size', 2.5, @isnumeric);
addParameter(p, 'output_format', 'mat', @ischar);
addParameter(p, 'save_models', false, @islogical);
addParameter(p, 'parallel', false, @islogical);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'ranknorm', false, @islogical);
addParameter(p, 'varnorm', false, @islogical);
addParameter(p, 'ico', 5, @isnumeric);
addParameter(p, 'pihat_file', [], @ischar);
addParameter(p, 'preg_file', [], @ischar);
addParameter(p, 'address_file', [], @ischar);

parse(p, varargin{:});
params = p.Results;

% Validate inputs
if ~ismember(lower(datatype), {'voxel', 'vertex', 'external', 'corrmat'})
    error('Invalid datatype: must be voxel, vertex, external, or corrmat');
end

if ~ismember(lower(params.classifier), {'svm', 'logistic', 'randomforest', 'ensemble'})
    error('Invalid classifier: must be svm, logistic, randomforest, or ensemble');
end

% Ensure cell arrays
if ~iscell(fname_design), fname_design = {fname_design}; end
if ~iscell(dirname_out), dirname_out = {dirname_out}; end
if ~iscell(target), target = {target}; end

if length(fname_design) ~= length(dirname_out)
    error('fname_design and dirname_out must have equal number of items');
end

if length(target) ~= length(fname_design)
    target = repmat(target, 1, length(fname_design));
end

% Initialize results structure
results = struct();
models = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND PROCESS IMAGING DATA

if params.verbose
    fprintf('Loading and processing imaging data...\n');
end

[ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat, preg, address] = ...
    FEMA_process_data(fstem_imaging, dirname_imaging, datatype, ...
    'ranknorm', params.ranknorm, 'varnorm', params.varnorm, 'ico', params.ico, ...
    'pihat_file', params.pihat_file, 'preg_file', params.preg_file, 'address_file', params.address_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ANALYSIS LOOP

for des = 1:length(fname_design)
    if params.verbose
        fprintf('Processing design matrix %d/%d: %s\n', des, length(fname_design), fname_design{des});
    end
    
    % Intersect with design matrix
    [X, iid, eid, fid, agevec, ymat_intersect, ~, colnames_model, pihatmat, PregID, HomeID] = ...
        FEMA_intersect_design(fname_design{des}, ymat, iid_concat, eid_concat, ...
        'pihat', pihat, 'preg', preg, 'address', address);
    
    % Find target variable
    colnum = find(strcmp(colnames_model, target{des}));
    if isempty(colnum)
        warning('Target %s not found in design matrix %s', target{des}, fname_design{des});
        continue;
    end
    
    % Prepare target variable
    Y = X(:, colnum) == 1;
    
    % Remove target column from features
    X_features = X(:, setdiff(1:size(X, 2), colnum));
    
    % Residualize imaging data by covariates
    if params.verbose
        fprintf('Residualizing imaging data...\n');
    end
    ymat_resid = ymat_intersect - X_features * pinv(X_features) * ymat_intersect;
    
    % Feature selection
    if params.verbose
        fprintf('Performing feature selection (%s)...\n', params.feature_selection);
    end
    
    [features, feature_info] = perform_feature_selection(ymat_resid, Y, params.feature_selection, params.n_components);
    
    % Age binning setup
    if params.age_binning
        age_bins = create_age_bins(agevec, params.age_bin_size);
        n_age_bins = length(age_bins);
    else
        age_bins = {1:length(agevec)};
        n_age_bins = 1;
    end
    
    % Initialize results for this design
    results(des).design_file = fname_design{des};
    results(des).target = target{des};
    results(des).classifier = params.classifier;
    results(des).cv_method = params.cv_method;
    results(des).n_folds = params.n_folds;
    results(des).n_repetitions = params.n_repetitions;
    results(des).feature_selection = params.feature_selection;
    results(des).n_components = params.n_components;
    results(des).age_binning = params.age_binning;
    results(des).age_bin_size = params.age_bin_size;
    
    % Cross-validation
    if params.verbose
        fprintf('Performing cross-validation (%d folds, %d repetitions)...\n', params.n_folds, params.n_repetitions);
    end
    
    [cv_results, cv_models] = perform_cross_validation(features, Y, fid, agevec, age_bins, params);
    
    % Store results
    results(des).cv_results = cv_results;
    if params.save_models
        models(des).cv_models = cv_models;
    end
    
    % Generate summary statistics
    results(des).summary = generate_summary_statistics(cv_results, age_bins);
    
    % Create visualizations
    if params.verbose
        fprintf('Creating visualizations...\n');
    end
    
    create_classification_plots(results(des), dirname_out{des}, fstem_imaging, datatype, target{des});
    
    % Save results
    save_results(results(des), dirname_out{des}, fstem_imaging, datatype, target{des}, params.output_format);
    
    if params.verbose
        fprintf('Completed design matrix %d/%d\n', des, length(fname_design));
    end
end

% Print final summary
if params.verbose
    print_final_summary(results);
end

logging('***FEMA_classify_enhanced completed***');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS

function [features, feature_info] = perform_feature_selection(data, labels, method, n_components)
    % Perform feature selection on the data
    
    switch lower(method)
        case 'svd'
            % SVD-based feature selection
            [U, S, V] = svd(data, 'econ');
            features = U(:, 1:min(n_components, size(U, 2)));
            feature_info.singular_values = diag(S);
            feature_info.explained_variance = cumsum(diag(S).^2) / sum(diag(S).^2);
            
        case 'lasso'
            % LASSO-based feature selection
            [B, FitInfo] = lasso(data, labels, 'CV', 5);
            features = data * B(:, FitInfo.IndexMinMSE);
            feature_info.lasso_coefficients = B(:, FitInfo.IndexMinMSE);
            feature_info.lasso_fitinfo = FitInfo;
            
        case 'mutual_info'
            % Mutual information-based feature selection
            mi_scores = zeros(size(data, 2), 1);
            for i = 1:size(data, 2)
                mi_scores(i) = mutual_information(data(:, i), labels);
            end
            [~, sorted_idx] = sort(mi_scores, 'descend');
            features = data(:, sorted_idx(1:min(n_components, length(sorted_idx))));
            feature_info.mi_scores = mi_scores;
            feature_info.selected_features = sorted_idx(1:min(n_components, length(sorted_idx)));
            
        case 'none'
            % No feature selection
            features = data;
            feature_info = struct();
            
        otherwise
            error('Unknown feature selection method: %s', method);
    end
end

function age_bins = create_age_bins(agevec, bin_size)
    % Create age bins for analysis
    
    age_min = floor(min(agevec));
    age_max = ceil(max(agevec));
    bin_edges = age_min:bin_size:age_max;
    
    age_bins = cell(length(bin_edges)-1, 1);
    for i = 1:length(bin_edges)-1
        age_bins{i} = find(agevec >= bin_edges(i) & agevec < bin_edges(i+1));
    end
end

function [cv_results, cv_models] = perform_cross_validation(features, labels, family_ids, agevec, age_bins, params)
    % Perform cross-validation with multiple algorithms
    
    cv_results = struct();
    cv_models = struct();
    
    % Create cross-validation indices
    cv_indices = create_cv_indices(family_ids, labels, params);
    
    % Initialize results arrays
    n_bins = length(age_bins);
    n_reps = params.n_repetitions;
    n_folds = params.n_folds;
    
    cv_results.auc = NaN(n_bins, n_reps);
    cv_results.accuracy = NaN(n_bins, n_reps);
    cv_results.sensitivity = NaN(n_bins, n_reps);
    cv_results.specificity = NaN(n_bins, n_reps);
    cv_results.f1_score = NaN(n_bins, n_reps);
    cv_results.predictions = cell(n_bins, n_reps);
    cv_results.probabilities = cell(n_bins, n_reps);
    
    % Perform cross-validation for each repetition
    for rep = 1:n_reps
        if params.verbose
            fprintf('  Repetition %d/%d\n', rep, n_reps);
        end
        
        % Get CV indices for this repetition
        cv_idx = cv_indices{rep};
        
        % Perform CV for each age bin
        for bin = 1:n_bins
            bin_indices = age_bins{bin};
            if isempty(bin_indices)
                continue;
            end
            
            [bin_results, bin_models] = perform_bin_cv(features(bin_indices, :), ...
                labels(bin_indices), cv_idx(bin_indices), params);
            
            % Store results
            cv_results.auc(bin, rep) = bin_results.auc;
            cv_results.accuracy(bin, rep) = bin_results.accuracy;
            cv_results.sensitivity(bin, rep) = bin_results.sensitivity;
            cv_results.specificity(bin, rep) = bin_results.specificity;
            cv_results.f1_score(bin, rep) = bin_results.f1_score;
            cv_results.predictions{bin, rep} = bin_results.predictions;
            cv_results.probabilities{bin, rep} = bin_results.probabilities;
            
            if params.save_models
                cv_models.bin_models{bin, rep} = bin_models;
            end
        end
    end
end

function cv_indices = create_cv_indices(family_ids, labels, params)
    % Create cross-validation indices based on method
    
    unique_families = unique(family_ids, 'stable');
    n_families = length(unique_families);
    n_folds = params.n_folds;
    n_reps = params.n_repetitions;
    
    cv_indices = cell(n_reps, 1);
    
    for rep = 1:n_reps
        switch lower(params.cv_method)
            case 'kfold'
                % Standard k-fold CV
                cv_indices{rep} = crossvalind('Kfold', length(family_ids), n_folds);
                
            case 'stratified'
                % Stratified k-fold CV
                cv_indices{rep} = crossvalind('StratifiedKfold', labels, n_folds);
                
            case 'family_holdout'
                % Family-based holdout (families in same fold)
                family_cv = crossvalind('Kfold', n_families, n_folds);
                cv_indices{rep} = zeros(length(family_ids), 1);
                
                for i = 1:n_families
                    family_mask = strcmp(family_ids, unique_families{i});
                    cv_indices{rep}(family_mask) = family_cv(i);
                end
                
            otherwise
                error('Unknown CV method: %s', params.cv_method);
        end
    end
end

function [results, models] = perform_bin_cv(features, labels, cv_idx, params)
    % Perform cross-validation for a single age bin
    
    n_folds = params.n_folds;
    predictions = zeros(length(labels), 1);
    probabilities = zeros(length(labels), 1);
    models = cell(n_folds, 1);
    
    for fold = 1:n_folds
        % Split data
        train_idx = cv_idx ~= fold;
        test_idx = cv_idx == fold;
        
        X_train = features(train_idx, :);
        y_train = labels(train_idx);
        X_test = features(test_idx, :);
        
        % Train model
        model = train_classifier(X_train, y_train, params);
        models{fold} = model;
        
        % Make predictions
        [pred, prob] = predict_classifier(model, X_test, params);
        predictions(test_idx) = pred;
        probabilities(test_idx) = prob;
    end
    
    % Calculate performance metrics
    results = calculate_performance_metrics(labels, predictions, probabilities);
    results.predictions = predictions;
    results.probabilities = probabilities;
end

function model = train_classifier(X, y, params)
    % Train classification model
    
    switch lower(params.classifier)
        case 'svm'
            if params.hyperparameter_tuning
                % Hyperparameter tuning for SVM
                model = fitcsvm(X, y, 'Standardize', true, 'KernelFunction', 'RBF', ...
                    'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'}, ...
                    'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                    'MaxObjectiveEvaluations', 30, 'ShowPlots', false));
            else
                model = fitcsvm(X, y, 'Standardize', true, 'KernelFunction', 'RBF', 'KernelScale', 'auto');
            end
            
        case 'logistic'
            if params.hyperparameter_tuning
                % Hyperparameter tuning for logistic regression
                model = fitclinear(X, y, 'Learner', 'logistic', ...
                    'OptimizeHyperparameters', {'Lambda'}, ...
                    'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
                    'MaxObjectiveEvaluations', 30, 'ShowPlots', false));
            else
                model = fitclinear(X, y, 'Learner', 'logistic');
            end
            
        case 'randomforest'
            if params.hyperparameter_tuning
                % Hyperparameter tuning for random forest
                model = TreeBagger(50, X, y, 'Method', 'classification', ...
                    'NumPredictorsToSample', 'all', 'MinLeafSize', 1);
            else
                model = TreeBagger(100, X, y, 'Method', 'classification');
            end
            
        case 'ensemble'
            % Ensemble of multiple classifiers
            svm_model = fitcsvm(X, y, 'Standardize', true, 'KernelFunction', 'RBF');
            log_model = fitclinear(X, y, 'Learner', 'logistic');
            rf_model = TreeBagger(50, X, y, 'Method', 'classification');
            
            model = struct('svm', svm_model, 'logistic', log_model, 'randomforest', rf_model);
            
        otherwise
            error('Unknown classifier: %s', params.classifier);
    end
end

function [predictions, probabilities] = predict_classifier(model, X, params)
    % Make predictions using trained model
    
    switch lower(params.classifier)
        case 'svm'
            [predictions, scores] = predict(model, X);
            probabilities = scores(:, 2); % Probability of positive class
            
        case 'logistic'
            [predictions, scores] = predict(model, X);
            probabilities = scores(:, 2);
            
        case 'randomforest'
            [predictions, scores] = predict(model, X);
            probabilities = scores(:, 2);
            
        case 'ensemble'
            % Ensemble prediction (majority vote)
            svm_pred = predict(model.svm, X);
            log_pred = predict(model.logistic, X);
            rf_pred = predict(model.randomforest, X);
            
            predictions = mode([svm_pred, log_pred, rf_pred], 2);
            
            % Average probabilities
            [~, svm_scores] = predict(model.svm, X);
            [~, log_scores] = predict(model.logistic, X);
            [~, rf_scores] = predict(model.randomforest, X);
            
            probabilities = mean([svm_scores(:, 2), log_scores(:, 2), rf_scores(:, 2)], 2);
            
        otherwise
            error('Unknown classifier: %s', params.classifier);
    end
end

function results = calculate_performance_metrics(labels, predictions, probabilities)
    % Calculate performance metrics
    
    % Confusion matrix
    TP = sum(predictions == 1 & labels == 1);
    TN = sum(predictions == 0 & labels == 0);
    FP = sum(predictions == 1 & labels == 0);
    FN = sum(predictions == 0 & labels == 1);
    
    % Basic metrics
    results.accuracy = (TP + TN) / (TP + TN + FP + FN);
    results.sensitivity = TP / (TP + FN); % Recall
    results.specificity = TN / (TN + FP);
    results.precision = TP / (TP + FP);
    
    % F1 score
    if results.precision + results.sensitivity > 0
        results.f1_score = 2 * (results.precision * results.sensitivity) / (results.precision + results.sensitivity);
    else
        results.f1_score = 0;
    end
    
    % AUC
    if length(unique(labels)) == 2
        [~, ~, ~, results.auc] = perfcurve(labels, probabilities, 1);
    else
        results.auc = NaN;
    end
    
    % Additional metrics
    results.balanced_accuracy = (results.sensitivity + results.specificity) / 2;
    results.matthews_correlation = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
end

function summary = generate_summary_statistics(cv_results, age_bins)
    % Generate summary statistics across cross-validation repetitions
    
    summary = struct();
    
    % Calculate mean and std across repetitions
    metrics = {'auc', 'accuracy', 'sensitivity', 'specificity', 'f1_score'};
    
    for i = 1:length(metrics)
        metric = metrics{i};
        if isfield(cv_results, metric)
            data = cv_results.(metric);
            summary.([metric '_mean']) = nanmean(data, 2);
            summary.([metric '_std']) = nanstd(data, 0, 2);
            summary.([metric '_median']) = nanmedian(data, 2);
            summary.([metric '_ci_95']) = 1.96 * nanstd(data, 0, 2) / sqrt(size(data, 2));
        end
    end
    
    % Overall performance (across all age bins)
    summary.overall_auc_mean = nanmean(cv_results.auc(:));
    summary.overall_auc_std = nanstd(cv_results.auc(:));
    summary.overall_accuracy_mean = nanmean(cv_results.accuracy(:));
    summary.overall_accuracy_std = nanstd(cv_results.accuracy(:));
end

function create_classification_plots(results, output_dir, fstem_imaging, datatype, target)
    % Create visualization plots
    
    % Ensure output directory exists
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % 1. Overall performance plot
    figure('Visible', 'off');
    plot_performance_across_age(results);
    saveas(gcf, fullfile(output_dir, sprintf('%s_%s_predicting_%s_performance.png', datatype, fstem_imaging, target)));
    close(gcf);
    
    % 2. ROC curves
    figure('Visible', 'off');
    plot_roc_curves(results);
    saveas(gcf, fullfile(output_dir, sprintf('%s_%s_predicting_%s_roc.png', datatype, fstem_imaging, target)));
    close(gcf);
    
    % 3. Confusion matrix heatmap
    figure('Visible', 'off');
    plot_confusion_matrix(results);
    saveas(gcf, fullfile(output_dir, sprintf('%s_%s_predicting_%s_confusion.png', datatype, fstem_imaging, target)));
    close(gcf);
end

function plot_performance_across_age(results)
    % Plot performance metrics across age bins
    
    if ~isfield(results, 'summary')
        return;
    end
    
    metrics = {'auc', 'accuracy', 'sensitivity', 'specificity'};
    colors = {'b', 'r', 'g', 'm'};
    
    for i = 1:length(metrics)
        metric = metrics{i};
        mean_field = [metric '_mean'];
        std_field = [metric '_std'];
        
        if isfield(results.summary, mean_field)
            mean_vals = results.summary.(mean_field);
            std_vals = results.summary.(std_field);
            
            x = 1:length(mean_vals);
            errorbar(x, mean_vals, std_vals, colors{i}, 'LineWidth', 2, 'DisplayName', upper(metric));
            hold on;
        end
    end
    
    xlabel('Age Bin');
    ylabel('Performance');
    title('Classification Performance Across Age Bins');
    legend('show');
    grid on;
end

function plot_roc_curves(results)
    % Plot ROC curves
    
    if ~isfield(results, 'cv_results') || ~isfield(results.cv_results, 'probabilities')
        return;
    end
    
    hold on;
    
    % Plot ROC for each age bin
    n_bins = size(results.cv_results.probabilities, 1);
    colors = lines(n_bins);
    
    for bin = 1:n_bins
        all_probs = [];
        all_labels = [];
        
        % Collect all probabilities and labels for this bin
        for rep = 1:size(results.cv_results.probabilities, 2)
            probs = results.cv_results.probabilities{bin, rep};
            if ~isempty(probs)
                all_probs = [all_probs; probs];
                all_labels = [all_labels; results.cv_results.predictions{bin, rep}];
            end
        end
        
        if ~isempty(all_probs)
            [X, Y, ~, AUC] = perfcurve(all_labels, all_probs, 1);
            plot(X, Y, 'Color', colors(bin, :), 'LineWidth', 2, ...
                'DisplayName', sprintf('Age Bin %d (AUC=%.3f)', bin, AUC));
        end
    end
    
    plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', 'Random');
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title('ROC Curves by Age Bin');
    legend('show');
    grid on;
end

function plot_confusion_matrix(results)
    % Plot confusion matrix heatmap
    
    if ~isfield(results, 'cv_results') || ~isfield(results.cv_results, 'predictions')
        return;
    end
    
    % Aggregate predictions and labels across all repetitions
    all_preds = [];
    all_labels = [];
    
    for bin = 1:size(results.cv_results.predictions, 1)
        for rep = 1:size(results.cv_results.predictions, 2)
            preds = results.cv_results.predictions{bin, rep};
            if ~isempty(preds)
                all_preds = [all_preds; preds];
                all_labels = [all_labels; ones(length(preds), 1)]; % Assuming binary classification
            end
        end
    end
    
    if ~isempty(all_preds)
        cm = confusionmat(all_labels, all_preds);
        confusionchart(cm, {'Negative', 'Positive'});
        title('Confusion Matrix (Aggregated)');
    end
end

function save_results(results, output_dir, fstem_imaging, datatype, target, format)
    % Save results in specified format
    
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    base_filename = sprintf('%s_%s_predicting_%s_results', datatype, fstem_imaging, target);
    
    switch lower(format)
        case 'mat'
            save(fullfile(output_dir, [base_filename '.mat']), 'results', '-v7.3');
            
        case 'json'
            json_str = jsonencode(results);
            fid = fopen(fullfile(output_dir, [base_filename '.json']), 'w');
            fprintf(fid, '%s', json_str);
            fclose(fid);
            
        case 'csv'
            % Save summary statistics as CSV
            if isfield(results, 'summary')
                summary_table = struct2table(results.summary);
                writetable(summary_table, fullfile(output_dir, [base_filename '_summary.csv']));
            end
            
        otherwise
            warning('Unknown output format: %s. Saving as .mat', format);
            save(fullfile(output_dir, [base_filename '.mat']), 'results', '-v7.3');
    end
end

function print_final_summary(results)
    % Print final summary to console
    
    fprintf('\n=== FEMA Classification Results Summary ===\n');
    
    for i = 1:length(results)
        fprintf('\nDesign %d: %s\n', i, results(i).target);
        fprintf('  Classifier: %s\n', results(i).classifier);
        fprintf('  Overall AUC: %.3f ± %.3f\n', ...
            results(i).summary.overall_auc_mean, results(i).summary.overall_auc_std);
        fprintf('  Overall Accuracy: %.3f ± %.3f\n', ...
            results(i).summary.overall_accuracy_mean, results(i).summary.overall_accuracy_std);
    end
    
    fprintf('\n=== End Summary ===\n');
end

function mi = mutual_information(x, y)
    % Calculate mutual information between two variables
    
    % Discretize continuous variables
    x_bins = 10;
    y_bins = 2; % Binary target
    
    x_edges = linspace(min(x), max(x), x_bins + 1);
    y_edges = [0, 0.5, 1];
    
    % Create joint histogram
    [counts, ~, ~] = histcounts2(x, y, x_edges, y_edges);
    
    % Calculate marginal distributions
    p_x = sum(counts, 2) / sum(counts(:));
    p_y = sum(counts, 1) / sum(counts(:));
    p_xy = counts / sum(counts(:));
    
    % Calculate mutual information
    mi = 0;
    for i = 1:size(p_xy, 1)
        for j = 1:size(p_xy, 2)
            if p_xy(i, j) > 0 && p_x(i) > 0 && p_y(j) > 0
                mi = mi + p_xy(i, j) * log2(p_xy(i, j) / (p_x(i) * p_y(j)));
            end
        end
    end
end 