%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEMA CLASSIFY ENHANCED DEMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script demonstrates the enhanced FEMA classification capabilities
% with multiple algorithms, feature selection methods, and cross-validation approaches

% This demo will work for ABCD INVESTIGATORS with access to abcd-sync.
% You need to have abcd-sync mirrored onto your server in the exact same directory structure.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP AND CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add cmig_tools to path if not already added
if ~exist('FEMA_classify_enhanced', 'file')
    fprintf('Adding cmig_tools to MATLAB path...\n');
    addpath(genpath('.')); % Adjust path as needed
end

% Configuration
fstem_imaging = 'thickness-sm16';  % Example: cortical thickness
dirname_imaging = '/path/to/abcd-sync';  % Update with your path
datatype = 'vertex';
target = {'diagnosis'};  % Example target variable

% Create output directories
output_base = './FEMA_classify_enhanced_demo_results';
if ~exist(output_base, 'dir')
    mkdir(output_base);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO 1: BASIC CLASSIFICATION WITH DIFFERENT ALGORITHMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== DEMO 1: Basic Classification with Different Algorithms ===\n');

% Test different classifiers
classifiers = {'svm', 'logistic', 'randomforest', 'ensemble'};
results_basic = cell(length(classifiers), 1);

for i = 1:length(classifiers)
    fprintf('\nTesting classifier: %s\n', classifiers{i});
    
    output_dir = fullfile(output_base, sprintf('basic_%s', classifiers{i}));
    
    try
        [results_basic{i}, ~] = FEMA_classify_enhanced(fstem_imaging, ...
            {'design_matrix.txt'}, {output_dir}, dirname_imaging, ...
            datatype, target, ...
            'classifier', classifiers{i}, ...
            'cv_method', 'stratified', ...
            'n_folds', 5, ...
            'n_repetitions', 5, ...
            'feature_selection', 'svd', ...
            'n_components', 100, ...
            'hyperparameter_tuning', false, ... % Disable for faster demo
            'age_binning', false, ... % Disable for basic demo
            'verbose', true);
        
        fprintf('Completed %s classifier\n', classifiers{i});
        
    catch ME
        fprintf('Error with %s classifier: %s\n', classifiers{i}, ME.message);
        results_basic{i} = [];
    end
end

% Compare results
fprintf('\n--- Basic Classification Results Comparison ---\n');
for i = 1:length(classifiers)
    if ~isempty(results_basic{i})
        fprintf('%s: AUC = %.3f ± %.3f, Accuracy = %.3f ± %.3f\n', ...
            classifiers{i}, ...
            results_basic{i}.summary.overall_auc_mean, ...
            results_basic{i}.summary.overall_auc_std, ...
            results_basic{i}.summary.overall_accuracy_mean, ...
            results_basic{i}.summary.overall_accuracy_std);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO 2: FEATURE SELECTION METHODS COMPARISON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== DEMO 2: Feature Selection Methods Comparison ===\n');

feature_methods = {'svd', 'lasso', 'mutual_info', 'none'};
results_features = cell(length(feature_methods), 1);

for i = 1:length(feature_methods)
    fprintf('\nTesting feature selection: %s\n', feature_methods{i});
    
    output_dir = fullfile(output_base, sprintf('features_%s', feature_methods{i}));
    
    try
        [results_features{i}, ~] = FEMA_classify_enhanced(fstem_imaging, ...
            {'design_matrix.txt'}, {output_dir}, dirname_imaging, ...
            datatype, target, ...
            'classifier', 'svm', ...
            'cv_method', 'kfold', ...
            'n_folds', 5, ...
            'n_repetitions', 5, ...
            'feature_selection', feature_methods{i}, ...
            'n_components', 100, ...
            'hyperparameter_tuning', false, ...
            'age_binning', false, ...
            'verbose', true);
        
        fprintf('Completed %s feature selection\n', feature_methods{i});
        
    catch ME
        fprintf('Error with %s feature selection: %s\n', feature_methods{i}, ME.message);
        results_features{i} = [];
    end
end

% Compare feature selection results
fprintf('\n--- Feature Selection Results Comparison ---\n');
for i = 1:length(feature_methods)
    if ~isempty(results_features{i})
        fprintf('%s: AUC = %.3f ± %.3f\n', ...
            feature_methods{i}, ...
            results_features{i}.summary.overall_auc_mean, ...
            results_features{i}.summary.overall_auc_std);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO 3: CROSS-VALIDATION METHODS COMPARISON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== DEMO 3: Cross-Validation Methods Comparison ===\n');

cv_methods = {'kfold', 'stratified', 'family_holdout'};
results_cv = cell(length(cv_methods), 1);

for i = 1:length(cv_methods)
    fprintf('\nTesting CV method: %s\n', cv_methods{i});
    
    output_dir = fullfile(output_base, sprintf('cv_%s', cv_methods{i}));
    
    try
        [results_cv{i}, ~] = FEMA_classify_enhanced(fstem_imaging, ...
            {'design_matrix.txt'}, {output_dir}, dirname_imaging, ...
            datatype, target, ...
            'classifier', 'svm', ...
            'cv_method', cv_methods{i}, ...
            'n_folds', 5, ...
            'n_repetitions', 5, ...
            'feature_selection', 'svd', ...
            'n_components', 100, ...
            'hyperparameter_tuning', false, ...
            'age_binning', false, ...
            'verbose', true);
        
        fprintf('Completed %s CV method\n', cv_methods{i});
        
    catch ME
        fprintf('Error with %s CV method: %s\n', cv_methods{i}, ME.message);
        results_cv{i} = [];
    end
end

% Compare CV results
fprintf('\n--- Cross-Validation Results Comparison ---\n');
for i = 1:length(cv_methods)
    if ~isempty(results_cv{i})
        fprintf('%s: AUC = %.3f ± %.3f\n', ...
            cv_methods{i}, ...
            results_cv{i}.summary.overall_auc_mean, ...
            results_cv{i}.summary.overall_auc_std);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO 4: AGE-BINNED ANALYSIS WITH HYPERPARAMETER TUNING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== DEMO 4: Age-Binned Analysis with Hyperparameter Tuning ===\n');

output_dir = fullfile(output_base, 'age_binned_advanced');

try
    [results_age, models_age] = FEMA_classify_enhanced(fstem_imaging, ...
        {'design_matrix.txt'}, {output_dir}, dirname_imaging, ...
        datatype, target, ...
        'classifier', 'ensemble', ...
        'cv_method', 'family_holdout', ...
        'n_folds', 5, ...
        'n_repetitions', 10, ...
        'feature_selection', 'lasso', ...
        'n_components', 50, ...
        'hyperparameter_tuning', true, ...
        'age_binning', true, ...
        'age_bin_size', 2.0, ...
        'output_format', 'json', ...
        'save_models', true, ...
        'verbose', true);
    
    fprintf('Completed age-binned analysis with hyperparameter tuning\n');
    
    % Display age-specific results
    if isfield(results_age, 'summary') && isfield(results_age.summary, 'auc_mean')
        fprintf('\n--- Age-Specific Results ---\n');
        for bin = 1:length(results_age.summary.auc_mean)
            fprintf('Age Bin %d: AUC = %.3f ± %.3f\n', ...
                bin, ...
                results_age.summary.auc_mean(bin), ...
                results_age.summary.auc_std(bin));
        end
    end
    
catch ME
    fprintf('Error with age-binned analysis: %s\n', ME.message);
    results_age = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO 5: COMPREHENSIVE COMPARISON AND SUMMARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== DEMO 5: Comprehensive Comparison and Summary ===\n');

% Create summary table
fprintf('\n--- Comprehensive Results Summary ---\n');
fprintf('%-20s %-15s %-15s %-15s\n', 'Method', 'AUC', 'Accuracy', 'F1 Score');
fprintf('%-20s %-15s %-15s %-15s\n', '--------------------', '---------------', '---------------', '---------------');

% Basic classifiers
for i = 1:length(classifiers)
    if ~isempty(results_basic{i})
        fprintf('%-20s %-15.3f %-15.3f %-15.3f\n', ...
            classifiers{i}, ...
            results_basic{i}.summary.overall_auc_mean, ...
            results_basic{i}.summary.overall_accuracy_mean, ...
            nanmean(results_basic{i}.summary.f1_score_mean));
    end
end

% Feature selection methods
for i = 1:length(feature_methods)
    if ~isempty(results_features{i})
        fprintf('%-20s %-15.3f %-15.3f %-15.3f\n', ...
            ['Feature_' feature_methods{i}], ...
            results_features{i}.summary.overall_auc_mean, ...
            results_features{i}.summary.overall_accuracy_mean, ...
            nanmean(results_features{i}.summary.f1_score_mean));
    end
end

% CV methods
for i = 1:length(cv_methods)
    if ~isempty(results_cv{i})
        fprintf('%-20s %-15.3f %-15.3f %-15.3f\n', ...
            ['CV_' cv_methods{i}], ...
            results_cv{i}.summary.overall_auc_mean, ...
            results_cv{i}.summary.overall_accuracy_mean, ...
            nanmean(results_cv{i}.summary.f1_score_mean));
    end
end

% Age-binned analysis
if ~isempty(results_age)
    fprintf('%-20s %-15.3f %-15.3f %-15.3f\n', ...
        'Age_Binned_Advanced', ...
        results_age.summary.overall_auc_mean, ...
        results_age.summary.overall_accuracy_mean, ...
        nanmean(results_age.summary.f1_score_mean));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO 6: VISUALIZATION AND PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== DEMO 6: Creating Summary Visualizations ===\n');

% Create comparison plots
figure('Position', [100, 100, 1200, 800]);

% Subplot 1: Classifier comparison
subplot(2, 3, 1);
classifier_names = {};
classifier_aucs = [];
for i = 1:length(classifiers)
    if ~isempty(results_basic{i})
        classifier_names{end+1} = classifiers{i};
        classifier_aucs(end+1) = results_basic{i}.summary.overall_auc_mean;
    end
end
bar(classifier_aucs);
set(gca, 'XTickLabel', classifier_names);
title('Classifier Performance Comparison');
ylabel('AUC');
xtickangle(45);

% Subplot 2: Feature selection comparison
subplot(2, 3, 2);
feature_names = {};
feature_aucs = [];
for i = 1:length(feature_methods)
    if ~isempty(results_features{i})
        feature_names{end+1} = feature_methods{i};
        feature_aucs(end+1) = results_features{i}.summary.overall_auc_mean;
    end
end
bar(feature_aucs);
set(gca, 'XTickLabel', feature_names);
title('Feature Selection Performance');
ylabel('AUC');
xtickangle(45);

% Subplot 3: CV method comparison
subplot(2, 3, 3);
cv_names = {};
cv_aucs = [];
for i = 1:length(cv_methods)
    if ~isempty(results_cv{i})
        cv_names{end+1} = cv_methods{i};
        cv_aucs(end+1) = results_cv{i}.summary.overall_auc_mean;
    end
end
bar(cv_aucs);
set(gca, 'XTickLabel', cv_names);
title('Cross-Validation Method Performance');
ylabel('AUC');
xtickangle(45);

% Subplot 4: Age-binned results (if available)
subplot(2, 3, 4);
if ~isempty(results_age) && isfield(results_age.summary, 'auc_mean')
    age_bins = 1:length(results_age.summary.auc_mean);
    errorbar(age_bins, results_age.summary.auc_mean, results_age.summary.auc_std, 'o-', 'LineWidth', 2);
    title('Age-Binned Classification Performance');
    xlabel('Age Bin');
    ylabel('AUC');
    grid on;
end

% Subplot 5: Performance metrics comparison
subplot(2, 3, 5);
if ~isempty(results_basic{1}) % Use first classifier as reference
    metrics = {'auc_mean', 'accuracy_mean', 'sensitivity_mean', 'specificity_mean'};
    metric_names = {'AUC', 'Accuracy', 'Sensitivity', 'Specificity'};
    metric_values = [];
    for i = 1:length(metrics)
        if isfield(results_basic{1}.summary, metrics{i})
            metric_values(i) = results_basic{1}.summary.(metrics{i});
        end
    end
    bar(metric_values);
    set(gca, 'XTickLabel', metric_names);
    title('Performance Metrics (SVM)');
    ylabel('Score');
    ylim([0, 1]);
end

% Subplot 6: Summary statistics
subplot(2, 3, 6);
text(0.1, 0.9, 'FEMA Classification Enhanced Demo', 'FontSize', 14, 'FontWeight', 'bold');
text(0.1, 0.8, sprintf('Total Tests: %d', length(classifiers) + length(feature_methods) + length(cv_methods) + 1), 'FontSize', 12);
text(0.1, 0.7, sprintf('Best AUC: %.3f', max([classifier_aucs, feature_aucs, cv_aucs])), 'FontSize', 12);
text(0.1, 0.6, 'Output saved to:', 'FontSize', 12);
text(0.1, 0.5, output_base, 'FontSize', 10);
axis off;

% Save the summary figure
saveas(gcf, fullfile(output_base, 'FEMA_classify_enhanced_demo_summary.png'));
fprintf('Summary visualization saved to: %s\n', fullfile(output_base, 'FEMA_classify_enhanced_demo_summary.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FINAL SUMMARY AND CLEANUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n=== DEMO COMPLETED ===\n');
fprintf('All results saved to: %s\n', output_base);
fprintf('Check the output directory for detailed results and visualizations.\n');
fprintf('Refer to FEMA_classify_enhanced_README.md for detailed documentation.\n');

% Display best performing method
all_aucs = [classifier_aucs, feature_aucs, cv_aucs];
all_names = [classifier_names, feature_names, cv_names];
[best_auc, best_idx] = max(all_aucs);

fprintf('\nBest performing method: %s (AUC = %.3f)\n', all_names{best_idx}, best_auc);

% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE. 