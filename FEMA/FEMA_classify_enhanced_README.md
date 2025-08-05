# FEMA Classification Enhanced

## Overview

`FEMA_classify_enhanced.m` is an improved version of the original `FEMA_classify.m` function that provides enhanced classification capabilities for neuroimaging data within the FEMA framework. This enhanced version addresses the limitations of the original function and adds numerous new features for robust machine learning classification.

## Key Improvements Over Original Function

### 1. **Fixed Dependencies**
- Added missing `plot_ROC()` function with comprehensive performance metrics
- Implemented `cmig_tools_cohensd()` for effect size calculations
- Fixed `nancorr()` function for correlation analysis with missing data

### 2. **Enhanced Algorithm Support**
- **Multiple Classifiers**: SVM, Logistic Regression, Random Forest, Ensemble methods
- **Advanced Feature Selection**: SVD, LASSO, Mutual Information, None
- **Hyperparameter Tuning**: Automatic optimization of model parameters
- **Cross-validation Methods**: K-fold, Stratified, Family-based holdout

### 3. **Improved Robustness**
- Comprehensive input validation
- Better error handling and informative messages
- Progress tracking and verbose output options
- Memory management and parallel processing support

### 4. **Enhanced Output and Visualization**
- Multiple output formats (MAT, JSON, CSV)
- Comprehensive performance metrics
- Automated visualization generation
- Detailed summary statistics

## Usage

### Basic Usage

```matlab
% Basic classification with default settings
[results, models] = FEMA_classify_enhanced('thickness-sm16', ...
    {'design_matrix.txt'}, {'output_dir'}, '/path/to/imaging', ...
    'vertex', {'diagnosis'});
```

### Advanced Usage with Custom Parameters

```matlab
% Advanced classification with custom parameters
[results, models] = FEMA_classify_enhanced('FA', ...
    {'design_matrix.txt'}, {'output_dir'}, '/path/to/imaging', ...
    'voxel', {'diagnosis'}, ...
    'classifier', 'ensemble', ...
    'cv_method', 'family_holdout', ...
    'n_folds', 5, ...
    'n_repetitions', 20, ...
    'feature_selection', 'lasso', ...
    'n_components', 100, ...
    'hyperparameter_tuning', true, ...
    'age_binning', true, ...
    'age_bin_size', 2.0, ...
    'output_format', 'json', ...
    'save_models', true, ...
    'parallel', true, ...
    'verbose', true);
```

## Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `fstem_imaging` | char | Name of imaging phenotype (e.g., 'thickness-sm16', 'FA') |
| `fname_design` | cell | Cell array of design matrix file paths |
| `dirname_out` | cell | Cell array of output directory paths |
| `dirname_imaging` | char | Path to imaging data directory |
| `datatype` | char | Data type: 'voxel', 'vertex', 'external', 'corrmat' |
| `target` | cell | Target variable names for classification |

### Optional Parameters

#### Classification Algorithm
| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `classifier` | 'svm' | 'svm', 'logistic', 'randomforest', 'ensemble' | Classification algorithm |
| `hyperparameter_tuning` | true | true/false | Whether to tune hyperparameters |

#### Cross-Validation
| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `cv_method` | 'kfold' | 'kfold', 'stratified', 'family_holdout' | Cross-validation method |
| `n_folds` | 10 | integer | Number of CV folds |
| `n_repetitions` | 10 | integer | Number of CV repetitions |

#### Feature Selection
| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `feature_selection` | 'svd' | 'none', 'svd', 'lasso', 'mutual_info' | Feature selection method |
| `n_components` | 250 | integer | Number of components/features |

#### Age Analysis
| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `age_binning` | true | true/false | Whether to analyze by age bins |
| `age_bin_size` | 2.5 | float | Age bin size in years |

#### Output and Processing
| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `output_format` | 'mat' | 'mat', 'json', 'csv' | Output file format |
| `save_models` | false | true/false | Whether to save trained models |
| `parallel` | false | true/false | Whether to use parallel processing |
| `verbose` | true | true/false | Whether to display progress |

#### Data Processing
| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `ranknorm` | false | true/false | Apply rank normalization |
| `varnorm` | false | true/false | Apply variance normalization |
| `ico` | 5 | integer | Icosahedron order for vertex data |

## Output Structure

### Results Structure
```matlab
results(design_number).design_file     % Design matrix file path
results(design_number).target          % Target variable name
results(design_number).classifier      % Used classifier
results(design_number).cv_method       % CV method used
results(design_number).n_folds         % Number of CV folds
results(design_number).n_repetitions   % Number of CV repetitions
results(design_number).feature_selection % Feature selection method
results(design_number).n_components    % Number of components
results(design_number).age_binning     % Whether age binning was used
results(design_number).age_bin_size    % Age bin size
results(design_number).cv_results      % Cross-validation results
results(design_number).summary         % Summary statistics
```

### CV Results Structure
```matlab
cv_results.auc              % AUC values [n_bins x n_reps]
cv_results.accuracy         % Accuracy values [n_bins x n_reps]
cv_results.sensitivity      % Sensitivity values [n_bins x n_reps]
cv_results.specificity      % Specificity values [n_bins x n_reps]
cv_results.f1_score         % F1 score values [n_bins x n_reps]
cv_results.predictions      % Predictions [n_bins x n_reps cell]
cv_results.probabilities    % Probabilities [n_bins x n_reps cell]
```

### Summary Structure
```matlab
summary.auc_mean            % Mean AUC across repetitions
summary.auc_std             % Standard deviation of AUC
summary.auc_median          % Median AUC
summary.auc_ci_95           % 95% confidence interval for AUC
summary.accuracy_mean       % Mean accuracy across repetitions
summary.accuracy_std        % Standard deviation of accuracy
% ... similar for other metrics
summary.overall_auc_mean    % Overall AUC across all age bins
summary.overall_auc_std     % Overall AUC standard deviation
```

## Performance Metrics

The enhanced function calculates comprehensive performance metrics:

- **AUC**: Area Under the ROC Curve
- **Accuracy**: Overall classification accuracy
- **Sensitivity/Recall**: True Positive Rate
- **Specificity**: True Negative Rate
- **Precision**: Positive Predictive Value
- **F1 Score**: Harmonic mean of precision and recall
- **Balanced Accuracy**: Average of sensitivity and specificity
- **Matthews Correlation**: Correlation coefficient for binary classification

## Visualization Outputs

The function automatically generates several visualization plots:

1. **Performance Across Age Bins**: Line plot showing performance metrics across age bins
2. **ROC Curves**: ROC curves for each age bin with AUC values
3. **Confusion Matrix**: Heatmap of confusion matrix for aggregated results

## Best Practices

### 1. **Data Preparation**
- Ensure your design matrix is properly formatted
- Check for missing values and handle appropriately
- Verify target variable coding (binary: 0/1)

### 2. **Feature Selection**
- **SVD**: Good for dimensionality reduction, preserves variance
- **LASSO**: Good for feature selection with regularization
- **Mutual Information**: Good for identifying relevant features
- **None**: Use when you want to use all features

### 3. **Cross-Validation**
- **K-fold**: Standard approach, good for most cases
- **Stratified**: Maintains class balance in folds
- **Family Holdout**: Important for family-based studies

### 4. **Classifier Selection**
- **SVM**: Good for high-dimensional data, robust
- **Logistic Regression**: Interpretable, good baseline
- **Random Forest**: Handles non-linear relationships well
- **Ensemble**: Combines multiple classifiers for better performance

### 5. **Hyperparameter Tuning**
- Enable for better performance (default: true)
- May increase computation time
- Particularly important for SVM and Random Forest

## Example Workflows

### 1. **Basic Diagnostic Classification**
```matlab
% Classify diagnostic groups using cortical thickness
[results, ~] = FEMA_classify_enhanced('thickness-sm16', ...
    {'design_diagnosis.txt'}, {'results_diagnosis'}, ...
    '/path/to/abcd-sync', 'vertex', {'diagnosis'}, ...
    'classifier', 'svm', ...
    'cv_method', 'stratified', ...
    'n_repetitions', 20);
```

### 2. **Age-Specific Analysis**
```matlab
% Analyze classification performance across age groups
[results, ~] = FEMA_classify_enhanced('FA', ...
    {'design_behavior.txt'}, {'results_behavior'}, ...
    '/path/to/abcd-sync', 'voxel', {'behavior_score'}, ...
    'age_binning', true, ...
    'age_bin_size', 1.0, ...
    'classifier', 'ensemble', ...
    'feature_selection', 'lasso');
```

### 3. **Feature Selection Comparison**
```matlab
% Compare different feature selection methods
methods = {'svd', 'lasso', 'mutual_info'};
for i = 1:length(methods)
    [results(i), ~] = FEMA_classify_enhanced('thickness-sm16', ...
        {'design.txt'}, {sprintf('results_%s', methods{i})}, ...
        '/path/to/abcd-sync', 'vertex', {'target'}, ...
        'feature_selection', methods{i}, ...
        'n_components', 100);
end
```

## Troubleshooting

### Common Issues

1. **"Target not found in design matrix"**
   - Check target variable name spelling
   - Verify target variable exists in design matrix

2. **"Invalid datatype"**
   - Ensure datatype is one of: 'voxel', 'vertex', 'external', 'corrmat'

3. **"Incorrect number of input arguments"**
   - Ensure all required parameters are provided
   - Check parameter order

4. **Memory issues with large datasets**
   - Reduce `n_components` for feature selection
   - Use `'single'` precision if available
   - Enable parallel processing

### Performance Optimization

1. **For large datasets**:
   - Use SVD feature selection with fewer components
   - Reduce number of CV repetitions
   - Disable hyperparameter tuning

2. **For better accuracy**:
   - Use ensemble classifier
   - Enable hyperparameter tuning
   - Increase number of CV repetitions
   - Use family-based cross-validation for family studies

## Dependencies

The enhanced function requires:
- MATLAB Statistics and Machine Learning Toolbox
- MATLAB Parallel Computing Toolbox (for parallel processing)
- FEMA core functions (`FEMA_process_data`, `FEMA_intersect_design`)

## Citation

When using this enhanced classification function, please cite:

1. The original FEMA paper:
   Parekh et al., (2024). FEMA: Fast and efficient mixed-effects algorithm for large sample whole-brain imaging data. Human Brain Mapping, 45(2), e26579.

2. This enhanced classification implementation (if appropriate for your work).

## Support

For issues, questions, or feature requests, please:
1. Check the troubleshooting section above
2. Review the example workflows
3. Open an issue on the GitHub repository
4. Contact the CMIG research group

## Version History

- **v2.0**: Enhanced version with multiple classifiers, improved CV, and comprehensive metrics
- **v1.0**: Original FEMA_classify.m function 