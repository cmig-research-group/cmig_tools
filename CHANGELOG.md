# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) starting from ``v1.0.0``.

## [Unreleased - 2024-07-04]

### Added
* ``FEMA_parse_PLINK`` function that calls either ``PlinkRead_binary2`` or ``PlinkRead_binary2_subj`` for reading genotyping data
* ``FEMA_compileTerms`` for compiling inverse of the V term for every bin
* ``FEMA_OLSResiduals`` for calculating OLS residuals, given a set of fixed effects
* ``FEMA_GLSResiduals`` for calculating GLS residuals, given a set of fixed effects and random effects information (output from ``FEMA_compileTerms``)
* ``divideSNPs`` for dividing SNP data into chunks or dividing based on chromosome

### Changed
* ``createBasisFunctions``:
    - fixed dropping columns
    - additionally saving settings and returning timing information
    - added SVD based modification of basis functions
    - default splineType is now ns, consistent with the R version
    - default column to drop is now the middle column, consistent with R version
    - added functionality for instance number (useful if called in parallel)
    - cleaned up documentation / explanation of settings
    - now allows character type input for toDrop
    - generates a warning if there is a mismatch in number of observations being read back into MATLAB
* ``caller_createBasis.R``: 
    - remove_col is set to "none"
    - added a comment explaining that the functionality for dropping columns is handled directly within MATLAB

### Deleted
* ``createBasisNS.R``
* ``FEMA_create_basisFunctions.m``: functionality merged with ``create_basisFunctions.m``

## [Unreleased - 2023-11-11]

### Added
* ``lsqnonneg_amd3``: modified version of ``lsqnonneg_amd2`` that avoids using pinv; should lead to performance improvement; additionally, residual and resnorm are only returned if user wants (should reduce memory use)
* ``doFEMA_tests``: a script that performs a series of ``FEMA_fit`` tests using synthetic data (including comparing estimates with ``fitlmematrix``

### Changed
* ``colvec``: switched to numel instead of prod (trivial performance improvement)
* ``rowvec``: switched to numel instead of prod (trivial performance improvement)
* ``FEMA_parse_family``: performance improvement (summary below)
  - performing a string comparison by creating a string version of ``freq_unique`` and ``famtypelist`` (instead of ``ivec = find(cellfun(@(x)isequal(x,freq_unique),famtypelist));``)
  - replaced repeated accessing of clusterinfo by converting to cell and then accessing contents of the cell
  - constructing Ss using rows, columns, and values instead of spalloc initially and then iteratively changing values inside - should lead to much better performance
* ``FEMA_fit``: combined version of Anders' ``FEMA_DEAP_fit`` and `GWAS_v2` branch ``FEMA_fit``
  - now includes IGLS
  - additionally returning reusableVars, if the user specifies ``returnReusable`` as ``true`` - a structure having some settings, residuals, and MSE which can be reused for FEMA-GWAS
  - slight changes to OLS estimation
  - returning non-permuted binvec (i.e., returns binvec_save)
  - additionally returning FamilyStruct - for use in DEAP
  - improved documentation of input and output variables
  - uses ``lsqnonneg_amd3`` instead of ``lsqnonneg_amd2`` - should improve performance
  - subvec1, subvec2, and indvec are updated to ``find(tril(S_sum))`` - should fix increased memory requirement
  - removed redundant fields between ``reusableVars`` and ``FamilyStruct`` - should reduce memory requirement
  - no longer explicitly saving ``ymat_hat_ols`` and ``ymat_res_ols`` as separate variables - directly using ``ymat_hat`` and ``ymat_res`` - should reduce memory requirement
  - merged GLS estimation from ``FEMA_sig2binseg_parfeval``
  - dropped ``reverse_cols`` and ``reverseinferenceflag``
  - roll back of the binning strategy to the older version
  - additionally handling cases where bin might be assigned as NaN
  - introduced a quick check of X and ymat to ensure no NaN or Inf
  - ``binvec_save`` is saved when ``permi==0``
  - more robust handling of singular or nearly singular matrix; no longer displaying warnings (use pinv if a warning is generated)
* ``FEMA_run_on_synthetic_data``: updated to match the output(s) from ``FEMA_fit``

### Deleted
* ``FEMA_sig2binseg_parfeval``: now merged with FEMA_fit (summary of changes done prior to merging below; these changes are now part of ``FEMA_fit``):
  - major code rewrite and formatting changes
  - inverses are calculated using backslash (only in rank deficient cases we either use lsqminnorm or fall back to pinv)
  - replaced repeated accessing of clusterinfo by converting to cell and then accessing contents of the cell
  - no longer computing XtWVsWt

## [Unreleased]

### Added (major options)
* _add your changes here_

### Added (minor options)
* ``FEMA_imageMosaic``: this function creates mosaic plots for volumetric statistics (fixed or random effects)
* ``FEMA_lookupVertices``: this function returns the DK40-based labels for vertex-wise statistics (fixed or random effects)
* Supporting functions for ``FEMA_imageMosaic``:
  - ``cmig_tools_utils/matlab/calc_rows_cols_subplot.m``: this function attempts to figure out an optimal combination of rows and columns for the mosaic plot
  - ``cmig_tools_utils/matlab/plotboxpos.m``: this function, from https://www.mathworks.com/matlabcentral/fileexchange/9615-plotboxpos, returns position for actual plotted data/image and not just the overall axes
* Added ``CHANGELOG.md``

### Changed
* ``cmig_tools_utils/r/makeDEAPdemos.R``: updated for the ABCD Study 5.0 data release

### Fixed

* _add your changes here_

### Removed

* _add your changes here_

## [2.2.1] - 2023-01-20

### Fixed

* amd update NIFTI functions by @dmysmith in 10eaa08
* ML CI calculation: allow nuisance parameters to vary by @dmysmith in b92e7b4
* amd update output to json for DEAP by @dmysmith in b1120ae

## [2.2.0] - 2022-11-10

### Added

* M,P,T,H random effects for MoBa
* add "home ID' (compatible with ABCD)
* helper function ``mvnpdfln.m``
* ``ABCD_DEMO_test.m``: script for testing the files in ``ABCD_DEMOS``

### Fixed

* anders edits to clean up ML estimation and random effects
* set minimum sig2vec to 0 in ``FEMA_fit.m``
* fix path to GRM file in ``FEMA_wrapper_demo.m``

## [v2.1.0] - 2022-10-24

### Added
* option to specify ML as estimation method

### Changed

* Enable logLikvec computation by default
* Update ``README.md``

## [v2.0.1] - initial public release
