# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) starting from ``v1.0.0``.

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
