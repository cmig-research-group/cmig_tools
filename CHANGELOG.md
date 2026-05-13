# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) starting from ``v1.0.0``.

## [4.0.0 - 2026-05-13]
### Added
* Major new functions added:
    - `FEMA/FEMA_makeDesign.m`: intended replacement for `cmig_tools_utils/r/makeDesign.R`
    - `FEMA/FEMA_read_GRM.m`: function for reading in genetic relationship matrix
    - `FEMA/FEMA_save.m`: multi-format saving function
    - `FEMA/writeGIfTI.m`: function that writes out statistical overlays as GIfTI files
    - `FEMA/writeGIfTI_underlay.m`: function to generate underlay and atlas files in GIfTI format
    - `FEMA/FEMA_propagateStats.m`: function that propagates ROI-level statistics to voxel/vertex-level
* To support the above functions, various supporting functions have been added:
    - `FEMA/FEMA_createInputJSON.m`
    - `FEMA/FEMA_loadDataDir.m`
    - `FEMA/FEMA_mergeArgs.m`
    - `FEMA/FEMA_parseInputs.m`
    - `FEMA/FEMA_parse_JSON.m`
    - `FEMA/FEMA_parse_config.m`
    - `FEMA/FEMA_readToml.m`
    - `cmig_tools_utils/matlab/numIcoVertices.m`
    - `cmig_tools_utils/matlab/read_fs_annot.m`
    - `cmig_tools_utils/matlab/dispNameStat.m`
    - `cmig_tools_utils/matlab/doTransformation.m`
    - `cmig_tools_utils/matlab/winsorize.m`
    - `cmig_tools_utils/matlab/local_mapEmbeddedToken.m`
    - `cmig_tools_utils/matlab/save_jsonencode.m`
    - `cmig_tools_utils/matlab/save_nonfiniteScalarTokens.m`
    - `cmig_tools_utils/matlab/save_numericToJsonCells.m`
    - `cmig_tools_utils/matlab/save_prepareStrictJson.m`
    - `cmig_tools_utils/matlab/toml_parse.m`
* Additional supporting files and functions for atlases:
    - `cmig_tools_utils/matlab/aparc2009s2tabroinames.m`
    - `cmig_tools_utils/matlab/aparc2tabroinames.m`
    - `cmig_tools_utils/matlab/aseg2tabroinames.m`
    - `cmig_tools_utils/matlab/atlas2tabNames.m`
    - `cmig_tools_utils/matlab/createDEAPlut.m`
    - `cmig_tools_utils/matlab/fiber2tabroinames.m`
    - `cmig_tools_utils/matlab/roi2ABCDAtlas.m`
    - `cmig_tools_utils/matlab/roi2FreesurferAtlas.m`
    - `support_files/aparc_a2009s_lut.csv`
    - `support_files/aparc_lut.csv`
    - `support_files/aparcaseg_lut.csv`
    - `support_files/aseg_lut.csv`
    - `support_files/fiber_lut.csv`
    - `support_files/roi2ABCDAtlasMaps.mat`
    - `support_files/roi2SurfaceAtlasMaps.mat`
* New recipes have been introduced:
    - `recipes/demo_run_unstructuredCovariance.m`
    - `recipes/run_GIfTI_tests.m`
* Miscellaneous supporting files
    - `FEMA/Makefile`
    - `FEMA/compile_stats.sh`
    - `cmig_tools_utils/matlab/tight_subplot.m`
    - `showVol/utils/jsystem.m`
* New files for DEAP parallelization are included (feature under development):
    - `DEAP/FEMA_DEAP_client.m`
    - `DEAP/FEMA_DEAP_gencache.m`
    - `DEAP/FEMA_DEAP_init.m`
    - `DEAP/FEMA_DEAP_server.m`
    - `DEAP/FEMA_DEAP_worker.m`
    - `DEAP/FEMA_DEAP_wrapper.m`
    - `DEAP/FEMA_DEAP_wrapper_submit.m`
    - `DEAP/GlobalConfig.m`
    - `DEAP/compile_DEAP.m`
    - `DEAP/run_client.sh`
    - `DEAP/run_gencache.sh`
    - `DEAP/start_workers.sh`
* A copy of the [GIfTI toolbox](https://github.com/gllmflndn/gifti) is included; license file added to `docs/licenses`

### Changed
* `FEMA/FEMA_fit.m`: major changes impacting backward compatiblity
    - Output changes:
        - replaced `reusableVars` with `unstructParams` and `residuals_GLS`
        - `unstructParams` is a structure with fields relevant to unstructured covariance; this is only output when `CovType` is `unstructured`
        - `residuals_GLS` is a vector or matrix having GLS residuals of `ymat`; this is only output when `returnResiduals` (optional input) is `true`
        - `info` is a new output variable containing settings and various timing information
    - Input changes:
        - `niter` is now an optional parameter (was previously mandatory; not backward compatible)
        - `SingleOrDouble` parameter replaced by `precision` (functionality unchanged; not backward compatible)
    - included checks to make sure that `X` is either single or double precision; if not, the data type is changed
    - warning generated if `ymat` has any constant value columns
    - all rank calculations use double precision
    - ensure that initialization follows `precision` parameter, if specified
    - calculating and returning timing for different sections (see `info.timing`)
    - removed iterative GLS code
    - code cleanup
* `FEMA/FEMA_WaldTest.m`: major changes impacting backward compatiblity
    - only `logp` values are now returned, consistent with `FEMA_fit` (not backward compatibile)
    - `logp` values are not truncated to realmin; it is possible to get `Inf` logp values (if p value was zero) (changed behaviour from previous versions)
    - performance update
* `FEMA/FEMA_fit_GWAS.m`: small but important changes that impact backward compatibility
    - `Wald_p` replaced by `Wald_logp`
    - `SingleOrDouble` parameter replaced by `precision` (functionality unchanged; not backward compatible)
    - using `nearestSPD_timeout` instead of `nearestSPD`; if `nearestSPD_timeout` fails to converge, `errFlag` is set to true
* `FEMA/FEMA_compileTerms.m`: `SingleOrDouble` parameter replaced by `precision` (functionality unchanged; not backward compatible)
* `FEMA/FEMA_gatherGWAS.m`: `Wald_p` replaced by `Wald_logp`
* `FEMA/FEMA_intersect_design.m`: changes impact backward compatibility
    - `contrats` and `colnames_model` are not output
    - `info` is an additional output
    - code cleanup
    - missingness information is saved in `info`
* `FEMA/FEMA_parse_contrastFile`: `contrastNames` is an additional output
* `FEMA/FEMA_parse_family.m`: for `F,E` model, `iid` is set to `fid`; for `S,E` model, `fid` is set to `iid`
* `FEMA_process_data.m`: major changes, including ones that impact backward compatibility
    - `colnames_imaging` is not output
    - `info` is an additional output
    - `wholeSampleTransform`, `fname_qc`, `qc_var` are new input parameters
    - `ranknorm` and `varnorm` are no longer input parameters
    - added validation functions for optional input parameters
    - code cleanup
    - included support for parquet files
    - major changes to handling `external` data type
    - included filtering on QC
    - recording and returning missingess via the `info` variable
    - optionally filter on iid and/or eid
* `FEMA_wrapper.m`: major changes, including changes that impact backward compatibility
    - `fpaths_out` is not output; changes to output parameters, consistent with output from `FEMA_fit`
    - `SingleOrDouble` parameter replaced by `precision` (functionality unchanged; not backward compatible)
    - added support for handling DEAP inputs: JSON config file, data file, and output directory
    - added support for non-DEAP inputs: TOML file
    - changes to optional input arguments to be compatible with other changes in 4.0.0
    - HBCD family ID file path hardcoded
    - added support for making design matrices by calling `FEMA_makeDesign`
    - logging timing and other information at different places using the `info` structure
    - changed call to `FEMA_fit` and gathering outputs from `FEMA_fit` to be compatible with the updates above
    - results savings has been delegated to `FEMA_save`
* `FEMA/FEMA_classify.m`: replaced the use of `crossvalind` with `cvparition`
* `FEMA/FEMA_run_on_synthetic_data.m`:
    - performance update to synthesising IDs
    - changes to calling `FEMA_fit` to be compatible with the changes above
* `FEMA/caller_FEMA_fit.m`:
    - `SingleOrDouble` parameter replaced by `precision` (functionality unchanged; not backward compatible)
    - `returnReusable` parameter replaced by `returnResiduals`
    - generates an error if no data is left after intersection of `X` and `ymat`
    - changed call to `FEMA_fit` and gathering outputs from `FEMA_fit` to be compatible with the updates above
    - results savings has been delegated to `FEMA_save`
* `FEMA/caller_FEMA_fit_GWAS.m`:
    - `doCoeffCovar` default changed to `true`
    - `SingleOrDouble` parameter replaced by `precision` (functionality unchanged; not backward compatible)
    - added `cleanUp` as an additional input argument for removing piece-wise GWAS statistics after the overall summary statistics has been generated
    - no longer needs settings file to be saved; uses `info` from the estimates file
    - updates to make code compatible with above changes
    - additionally gathering summary statistics into a single file
* `FEMA/caller_createBasisFunctions.m`: 
    - bug fix: replace incorrect use of `strreplace` with `replace`
    - bug fix handling of `knots` not being specified
* `FEMA/logging.m`: code clean up; using `datetime`
* `FEMA/niftiwrite_amd.m`:
    - replaced `symmetric_cLim` input argument with `cLim_range`
    - calculates colour limits omitting missing
    - if input `cLim_range` is scalar, symmetric colour range is written
* `FEMA/FEMA_info.m`: updated version information
* `cmig_tools_utils/matlab/PrintMemoryUsage.m`: code clean up
* `cmig_tools_utils/matlab/SwEfit_amd2.m`: fixed mismatched function name
* `cmig_tools_utils/matlab/blueblackred.m`: fixed mismatched function name
* `cmig_tools_utils/matlab/colvec.m`: using `numel` instead of `prod`
* `cmig_tools_utils/matlab/createBasisFunctions.m`: `knots` can be specified as `quartiles`
* `cmig_tools_utils/matlab/sub2ind_amd.m`: fixed mismatched function name
* `cmig_tools_utils/matlab/subsample_volume.m`: defaulting to 2, if input unspecified
* `recipes/demo_doGWAS_Long_splines.m`: edits for compatibility with above changes; added demo for using `FEMA_convert_splines`
* corresponding `docs/help` files updated

### Deleted
* `cmig_tools_utils/matlab/read_annotations.m`
* `cmig_tools_utils/matlab/mmil_args2parms.m`
* `cmig_tools_utils/matlab/mmil_check_nargs.m`
* `cmig_tools_utils/matlab/mmil_cmap_blueblackred.m`
* `cmig_tools_utils/matlab/mmil_colvec.m`
* `cmig_tools_utils/matlab/mmil_dilate_mask.m`
* `cmig_tools_utils/matlab/mmil_parms2args.m`
* `cmig_tools_utils/matlab/mmil_smooth3d.m`
* `cmig_tools_utils/matlab/plotboxpos.m`
* `cmig_tools_utils/matlab/struct2args.m`
* `docs/FEMAv2.0_UserGuide.pdf`
* `docs/licenses/LICENSE_plotboxpos.txt`
* `admin/ABCD_DEMO_test.m`
* corresponding `docs/help` files removed

## [3.0.0 - 2025-09-17]
### Changed
* `FEMA_wrapper.m`: fixed parameter `pihat` to `GRM`
* `recipes/demo_doGWAS_Long_splines.m`: added `outPrefix` as additional input to `divideSNPs`
* `FEMA_info.m`: Updated version information to include date
* `caller_FEMA_fit.m`: patched `RandomEffects` parameter check for deployed version

## [3.0.0 - 2025-09-11]
### Added
* `FEMA/FEMA_info.m`: function that returns FEMA version and citation information
* `docs/citationInfo.txt`: FEMA citation information
* `docs/help/FEMA_info.txt`: help for `FEMA_info`
* Included a Linux and MacOS compiled version of the FEMA app

### Changed
* `FEMA/caller_FEMA_fit_GWAS.m`: handling logical values which are input as character type
* Included displaying version information for:
    - `caller_FEMA_fit`
    - `caller_FEMA_fit_GWAS`
    - `FEMA_classify`
    - `FEMA_fit`
    - `FEMA_wrapper`
* Included version and citation calls in `caller_FEMA`
* Labeling version as 3.0.0

## Deleted
* Linux and MacOS compiled version of the FEMA app (release binaries now)

## [Unreleased - 2025-08-25]
### Changed
* `FEMA/caller_FEMA.m`:
    - passing unmatched arguments to `caller_FEMA_fit_GWAS`
    - additional name possibility for the GWAS module
* `FEMA/caller_FEMA_fit.m`:
    - fixed incorrect specification for binary file for GRM
    - updated requirement to have `uqObservations` available as part of GRM
    - if input is a binary file, allow same basename `.mat` file to be present which specifies the ID list
* `FEMA/caller_FEMA_fit_GWAS.m`:
    - fixed typo in accessing basis function variable
    - initializing `contrasts` and `hypValues` to empty when `file_contrast` is not passed
    - fixed checking correct variables for design matrix
* Removed unnecessary files from `docs/help`

## [Unreleased - 2025-08-25]
### Added
* `FEMA/FEMA_exportDoc.m`: function that exports headers as text file
* `FEMA/caller_FEMA.m`: main caller function for compiled version of FEMA
* `FEMA/caller_FEMA_fit.m`: caller function for compiled version of `FEMA_fit`
* `FEMA/caller_FEMA_fit_GWAS.m`: caller function for compiled version of `FEMA_fit_GWAS`
* `FEMA/caller_createBasisFunctions.m`: caller function for compiled version of `createBasisFunctions`
* Help files created by `FEMA_exportDoc` to be used by compiled version of FEMA
* `FEMA/FEMA_parse_contrastFile.m`: function that reads in a contrast file and constructs zero-padded univarate and multivariate contrasts

### Changed
* `FEMA/FEMA_WaldTest.m`: additionally returns `LB_hat` and `LB_SE`
* `FEMA/FEMA_fit.m`: additionally record time taken to create and delete a parallel pool
* `cmig_tools_utils/matlab/divideSNPs.m`: added option to pass in an output prefix
* `FEMA/FEMA_gatherGWAS.m`: bug fix #34 via PR#37
* replaced occurrence of pihat (and related words) in several functions to GRM
* `cmig_tools_utils/matlab/PlinkRead_binary2_subj.m`: added functionality to read SNPs in chunks
* `FEMA/niftiwrite_amd.m`: saving the min and max values in the NIfTI header

## [Unreleased - 2025-08-13]
### Added
* `showVol/utils/rcs2rcs.m`: function to convert rcs coordinates of one image to another

## [Unreleased - 2025-08-01]
### Added
* multiple helper functions for `save_showVol_images`:
    - `showVol/utils/anatomyDrawRoiOutline_save_images.m`
    - `showVol/utils/anatomyAddRoiOverlay_save_images.m`

### Changed
* `‎showVol/utils/save_showVol_images.m`: added option for ROI outlining and taking screenshots

## [Unreleased - 2025-07-30]
### Changed
* `‎showVol/showVol.m`: added option for selecting orientation for initial view

## [Unreleased - 2025-07-13]
### Changed
* `‎showSurf/showSurfPlot`: allow controlling for subplot number and dimension
* renamed `ABCD_DEMOS` to `recipes`

### Added
* recipe `makeGRM`
* recipe `doGWAS_Long_splines`

## [Unreleased - 2025-07-11]
### Changed
* `showSurf/showSurfPlot`: added additional view options for surface results
### Added
* `save_showVol_images`: automatically save showVol images as screenshots
* `cmig_tools_utils/r/makeDesign_DEAP.R`: DEAP compatible version of `makeDesign`

## [Unreleased - 2025-07-09]
### Changed
* `FEMA/FEMA_process_data`:
    - robust handling of path to `SurfView_surfs.mat`
    - only loading `icsurfs` from `SurfView_surfs.mat`
    - minor formatting fixes
* `FEMA/FEMA_wrapper`:
    - robust handling of path to `SurfView_surfs.mat`
    - only loading `icsurfs` from `SurfView_surfs.mat`

## [Unreleased - 2025-07-07]
### Changed
* `FEMA/FEMA_gatherGWAS.m`: bug fix incorrect dimensionality when single phenotype, multiple coefficients [issue #34]

## [Unreleased - 2025-06-24]
### Changed
* `cmig_tools_utils/matlab/lsqnonneg_amd3.m`: reverting to using `pinv` for calculating inverse

## [Unreleased - 2025-06-17]
### Added
* `cmig_tools_utils/r/fema_env.yml`: configuration file to create conda environment to run `makeDesign.R`

## [Unreleased - 2025-06-14]
### Changed
* `FEMA/FEMA_fit.m`: bug fix OLS inversion when using `pinv` [issue #33]

## [Unreleased - 2025-05-22]
### Changed
* `FEMA/FEMA_process_data.m`: added `participant_id` and `session_id` for external analyses

## [Unreleased - 2025-05-17]
### Changed
* `cmig_tools_utils/r/makeDesign.R`: added a `study` parameter

## [Unreleased - 2025-05-02]
### Changed
* `FEMA/FEMA_classify.m`: added option to skip bin by age
* `FEMA/FEMA_process_data.m`: added option to vary threshold for correlation to atlas
* `FEMA/FEMA_wrapper.m`: added option to vary threshold for correlation to atlas

## [Unreleased - 2025-04-19]
### Changed
* Updated `README`

## [Unreleased - 2025-04-18]
### Added
* Added figures to `docs/pictures/showVolWiki/` and `docs/pictures/showSurfWiki/`
* Added `showVol` and `showSurf` wiki
* Added `.gitkeep`

## [Unreleased - 2025-04-01]
### Changed
* `cmig_tools_utils/r/makeDesign.R`:
    - added `filtervar` option
    - documentation update

## [Unreleased - 2025-03-27]
### Changed
* `cmig_tools_utils/matlab/plotManhattan.m`:
    - added option to filter p values
    - added addtional `stark` style option
    - added primitive support to highlight certain SNPs
    - additionally returning various components of the plot for selective rasterization
    - setting default plotting point size to 30
    - updated min and max range of y axis
    - detecting if the input values are p values instead of -log10 p values
* `cmig_tools_utils/matlab/plotMiami.m`:
    - added option to filter p values
    - added additional styles: 'mono2', 'diverge-diff2', 'stark'
    - allow user to control tick marks
    - added primitive support to highlight certain SNPs
    - additionally returning various components of the plot for selective rasterization
    - axis colours are handled separately
    - setting default plotting point size to 30
    - updated min and max range of y axis
    - detecting if the input values are p values instead of -log10 p values

## [Unreleased - 2025-03-13]
### Changed
* `cmig_tools_utils/matlab/extract_roi_val.m`: added option to extract ROI values from 2D matrix

## [Unreleased - 2025-03-07]
### Changed
* `FEMA_fit.m`: ensuring `contrasts` is a matrix if reading from a file

## [Unreleased - 2025-03-06]
### Added
* `FEMA_classify.m`: for classification of binary outcome variables

## [Unreleased - 2025-03-05]
### Changed
* `cmig_tools_utils/r/makeDesign.R`:
    - fixed typo in family ID variable
    - installing libraries if they are not found

## [Unreleased - 2025-03-04]
### Changed
* `FEMA_fit.m`: 
    - Added missing parameter `useLSQ` for FEMA_unstructuredGLS.m`
    - Only returning `visitnum` if `returnReusable` is `true`
* `FEMA_unstructuredGLS`: added missing parameter `useLSQ`
* `showVol/utils/expandVol.m`: fixed stray character

## [Unreleased - 2025-03-03]
### Changed
* `showSurf/showSurfPlot.m`: minor edit in function documentation

## [Unreleased - 2025-03-01]
### Changed
* `FEMA_fit.m`:
    - Added support for parallel processing for implementing GLS solution in case of unstructured covariance
    - Explicitly compiling `allR`, `allC`, and `allSz` once prior to performing GLS
    - Moved the GLS estimation for unstructured covariance to `FEMA_unstructuredGLS.m`
* `cmig_tools_utils/matlab/plotManhattan.m`:
    - Fixed outdated documentation
    - Can optionally take custom GWAS threshold
    - Introduced five different color styles
    - Output figure size is 12 cm tall and 18 cm wide (better aligned with publication style)
    - Can optionally plot separate colors for significant and non-significant SNPs
    - Better alignment of chromosome numbers
    - Handling NaN / Inf / complex valued p values: these are removed
* `cmig_tools_utils/matlab/plotQQ.m`:
    - Fixed outdated documentation
    - Optionally user can turn off axis labels
    - Additionally returning axis and line handles
    - Handling NaN / Inf / complex valued p values: these are removed
    - Automatically converts log p values to p values
    - Making sure y ticks are not too crowded
    - Handling case when a figure handle is passed instead of axis handle

### Added
* `FEMA_unstructuredGLS.m`: performs GLS estimation for unstructured covariance for a given phenotype
* `FEMA_convert_splines.m`: generic utility function that computes weighted combinations of splines and its derivatives
* `cmig_tools_utils/matlab/plotMiami.m`: utility function to create Miami plots (mirrored Manhattan plots)

## [Unreleased - 2025-02-27]
### Added
* `CONTRIBUTING.md`
* Issue templates
* The Wiki is live!
* Discussion board is live!

### Changed
* `FEMA_fit_GWAS.m`:
    - `genoMat` and `genStruct` are merged into `genoMat`
    - additionally calculates `coeffCovar` if user wants (default: false)
    - `coeffCovar` is only computed if basis functions are specified
    - default omnibus test is `eye`
    - default `roundOff` is `false`
    - Updated code to GLS solution for standard GWAS
    - Keeping track of any complex valued standard errors being produced
* `cmig_tools_utils/r/makeDesign.R`:
    - Included option for different family IDs
    - Making `ab_g_stc__design_id__fam` as the default family ID

## [Unreleased - 2025-02-26]
### Changed
* `cmig_tools_utils/r/makeDesign.R`: update to handle case of no categorical variable

## [Unreleased - 2025-02-24]
### Changed
* `FEMA_fit.m`: introducing unstructured covariance
    - Added support for `unstructured` `CovType`
    - Always including `E` as `RandomEffects`
    - Always ensuring `E` is the last entry in `RandomEffects`
    - Using the full matrix for `S_sum` instead of lower triangle (revert change from 2023-11-11)
    - No longer saving OLS residuals as part of `reusableVars`
    - For method of moments, using the unbiased estimator `sig2tvec` for normalizing `LHS` instead of `mean(ymat_res.^2,1)`
    - Unstructured covariances are estimated as visit * visit * RFX * y covariance matrices
    - Additional outputs in `reusableVars`
    - Binning is disabled for `unstructed` `CovType`
    - Added option for disabling binning by setting `nbins==0`
    - No longer saving binning info as part of `reusableVars` (already output as `binvec_save`)
    - Saving the time taken for entire `FEMA_fit` as `time` variable in `reusableVars`
    - Displaying correct time taken at the end of `FEMA_fit`
    - Some documentation update; fixed typo
* `FEMA_run_on_synthetic_data.m`: formatting and documentation update
* `FEMA_compileTerms.m`:
    - Updated to handle unstructured covariance
    - No longer requires `nfamtypes` as an input
    - Performance tweaks
    - Consistent handling of calculating the inverse of `V` (similar to `FEMA_fit`)

### Added
* `cmig_tools_utils/matlab/nearestSPD_timeout.m`: 
    - Modified version of `nearestSPD` which exits after a specified number of iteration
    - Additionally outputs convergance status
* `FEMA_gatherGWAS.m`: function that aggregates GWAS summary statistics across chunks
* `cmig_tools_utils/matlab/plotQQ.m`: utility function to create Q-Q plots
* `cmig_tools_utils/matlab/plotManhattan.m`: utility function to create Manhattan plots

## [Unreleased - 2025-02-19]
### Changed
* `cmig_tools_utils/matlab/Mvxl2lph_atlas.m`:
    - Atlas is a required input
    - Support for 6.0 atlas
    - Consistent atlas names with release prefix
* `cmig_tools_utils/matlab/atlas_T1.m`:
    - Atlas is a required input
    - `mask` defaults to `true`
    - Consistent atlas names with release prefix
* `cmig_tools_utils/matlab/maskBrain.m`:
    - Consistent atlas names with release prefix
* `showVol/ABCD_anatomy_visualization_example.m`:
    - Defaults to 6.0 atlas
    - Examples for 5.0 and 6.0 atlases
* `showVol/showVol.m`:
    - Defaults to 6.0 atlas
    - Tweaks for consistency across atlas names
* `showVol/utils/loadPrerenderedIfNeeded.m`:
    - Atlas is a required input
    - Consistent atlas names with release prefix

### Added
* `showVol/utils/compressVol.m`: compression utility for only storing non-zero values
* `showVol/utils/expandVol.m`: expansion counterpart of `compressVol`
* `showVol/utils/parseAtlasVersion.m`: utility function for parsing atlas names
* `showVol/utils/prepareAtlases__5060_ABCD3_cor10_orig.m`: utility for preparing atlas file(s) and 6.0 support
* `showVol/utils/showVolAtlasFile.m`: returns default atlas
* `showVol/utils/validateAtlasVersion.m`: utility function for validating atlas names

### Merged
* `prepareAtlases_ABCD3_cor10` with `showVol/utils/prepareAtlases__5060_ABCD3_cor10.m`

### Deleted
* `showVol/utils/showVol.fig_orig`
* `showVol/utils/showVol_prerend.fig`

## [Pre-3.0 release - 2025-02-15]
* Tagging branch and preparing for FEMA release 3.0

## [Unreleased - 2025-02-07]
### Changed
* `FEMA_intersect_design`:
    - Added `demean` option for design matrix
    - Handling of different event ID options
* `FEMA_process_data`:
    - deprecated `dirname_tabulated` option
    - Streamlined workflow for ABCD releases
* `FEMA_wrapper`:
    - deprecated `dirname_tabulated` option
    - Added `demean` option for design matrix

## [Unreleased - 2025-02-05]
### Changed
* `makeDesign.R`: updated visit age variable name for delta to match 6.0 naming

## [Unreleased - 2025-01-29]
* Fixed typo in documentation of `showSurfPlot.m`

## [Unreleased - 2025-01-23]
### Added
* `showSurfPlot`
### Changed
* `FEMA_convert_splinesurf`: Bug fixing incorrect varargout

## [Unreleased - 2025-01-15]
### Changed
* `FEMA_convert_splinesurf`: Bug fix and ensuring common variable naming as createBasisFunctions

## [Unreleased - 2025-01-14]
### Added
* `FEMA_convert_splinesurf`

## [Unreleased - 2025-01-13]
### Changed
* `FEMA_fit`: 
    - Minor bug fixes to ensure all variables are defined when `niter` is `0`
    - Additionally computing `coeffCovar` for OLS solution
* `FEMA_run_on_synthetic_data`: additionally have `coeffCovar` as an output variable
* `FEMA_wrapper`: ensuring `coeffCovar` is converted to full dimensions of mask 

## [Unreleased - 2025-01-09]
### Changed
* `FEMA_convert_splinevols`: minor bug fix in indexing
* Major updates to `extract_roi_val` with improved documentation

## [Unreleased - 2025-01-08]
### Changed
* `makeDesign.R`: updates to support ABCD 6.0 variable names

## [Unreleased - 2025-01-07]
### Changed
* `createBasisFunctions`:
    - Removed redundant `dfFlag` output in settings
    - Some settings are only output when `doR` is `true`
    - Removed checking the number of outputs
    - Fixed recording time taken for calculating rank
    - Explicitly output `Xvars` and removed from `settings`
    - Ensuring `Xvars` is a column vector; appropriate fix when calculating number of columns to drop
* `FEMA_convert_splinevols`:
    - Consistent variable names with `createBasisFunctions`
    - Avoiding growing arrays within loop

## [Unreleased - 2024-12-18]
### Changed
* Introduced `varnorm` flag which can be used to perform variance normalization (zero mean and unit standard deviation); appropriate changes in `FEMA_process_data` and `FEMA_wrapper`

## [Unreleased - 2024-11-28]
### Changed
* `FEMA_WaldTest` displays a warning if p values are smaller than 2.2251e-308 and truncates the log10 values to 2.2251e-308 (instead of returning Inf)

## [Unreleased - 2024-11-26]
### Changed
* `FEMA_WaldTest` additionally returns log 10 p values

## [Unreleased - 2024-11-08]
### Changed
* `createBasisFunctions`:
    - documentation update
	- added support for creating `nsk` splines in MATLAB
	- renamed `nsk` functionalty to `nsk-R`
	- updated methodology where linearly spaced vector of values are used for creating bases
	- replaced `Xvars` with `Xpowers` which are powers of the linearly spaced vector of values
	- user can control min, max, and number of linearly spaced values by providing `minMax`
	- renamed `default` `method` to `raw` to avoid confusion
	- removed `ageSubset` and `addConst` features
	- removed the option of providing degrees of freedom: removed `dfFlag`
	- default spline is nsk created in MATLAB
	- renamed `age` to `valvec`
	- using `svd` instead of `orth`
* `createBasis.R`:
	- documentation update, now consistent with MATLAB version
	- default `splineType` is `nsk`
	- removed support for `df` specification
	- dropping of columns functionality removed
	- deprecated `get_basis_values` and `extract_basis` functions
	- overall consistent with `createBasisFunctions` - same set of changes incorporated
* `caller_createBasis.R`:
	- updated to be consistent with changes to `createBasis.R` and `createBasisFunctions`
	- reordered inputs
	- removed `df` flag
	
## [Unreleased - 2024-10-31]
### Changed
* `README.MD`: fixed broken links

## [Unreleased - 2024-10-30]
### Added
* Guide for non-ABCD data

### Changed
* `makeDesign.R`: removed duplicate age column

## [Unreleased - 2024-10-29]
### Changed
* `FEMA_intersect_design`, `FEMA_process_data`, and `makeDesign.R`: support for ABCD 6.0

## [Unreleased - 2024-09-30]
### Changed
* ``FEMA_fit``: `niter` defaults to 1

## [Unreleased - 2024-08-21]
### Changed
* ``createBasisFunctions``:
    - removed ``regress`` and ``demean`` methods
    - removed functionality of dropping basis functions
    - can regress out the effect of multiple variables from the basis functions
    - updated ``svd`` functionality with scaling of the orthonormal basis
    - ``addConst`` is ``false`` by default

## [Unreleased - 2024-07-29]
### Changed
* ``FEMA_fit_GWAS``:
    - fixed standard errors in OLS estimation by calculating SNP-wise mean squared errors - now almost exactly matches a call to standard linear model
    - fixed standard errors in GLS estimation by calculating SNP-wise mean squared errors - similar to OLS solution
    - dropped requirement for ``sig2tvec`` as MSE is being calculated for every SNP
    - fixed incorrect variable name check when calculating size of variables to be saved
    - re-written parts of OLS estimation - no requirement of bins for OLS
    - added description of Wald statistics and p values as output

## [Unreleased - 2024-07-24]
### Changed
* ``FEMA_parse_PLINK``:
    - returns ``genInfo`` with ``locIID``, ``locSNPs``, and ``numSubjs`` as fields
    - accepts ``genInfo`` as an input and skips checking of data
    - output ``benchmark`` renamed to ``tInfo``
    - output ``basePair`` renamed to ``BP``
    - additionally calculates overall timing information
* ``divideSNPs``: 
    - documentation update
    - ``bFile`` is mandatory input now; the file is not parsed if ``Chr``, ``SNPID``, ``BP``, and ``genInfo`` are provided as additional inputs
    - additional fields are now output for each part: ``fname``, ``Locs``, ``Chr``, ``SNPID``, ``BP``, ``outName``, and ``genInfo``

## [Unreleased - 2024-07-11]
### Added
* ``FEMA_fit_GWAS`` function for performing GWAS using FEMA

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

## Miscellaneous
* Included license information for various external code

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
