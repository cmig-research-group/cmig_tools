# Fast Efficient Mixed effects Analysis (FEMA)

This package includes code to run mass univariate linear mixed effects analysis in a fast and efficient manner.  This can be done at the whole brain level using vertexwise and voxelwise data, and connectivity matrices.  FEMA includes an external data option (.txt format) which can include columms of imaging or behavioral data.  The code includes an ABCD Study specific wrapper ([`FEMA_wrapper.m`](FEMA/FEMA_wrapper.m)), which uses the concatenated imaging files saved on abcd-sync.  The internal function(s) within the wrapper can be used to analyse any data as long as its inputs are in the correct format within MATLAB.  Specifically, ([`FEMA_fit.m`](FEMA/FEMA_fit.m)) is the core function that performs model fitting (see [`FEMA_Guide.md`](Guide_FEMA.md) for more details on how to format the data and use FEMA_fit outside of the wrapper). The code was developed and tested using MATLAB 2020a / 2023a.

Technical details on how FEMA works can be found [here](https://doi.org/10.1002/hbm.26579) and [here](https://doi.org/10.1101/2025.05.09.653146)

**Reference**: 

Parekh, P., Fan, C. C., Frei, O., Palmer, C. E., Smith, D. M., Makowski, C., Iversen, J. R., Pecheva, D., Holland, D., Loughnan, R., Nedelec, P., Thompson, W. K., Hagler, D. J. Jr, Andreassen, O. A., Jernigan, T. L., Nichols, T. E., & Dale, A. M. (2024). FEMA: Fast and efficient mixed-effects algorithm for large sample whole-brain imaging data. Human Brain Mapping, 45(2), e26579. https://doi.org/10.1002/hbm.26579 

If you use splines, unstructured covariance, or the GWAS module, please additionally cite:

Parekh, P., Parker, N., Pecheva, D., Frei, E., Vaudel, M., Smith, D.M., Rigby, A., Jahołkowski, P., Sønderby, I.E., Birkenæs, V., Bakken, N.R., Fan, C.C., Makowski, C., Kopal, J., Loughnan, R., Hagler, D.J., Meer, D. van der, Johansson, S., Njølstad, P.R., Jernigan, T.L., Thompson, W.K., Frei, O., Shadrin, A.A., Nichols, T.E., Andreassen, O.A., Dale, A.M., 2025. FEMA-Long: Modeling unstructured covariances for discovery of time-dependent effects in large-scale longitudinal datasets. bioRxiv. https://doi.org/10.1101/2025.05.09.653146


Please see the see the tutorials in the [`CMIG Tools Wiki`](https://github.com/cmig-research-group/cmig_tools/wiki) for more details on how to use FEMA.

* [Installation](#installation)
  + [Data](#data)
  + [Code dependencies](#code-dependencies)
  + [Local setup](#local-setup)
* [Usage](#usage)
* [Reference](#reference)
* [Copyright](#copyright)
* [Contributing](#contributing)

## Installation
### Data
Whole-brain or surface-based imaging data should already by pre-processed and concatenated across observations in the same atlas space into a 2D matrix in MATLAB that can be input to `FEMA_fit` for analysis (`ymat`).  These data (`ymat`) should be *e.g.* nobservations x nvoxels.

For ABCD Investigators - To run FEMA using the ABCD specific `FEMA_wrapper` we advise that you mirror `abcd-sync` on your own server maintaining the data structure.  You can then provide `FEMA_wrapper` with the appropriate paths to run your analysis.  Abcd-sync is still being developed, so make sure to sync regularly in order to get access to all of the latest imaging data.  Not all imaging modalities are available yet but will be uploaded in the coming weeks.

### Code dependencies
Necessary directories to add to MATLAB path, all part of this repository:
- [`FEMA`](FEMA)
- [`cmig_tools_utils/matlab`](cmig_tools_utils/matlab)
- [`PALM`](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM) (only required if you plan to use PALM for permutations)

To visualize results:
- for voxel data, use [`showVol`](showVol)
- for vertex data, use [`showSurf`](showSurf)

### Local setup
If you want to run the [`FEMA_wrapper_demo.m`](ABCD_DEMOS/FEMA_wrapper_demo.m), run [`abcdConfig.m`](cmig_tools_utils/matlab/abcdConfig.m) to setup your local computer to find the path to the code and data directories (after they have been downloaded).


## Usage
For instructions on how to create design matrices for FEMA, run FEMA analyses and visualise the results see the tutorials in the [`CMIG Tools Wiki`](https://github.com/cmig-research-group/cmig_tools/wiki). 

## Reference
Parekh, P., Fan, C. C., Frei, O., Palmer, C. E., Smith, D. M., Makowski, C., Iversen, J. R., Pecheva, D., Holland, D., Loughnan, R., Nedelec, P., Thompson, W. K., Hagler, D. J. Jr, Andreassen, O. A., Jernigan, T. L., Nichols, T. E., & Dale, A. M. (2024). FEMA: Fast and efficient mixed-effects algorithm for large sample whole-brain imaging data. Human Brain Mapping, 45(2), e26579. https://doi.org/10.1002/hbm.26579

Parekh, P., Parker, N., Pecheva, D., Frei, E., Vaudel, M., Smith, D.M., Rigby, A., Jahołkowski, P., Sønderby, I.E., Birkenæs, V., Bakken, N.R., Fan, C.C., Makowski, C., Kopal, J., Loughnan, R., Hagler, D.J., Meer, D. van der, Johansson, S., Njølstad, P.R., Jernigan, T.L., Thompson, W.K., Frei, O., Shadrin, A.A., Nichols, T.E., Andreassen, O.A., Dale, A.M., 2025. FEMA-Long: Modeling unstructured covariances for discovery of time-dependent effects in large-scale longitudinal datasets. bioRxiv. https://doi.org/10.1101/2025.05.09.653146

## Contributing
Contributions to this work are welcome, thank you for your time and interest! 
- found a bug? thought of an enhancement? submit a detailed issue for review
- if you want to make changes yourself, you can fork this repo, make changes on your forked version, test the changes, and then submit a pull request
