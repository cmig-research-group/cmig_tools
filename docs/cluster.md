# Running FEMA on a cluster

Why run on a cluster?
- speed up number of permutations
- run several analyses at once

## Getting started
1. modify [`submit_cluster_do`](../FEMA/submit_cluster_do) to include any setup required for matlab to load on the compute node, and the path to your local PALM setup, *e.g.*:
    ```bash
    module load matlab/2021a
    export MATLABPATH=$MATLABPATH:/path/to/your/local/PALM
    ```
2. on a submission host, run `python3 submit_cluster.py --required_arguments`. Only the base install of Python 3.6 or higher is required, no special virtual environment is needed.

## Supported arguments
Run `python3 submit_cluster.py -h` for a complete list of supported inputs.

Here is a an example call to get started:
```bash
python3 submit_cluster.py -dm /path/to/design_matrix.txt -dt vertex -mod smri -img area -cols "[11,12]" -o /path/to/outputdir
```

## Output directory structure

```bash
# taking the example from above
output_dir
└── nperms-XX_rel-4.0_dt-vertex_mode-smri_img-area
    ├── logs
    │   ├── cluster_submission.log
    │   ├── cluster_gather.oXX
    │   ├── FEMA_fit.oJOBID
    │   ├── FEMA_fit.oJOBID
        ...
    ├── tmp
    │   ├── 1
    │   │   └── permutation-type_nperms / FEMA_wrapper_output_vertex_area_icX_smXX.mat
    │   ├── 2
    │   ├── 3
    │   ...
    │   └── njobs
    └── permutation-type_nperms
        ├── FEMA_wrapper_output_vertex_area_icX_smXX.mat
        ├── FEMA_perm_significance_z_vertex_area_icX_smXX.mat
        ├── FEMA_perm_significance_tfce_vertex_area_icX_smXX.mat
        └── FEMA_mostest_results_vertex_area_icX_smXX.mat
```