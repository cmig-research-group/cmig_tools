#!/bin/bash

MCRROOT="/usr/pubsw/packages/matlab/R2023a"

designid="des000009"
#fstem_imaging="JA"
fstem_imaging="RND"
datatype="voxel"

#abcd_data_path="/space/syn65/1/data/abcd-sync/4.0"
abcd_data_path="/space/syn65/1/data/abcd-sync/5.0"
#fname_design="$abcd_data_path/support_files/design_matrices/ABCD_rel4.0_cross_desmat_PCs_SES_motion_interview_age.txt"
fname_design="~/design.txt"
#dirname_tabulated="$abcd_data_path/tabulated/released"
dirname_tabulated="$abcd_data_path/tabulated/img"
#dirname_imaging="$abcd_data_path/imaging_concat/voxelwise/ABCD2_cor10/smri"
dirname_imaging="$abcd_data_path/imaging_concat/voxelwise/dmri"
dirname_out='./data_v5.0/output'
dirname_cache='./data_v5.0/FEMA_cache'
dirname_jobs='./data_v5.0/jobs'
nfrac="10"

mkdir -p $dirname_jobs

echo "Running run_FEMA_DEAP_wrapper.sh $MCRROOT $fstem_imaging $fname_design $dirname_out $dirname_tabulated $dirname_imaging $datatype $dirname_cache $dirname_jobs $designid $nfrac"

./FEMA_DEAP_wrapper/run_FEMA_DEAP_wrapper.sh $MCRROOT $fstem_imaging $fname_design $dirname_out $dirname_tabulated $dirname_imaging $datatype $dirname_cache $dirname_jobs $designid $nfrac
