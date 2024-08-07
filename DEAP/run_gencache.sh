#!/bin/bash

#abcd_root='/abcd/abcd-sync'
abcd_root="/space/syn65/1/data/abcd-sync"

MCRROOT="/usr/pubsw/packages/matlab/R2023a"

#fname_design="$abcd_root/4.0/support_files/design_matrices/ABCD_rel4.0_cross_desmat_PCs_SES_motion_interview_age.txt"
#dirname_tabulated="$abcd_root/4.0/tabulated/released"
#dirname_imaging="$abcd_root/4.0/imaging_concat/voxelwise/ABCD2_cor10"
fname_design="~/designmat_v5.txt"
dirname_tabulated="$abcd_root/5.0/tabulated/img"
dirname_imaging="$abcd_root/5.0/imaging_concat/voxelwise"
dirname_out='./data_v5.0/output'
dirname_cache='./data_v5.0/FEMA_cache'
nfrac=10

./FEMA_DEAP_gencache/run_FEMA_DEAP_gencache.sh $MCRROOT $nfrac $dirname_cache $dirname_tabulated $dirname_imaging $dirname_out $fname_design
