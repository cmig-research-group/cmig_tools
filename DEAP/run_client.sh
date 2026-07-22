#!/bin/bash

# All configurable settings go here
MCRROOT="/home/pparekh/matlabruntime_R2024b_Update_6/R2024b"
workDir="/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel"
abcd_root="/space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed"
designid="des00001"
fstem_imaging="thk_sm16"
datatype="vertex"
fname_design="$workDir/test_designMatrix.txt"
dirname_imaging="$abcd_root/6_0/concatenated/imaging/vertexwise/smri"
dirname_out="$workDir/output"
dirname_cache="$workDir/FEMA_cache"
dirname_jobs="$workDir/jobs"
nfrac="10"
colsinterest="1:20"

mkdir -p $dirname_jobs

echo "Running run_FEMA_DEAP_wrapper.sh $MCRROOT $fstem_imaging $fname_design $dirname_out $dirname_imaging $datatype $dirname_cache $dirname_jobs $designid $nfrac"

sh /space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/compiledCode/run_FEMA_DEAP_wrapper.sh $MCRROOT $fstem_imaging $fname_design $dirname_out $dirname_imaging $datatype $dirname_cache $dirname_jobs $designid $nfrac colsinterest $colsinterest