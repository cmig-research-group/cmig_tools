#!/bin/bash

MCRROOT="/home/pparekh/matlabruntime_R2024b_Update_6/R2024b"
abcd_root="/space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed/"
workDir="/space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel"
fname_design="$workDir/test_designMatrix.txt"
dirname_imaging="$abcd_root/6_0/concatenated/imaging/vertexwise"
dirname_out="$workDir/output"
dirname_cache="$workDir/FEMA_cache"
nfrac="10"

sh /space/ceph/1/ABCD/users/parekh/2026-01-02_DEAP_Parallel/compiledCode/run_FEMA_DEAP_gencache.sh $MCRROOT $nfrac $dirname_cache $dirname_imaging $dirname_out $fname_design