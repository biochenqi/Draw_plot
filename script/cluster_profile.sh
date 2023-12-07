#!/bin/bash
# **********************************************************
# * Author        : chenqi
# * Email         : chenqi@gooalgene.com
# * Create time   : 2023-11-10 10:50
# * Last modified : 2023-11-10 10:50
# * Filename      : cluster_profile.sh
# * Description   : 
# **********************************************************
if [ ! -n "$1" ] || [ $1 == "-h" ] || [ $1 == "--help" ]; then
    echo "Usage"
    echo "  sh $0 <file:gene> <prefix> <outdir>"
else
echo start at time `date +%F'  '%H:%M:%S`
/haplox/haprs/zhaoshuhua/software/miniconda3/envs/Infiltration/bin/Rscript /haplox/haprs/chenqi/project/DSP/script/clusterProfiler.test.R -f $1 -n $2 -o $3 -a true -s org.Hs.eg.db,hsa,human -g 1 -t SYMBOL -d  KEGG -C 0.05 
echo finish at time `date +%F'  '%H:%M:%S`
fi

