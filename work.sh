#!/bin/bash
# **********************************************************
# * Author        : chenqi
# * Email         : chenqi@gooalgene.com
# * Create time   : 2023-12-06 10:31
# * Last modified : 2023-12-06 10:31
# * Filename      : work.sh
# * Description   : 
# **********************************************************
python=/haplox/haprs/chenqi/software/anaconda3/bin/python
#通过机器学习模型筛选重要的特征
$python script/ml_feature_select.py -i 00.data/q_norm.txt -m 01.Machine_learning/p1im_n1im/new_meta.csv -o 01.Machine_learning/p1im_n1im/
pwd=`pwd`
###对特征基因进行功能富集分析
ls 01.Machine_learning/p1im_n1im/*/feature.txt |while read line;
do
    prefix=`echo ${line}|awk -F '/' '{print $(NF-1)}' `
    echo $line,$prefix
    sh script/cluster_profile.sh ${pwd}/${line} ${prefix} 01.Machine_learning/p1im_n1im/${prefix}/cluster_profile/
done
#通过SpatialDecon进行细胞类型鉴定
/haplox/haprs/zhaoshuhua/software/miniconda3/envs/Infiltration/bin/Rscript script/spatialDecon_ana.r 00.data/q_norm.txt 00.data/group.csv 02.Cell_type/SpatialDecon/
#通过细胞类型鉴定结果绘制图像
$python script/draw_bar_heatmap.py 02.Cell_type/SpatialDecon/cell_type.xls 00.data/group.csv 02.Cell_type/SpatialDecon/png_result/