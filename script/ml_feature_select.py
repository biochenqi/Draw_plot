#!/usr/bin/python
# -*- coding: UTF-8 -*-

# **********************************************************
# * Author        : chenqi
# * Email         : chenqi@novogene.com
# * Create time   : 2023-11-03 14:09
# * Last modified : 2023-11-03 14:09
# * Filename      : ml.py
# * Description   : 
# **********************************************************
from sklearn.ensemble import RandomForestClassifier
import argparse,warnings,sys,os
warnings.filterwarnings('ignore')
#数据处理及分类
import pandas as pd
import numpy as np
from sklearn.svm import LinearSVC,SVC
from sklearn.feature_selection import RFE
from sklearn.metrics import roc_curve,auc,roc_auc_score
from sklearn.model_selection import LeaveOneOut
import gseapy as gp
from sklearn.utils import resample
from gseapy import barplot, dotplot
import matplotlib.pyplot as plt
import seaborn as sns 
from glob import glob

usages = """ 
Author : chenqi
Email  : chenqi@haplox.com
Date   : 2023/11/01
Version: v1.0
Description:
    通过随机森林模型筛选特征
Example:
    python3 %s -i <matrix> -g <group> 

Example path:

Example cmd:
    python3 %s -i /haplox/haprs/danny/project/16s/20230428_HGC20230320002-0001_16s/hapyun_output/qiime2_asv_annot_node_8/asv/otu_norm_tax.xls  -g /haplox/haprs/danny/project/16s/20230428_HGC20230320002-0001_16s/meta.tsv

"""%(sys.argv[0],sys.argv[0])

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def getopt():
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter, description=usages)
    parser.add_argument('-o','--output',help = 'output dir',type = str,default='.',dest='outdir')
    # parser.add_argument('-p','--prefix',help = 'prefix of output file',type = str,default='result',dest='prefix')
    parser.add_argument('-i','--input',help='input matrix file',type=str,dest='input',required=True)
    parser.add_argument('-m','--meta',help='input meta file',type=str,dest='meta',required=True)
    parser.add_argument('-c','--circle',help='num of circle ',type=int,dest='circle',default=20)
    parser.add_argument('--feature_left',help='final left num of feature in each circle ',type=int,dest='feature_left',default=50)
    parser.add_argument('--feature_step',help='each time drop num of feature in each circle ',type=int,dest='feature_step',default=500)
    parser.add_argument('--model',help='chose model to select feature', dest='model', choices=['svc','rf'],default='rf')
    args = parser.parse_args()
    if not args.input:
        print('matrix file must be given!!')
        sys.exit(0)
    elif not args.meta:
        print('group file must be given!!')
    return args

args = getopt()

#绘制roc曲线图
def bootstrap_auc(y_true, y_pred, n_bootstraps=1000):
    bootstrap_auc_scores = []
    for i in range(n_bootstraps):
        # 创建一个自助样本
        boot_y_true, boot_y_score = resample(y_true, y_pred)
        
        # 计算并存储AUC分数
        if len(np.unique(boot_y_true)) > 1:
            # 计算并存储AUC分数
            auc = roc_auc_score(boot_y_true, boot_y_score)
            bootstrap_auc_scores.append(auc)

    # 计算95%置信区间
    lower_bound = np.percentile(bootstrap_auc_scores, 2.5)
    upper_bound = np.percentile(bootstrap_auc_scores, 97.5)
    return lower_bound, upper_bound

def draw_roc(y,y_score,y_predict,out_dir):
    fpr,tpr,threshold = roc_curve(y,y_score)
    lower_bound, upper_bound = bootstrap_auc(y,y_predict)
    roc_auc = auc(fpr,tpr)
    print(caculate_acc_sens_spc(y,y_predict),roc_auc)
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f(CI:%.4f-%.4f))' % (roc_auc,lower_bound,upper_bound))
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic example')
    plt.legend(loc="lower right")
    plt.savefig('%s/roc.png'%out_dir,bbox_inches = 'tight')

def one_eva(x,y,outdir):
    model = LinearSVC(dual=True)
    loo = LeaveOneOut()
    y_predict,y_score = [],[]
    for train,test in loo.split(x):
        model.fit(x.iloc[train,:],y.iloc[train])
        y_predict.append(model.predict(x.iloc[test,:]))
        y_score.append(model.decision_function(x.iloc[test,:]))
    print("leave one out method:",x.shape)
    draw_roc(y,y_score,y_predict,outdir)

#计算准确性、特异性、灵敏度
def caculate_acc_sens_spc(y,y_p):
    tn, fp, fn, tp = 0, 0, 0, 0
    for i in range(len(y)):
        s, s_p = y[i], y_p[i]
        if s == 0 and s_p == 0:
            tn += 1
        elif s == 0 and s_p == 1:
            fp += 1
        elif s == 1 and s_p ==0:
            fn += 1
        elif s == 1 and s_p == 1:
            tp += 1
    # print(tn, fp, fn, tp)
    #           accuracy            sensitive       Specificity
    return (tn+tp)/(fn+fp+tn+tp), tp / (tp + fn), tn / (tn + fp)

def gini_get(x,y,outdir):
    model = RandomForestClassifier()
    model.fit(x,y)
    imp = model.feature_importances_
    x_out = pd.DataFrame({'gini':imp},index=x.columns)
    x_out = x_out.sort_values('gini',ascending=False)
    x_out.to_csv('%s/feature.txt'%outdir,sep='\t')
    return x_out

def diff_group(x,y):
    x_normal = x.loc[y==0,:].mean()
    x_tumor = x.loc[y==1,:].mean()
    return x_normal/x_tumor

def feature_chose(x,y,outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    x = x.transpose()
    columns = x.columns.values
    # y['meta'] = [0 if 'N' in i else 1 for i in y['group'].values]
    model = RandomForestClassifier() if args.model=='rf' else SVC(kernel="linear")
    feature_all = []
    for i in range(args.circle):
        selector = RFE(model,n_features_to_select=args.feature_left,step=args.feature_step)
        selector = selector.fit(x,y['meta'])
        feature = list(columns[np.where(selector.support_==True)[0]])
        if not feature_all:
            feature_all = feature
        else:
            feature_all = list(set(feature_all)&set(feature))
    df_bubble = pd.DataFrame()
    x = x.loc[:,feature_all]
    one_eva(x,y['meta'],outdir)
    df_bubble['gini'] = gini_get(x,y['meta'],outdir)
    df_bubble['time'] = diff_group(x,y['meta'])
    draw_pic(df_bubble,outdir)
    # try:
    #     enr = gp.enrichr(gene_list=feature_all, # or "./tests/data/gene_list.txt",
    #                 gene_sets=['MSigDB_Hallmark_2020','KEGG_2021_Human','GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023'],
    #                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
    #                 outdir=outdir, # don't write to disk
    #                 verbose=False #取消显示Log日志
    #                 )
    #     draw_plot(enr,outdir)
    # except:
    #     ##删除掉所有的log日志文件
    #     for i in glob('%s/*log'%outdir):
    #         os.remove(i)

#绘制功能富集图
def draw_plot(enr,outdir):
    ax = dotplot(enr.results,
              column="Adjusted P-value",
              x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
              size=10,
              top_term=5,
              figsize=(6,12),
              title = "KEGG",
              xticklabels_rot=45, # rotate xtick labels
              show_ring=True, # set to False to revmove outer ring
              marker='o',
              ofname='%s/dot.pdf'%outdir
             )
    ax = barplot(enr.results,
              column="Adjusted P-value",
              group='Gene_set', # set group, so you could do a multi-sample/library comparsion
              size=10,
              top_term=5,
              figsize=(6,12),
              color=['darkred', 'darkblue','yellow','green','purple'], # set colors for group
              ofname='%s/bar.pdf'%outdir
             )

#绘制气泡图
def draw_pic(df_mda,outdir):
    plt.figure(figsize=(6,10))
    g = sns.scatterplot( 
                    x=df_mda['gini'], 
                    y=df_mda.index, 
                    size=df_mda['time'], 
                    # hue = df_mda['Group'],
                    # hue_order=list(dict_y.keys()),
                    alpha=1,
                    zorder=10)
    g.set(xlabel="Mean Decrease Gini", ylabel="Most important feature")
    # 获取坐标轴对象
    ax = plt.gca()
    ax.grid(True)
    # 取消上边框和右边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend(bbox_to_anchor=(1.4, 1), loc=1, borderaxespad=0.)
    plt.savefig('%s/bubble.png'%outdir,bbox_inches = 'tight', dpi =300 )
    plt.savefig('%s/bubble.pdf'%outdir,bbox_inches = 'tight')

df = pd.read_table(args.input,index_col=0)

y_info = pd.read_table(args.meta,index_col=0)

df_info = df.loc[:,y_info.index]

y_cd45 = y_info[y_info['type']=='CD45']
df_cd45 = df.loc[:,y_cd45.index]
y_panck = y_info[y_info['type']=='panck']
df_panck = df.loc[:,y_panck.index]
feature_chose(df_cd45,y_cd45,'%s/cd_45'%args.outdir)
feature_chose(df_panck,y_panck,'%s/panck'%args.outdir)
feature_chose(df_info,y_info,'%s/total'%args.outdir)

