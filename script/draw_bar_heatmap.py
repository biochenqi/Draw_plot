import sys
import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import combinations
from scipy.stats import ttest_ind
from statannotations.Annotator import Annotator
from statannot import add_stat_annotation
#########################################后期需要改进的地方：将输入参数化并且可以提供分组类别，例如输入 --group N_1_IM,P_1_IM,P_2_CA 就只对这三个类别进行箱线图的比较
df = pd.read_table(sys.argv[1],sep='\t',index_col=0) 
group = pd.read_table(sys.argv[2],sep='\t',index_col=0)
#对group所有种类设置颜色类型
uniq_g = group['group'].unique()
colors = np.random.choice(list(mcolors.CSS4_COLORS.keys()), len(uniq_g ), replace=False)
color_dict = dict(zip(uniq_g, colors))
# df = df.drop('B lineage')#用于去除mcp结果中的B细胞类型  
out_dir = sys.argv[3] if len(sys.argv) ==4 else '.'
#######热土需要更改 样本聚类，基因不需要聚类############################已经更改
def draw_heatmap(df,group,prefix=""):
    # df = df.applymap(np.log10) ###目前仅用于MCPcounter的结果
    ##按照group对矩阵进行排序
    group = group.sort_values('group')
    df = df.reindex(columns = group.index.to_list())
    ##重新对df的index进行命名
    group['name'] = group.index.astype('str')+"-"+group['group2'].astype('str') + "-" + group['patient'].astype('str')
    df.columns = df.columns.map(group['name'])
    group.set_index('name',inplace=True)
    #重置画布
    plt.figure()
    col_name = pd.DataFrame({'ROIs':[color_dict[group.loc[i,'group']] for i in df.columns]},index=df.columns)
    ###当使用col_cluster=True聚类算法之后列的顺序就会被打乱
    g = sns.clustermap(df,col_cluster=True,figsize=(10, 7),
                    col_colors=col_name,
                    #设置bar的刻度以及将竖直的转换成水平的bar                                   yticklabels=True用来显示所有的行
                    # cbar_kws=dict(ticks=[0, 0.25, 0.50, 0.75,1],orientation='horizontal'),
                    cmap="Reds", vmin=0, vmax=df.max().max(),
                    yticklabels=True,
                    row_cluster=False)
    
    group_order = group['group']
    group_order = group_order.groupby(group_order).size()
    for i in group_order.index.to_list():
        g.ax_col_dendrogram.bar(1,0,color=color_dict[i],label=i,linewidth=0)
    g.ax_col_dendrogram.legend(loc='lower right', bbox_to_anchor=(1.6, -4), ncol=1,title='ROIs')
    #去除掉热图的x坐标轴标签
    # g.ax_heatmap.set_xticklabels([])
    # g.ax_heatmap.set_xlabel('sample')
    #去除掉热图的y坐标轴标签
    # g.ax_heatmap.set_yticklabels([])
    # g.ax_heatmap.set_ylabel('feature')
    #去掉热图的坐标轴的轴线
    g.ax_heatmap.tick_params(axis='both', which='both', length=0)
    #去除掉col_colors的轴线
    g.ax_col_colors.tick_params(axis='both', which='both', length=0)
    #设置col_colors标签的位置
    g.ax_col_colors.set_yticks([1.5, 1.5])
    #设置col_colors标签的名称
    # g.ax_col_colors.set_yticklabels(['ROIs'])
    # g.ax_col_colors.set_ylabel('Log10')
    #设置bar的位置大小
    # g.ax_cbar.set_position([0.32, 0.83, g.ax_row_dendrogram.get_position().width, 0.02])
    #去除掉bar的轴线
    # g.ax_cbar.tick_params(axis='both', which='both', length=0)
    g.ax_cbar.set_title('Log10')
    #文件保存
    plt.savefig('%s/%sHeatmap.png'%(out_dir,prefix),dpi=100,bbox_inches='tight')
    plt.savefig('%s/%sHeatmap.pdf'%(out_dir,prefix),dpi=100,bbox_inches='tight')

##需要给每个图片加一个title用于区分 CD45和panck########################已经更改
def draw_bowplot(df,list_type,prefix):
    #固定颜色
    dict_color = {'N_1_IM':'g','P_1_IM':'b','P_2_CA':'r'}
    colors = [dict_color[i] for i in list_type]
    df = df[df['group2'].isin(list_type)]
    #重置画布
    plt.figure(figsize=(15,8))
    sns.set_theme(style="ticks", palette="pastel")

    ax = sns.boxplot(x="cell_type", y="value",
            hue="group2", 
            palette=colors,
            hue_order = list_type,
            data=df)
    sns.despine(offset=10, trim=True)
    #设置X和Y轴label
    plt.xlabel('Cell type')
    plt.ylabel('Relative Abundance')
    plt.xticks(rotation=90)
    x_labels = df['cell_type'].unique()
    hue_labels = df['group2'].unique()
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    for i, cell_type in enumerate(x_labels):
        group1 = df[(df['cell_type'] == cell_type) & (df['group2'] == hue_labels[0])]['value']
        group2 = df[(df['cell_type'] == cell_type) & (df['group2'] == hue_labels[1])]['value']
        t, p = ttest_ind(group1, group2)
        
        #添加注释
        y = max(max(group1), max(group2))
        plt.plot([i-0.2, i+0.2], [y+3, y+3], color='k')
        plt.text(i-0.2, y+4, '{:.3f}'.format(p),fontsize=6,color='red' if p<0.05 else 'black')
        # 添加注释到图形中
        # if p < 0.05:
        #     ax.annotate(f"{p:.3f}", xy=(i, 0), xycoords='data', ha='center',fontsize=6,color='red')
        # else:
        #     ax.annotate(f"{p:.3f}", xy=(i, 0), xycoords='data', ha='center',fontsize=6,color='black')
    title_name = prefix.split('_')[0]
    title_name = title_name.capitalize() if title_name == "panck" else title_name
    plt.title(title_name)
    plt.savefig('%s/%s_boxplot.png'%(out_dir,prefix),dpi=300,bbox_inches='tight')
    plt.savefig('%s/%s_boxplot.pdf'%(out_dir,prefix),dpi=300,bbox_inches='tight')

##此处通过注释标记表明显著差异情况 其中 'star': 只显示星号，表示 p 值的大小 (****: p <= 0.0001, ***: p <= 0.001, **: p <= 0.01, *: p <= 0.05)。
def draw_bowplot_alltype(df,prefix):
    #重置画布
    plt.figure(figsize=(18,8))
    sns.set_theme(style="ticks", palette="pastel")
    hue_order=['N_1_IM','P_1_IM','P_2_CA']
    pairs = []
    for i in combinations(hue_order,2):
        i = list(i)
        pairs = pairs + [((clar, i[0]), (clar, i[1]))
                            for clar in df['cell_type'].unique()]
    ax = sns.boxplot(data=df, x="cell_type", y="value", hue="group2", palette=["g", "b",'r'],hue_order = hue_order)
    add_stat_annotation(ax, data=df, x="cell_type", y="value", hue="group2",
                    box_pairs=pairs, 
                    test='t-test_ind', text_format='star',
                    # test="Mann-Whitney", text_format='star',
                    loc='inside', verbose=2,
                    comparisons_correction=None, #此处不进行任何矫正，默认为bonferroni矫正
                    hue_order = hue_order)
    sns.despine(offset=10, trim=True)
    #设置X和Y轴label
    plt.xlabel('Cell type')
    plt.ylabel('Relative Abundance')
    plt.xticks(rotation=90)
    prefix = prefix.capitalize() if prefix == "panck" else prefix
    plt.title(prefix)
    plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
    plt.savefig('%s/%s_boxplot_all.png'%(out_dir,prefix),dpi=300,bbox_inches='tight')
    plt.savefig('%s/%s_boxplot_all.pdf'%(out_dir,prefix),dpi=300,bbox_inches='tight')

def draw_boxplot_all(df,group,types):
    list_info = ['N_1_IM','P_1_IM','P_2_CA']
    groups = group[group['type']==types]
    groups = groups[groups['group2'].isin(list_info)]
    dfs = df.loc[:,groups.index]
    dfs = dfs.transpose()
    dfs = dfs.reset_index().melt(id_vars='tags',var_name='cell_type',value_name='value')
    dfs.set_index('tags',inplace=True)
    # 在tag上合并dfs和groups 将group2复制到dfs
    dfs= dfs.merge(groups, left_index=True, right_index=True)
    for i in combinations(list_info,2):
        i = list(i)
        draw_bowplot(dfs,i,types+"_" + '_'.join(i))
    draw_bowplot_alltype(dfs,types)

def survey(results, category_names,ax):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    # category_colors = plt.colormaps['RdYlGn'](
    #     np.linspace(0.15, 0.85, data.shape[1]))
    category_colors = color_dict.values()
    # fig, ax = plt.subplots(figsize=(9.2, 5))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.01,
                        label=colname, color=color)

        # r, g, b, _ = color
        # text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
        # ax.bar_label(rects, label_type='center', color=text_color)
    ax.legend(ncols=5, bbox_to_anchor=(0, 1),
              loc='lower left', fontsize='small')

    return  ax

def draw_bar(df,group):
    ##重新对df的index进行命名
    group['name'] = group.index.astype('str')+"-"+group['group2'].astype('str') + "-" + group['patient'].astype('str')
    group = group.sort_values('group')
    df.columns = df.columns.map(group['name'])
    df = df.transpose()
    df = df.reindex(group['name'])  #将df index的顺序变得和group的 name一致
    df = df.div(df.sum(axis=1),axis=0)
    ##按照原本的顺序排列计算每个组名的数量
    group_order = group['group']
    group_order = group_order.groupby(group_order).size()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [1, 24]})
    plt.subplots_adjust(hspace = 0.05)
    # 在上方的子图中画出分组信息
    # group['group2'].value_counts().plot(kind='bar', ax=ax1)
    # ax1.set_xticklabels([])  # Remove x-tick labels
    # ax1.set_xlabel('')
    ax = survey({'ROIs':list(group_order.values)}, group_order.index.to_list(),ax1)
    df.plot(kind='bar',stacked=True,ax=ax2)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('')
    ax2.set_ylim([0, 1.0])
    plt.xticks(fontsize=3)
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(reversed(handles), reversed(labels),bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.,fontsize=6,title='Cell_type')
    plt.savefig('%s/cell_ratio.png'%out_dir,dpi=300,bbox_inches='tight')
    plt.savefig('%s/cell_ratio.pdf'%out_dir,dpi=300,bbox_inches='tight')

##分离CD45和panck的结果
def split_info(df,group,types):
    groups = group[group['type']==types]
    dfs = df.loc[:,groups.index]
    return dfs,groups

draw_heatmap(df,group)
dfs,groups = split_info(df,group,"CD45")
draw_heatmap(dfs,groups,'CD45')
dfs,groups = split_info(df,group,"panck")
draw_heatmap(dfs,groups,'Panck')

# draw_boxplot_all(df,group,'CD45')
# draw_boxplot_all(df,group,'panck')
# draw_bar(df,group)