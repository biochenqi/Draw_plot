library(SpatialDecon)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

df <- read.table(args[1],header=TRUE)

group <- read.table(args[2],header=TRUE)

boao_dea <- function(file){
  boao <- read.table(file,header=TRUE)
  rownames(boao) <-boao$Gene
  boao <- boao[,-1]
  colnames(boao)= gsub("^PANCK","panck",colnames(boao),ignore.case=TRUE)
}
deal_matrix <- function(df,group){
##只分析CD45的tumor数据,现在分析所有的样本
# group <- group[group['type']=='CD45' & group['meta']==1,]
# group <- group[group['type']=='CD45',]
df <- df[,colnames(df)%in% group$tags]
df <- as.matrix(df)
new_row <- data.frame(matrix(data = 1, nrow = 1, ncol = ncol(df)))
rownames(new_row) <- "NegProbe"##增加一行NegProbe用于后续生成背景矩阵
colnames(new_row) <- colnames(df)
df <- rbind(df, new_row)

df<- as.matrix(df)

bg <- derive_GeoMx_background(
       norm = df,
       probepool = rep(1, nrow(df)),
       negnames = "NegProbe"
     )
##运行计算spatialDecon用来计算各种细胞类型的丰度情况
res0 <- spatialdecon(
       norm = df,
       bg = bg,
       X = safeTME
     )
return(res0)
}

res = deal_matrix(df,group)
write.table(res$beta,paste0(args[3],"/cell_type.xls"),quote = FALSE, sep = "\t")


draw_bar_col <- function(matrix,prefix){
    df_long <- matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("cell_type") %>%
  gather(key = "sample", value = "count", -cell_type)

# 计算每个样本中每种细胞类型的相对丰度比例
df_long <- df_long %>%
  group_by(sample) %>%
  mutate(relative_count = count / sum(count) * 100)

df_long <- df_long %>%
  left_join(group, by = c("sample" = "tags")) %>%
  rename(sample_type = group)

p2 = ggplot(df_long , aes(x=factor(sample,levels=unique(newrank$sample)),y=1,fill=sample_type))+geom_col(width=1)+theme_bw()+theme(legend.position="none",axis.text.x=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),panel.border=element_blank())

}

draw_bar <- function(matrix,prefix){
    data <- matrix

# 转置数据框并将其转化为长格式
data_long <- melt(t(data))

# 重命名列
colnames(data_long) <- c("Sample", "CellType", "Count")

# 计算每个样本的细胞类型的百分比
data_long$Percentage <- with(data_long, ave(Count, Sample, FUN = function(x) 100*x/sum(x)))

# 制作堆叠的百分比柱状图
p <- ggplot(data_long, aes(x = Sample, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y = "Percentage", fill = "Cell Type", title = "Cell type proportion per sample")
ggsave(paste0(args[3],"/",prefix,"_bar.png"),p,width=7,height=10,dpi=400)
}

draw_heatmap <- function(a,prefix,size,name) {
write.table(a,paste0(args[3],"/",prefix,".csv"),quote = FALSE, sep = "\t")
data_long <- melt(a)
p <- ggplot(data_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12,angle = 90),axis.text.y = element_text(size=size))+
  labs(x = "", y = "", fill = name)
ggsave(paste0(args[3],"/",prefix,".png"),p,width=7,height=10,dpi=400)
}

draw_residual <- function(a,prefix,size,name) {
write.table(a,paste0(args[3],"/",prefix,".csv"),quote = FALSE, sep = "\t")
data_long <- melt(a)
p <- ggplot(data_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0, limit = c(-3, 3), space = "Lab", 
                       na.value = "grey50", guide = "colourbar") +
  theme_minimal() +
  theme(axis.text.x = element_text(size=12,angle = 90),axis.text.y = element_text(size=size))+
  labs(x = "", y = "", fill = name)
ggsave(paste0(args[3],"/",prefix,".png"),p,width=7,height=10,dpi=400)
}

#绘制柱状图占比
# draw_bar(res0$beta,"cell_abundence")
##绘制热图
# draw_heatmap(res0$beta,"cell_abundence",12,"Abundence")
# draw_heatmap(res0$X,"gene_cell_abundence",3,"Abundence")
# draw_residual(res0$resids,"gene_sample_abundence",3,"Residual")
