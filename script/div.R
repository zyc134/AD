rm(list = ls())
library(Maaslin2)
library(MicrobiotaProcess)
library(tidyverse)
library(compositions)
library(pheatmap)
library(phyloseq)
library(ggalluvial)
library(dplyr)
library(qvalue)
library(gridExtra)
library(ggsignif)
library(ggpubr)
library(data.table)
library(ggplot2)
library(ggExtra)
library(ggtext)
library(vegan)
library(caret)
library(data.table)
library(ggtext)
library(ggbreak)
library(phyloseq)
library(ggstar)
library(ggnewscale)
abundance_D<-read.csv('D:/microbiome/AD/data/AD/genus/abundance_genus.csv',row.names = 1)
meta_AD<-read.csv('D:/microbiome/AD/data/meta_AD_cn.csv')
meta_AD$group<-'AD'
abundance_D<-abundance_D[,colnames(abundance_D)%in%meta_AD$Run.ID]

abundance_H<-read.csv('D:/microbiome/AD/data/Healthy/genus/abundance_genus.csv',row.names = 1)
meta_H<-read.csv('D:/microbiome/AD/data/meta_h_cn.csv')
meta_H$group<-'Healthy'
colnames(meta_H)[1]<-colnames(meta_AD)[1]
abundance_H<-abundance_H[,colnames(abundance_H)%in%meta_H$Run.ID]
genus<-intersect(rownames(abundance_D),rownames(abundance_H))
abundance_D<-abundance_D[genus,]
abundance_H<-abundance_H[genus,]
abundance<-cbind(abundance_D,abundance_H)
#write.csv(abundance,'D:/microbiome/AD/data/abundance_cn.csv')
abundance<-read.csv('D:/microbiome/AD/data/abundance_cn.csv',row.names = 1)
meta<-rbind(meta_AD,meta_H)
rownames(meta)<-meta$Run.ID
meta2<-meta[,6,drop=F]
rownames(meta2)<-meta$Run.ID
OTU<-phyloseq::otu_table(abundance,taxa_are_rows = T)
meta1<-phyloseq::sample_data(meta2)
physeq<-phyloseq(OTU,meta1)
mpse<-as.MPSE(physeq)
mpse@assays@data@listData[["Abundance"]]<-mpse@assays@data@listData[["Abundance"]]*1e10
mpse %<>% 
  mp_decostand(.abundance=Abundance)
cols<-c('AD'='#E7B800', 'Healthy'='#00AFBB')
mpse %<>% mp_cal_alpha(.abundance = Abundance,force = T)

alpha_res<-mpse@colData@listData%>%as.data.frame()
alpha_res$ACE<-NULL

# 1. 宽转长，并固定因子顺序
df_long <- alpha_res %>%
  pivot_longer(
    cols = c(Observe, Shannon, Chao1, Simpson, Pielou),
    names_to = "Index",
    values_to = "Value"
  ) %>%
  mutate(
    Index = factor(Index, levels = c("Observe", "Shannon", "Chao1", "Simpson", "Pielou")),
    group = factor(group, levels = c("AD", "Healthy"))
  )

# 2. 按不同指数进行 Wilcoxon 检验，并计算多重检验校正后的 p.adj (默认 Holm 或可改为 BH/FDR)
stat_test <- df_long %>%
  group_by(Index) %>%
  wilcox_test(Value ~ group) %>%
  adjust_pvalue(method = "BH") %>% # 这里使用 BH (FDR) 校正，也可改为 "holm", "bonferroni" 等
  add_significance("p.adj") %>%
  add_xy_position(x = "group", scales = "free_y") %>%
  # 自定义 p.adj 的显示格式（小于 2.22e-16 时按原图显示，否则保留合适小数/科学计数法）
  mutate(
    p_label = ifelse(p.adj < 2.22e-16, "P.adj < 2.22e-16", signif(p.adj, 2))
  )
# 自定义配色（还原原图的黄褐色与青蓝色）
my_colors <- c("AD" = "#E6B800", "Healthy" = "#00A6B4")

alpha <- ggplot(df_long, aes(x = group, y = Value)) +
  # 1. 左侧：半小提琴图 (side = "l")
  geom_half_violin(
    aes(fill = group),
    side = "l",
    position = position_nudge(x = -0.05), # 略微向左移动，留出中间空隙
    color = NA,                           # 无边框
    trim = FALSE,
    alpha = 0.9
  ) +
  # 2. 右侧：带黑色细边框的抖动散点
  geom_jitter(
    aes(fill = group),
    shape = 21,                           # 支持内部填充和外部描边的点形状
    color = "black",                      # 黑色边框
    stroke = 0.4,                         # 边框粗细
    size = 1.6,
    alpha = 0.85,
    position = position_jitterdodge(
      jitter.width = 0.22,                # 散点水平分布宽度
      dodge.width = 0.1
    )
  ) +
  # 3. 右侧：同色系箱线图（叠加在散点上方）
  geom_boxplot(
    aes(color = group),
    width = 0.12,
    position = position_nudge(x = 0.12),  # 向右偏移，与散点对齐
    fill = NA,                            # 透明填充，露出底层散点
    outlier.shape = NA,                   # 隐藏离群点（避免与散点重复）
    linewidth = 0.6
  ) +
  # 4. 添加 p.adj 显著性括号
  stat_pvalue_manual(
    stat_test,
    label = "p_label",                    # 使用我们格式化好的 p.adj 标签
    tip.length = 0.02,                    # 括号两端竖线长度
    bracket.size = 0.5,
    size = 3.8
  ) +
  # 5. 分面与坐标轴设置
  facet_wrap(~ Index, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) + # 顶部留出空间放 P 值括号
  labs(x = NULL, y = "Alpha Index Value") +
  # 6. 主题美化（还原灰色分面标题栏与白底边框）
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),                                      # 去除网格线
    strip.background = element_rect(fill = "#BDBDBD", color = NA),     # 灰色分面标题背景
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "grey30", size = 10),
    axis.text.y = element_text(color = "grey30", size = 9),
    axis.title.y = element_text(size = 12, color = "black"),
    panel.spacing = unit(0.6, "lines")                                 # 分面之间的间距
  )

print(alpha)





mpse %>% 
  mp_anosim(.abundance=Abundance, .group=group, action="get")
mpse %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
#dist <- mpse %>% mp_plot_dist(.distmethod = bray, .group = group, group.test=TRUE, textsize=2)+
#  ggtitle('')+
#  theme(plot.title = element_text(face = 'bold'))
#dist 
#library(patchwork)
#combined_plot <- (alpha | dist) + 
#  plot_layout(widths = c(7, 3))
#combined_plot
pdf('D:/microbiome/AD/1118/alpha.pdf',width = 8,height = 6)
alpha
dev.off()
##########pcoa
dist_bray <- phyloseq::distance(physeq, method = "bray")
pcoa <- ordinate(physeq, method = "PCoA", distance = dist_bray)
pcoa_df <- data.frame(pcoa$vectors)
pcoa_df$SampleID <- rownames(pcoa_df)
sample_data<-cbind(rownames(meta),meta$group)%>%as.data.frame()
colnames(sample_data)<-c('SampleID','disease')
rownames(sample_data)<-sample_data$SampleID
pcoa_df <- merge(pcoa_df, sample_data, by.x = "SampleID", by.y = "SampleID")
pcoa_plot<-ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = disease)) +
  geom_point(alpha = 0.4) +
  labs(
    x = paste0("PCOA1 (", round(pcoa[["values"]][["Eigenvalues"]][1], 2), "%)"),
    y = paste0("PCOA2 (", round(pcoa[["values"]][["Eigenvalues"]][2], 2), "%)"),
    color = "label"
  ) +
  theme_minimal()+
  theme(legend.position = 'none')+
  stat_ellipse(aes(color = disease),type = "t", level = 0.95)+
  scale_color_manual(values=c('AD'='#E7B800', 'Healthy'='#00AFBB'))+
  annotate("text", x = -0.3, y = -0.3, label = "ANOSIM R: 0.467 
           P: 0.001")+
  ggtitle('')+
  theme(plot.title = element_text(face = 'bold'))
pcoa_plot
pcoa_plot<-ggMarginal(pcoa_plot, type = "boxplot", margins = "both", size = 5, groupColour = TRUE, groupFill = TRUE)
pdf('D:/microbiome/AD/1118/div.pdf',width = 8,height = 10)
alpha/pcoa_plot
dev.off()
pdf('D:/microbiome/AD/1118/pcoa.pdf')
pcoa_plot
dev.off()
##########
tree<-readRDS('D:/microbiome/AD/genus_tree.rds')
MicrobiotaProcess::otutree()
mpse %<>%
  mp_diff_analysis(
    .abundance = Abundance,
    .group = group,
    force = T,
    first.test.alpha = 1e-5
  )
old_theme <- theme_get()
theme_update(axis.text.y = element_text(face = "italic"))

f.box <- mpse %>%
  mp_plot_diff_boxplot(.group = group) %>%
  set_diff_boxplot_color(
    values = c("#E7B800", "#00AFBB"),
    guide = guide_legend(title = NULL)
  )

print(f.box)

# 恢复默认主题，避免影响后续其他图

pdf('D:/microbiome/AD/1118/diff_box.pdf',width = 6,height = 10)
print(f.box)
dev.off()
theme_set(old_theme)



res<-cbind(mpse@NAMES,mpse@elementMetadata@listData[["LDAupper"]],
           mpse@elementMetadata@listData[["LDAmean"]],
           mpse@elementMetadata@listData[["LDAlower"]],
           mpse@elementMetadata@listData[["Sign_group"]],
           mpse@elementMetadata@listData[["pvalue"]],
           mpse@elementMetadata@listData[["fdr"]])%>%as.data.frame()
colnames(res)<-c('genus','LDAupper','LDAmean','LDAlower','Sign_group','pvalue','fdr')
res<-na.omit(res)
res[,c(2,3,4,6,7)]<-as.numeric(res[,c(2,3,4,6,7)])
write.csv(res,'D:/microbiome/AD/1118/diff_genus.csv')

grid.arrange(arrangeGrob(arrangeGrob(alpha,pcoa_plot,ncol=1),final_plot,ncol = 2))
########
