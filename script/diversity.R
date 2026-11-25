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
abundance<-read.csv('D:/microbiome/AD/data/abundance_all_filtered.csv',row.names = 1)
meta<-read.csv('D:/microbiome/AD/data/meta_958.csv',row.names = 1)
OTU<-phyloseq::otu_table(abundance,taxa_are_rows = T)
meta1<-phyloseq::sample_data(meta)
physeq<-phyloseq(OTU,meta1)
mpse<-as.MPSE(physeq)
mpse@assays@data@listData[["Abundance"]]<-mpse@assays@data@listData[["Abundance"]]*1e10
mpse %<>% 
  mp_decostand(.abundance=Abundance)
cols<-c('AD'='#E7B800', 'Healthy'='#00AFBB')
mpse %<>% mp_cal_alpha(.abundance = Abundance,force = T)
alpha <- mpse %>% 
  mp_plot_alpha(
    .alpha = c(Observe, Shannon,Chao1,Simpson, Pielou),
    .group = group,
  ) +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  theme(
    legend.position="none",
    strip.background = element_rect(colour=NA, fill="grey")
  )+
  ggtitle('')+
  theme(plot.title = element_text(face = 'bold'))
alpha

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
pdf('D:/microbiome/AD/alpha.pdf',width = 8,height = 6)
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
  annotate("text", x = -0.3, y = -0.3, label = "ANOSIM R: 0.592 
           P: 0.001")+
  ggtitle('')+
  theme(plot.title = element_text(face = 'bold'))
pcoa_plot
pcoa_plot<-ggMarginal(pcoa_plot, type = "boxplot", margins = "both", size = 5, groupColour = TRUE, groupFill = TRUE)
pdf('D:/microbiome/AD/div.pdf',width = 8,height = 10)
alpha/pcoa_plot
dev.off()
pdf('D:/microbiome/AD/pcoa.pdf')
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
f.box <- mpse %>%
  mp_plot_diff_boxplot(
    .group = group,
  ) %>%
  set_diff_boxplot_color(
    values = c('#E7B800', '#00AFBB'),
    guide = guide_legend(title=NULL)
  )
pdf('D:/microbiome/AD/diff_box.pdf',width = 6,height = 10)
f.box
dev.off()
res<-cbind(mpse@NAMES,mpse@elementMetadata@listData[["LDAupper"]],
           mpse@elementMetadata@listData[["LDAmean"]],
           mpse@elementMetadata@listData[["LDAlower"]],
           mpse@elementMetadata@listData[["Sign_group"]],
           mpse@elementMetadata@listData[["pvalue"]],
           mpse@elementMetadata@listData[["fdr"]])%>%as.data.frame()
colnames(res)<-c('genus','LDAupper','LDAmean','LDAlower','Sign_group','pvalue','fdr')
write.csv(res,'D:/microbiome/AD/diff_genus.csv')
(alpha+pcoa_plot)/f.box
grid.arrange(arrangeGrob(arrangeGrob(alpha,pcoa_plot,nrow=1),f.box,nrow = 2))
########
