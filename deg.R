rm(list = ls())
library(dplyr)
library(tidyverse)
library(ggrepel)
library(UpSetR)
library(igraph)
library(RColorBrewer)
library(ggraph)
library(tidygraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(DESeq2)
library(tibble)
fba_res<-read.csv('D:/microbiome/AD/res_grow_FBA_agora2_CN/exchanges.csv')
meta<-read.csv('D:/microbiome/AD/data/meta_958.csv',row.names = 1)
fba_res<-fba_res[grep('_m',fba_res$metabolite),]
fba_res <- fba_res %>%
  pivot_wider(
    names_from = sample_id,
    values_from = flux,
    id_cols = metabolite,
    values_fill = 0
  )%>%column_to_rownames('metabolite')
meta<-meta[colnames(fba_res),,drop=F]
table(meta$group)
rownames(fba_res)<-gsub('_m','',rownames(fba_res))
if (!identical(colnames(fba_res), rownames(meta)[match(colnames(fba_res), rownames(meta))])) {
  stop("样本ID在代谢物数据和分组数据之间不匹配或顺序不一致。")
}

AD<-rownames(meta)[meta$group=='AD']
Healthy<-rownames(meta)[meta$group=='Healthy']



im<-fba_res
im[im<0]<-0
im<-abs(im)
res_im<-matrix(0,nrow(im),4)%>%as.data.frame()
colnames(res_im)<-c('mets','fc','pval','padj')
res_im$mets<-rownames(im)
for(i in 1:nrow(im)){
  print(i)
  res_im$fc[i]<-mean((im[i,match(AD,colnames(im))]%>%as.matrix()%>%as.numeric()))/mean((im[i,match(Healthy,colnames(im))]%>%as.matrix()%>%as.numeric()))
  res_im$pval[i]<-wilcox.test(im[i,match(AD,colnames(im))]%>%as.matrix()%>%as.numeric(),
                              im[i,match(Healthy,colnames(im))]%>%as.matrix()%>%as.numeric())$p.val
}
res_im$padj<-p.adjust(res_im$pval,method = 'fdr')
res_im<-na.omit(res_im)
res_im$label<-ifelse(res_im$padj<0.05,'yes','no')

ex<-fba_res
ex[ex>0]<-0
ex<-abs(ex)
res_ex<-matrix(0,nrow(ex),4)%>%as.data.frame()
colnames(res_ex)<-c('mets','fc','pval','padj')
res_ex$mets<-rownames(ex)
for(i in 1:nrow(ex)){
  print(i)
  res_ex$fc[i]<-mean((ex[i,match(AD,colnames(ex))]%>%as.matrix()%>%as.numeric()))/mean((ex[i,match(Healthy,colnames(ex))]%>%as.matrix()%>%as.numeric()))
  res_ex$pval[i]<-wilcox.test(ex[i,match(AD,colnames(ex))]%>%as.matrix()%>%as.numeric(),
                              ex[i,match(Healthy,colnames(ex))]%>%as.matrix()%>%as.numeric())$p.val
}
res_ex$padj<-p.adjust(res_ex$pval,method = 'fdr')
res_ex<-na.omit(res_ex)
res_ex$label<-ifelse(res_ex$padj<0.05,'yes','no')

sig_mets<-c(res_ex$mets[res_ex$label=='yes'],res_im$mets[res_im$label=='yes'])%>%unique

#####
library(data.table)
library(tidyr)
library(dplyr)
reaction<-fread('D:/microbiome/micom/202404/paper/data/ReactionDatabase.txt')%>%as.data.frame()
pathway<-unique(reaction$Subsystem)
pathway_data<-list()
a<-c()
for(i in 1:length(pathway)){
  if(pathway[i]!=''){
    temp<-reaction[reaction$Subsystem%in%pathway[i],]
    metabolite_temp<-strsplit(temp$Formula, split = "\\+|->|<=>")%>%unlist()
    metabolite_temp <- gsub("\\s*\\[[^\\]]*\\]\\s*|\\s+\\d+\\s+|\\s+", "", metabolite_temp)%>%unique%>%list
    names(metabolite_temp)<-pathway[i]
    pathway_data<-c(pathway_data,metabolite_temp)
  }
}
#saveRDS(pathway_data,'D:/lsr/202505/pathway_data.rds')


enrichment_analysis <- function(pathway_data, diff_metabolites) {
  # 初始化结果数据框
  result <- data.frame(
    pathway = character(),
    p_adj = numeric(),
    count = integer(),
    metabolites = character(),
    stringsAsFactors = FALSE
  )
  
  # 总代谢物数（用于超几何检验）
  total_metabolites <- unique(unlist(pathway_data))
  N <- length(total_metabolites)
  
  # 差异代谢物数
  n <- length(diff_metabolites)
  
  # 遍历每个通路
  for (pathway in names(pathway_data)) {
    # 通路中的代谢物
    pathway_metabolites <- pathway_data[[pathway]]
    K <- length(pathway_metabolites)
    
    # 差异代谢物中属于该通路的
    overlap <- intersect(diff_metabolites, pathway_metabolites)
    k <- length(overlap)
    
    if (k > 0) {
      # 超几何检验 p-value
      p_value <- phyper(
        q = k - 1,  # P(X >= k)
        m = K,
        n = N - K,
        k = n,
        lower.tail = FALSE
      )
      
      # 保留前10个代谢物
      if (k > 10) {
        overlap <- overlap[1:10]
      }
      
      # 添加到结果
      result <- rbind(result, data.frame(
        pathway = pathway,
        p_adj = p_value,  # 这里可以进一步做多重检验校正
        count = k,
        metabolites = paste(overlap, collapse = "/"),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # 对 p-value 进行多重检验校正（BH方法）
  result$p_adj <- p.adjust(result$p_adj, method = "BH")
  
  # 按 p_adj 排序
  result <- result[order(result$p_adj), ]
  
  return(result)
}
flux<-read.csv("D:/microbiome/AD/res_grow_pFBA_AD_agora2_wd/exchanges.csv",header = T)
#flux<-flux[flux$sample_id%in%selected$x,]
flux$metabolite<-gsub("\\[e\\]","_e",flux$metabolite)
doubleml<-read.csv('D:/microbiome/AD/doubleml.csv')
doubleml$padj<-p.adjust(doubleml$pvalue,method = 'BH')
#sig1<-doubleml[doubleml$padj<0.01,]
sig<-doubleml[doubleml$RV>0.15&doubleml$padj<0.01,]
table(sig$ATE>0)
sig$metabolite<-gsub('_m','',sig$metabolite)
#diff<-sig$metabolite#diff是差异代谢物的名字
#diff<-res_im[res_im$padj<0.05,]
diff<-sig$metabolite[sig$metabolite%in%res_ex$mets]#diff是差异代谢物的名字
diff<-gsub("\\[.*?\\]", "", diff)#删除[]和[]中的内容
use_pathway <- enrichment_analysis(pathway_data, diff)
colnames(use_pathway)<-c('Description','p.adjust','Count','geneID')
use_pathway$Category<-reaction$Subsystem_general[match(use_pathway$Description,reaction$Subsystem)]
use_pathway<-use_pathway[order(use_pathway$p.adjust,decreasing = F),]
use_pathway<-use_pathway[1:5,]


library(gground)
library(ggprism)
width <- 0.3#最左侧色块的宽度、色块与圆点、圆点与y轴之间的宽度
# x 轴长度
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
# 左侧分类标签数据
rect.data <- group_by(use_pathway, Category) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

use_pathway$Description<-factor(use_pathway$Description,levels = use_pathway$Description[order(use_pathway$Category)])

p_im<-use_pathway %>%
  ggplot(aes(-log10(p.adjust), y = Description, fill = Category)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8 #width调整矩形宽度
  ) +
  geom_text(
    aes(x = 0.05,label = Description),
    hjust = 0, size = 5
  ) +
  geom_text(
    aes(x = 0.1,label = geneID, colour = Category), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -width, label = Count)
  ) +
  scale_size_continuous(name = 'Count', 
                        range = c(5, 16),
                        labels = scales::number_format(accuracy = 1)) +
  # 分类标签
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = Category),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +
  #geom_text(
  #  aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = Category),
  #  data = rect.data,
  #  inherit.aes = FALSE
  #) + #左侧色块中的文字
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  labs(y = NULL) +
  #scale_fill_manual(name = 'Category', values = pal) +
  #scale_colour_manual(values = pal) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 2), 
    expand = expansion(c(0, 0))
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text()
  )
p_im

diff<-sig$metabolite[sig$metabolite%in%res_im$mets]#diff是差异代谢物的名字
diff<-gsub("\\[.*?\\]", "", diff)
use_pathway <- enrichment_analysis(pathway_data, diff)
colnames(use_pathway)<-c('Description','p.adjust','Count','geneID')
use_pathway$Category<-reaction$Subsystem_general[match(use_pathway$Description,reaction$Subsystem)]
use_pathway<-use_pathway[order(use_pathway$p.adjust,decreasing = F),]
use_pathway<-use_pathway[1:5,]


library(gground)
library(ggprism)
# x 轴长度
width <- 0.3
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
# 左侧分类标签数据
rect.data <- group_by(use_pathway, Category) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

use_pathway$Description<-factor(use_pathway$Description,levels = use_pathway$Description[order(use_pathway$Category)])

p_ex<-use_pathway %>%
  ggplot(aes(-log10(p.adjust), y = Description, fill = Category)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8 #width调整矩形宽度
  ) +
  geom_text(
    aes(x = 0.05,label = Description),
    hjust = 0, size = 5
  ) +
  geom_text(
    aes(x = 0.1,label = geneID, colour = Category), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -width, label = Count)
  ) +
  scale_size_continuous(name = 'Count', 
                        range = c(5, 16),
                        labels = scales::number_format(accuracy = 1)) +
  # 分类标签
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = Category),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +
  #geom_text(
  #  aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = Category),
  #  data = rect.data,
  #  inherit.aes = FALSE
  #) + #左侧色块中的文字
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  labs(y = NULL) +
  #scale_fill_manual(name = 'Category', values = pal) +
  #scale_colour_manual(values = pal) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 2), 
    expand = expansion(c(0, 0))
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text()
  )
p_ex
pdf('D:/microbiome/AD/pathway.pdf',width = 12,height = 9)
p_im/p_ex
dev.off()
#####
flux<-read.csv("D:/microbiome/AD/res_grow_pFBA_AD_agora2_wd/exchanges.csv",header = T)
#flux<-flux[flux$sample_id%in%selected$x,]
flux$metabolite<-gsub("\\[e\\]","_e",flux$metabolite)
doubleml<-read.csv('D:/microbiome/AD/doubleml.csv')
doubleml$padj<-p.adjust(doubleml$pvalue,method = 'BH')
#sig1<-doubleml[doubleml$padj<0.01,]
sig<-doubleml[doubleml$RV>0.15&doubleml$padj<0.01,]
table(sig$ATE>0)
sig$metabolite<-gsub('_m','',sig$metabolite)
diff<-sig$metabolite#diff是差异代谢物的名字
diff<-gsub("\\[.*?\\]", "", diff)#删除[]和[]中的内容
use_pathway <- enrichment_analysis(pathway_data, diff)
colnames(use_pathway)<-c('Description','p.adjust','Count','geneID')
use_pathway$Category<-reaction$Subsystem_general[match(use_pathway$Description,reaction$Subsystem)]
use_pathway<-use_pathway[order(use_pathway$p.adjust,decreasing = F),]
use_pathway<-use_pathway[1:10,]


library(gground)
library(ggprism)
# x 轴长度
width <- 0.6
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
# 左侧分类标签数据
rect.data <- group_by(use_pathway, Category) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

use_pathway$Description<-factor(use_pathway$Description,levels = use_pathway$Description[order(use_pathway$Category)])

p_sig<-use_pathway %>%
  ggplot(aes(-log10(p.adjust), y = Description, fill = Category)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8 #width调整矩形宽度
  ) +
  geom_text(
    aes(x = 0.05,label = Description),
    hjust = 0, size = 5
  ) +
  geom_text(
    aes(x = 0.1,label = geneID, colour = Category), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = -width, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -width, label = Count)
  ) +
  scale_size_continuous(name = 'Count', 
                        range = c(5, 16),
                        labels = scales::number_format(accuracy = 1)) +
  # 分类标签
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = Category),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +
  #geom_text(
  #  aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = Category),
  #  data = rect.data,
  #  inherit.aes = FALSE
  #) + #左侧色块中的文字
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  labs(y = NULL) +
  #scale_fill_manual(name = 'Category', values = pal) +
  #scale_colour_manual(values = pal) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 2), 
    expand = expansion(c(0, 0))
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text()
  )
p_sig
pdf('D:/microbiome/AD/pathway_imp_met.pdf',width = 12,height = 6)
p_sig
dev.off()
########
mets_data<-fread('D:/microbiome/micom/202404/paper/data/MetaboliteDatabase.txt')%>%as.data.frame()
mets_data<-mets_data[match(sig_mets,mets_data$V1),]

mets_data<-mets_data[mets_data$V7!='',]
url1<-'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
id<-paste(mets_data$V7%>%unique(), collapse = ",")
url2<-'/property/CanonicalSMILES/CSV'
url<-paste0(url1,id,url2)
url

stp<-read.csv('D:/microbiome/AD/SwissTargetPrediction.csv')
stp<-stp[stp$Probability>0,]
gene_stp<-unique(stp$Common_name)
gene_stp<-unlist(strsplit(gene_stp, split = "\\s+"))

gutmgene<-read.csv('D:/microbiome/AD/data/gutmgene.csv')
gutmgene$Metabolite..ID.<-sub("\\([^()]*\\)$", "", gutmgene$Metabolite..ID.)
gutmgene<-gutmgene[gutmgene$Metabolite..ID.%in%mets_data$V2,]

colnames(gutmgene)
edges_df <- gutmgene[, c("Metabolite..ID.", "Gene..ID.", "Alteration", "Evidence.Number")]
colnames(edges_df) <- c("Source", "Target", "Alteration", "Evidence.Number")
g <- graph_from_data_frame(d = edges_df, directed = TRUE)
min_width <- 1
max_width <- 3
E(g)$width <- min_width + (E(g)$Evidence.Number - min(E(g)$Evidence.Number)) / 
  (max(E(g)$Evidence.Number) - min(E(g)$Evidence.Number)) * (max_width - min_width)

# 如果 Evidence.Number 全是相同的，则 width 也会全是 NA，需要处理
if(is.na(min(E(g)$width))) {
  E(g)$width <- 1 # 如果所有 Evidence.Number 都相同，则设置为默认值1
}

alteration_levels <- unique(E(g)$Alteration)


# 2. 为每个级别生成颜色
# 选择一个好的调色板，例如 "Set1"
num_colors <- length(alteration_levels)
edge_color_palette <- brewer.pal(n = max(3, num_colors), name = "Set2")[1:num_colors] # 至少需要3种颜色

# 创建一个从 Alteration 值到颜色的映射
color_map <- setNames(edge_color_palette, alteration_levels)

g_tidy <- as_tbl_graph(g)

# 获取 Alteration 的唯一值数量，用于选择颜色
num_alterations <- length(unique(E(g)$Alteration))
g_tidy <- g_tidy %>%
  activate(nodes) %>%
  mutate(type = ifelse(name %in% gutmgene$`Metabolite..ID.`, 
                       "Metabolite", 
                       "Gene"))

ggraph(g_tidy, layout = 'nicely') +
  
  # 绘制边
  geom_edge_fan(aes(
    color = Alteration, 
    width = width
  ), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(3, 'mm')) + 
  
  # 绘制节点，根据 type 设置颜色
  geom_node_point(aes(color = type), size = 5) + 
  
  # 绘制标签
  geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 20) + 
  
  # 设置颜色和主题
  scale_edge_color_manual(values = color_map, name = "Alteration Type") +
  scale_color_manual(values = c("Gene" = "skyblue", "Metabolite" = "orange"), name = "Node Type") + # 节点颜色手动设置
  scale_edge_width(range = c(min_width, max_width), name = "Evidence Number (Scaled)") +
  
  theme_graph() +
  guides(edge_width = "none") +
  labs(title = "")
graph_gutm<-ggraph(g_tidy, layout = 'nicely') +
  
  # 绘制边
  geom_edge_fan(aes(
    color = Alteration, 
    width = width
  ), arrow = arrow(length = unit(2, 'mm')), end_cap = circle(3, 'mm')) + 
  
  # 绘制节点，根据 type 设置颜色
  geom_node_point(aes(color = type), size = 5) + 
  
  # 绘制标签
  geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 20) + 
  
  # 设置颜色和主题
  scale_edge_color_manual(values = color_map, name = "Alteration Type") +
  scale_color_manual(values = c("Gene" = "skyblue", "Metabolite" = "orange"), name = "Node Type") + # 节点颜色手动设置
  scale_edge_width(range = c(min_width, max_width), name = "Evidence Number (Scaled)") +
  
  theme_graph() +
  guides(edge_width = "none") +
  labs(title = "")


gene_gut<-unique(gutmgene$Gene..ID.)
gene_gut <- sub(" \\(.*", "", gene_gut)

genecard<-read.csv('D:/microbiome/AD/data//GeneCards.csv')
genecard_gene<-unique(genecard$Gene.Symbol)


library(data.table)
AD<-fread('D:/microbiome/AD/data/all_gene_disease_associations.txt')%>%as.data.frame()
#AD<-AD[AD$diseaseId=='umls:C0002395',]
AD<-AD[AD$diseaseId=='umls:C0002395'&AD$score>0,]
AD_gene<-unique(AD$geneSymbol)
#intersect(AD_gene,target)


#genecard<-read.csv('D:/microbiome/AD/data/CNP0001162/GeneCards-AD.csv')
#genecard_gene<-unique(genecard$Gene.Symbol)

#omim<-read.csv('D:/microbiome/AD/data/CNP0001162/omim.csv')
#omim_gene<-unique(omim$OMIM)
gene_list<-list(gene_stp,gene_gut,AD_gene,genecard_gene)
target_AD<-Reduce(intersect, gene_list)

listInput <- list(#OMIM = symbol,
                  DisGeNET = AD_gene,
                  #GeneCards = genecard_gene,
                  #STP_target = gene_stp,
                  gutMGene_target = gene_gut)
upset(fromList(listInput), 
      sets = c("DisGeNET",
               #'STP_target',
               'gutMGene_target'), 
      order.by = "freq", 
      main.bar.color = "steelblue",
      sets.bar.color = "grey40")

upset_plot<-upset(fromList(listInput), 
                  sets = c("DisGeNET",'STP_target','gutMGene_target'), 
                  order.by = "freq", 
                  main.bar.color = "steelblue",
                  sets.bar.color = "grey40")

#############
rm(list = ls())
fba_res<-read.csv('D:/microbiome/AD/res_grow_FBA_agora2_CN/exchanges.csv')
meta<-read.csv('D:/microbiome/AD/data/meta_958.csv',row.names = 1)
#fba_res<-fba_res[grep('_m',fba_res$metabolite),]
fba_res$label<-paste0(fba_res$taxon,'_',fba_res$metabolite)
dup_rows <- fba_res %>%
  group_by(label, sample_id) %>%
  filter(n() > 1) %>%
  arrange(label, sample_id)
fba_res<-fba_res[!fba_res$label%in%dup_rows$label,]
setDT(fba_res)
fba_wide <- dcast(fba_res, label ~ sample_id, value.var = "flux")
fba_but <- fba_wide[grep('_but\\[e\\]',fba_wide$label),]%>%as.data.frame()
fba_but[is.na(fba_but)]<-0
fba_but<-fba_but%>%column_to_rownames('label')
meta<-meta[colnames(fba_but),,drop=F]
table(meta$group)
#rownames(fba_but)<-gsub('_m','',rownames(fba_but))
if (!identical(colnames(fba_but), rownames(meta)[match(colnames(fba_but), rownames(meta))])) {
  stop("样本ID在代谢物数据和分组数据之间不匹配或顺序不一致。")
}

AD<-rownames(meta)[meta$group=='AD']
Healthy<-rownames(meta)[meta$group=='Healthy']



im<-fba_but
im[im<0]<-0
im<-abs(im)
res_im<-matrix(0,nrow(im),4)%>%as.data.frame()
colnames(res_im)<-c('mets','fc','pval','padj')
res_im$mets<-rownames(im)
for(i in 1:nrow(im)){
  print(i)
  res_im$fc[i]<-mean((im[i,match(AD,colnames(im))]%>%as.matrix()%>%as.numeric()))/mean((im[i,match(Healthy,colnames(im))]%>%as.matrix()%>%as.numeric()))
  res_im$pval[i]<-wilcox.test(im[i,match(AD,colnames(im))]%>%as.matrix()%>%as.numeric(),
                              im[i,match(Healthy,colnames(im))]%>%as.matrix()%>%as.numeric())$p.val
}
res_im$padj<-p.adjust(res_im$pval,method = 'fdr')
res_im<-na.omit(res_im)
res_im$label<-ifelse(res_im$padj<0.05,'yes','no')
res_im$taxon<-gsub('_but\\[e\\]','',res_im$mets)


ex<-fba_but
ex[ex>0]<-0
ex<-abs(ex)
res_ex<-matrix(0,nrow(ex),4)%>%as.data.frame()
colnames(res_ex)<-c('mets','fc','pval','padj')
res_ex$mets<-rownames(ex)
for(i in 1:nrow(ex)){
  print(i)
  res_ex$fc[i]<-mean((ex[i,match(AD,colnames(ex))]%>%as.matrix()%>%as.numeric()))/mean((ex[i,match(Healthy,colnames(ex))]%>%as.matrix()%>%as.numeric()))
  res_ex$pval[i]<-wilcox.test(ex[i,match(AD,colnames(ex))]%>%as.matrix()%>%as.numeric(),
                              ex[i,match(Healthy,colnames(ex))]%>%as.matrix()%>%as.numeric())$p.val
}
res_ex$padj<-p.adjust(res_ex$pval,method = 'fdr')
res_ex<-na.omit(res_ex)
res_ex$label<-ifelse(res_ex$padj<0.05,'yes','no')
res_ex$taxon<-gsub('_but\\[e\\]','',res_ex$mets)

tax_im<-res_im$taxon
tax_ex<-res_ex$taxon
intersect(tax_im,tax_ex)
tax_im<-setdiff(tax_im,intersect(tax_im,tax_ex))
tax_ex<-setdiff(tax_ex,intersect(tax_im,tax_ex))

diff_genus<-read.csv('D:/microbiome/AD/diff_genus.csv',row.names = 1)
diff_genus<-diff_genus[diff_genus$fdr<0.05,]
tax_im<-setdiff(tax_im,diff_genus$genus)
tax_ex<-setdiff(tax_ex,diff_genus$genus)

##############




##############

















###########
#GSE276756<-read.csv('D:/microbiome/AD/data/GSE276756.csv',row.names = 1)
#
#symbol <- bitr(rownames(GSE276756),
#               fromType = "ENSEMBL",
#               toType   = "SYMBOL",
#               OrgDb    = org.Hs.eg.db
#)
#GSE276756$symbol<-NA
#GSE276756$symbol[match(symbol$ENSEMBL,rownames(GSE276756))]<-symbol$SYMBOL
#GSE276756<-na.omit(GSE276756)
#setDT(GSE276756)
#GSE276756 <- GSE276756[, lapply(.SD, function(x) round(sum(x, na.rm = TRUE))), 
#             by = 'symbol', 
#             .SDcols = which(sapply(GSE276756, is.numeric))]
#GSE276756<-as.data.frame(GSE276756)
#GSE276756<-GSE276756%>%column_to_rownames('symbol')
#meta_276576<-colnames(GSE276756)%>%as.data.frame()
#meta_276576$group<-meta_276576$.
#meta_276576<-meta_276576%>%column_to_rownames('.')
#meta_276576$group[grep('AD',meta_276576$group)]<-'AD'
#meta_276576$group[grep('HC',meta_276576$group)]<-'HC'
#
#meta_276576$group<-as.factor(meta_276576$group)
#dds <- DESeqDataSetFromMatrix(countData = GSE276756,colData = meta_276576,design = ~group)
#dds <- DESeq(dds, parallel = FALSE)
#res <- results(dds, contrast = c('group', 'HC', 'AD'), pAdjustMethod = 'fdr', alpha = 0.05)
#summary(res)
#res1 <- as.data.frame(res[order(res$padj), ])
#DEG<-na.omit(res1)
#deseq_res<-DEG[DEG$padj<0.05&abs(DEG$log2FoldChange)>1,]
#significant <- res1$padj < 0.05 & abs(res1$log2FoldChange)>1
#res1$color_val <- -log10(res1$padj)
#res1$color_val[!significant] <- NA 
#colnames(res1)[6]<-'adj.P.Val'
#
#P_threshold <- 0.05
#Y_cutoff <- -log10(P_threshold) # Y_cutoff 约等于 1.301
#
## 设置 logFC 的阈值
#LogFC_cutoff <- 1
#ggplot(data = res1,
#       mapping = aes(x = log2FoldChange, y = -log10(adj.P.Val))) +
#  
#  # 1. 添加水平辅助线 (P-value 阈值)
#  geom_hline(yintercept = Y_cutoff, 
#             linetype = "dashed", # 使用虚线
#             color = "black",     # 颜色
#             linewidth = 0.5) +   # 线的粗细
#  
#  # 2. 添加垂直辅助线 (logFC 阈值)
#  geom_vline(xintercept = c(-LogFC_cutoff, LogFC_cutoff), 
#             linetype = "dashed", 
#             color = "black", 
#             linewidth = 0.5) +
#  
#  geom_point(aes(color = color_val,
#                 size = -log10(adj.P.Val)),
#             alpha = 1) +
#  theme_bw() +
#  scale_color_gradientn(
#    colours = rev(brewer.pal(11, "RdYlBu")),
#    na.value = "gray",
#    name = "-log10(adj.P.Val)"
#  ) +
#  scale_size_continuous(range = c(0.3, 3)) +
#  xlim(c(-5,5))+
#  xlab("logFC") + ylab("-log10(adj.P.Val)") +
#  guides(color = guide_colorbar(title = "-log10(adj.P.Val)"))
#
##########
#GSE261050<-read.csv('D:/microbiome/AD/data/GSE261050.csv',row.names = 1)
#meta_261050<-read.csv('D:/microbiome/AD/data/GSE261050_metadata.csv',row.names = 1)
#
#symbol <- bitr(rownames(GSE261050),
#               fromType = "ENSEMBL",
#               toType   = "SYMBOL",
#               OrgDb    = org.Hs.eg.db
#)
#GSE261050$symbol<-NA
#GSE261050$symbol[match(symbol$ENSEMBL,rownames(GSE261050))]<-symbol$SYMBOL
#GSE261050<-na.omit(GSE261050)
#setDT(GSE261050)
#GSE261050 <- GSE261050[, lapply(.SD, function(x) round(sum(x, na.rm = TRUE))), 
#                       by = 'symbol', 
#                       .SDcols = which(sapply(GSE261050, is.numeric))]
#GSE261050<-as.data.frame(GSE261050)
#GSE261050<-GSE261050%>%column_to_rownames('symbol')
#meta_261050$Group<-factor(meta_261050$Group,levels = c('control','AD'))
#dds <- DESeqDataSetFromMatrix(countData = GSE261050,colData = meta_261050,design = ~Group)
#dds <- DESeq(dds, parallel = FALSE)
#res <- results(dds, contrast = c('Group', 'control', 'AD'), pAdjustMethod = 'fdr', alpha = 0.05)
#summary(res)
#res1 <- as.data.frame(res[order(res$padj), ])
#DEG<-na.omit(res1)
#deseq_res<-DEG[DEG$padj<0.05&abs(DEG$log2FoldChange)>1,]
#significant <- res1$padj < 0.05 & abs(DEG$log2FoldChange)>1
#res1$color_val <- -log10(res1$padj)
#res1$color_val[!significant] <- NA 
#colnames(res1)[6]<-'adj.P.Val'
#
#P_threshold <- 0.05
#Y_cutoff <- -log10(P_threshold) # Y_cutoff 约等于 1.301
#
## 设置 logFC 的阈值
#LogFC_cutoff <- 1
#ggplot(data = res1,
#       mapping = aes(x = log2FoldChange, y = -log10(adj.P.Val))) +
#  
#  # 1. 添加水平辅助线 (P-value 阈值)
#  geom_hline(yintercept = Y_cutoff, 
#             linetype = "dashed", # 使用虚线
#             color = "black",     # 颜色
#             linewidth = 0.5) +   # 线的粗细
#  
#  # 2. 添加垂直辅助线 (logFC 阈值)
#  geom_vline(xintercept = c(-LogFC_cutoff, LogFC_cutoff), 
#             linetype = "dashed", 
#             color = "black", 
#             linewidth = 0.5) +
#  
#  geom_point(aes(color = color_val,
#                 size = -log10(adj.P.Val)),
#             alpha = 1) +
#  theme_bw() +
#  scale_color_gradientn(
#    colours = rev(brewer.pal(11, "RdYlBu")),
#    na.value = "gray",
#    name = "-log10(adj.P.Val)"
#  ) +
#  scale_size_continuous(range = c(0.3, 3)) +
#  xlim(c(-5,5))+
#  xlab("logFC") + ylab("-log10(adj.P.Val)") +
#  guides(color = guide_colorbar(title = "-log10(adj.P.Val)"))
##########
#GSE249477<-read.csv('D:/microbiome/AD/data/GSE249477.csv',row.names = 1,quote = '')
#meta_249477<-read.csv('D:/microbiome/AD/data/meta_249477.csv',row.names = 1,header = F)
#colnames(meta_249477)<-'group'
#meta_249477$id<-sub(".*\\[([^][]*)\\]$", "\\1", meta_249477$group)
#meta_249477$id<-gsub(' ','',meta_249477$id)
#rownames(meta_249477)<-meta_249477$id
#meta_249477$id<-NULL
#GSE249477<-GSE249477[,grep('counts',colnames(GSE249477))]
#meta_249477$group[grep('control',meta_249477$group)]<-'control'
#meta_249477$group[grep('mild',meta_249477$group)]<-'MCI'
#meta_249477$group[grep('Alzheimer',meta_249477$group)]<-'AD'
#meta_249477<-meta_249477[meta_249477$group!='MCI',,drop=F]
#
#colnames(GSE249477)<-substr(colnames(GSE249477),1,10)
#GSE249477<-GSE249477[,rownames(meta_249477)]
#meta_249477$group<-factor(meta_249477$group,levels = c('control','AD'))
#GSE249477[is.na(GSE249477)]<-0
#GSE249477<-round(GSE249477)
#dds <- DESeqDataSetFromMatrix(countData = GSE249477,colData = meta_249477,design = ~group)
#dds <- DESeq(dds, parallel = FALSE)
#res <- results(dds, contrast = c('group', 'control', 'AD'), pAdjustMethod = 'fdr', alpha = 0.05)
#summary(res)
#res1 <- as.data.frame(res[order(res$padj), ])
#DEG<-na.omit(res1)
#deseq_res<-DEG[DEG$pval<0.05&abs(DEG$log2FoldChange)>1,]
#significant <- res1$pval < 0.05 & abs(res1$log2FoldChange)>1
#res1$color_val <- -log10(res1$pval)
#res1$color_val[!significant] <- NA 
#colnames(res1)[6]<-'adj.P.Val'
#
#P_threshold <- 0.05
#Y_cutoff <- -log10(P_threshold) # Y_cutoff 约等于 1.301
#
## 设置 logFC 的阈值
#LogFC_cutoff <- 1
#ggplot(data = res1,
#       mapping = aes(x = log2FoldChange, y = -log10(pvalue))) +
#  
#  # 1. 添加水平辅助线 (P-value 阈值)
#  geom_hline(yintercept = Y_cutoff, 
#             linetype = "dashed", # 使用虚线
#             color = "black",     # 颜色
#             linewidth = 0.5) +   # 线的粗细
#  
#  # 2. 添加垂直辅助线 (logFC 阈值)
#  geom_vline(xintercept = c(-LogFC_cutoff, LogFC_cutoff), 
#             linetype = "dashed", 
#             color = "black", 
#             linewidth = 0.5) +
#  
#  geom_point(aes(color = color_val,
#                 size = -log10(pvalue)),
#             alpha = 1) +
#  theme_bw() +
#  scale_color_gradientn(
#    colours = rev(brewer.pal(11, "RdYlBu")),
#    na.value = "gray",
#    name = "-log10(pvalue)"
#  ) +
#  scale_size_continuous(range = c(0.3, 3)) +
#  xlim(c(-5,5))+
#  xlab("logFC") + ylab("-log10(pvalue)") +
#  guides(color = guide_colorbar(title = "-log10(pvalue)"))
###################
GSE276756<-read.csv('D:/microbiome/AD/data/GSE276756.csv',row.names = 1)

symbol <- bitr(rownames(GSE276756),
               fromType = "ENSEMBL",
               toType   = "SYMBOL",
               OrgDb    = org.Hs.eg.db
)
GSE276756$symbol<-NA
GSE276756$symbol[match(symbol$ENSEMBL,rownames(GSE276756))]<-symbol$SYMBOL
GSE276756<-na.omit(GSE276756)
setDT(GSE276756)
GSE276756 <- GSE276756[, lapply(.SD, function(x) round(sum(x, na.rm = TRUE))), 
                       by = 'symbol', 
                       .SDcols = which(sapply(GSE276756, is.numeric))]
GSE276756<-as.data.frame(GSE276756)
GSE276756<-GSE276756%>%column_to_rownames('symbol')
meta_276576<-colnames(GSE276756)%>%as.data.frame()
meta_276576$group<-meta_276576$.
meta_276576<-meta_276576%>%column_to_rownames('.')
meta_276576$group[grep('AD',meta_276576$group)]<-'AD'
meta_276576$group[grep('HC',meta_276576$group)]<-'HC'

GSE261050<-read.csv('D:/microbiome/AD/data/GSE261050.csv',row.names = 1)
meta_261050<-read.csv('D:/microbiome/AD/data/GSE261050_metadata.csv',row.names = 1)

symbol <- bitr(rownames(GSE261050),
               fromType = "ENSEMBL",
               toType   = "SYMBOL",
               OrgDb    = org.Hs.eg.db
)
GSE261050$symbol<-NA
GSE261050$symbol[match(symbol$ENSEMBL,rownames(GSE261050))]<-symbol$SYMBOL
GSE261050<-na.omit(GSE261050)
setDT(GSE261050)
GSE261050 <- GSE261050[, lapply(.SD, function(x) round(sum(x, na.rm = TRUE))), 
                       by = 'symbol', 
                       .SDcols = which(sapply(GSE261050, is.numeric))]
GSE261050<-as.data.frame(GSE261050)
GSE261050<-GSE261050%>%column_to_rownames('symbol')

exp<-merge(GSE276756,GSE261050,by = "row.names",all = TRUE)
exp[is.na(exp)]<-0
rownames(exp)<-exp$Row.names
exp$Row.names<-NULL

meta_261050<-meta_261050[,2,drop=F]
colnames(meta_261050)<-'group'
meta<-rbind(meta_261050,meta_276576)
meta$group[meta$group=='HC']<-'control'
meta$group<-factor(meta$group,levels = c('AD','control'))

exp<-exp[,rownames(meta)]

dds <- DESeqDataSetFromMatrix(countData = exp,colData = meta,design = ~group)
dds <- DESeq(dds, parallel = FALSE)
res <- results(dds, contrast = c('group', 'AD', 'control'), pAdjustMethod = 'fdr', alpha = 0.05)
summary(res)
res1 <- as.data.frame(res[order(res$padj), ])
DEG<-na.omit(res1)
deseq_res<-DEG[DEG$padj<0.05,]
significant <- res1$padj < 0.05
res1$color_val <- -log10(res1$padj)
res1$color_val[!significant] <- NA 
colnames(res1)[6]<-'adj.P.Val'

P_threshold <- 0.05
Y_cutoff <- -log10(P_threshold) # Y_cutoff 约等于 1.301
res1$Gene<-rownames(res1)
res1$label_gene <- ifelse(res1$adj.P.Val < P_threshold, res1$Gene, NA)
target_gene <- "PPARG"
res1$is_PPARG <- (res1$Gene %in% intersect(gene_gut,AD_gene))&res1$adj.P.Val<0.05
#res1['PPARG',][8]<-'PPARG(P.adj = 0.004)'

deg<-ggplot(data = res1,
            mapping = aes(x = log2FoldChange, y = -log10(adj.P.Val))) +
  
  # 1. 添加辅助线
  geom_hline(yintercept = Y_cutoff, 
             linetype = "dashed", 
             color = "black",     
             linewidth = 0.5) +   
  
  # 2. 绘制所有非 PPARG 的点 (使用 color_val 渐变)
  geom_point(data = subset(res1, !is_PPARG), # 排除 PPARG
             aes(color = color_val,
                 size = -log10(adj.P.Val)),
             alpha = 0.7) +
  
  # 3. 绘制 PPARG 点 (突出显示)
  geom_point(data = subset(res1, is_PPARG), # 只选择 PPARG
             color = "red",              # 设置突出颜色为紫色
             #shape = 18,                    # 设置突出形状 (例如菱形)
             size = 2,                      # 增大点的大小
             alpha = 1) +
  
  # 4. 标记显著基因的名称 (使用 ggrepel 避免重叠)
  #geom_text_repel(data = subset(res1, !is.na(label_gene)), # 只标记显著基因
  #                aes(label = label_gene),
  #                size = 3.5,
  #                box.padding = 0.5, # 调整标签与点的距离
  #                point.padding = 0.5,
  #                min.segment.length = 0, # 确保始终有引导线
  #                segment.color = 'grey50') +
  geom_text_repel(data = subset(res1, is_PPARG), 
                  aes(label = Gene), 
                  color = "red", # 标签颜色与点颜色一致
                  size = 4, 
                  box.padding = 1,
                  point.padding = 1,
                  fontface = "bold") +
  # 5. 主题和标尺设置
  theme_bw() +
  scale_color_gradientn(
    colours = rev(brewer.pal(11, "RdYlBu")),
    na.value = "gray", # 非显著基因的点（NA）使用灰色
    name = "-log10(adj.P.Val)"
  ) +
  scale_size_continuous(range = c(0.3, 3)) +
  xlim(c(-5,5))+
  xlab("logFC") + ylab("-log10(adj.P.Val)") +
  guides(color = guide_colorbar(title = "-log10(adj.P.Val)"))


(graph_gutm+upset_plot)/deg
grid.arrange(arrangeGrob(arrangeGrob(graph_gutm,upset_plot,nrow=1),deg))

pdf('../graph_gutm.pdf',width = 9,height = 8)
graph_gutm
dev.off()

pdf('../upsetplot.pdf',width = 5,height = 4)
upset_plot
dev.off()



pdf('../deg.pdf',width = 4,height = 4,)
ggplot(data = res1,
       mapping = aes(x = log2FoldChange, y = -log10(adj.P.Val))) +
  
  # 1. 添加辅助线
  geom_hline(yintercept = Y_cutoff, 
             linetype = "dashed", 
             color = "black",     
             linewidth = 0.5) +   
  
  # 2. 绘制所有非 PPARG 的点 (使用 color_val 渐变)
  geom_point(data = subset(res1, !is_PPARG), # 排除 PPARG
             aes(color = color_val,
                 size = -log10(adj.P.Val)),
             alpha = 0.7) +
  
  # 3. 绘制 PPARG 点 (突出显示)
  geom_point(data = subset(res1, is_PPARG), # 只选择 PPARG
             color = "red",              # 设置突出颜色为紫色
             #shape = 18,                    # 设置突出形状 (例如菱形)
             size = 2,                      # 增大点的大小
             alpha = 1) +
  
  # 4. 标记显著基因的名称 (使用 ggrepel 避免重叠)
  #geom_text_repel(data = subset(res1, !is.na(label_gene)), # 只标记显著基因
  #                aes(label = label_gene),
  #                size = 3.5,
  #                box.padding = 0.5, # 调整标签与点的距离
  #                point.padding = 0.5,
  #                min.segment.length = 0, # 确保始终有引导线
  #                segment.color = 'grey50') +
  geom_text_repel(data = subset(res1, is_PPARG), 
                  aes(label = Gene), 
                  color = "red", # 标签颜色与点颜色一致
                  size = 4, 
                  box.padding = 1,
                  point.padding = 1,
                  fontface = "bold") +
  # 5. 主题和标尺设置
  theme_bw() +
  scale_color_gradientn(
    colours = rev(brewer.pal(11, "RdYlBu")),
    na.value = "gray", # 非显著基因的点（NA）使用灰色
    name = "-log10(adj.P.Val)"
  ) +
  scale_size_continuous(range = c(0.3, 3)) +
  xlim(c(-5,5))+
  xlab("logFC") + ylab("-log10(adj.P.Val)") +
  guides(color = guide_colorbar(title = "-log10(adj.P.Val)"))
dev.off()
############

