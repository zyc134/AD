rm(list = ls())
gc()
# ==============================================================================
# 1. 库加载 (Load Libraries)
# ==============================================================================
library(NetCoMi)        # 核心：用于微生物共现网络构建和比较
library(tidyverse)      # 集合：包括ggplot2, dplyr等，用于数据处理和可视化
library(tidyr)
library(reshape2)       # 用于数据框的重塑（melt/cast）
library(ggraph)         # 用于高级网络可视化
library(igraph)         # 用于图论和网络分析
library(colorspace)
library(RColorBrewer)
library(graphlayouts)
library(ggvenn)
library(data.table)
library(fgsea)          # 快速GSEA（基因集富集分析），这里用于物种集富集分析
library(grDevices)
library(ggheatmap)
library(gggsea)
library(ggrepel)        # 用于避免绘图中文本标签重叠
library(showtext)       # 用于在R中显示字体
library(ggforce)
library(clusterProfiler)
library(enrichplot)
library(GseaVis)
library(aplot)
library(patchwork)      # 用于组合多个ggplot图
library(gridExtra)      # 用于组合多个图（包括非ggplot对象）
library(cowplot)
#library(curatedMetagenomicData) # 可能用于验证数据或获取元数据
library(ggpubr)         # 用于添加统计比较结果等
showtext_auto()         # 自动开启字体显示功能

# ==============================================================================
# 2. 数据加载与预处理 (Data Loading and Preprocessing)
# ==============================================================================
# 加载属水平的物种丰度数据
#abundance_D<-read.csv('D:/microbiome/AD/data/AD/genus/abundance_genus.csv',row.names = 1)
#colSums(abundance_D)%>%sort()%>%plot
#colSums(abundance_D)%>%hist
#temp<-colSums(abundance_D)
##for(i in 1:ncol(abundance_D)){
##  abundance_D[,i]<-abundance_D[,i]/temp[i]
##}
##files<-list.files('D:/microbiome/micom/agora103_genus/890088c4-1a67-498b-8a54-3eb31ee06f63/data')
##for(i in 1:length(files)){
##  files[i]<-substr(files[i],1,nchar(files[i])-5)
##}
##abundance_D1<-abundance_D[intersect(rownames(abundance_D),files),]
##colSums(abundance_D1)%>%hist
##abundance_D1<-abundance_D1[,colSums(abundance_D1)>0.8]
##abundance_D<-abundance_D[,colnames(abundance_D1)]
##temp_D<-abundance_D1
##temp_D[temp_D!=0]<-1
##colSums(temp_D)%>%hist
#abundance_H<-fread('D:/microbiome/AD/data/Healthy/genus/abundance_genus.csv')%>%as.data.frame()%>%column_to_rownames('V1')
#colSums(abundance_H)%>%sort()%>%plot
#colSums(abundance_H)%>%hist
#abundance_H<-abundance_H[,colSums(abundance_H)<=100]
##abundance_H1<-abundance_H[intersect(rownames(abundance_H),files),]
##colSums(abundance_H1)%>%hist()
#abundance_H<-abundance_H[,colSums(abundance_H)>90&colSums(abundance_H)<97]
#abundance_H<-abundance_H[which(rowSums(abundance_H)!=0),]
#
#
##temp_H<-abundance_H
##temp_H[temp_H!=0]<-1
##colSums(temp_H)%>%hist
##abundance_H<-abundance_H[,colnames(abundance_H1)]
#data<-merge(abundance_D,abundance_H,by = "row.names",all = TRUE)
#data[is.na(data)]<-0
#rownames(data)<-data$Row.names
#data$Row.names<-NULL
#ab<-data
# 加载样本元数据
ab_genus<-read.csv('D:/microbiome/AD/data/abundance_cn.csv',row.names = 1)
meta<-read.csv('D:/microbiome/AD/1118/meta.csv',row.names = 1)
table(meta$group)
# 加载筛选后的样本列表

#selected<-read.csv("data/res_grow_wd_pFBA/selected_sample.csv")

# 筛选元数据，只保留被选中的样本
#meta<-meta[meta$Sample.Name%in%selected$x,]
# 格式化样本ID，将"-"替换为"."，以匹配丰度矩阵的列名
meta$sample_id<-rownames(meta)

# 确保丰度矩阵的列顺序与元数据中的样本顺序一致
#ab_genus<-ab
ab_genus<-ab_genus[,meta$sample_id]

# 过滤物种：保留在至少 10% 样本中丰度不为零的物种
#keep<-apply(ab_genus, 1, function(x){
#  mean(x!=0)>0.1
#})
#ab_genus<-ab_genus[keep,]

# 定义样本分组：根据WHO标准，将收缩压(Systolic_BP)或舒张压(diastolic_bp)升高定义为高血压
#meta$group<-ifelse(meta$diastolic_bp>=90|meta$Systolic_BP>=140,"AD","Healthy") ##WHO


# ==============================================================================
# 3. 网络构建 (Network Construction - NetCoMi)
# ==============================================================================

# 第一次 netConstruct 调用：用于计算稀疏化阈值 (thresh)
# 目的：使用 SparCC 方法构建两个组（高血压 vs 健康）的网络
net <- netConstruct(data = t(ab_genus)[meta$group=="AD",],data2 = t(ab_genus)[meta$group=="Healthy",], 
                    filtTax = "none",
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta$group=="Healthy")), # 确保两组样本数相等
                    measure = "sparcc", # 使用 SparCC 测量关联
                    normMethod = "clr", # 使用中心对数比 (CLR) 转换进行归一化
                    zeroMethod = "alrEM", # 使用 ALR-EM 方法处理零值
                    sparsMethod = "threshold",
                    thresh = 0, # 初始阈值设为0，计算所有相关性
                    dissFunc = "signed", # 使用带符号的距离/关联矩阵
                    verbose = 3,
                    seed = 123)

# 计算稀疏化阈值 (Threshold Calculation)
# 公式：网络中非对角线绝对关联值的总和 / (2 * 节点数 * (节点数 - 1))
# 目标：获得一个平均关联度
thresh=(sum(abs(net[["assoMat1"]])+abs(net[["assoMat2"]]))-2*nrow(net[["assoMat2"]]))/(2*nrow(net[["assoMat2"]])*(nrow(net[["assoMat2"]])-1))

# 第二次 netConstruct 调用：使用计算出的阈值进行最终网络构建
net <- netConstruct(data = t(ab_genus)[meta$group=="AD",],data2 = t(ab_genus)[meta$group=="Healthy",], 
                    filtTax = "none",
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta$group=="Healthy")),
                    measure = "sparcc",
                    normMethod = "clr", 
                    zeroMethod = "alrEM",
                    sparsMethod = "threshold",
                    thresh = thresh, # 应用计算出的稀疏化阈值
                    dissFunc = "signed",
                    verbose = 3,
                    seed = 123456)

# ==============================================================================
# 4. 网络拓扑分析与比较 (Network Topology Analysis and Comparison)
# ==============================================================================
# 分析网络拓扑性质 (Props)
props <- netAnalyze(net, 
                    centrLCC = FALSE,
                    avDissIgnoreInf = TRUE,
                    sPathNorm = FALSE,
                    clustMethod = "cluster_fast_greedy", # 使用快速贪婪聚类发现模块（公会）
                    hubPar = c("degree", "eigenvector"),
                    hubQuant = 0.9, # 识别度或特征向量中心性在前90%的节点为Hub
                    lnormFit = TRUE,
                    normDeg = FALSE,
                    normBetw = FALSE,
                    normClose = FALSE,
                    normEigen = FALSE)

# 比较两组网络的拓扑性质 (Comp)
comp <- netCompare(props, 
                   permTest = FALSE, # 这里不执行置换检验，因为作者可能想在后续 diffnet 中进行
                   verbose = FALSE,
                   seed = 123456)

# 汇总比较结果
summ_res<-summary(comp, 
                  groupNames = c("AD","Healthy"),pAdjust=T,
                  showCentr="all",
                  numbNodes = 200)
summ_res
# 差异网络分析 (DiffNet)
# 目的：识别两组间关联强度有显著差异的边
net<-readRDS('D:/microbiome/AD/1118/net.rds')
diff <- diffnet(net,discordThresh=0.8,alpha = 0.05,lfdrThresh = 0.2,
                diffMethod = "permute", nPerm = 100,cores = 4, # 使用置换检验（1000次）
                storeCountsPerm = TRUE, 
                fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
                storeAssoPerm = TRUE,
                fileStoreAssoPerm = "assoPerm",
                adjust = "none",seed = 123)
saveRDS(diff,"D:/microbiome/AD/1118/specie_network_res.rds") # 保存结果
diff<-readRDS("D:/microbiome/AD/specie_network_res_958_filtered.rds") # 读取结果

# 绘制差异网络可视化图 (DiffNet Visualization)
plot(diff,
     cexNodes = 1, 
     cexLegend = 0.8,
     cexTitle = 2,
     mar = c(10,10,10,20),layout=NULL,
     legendGroupnames = c("AD", "Healthy"),
     legendPos = c(1.1,0.8))

# 再次分析网络拓扑（可能用于绘图，这次不进行聚类，仅计算度）
props_pears <- netAnalyze(net, 
                          clustMethod = "cluster_fast_greedy",
                          weightDeg = TRUE,
                          normDeg = FALSE,
                          gcmHeat = FALSE)

# 绘制差异关联节点图 (Differential Association Nodes Plot)
# 筛选在差异网络中至少有一条差异边的节点
diffmat_sums <- rowSums(diff$diffAdjustMat)
diff_asso_names <- names(diffmat_sums[diffmat_sums > 0])
pdf("output/comp_network.pdf",width = 12,height = 6)
plot(props_pears, 
     nodeFilter = "names",
     nodeFilterPar = diff_asso_names, # 只绘制有差异关联的节点
     nodeColor = "gray",
     highlightHubs = FALSE,
     sameLayout = TRUE, 
     layoutGroup = "union",
     rmSingles = FALSE, 
     nodeSize = "clr",
     edgeTranspHigh = 20,
     labelScale = FALSE,
     cexNodes = 1, 
     cexLabels = 0.8,
     cexTitle = 0.8,
     groupNames = c("AD", "Healthy"),
     hubBorderCol  = "gray40")
dev.off()


# ==============================================================================
# 5. 提取和处理差异网络边 (Extracting Differential Edges)
# ==============================================================================

# 计算两组关联矩阵的差异
diffedge<-diff$assoMat1-diff$assoMat2 
diffpadj<-diff[["pAdjustMat"]]
# 将上三角矩阵（重复数据）的p值设为1
diffpadj[upper.tri(diffpadj)]<-1

# 筛选显著差异的边（p.adj < 0.05）
diffpadj[diffpadj>0.05]<-1
# 将矩阵转换为长格式
diffedge<-reshape2::melt(diffedge)
colnames(diffedge)<-c("m1","m2","diff")
diffpadj<-reshape2::melt(diffpadj)
colnames(diffpadj)<-c("m1","m2","p")
# 筛选显著的差异边
diffpadj<-diffpadj[diffpadj$p<0.05&diffpadj$m1!=diffpadj$m2,]
diffedge<-diffedge[paste(diffedge$m1,diffedge$m2,sep = ",")%in%paste(diffpadj$m1,diffpadj$m2,sep = ","),]
# 最终的差异网络边数据框
difnetwork<-data.frame(Var1=diffedge$m1,Var2=diffedge$m2,cor=diffedge$diff)

# 提取差异网络中的所有节点
nodes<-unique(c(difnetwork$Var1,difnetwork$Var2))

# 创建节点属性数据框
nodeatr<-data.frame(nodes = nodes,nodetype="difnetwork")
nodeatr<-nodeatr[nodeatr$nodes%in%unique(c(difnetwork$Var1,difnetwork$Var2)),]




# 加载由 DoubleML 确定的重要代谢物，并格式化名称
#importance<-read.csv("output/doubleml_res.csv")
#importance<-gsub("^m_","",importance$X)
#importance<-gsub("_m$","_e",importance) # 可能是指代谢物交换 flux

# 加载 FBA 预测的交换通量数据
flux<-read.csv("D:/microbiome/AD/1118/res_grow/exchanges.csv",header = T)
#flux<-flux[flux$sample_id%in%selected$x,]
flux$metabolite<-gsub("\\[e\\]","_e",flux$metabolite)
doubleml<-read.csv('D:/microbiome/AD/1118/doubleml1.csv')
doubleml$padj<-p.adjust(doubleml$pvalue,method = 'BH')
sig<-doubleml[doubleml$RV>0.15&doubleml$padj<0.01,]
sig$metabolite<-gsub('_m','',sig$metabolite)
# 筛选与重要代谢物相关的通量数据
flux$metabolite<-gsub('_e','',flux$metabolite)
flux$metabolite<-gsub('_m','',flux$metabolite)

flux<-flux[(flux$metabolite)%in%(sig$metabolite),c(1,2,7)]


# 识别与重要代谢物通量相关的物种
imp_species<-unique(flux$taxon)
imp_species<-intersect(imp_species,rownames(ab_genus))
# 更新节点属性：区分“重要物种”和“相关物种”
spedif_nodes<-nodeatr
spedif_nodes$nodetype="Related"
spedif_nodes$nodetype[spedif_nodes$nodes%in%imp_species]="Importance"

# 构建差异网络图对象 (igraph)
dif_graph<-graph_from_data_frame(d=difnetwork,vertices = spedif_nodes,directed = F)

# 使用 Louvain 算法进行社区检测（聚类），识别“公会”
set.seed(444)
cluster<-cluster_louvain(dif_graph)
cluster<-membership(cluster)
table(cluster) # 查看公会大小

# 为节点添加公会（cluster）属性
V(dif_graph)$cluster<-paste("cluster",cluster,sep = "_")
# 特别标记 cluster 3 为 "FERM" (可能是作者的研究重点)
tmp<-V(dif_graph)$cluster
tmp[tmp=="cluster_3"]<-"cluster_3"
V(dif_graph)$cluster<-factor(tmp,levels = unique(tmp))


# ==============================================================================
# 6. 差异网络可视化 (Differential Network Visualization - ggraph)
# ==============================================================================

# 设置网络布局（使用 igraph 的 backbone 布局，并进行手动调整）
# 以下代码通过迭代为每个 cluster 生成一个布局，使相似的 cluster 节点靠近
# (代码较为复杂，此处不逐行注释布局细节，总体目标是创建自定义的聚类布局)
g <- sample_islands(length(table(cluster)), max(table(cluster)), 0.5, length(table(cluster)))
g <- igraph::simplify(g)
set.seed(123)
bb <- layout_as_backbone(g, keep = 0.2)
pos_tmp<-bb$xy
pos<-data.frame()
nodelist<-V(dif_graph)$cluster%>%unique()
spedif_nodes$cluster<-V(dif_graph)$cluster
for(i in 1:length(table(cluster))){
  set.seed(123)
  x.pos<-sample((-7.5+i):(-6.5+i))
  tmp<-pos_tmp[((max(table(cluster))*(i-1))+1):(max(table(cluster))*i),]
  tmp<-tmp[sample(1:max(table(cluster)),sum(spedif_nodes$cluster==nodelist[i])),,drop=F]
  pos<-rbind(pos,tmp)
}
rownames(pos)<-names(V(dif_graph))[order(V(dif_graph)$cluster)]
pos<-pos[names(V(dif_graph)),]

# 使用 ggraph 绘制差异网络图 (图A)
difgraph<-ggraph(dif_graph,
                 layout = "manual",
                 x = pos[, 1],
                 y = pos[, 2]) +
  geom_edge_link0(aes(col = cor), width = 0.2) +
  geom_node_point(aes(fill = cluster), shape = 21, size = 3) +
  geom_mark_hull( 
    aes(x, y, group = cluster, fill = cluster),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  scale_color_manual(values = c("cluster_1"="orange3","cluster_2"="orchid4","cluster_3"="red","cluster_4"="orangered4",
                                "cluster_5"="lightsteelblue4","cluster_6"="linen","cluster_7"="orange","cluster_8"="royalblue4",
                                "cluster_9"="skyblue4","cluster_10"="thistle","cluster_11"="turquoise4","cluster_12"="palegreen4","cluster_13"="mediumseagreen"),
                     labels=c(paste("guild",1:2,sep="_"),"guild_3",paste("guild",4:13,sep="_"))) +
  scale_fill_manual(values = c("cluster_1"="orange3","cluster_2"="orchid4","cluster_3"="red","cluster_4"="orangered4",
                               "cluster_5"="lightsteelblue4","cluster_6"="linen","cluster_7"="orange","cluster_8"="royalblue4",
                               "cluster_9"="skyblue4","cluster_10"="thistle","cluster_11"="turquoise4","cluster_12"="palegreen4","cluster_13"="mediumseagreen"),
                    labels=c(paste("guild",1:2,sep="_"),"guild_3",paste("guild",4:13,sep="_"))) +
  scale_edge_colour_gradientn(name = "difcor",
                              colors = c("aquamarine3","grey90","mediumpurple"),
                              #limits = c(0, 0.28),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "none")+
  guides(fill = guide_legend(ncol = 1))+
  theme_graph(base_family = "sans")+theme(legend.box = 'horizontal',
                                          legend.box.just = 'top')+ggtitle("")+
  theme(plot.title = element_text(hjust = -0.02,vjust = 4.5),title = element_text(size = 16,face="bold"))

difgraph

pdf('D:/microbiome/AD/1118/difgraph.pdf',width = 6,height = 6)
difgraph
dev.off()
# ==============================================================================
# 7. GSEA 分析和代谢关联热图 (GSEA and Metabolic Association)
# ==============================================================================

# 重新加载和处理重要代谢物列表
#importance<-read.csv("output/doubleml_res.csv")
#importance<-gsub("^m_","",unique(importance$X))
#importance<-gsub("_m$","_e",importance)
## 重新加载和处理通量数据
#flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv",header = T)
#meta<-read.csv("data/test_cohort/test_16s/SraRunTable.csv",row.names = 1)
#selected<-read.csv("data/res_grow_wd_pFBA/selected_sample.csv")
#meta<-meta[meta$Sample.Name%in%selected$x,]
#flux<-flux[flux$sample_id%in%meta$Sample.Name,]
#flux$metabolite<-gsub("\\[e\\]","_e",flux$metabolite)
#flux<-flux[flux$metabolite%in%importance,c(1,2,7)]
# 计算物种与重要代谢物通量（被检测到的样本数）的关联，作为排序列表 (Ranked List)
list<-sort(table(flux$taxon),decreasing = T) 

# GSEA 准备：将 Louvain 聚类结果作为基因集 (Species Set)
species_cluster<-cluster
species_cluster<-as.data.frame(species_cluster)
species_cluster$spcies<-rownames(species_cluster)
species_cluster$x<-paste("cluster",species_cluster$x,sep = "_")
species_cluster<-split(as.data.table(species_cluster),by=c("x"))
species_cluster<-lapply(species_cluster,function(x){x[[2]]})
# 执行 GSEA
es<-fgsea(species_cluster, list,nperm=1000)
# 提取 GSEA 结果中最显著的 cluster (pval 最小) 的节点
myRankedlist<-as.numeric(list)
names(myRankedlist)<-names(list)
es
cluster_1nodes<-intersect(names(cluster)[cluster==which.min(es$pval)],names(myRankedlist))
write.csv(cluster_1nodes,"D:/microbiome/AD/1118/cluster_1.csv",row.names = F) # 保存显著 cluster 的节点列表

# 数据准备：用于绘制物种-代谢物关联热图和物种丰度柱状图

# 重新加载样本元数据并分组
#meta<-read.table("data/test_cohort/test_16s/sampleda.txt",row.names = 1)
#meta<-meta[meta$sample%in%selected$x,]
flux_h<-flux[flux$sample_id%in%meta$sample[meta$group=="Healthy"],]
flux_d<-flux[flux$sample_id%in%meta$sample[meta$group=="AD"],]
# 计算每个物种-代谢物对在健康组和高血压组中被检测到的样本数
flux_h$sample_id<-1
flux_d$sample_id<-1
aggreate_flux_h<-aggregate(sample_id~taxon+metabolite,data=flux_h,sum)
aggreate_flux_d<-aggregate(sample_id~taxon+metabolite,data=flux_d,sum)
aggreate_flux_h$metabolite<-paste(aggreate_flux_h$metabolite,"Healthy",sep = "_")
aggreate_flux_d$metabolite<-paste(aggreate_flux_d$metabolite,"AD",sep = "_")
aggreate_flux<-rbind(aggreate_flux_h,aggreate_flux_d)

agg_mat<-as.data.frame(pivot_wider(aggreate_flux,names_from = taxon,values_from = sample_id))
agg_mat[is.na(agg_mat)]<-0
rownames(agg_mat)<-agg_mat$metabolite
agg_mat<-agg_mat[,-1]
agg_mat<-agg_mat[,order(colSums(agg_mat),decreasing = T)]
agg_mat<-agg_mat[order(rowSums(agg_mat)),]
rownames(agg_mat)<-gsub("_e","",rownames(agg_mat))

# agg_mat1：用于绘制热图的数据
agg_mat1 <- reshape2::melt(as.matrix(agg_mat)) 
colnames(agg_mat1)[3]<-"Nsample"

# 绘制代谢物-物种关联热图 (gheat)
custom_palette <- colorRampPalette(c("grey94", "#FEE08B", "#FDAE61", "orangered3"))
set.seed(444)
cluster<-cluster_louvain(dif_graph)
#cluster<-cluster_walktrap(dif_graph)
subg_node<-names(membership(cluster))[membership(cluster)==1]
subg_node<-intersect(subg_node,colnames(agg_mat))

selected_s<-subg_node
label_data <- agg_mat1 %>% 
  filter(Var2 %in% selected_s) %>% 
  distinct(Var2)
gheat<-ggplot(agg_mat1, aes(x = Var2, y = Var1)) +  
  geom_tile(aes(fill = Nsample), colour = "white") +  
  scale_fill_gradientn(name = "Detected\nsample count",colors = custom_palette(100)) +  
  theme_minimal() +  
  labs(y = "Important Features", x = "")+
  theme(axis.text.x = element_blank(),#axis.text.x = element_text(angle = 45,hjust = 1,size = 4),
        axis.title.y = element_text(face = "bold",size = 12),
        axis.text.y = element_text(size=9,face = "bold"),
        legend.key.size = unit(0.3, "cm")) +coord_cartesian(clip = "off") +
  theme(plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))+
  geom_text_repel(
    data = label_data,
    aes(x = Var2, y = 0.5, label = Var2),
    nudge_y = -3,  # 调整标签位置
    direction = "x",
    angle = 45,color="black",vjust=0.1,facet="bold",
    hjust = 0.6,size=3.5,box.padding = 0.1,parse = T,segment.alpha=0.2
  )+
  expand_limits(y = c(-6, 1.5))+
  #scale_y_discrete(breaks = selected_s,labels=selected_s) + 
  theme(plot.margin = ggplot2::margin(t=-2,b=1,r=3),
        title = element_text(size = 12,face="bold") )
gheat
# 绘制不同组的物种总和样本数柱状图 (bar_species)
# ... (柱状图绘制代码) ...
aggreate_flux_h<-reshape2::melt(as.matrix(agg_mat[grep("Healthy",rownames(agg_mat)),]))
colnames(aggreate_flux_h)<-c("metabolites","taxon","count")
aggreate_flux_h<-aggregate(count~taxon,data=aggreate_flux_h,sum)
aggreate_flux_d<-reshape2::melt(as.matrix(agg_mat[grep("AD",rownames(agg_mat)),]))
colnames(aggreate_flux_d)<-c("metabolites","taxon","count")
aggreate_flux_d<-aggregate(count~taxon,data=aggreate_flux_d,sum)

aggreate_flux_h$group<-"Healthy"
aggreate_flux_d$group<-"AD"
species<-rbind(aggreate_flux_h,aggreate_flux_d)
colnames(species)[2]<-"counts"
species$taxon<-factor(species$taxon,levels = colnames(agg_mat))
bar_species<-ggplot()+
  geom_bar(data = species[species$group=="Healthy",],aes(taxon,counts,fill=group),
           stat = "identity", color=NA,alpha=.6)+
  geom_bar(data = species[species$group=="AD",],aes(taxon,counts,fill=group),
           stat = "identity", color=NA,alpha=.6)+
  scale_fill_manual(values = c('#E7B800','#00AFBB'))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        # axis.text.x = element_text(face = "italic", color="black"),
        #axis.line = element_line(colour = "black"),
        plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm"),
        axis.text.y = element_text(size=10, color="black")) +
  # scale_x_discrete(position = "top")+
  theme(axis.title.x = element_blank(),  
        title = element_text(size = 12,face = "bold"),
        legend.key.size = unit(0.5, "cm"),legend.text = element_text(size=10))+
  #ylim(c(0,900))+
  #geom_text(aes(label = freq),vjust = -0.5, size = 4.5, color = "black")+
  xlab("")+ylab("Sum of samples")

bar_species



# GSEA 富集图 (gseap)
names(species_cluster)[1]<-"Cluster_1" # 重命名 cluster 3 为 FERM
# 准备 GSEA 输入格式
term2gene <- stack(species_cluster) 
colnames(term2gene) <- c("feature", "set")
term2gene<-term2gene[,c(2,1)]
keep<-names(table(term2gene$set))[table(term2gene$set)>=5]
term2gene$set<-as.character(term2gene$set)
term2gene<-term2gene[term2gene$set%in%keep,]
term2gene$set<-factor(term2gene$set,levels = unique(term2gene$set))
gsea_res <- GSEA(myRankedlist, 
                 TERM2GENE = term2gene,
                 minGSSize = 0,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "fdr",
                 seed = 123
)
gsea_df<-as.data.frame(gsea_res@result)




#write.table(gsea_df,"output/gsea_df.txt",row.names = F,quote = F,sep = "\t")
source("D:/microbiome/AD/script/gseaplot3.R") # 依赖一个自定义的绘图函数
library(ggfun)
library(grid)
gseap<-gseaplot3(gsea_res,geneSetID = "Cluster_1",base_size=10,pvalue_table = T)
gseap
# 组合 GSEA 和代谢关联图（图B）
########
ATE<-read.csv("D:/microbiome/AD/1118/doubleml1.csv")
ATE$padj<-p.adjust(ATE$pvalue,method = "BH")
sign<-sign(ATE[,c("ATE","CI2.5","CI97.5","theta.lower","theta.upper")])
#keep<-apply(sign, 1, function(x){
#  length(unique(x))==1
#})
sig<-ATE[ATE$RV>0.15&ATE$padj<0.01,]
sig$color<-ifelse(sig$ATE>0,"H","A")
sig$metabolite<-factor(sig$metabolite,levels = sig$metabolite[order(sig$ATE)])
sig[order(sig$ATE),]->sig

#sig$class<-c("Health","Health","Health","Health","ASD","ASD","ASD","ASD")


sig$x=ifelse(sig$ATE<0,0.05,min(sig$CI2.5))
sig$padj1<-ifelse(sig$padj<0.001,"<0.001",round(sig$padj,3))
sig$label<-paste(sprintf("%.3f", sig$RV),sig$padj1,sep = "      ")
sig<-sig[order(sig$padj,decreasing = F)[1:15],]
colnames(sig)
p1=ggplot(data = sig, mapping = aes(y = reorder(metabolite, ATE))) +
  geom_text(aes(x = 0, label = label), hjust = 0) +
  xlim(-0.01, 0.05) +
  xlab("") +
  ylab("") +
  
  # 添加 RV 标签
  #annotate("text", 
  #         x = 0.004, # 调整 RV 的水平位置 (x轴)
  #         y = max(as.numeric(factor(sig$metabolite))) + 0.5, # 调整 RV 的垂直位置 (y轴) - 放在最上方
  #         label = "RV", 
  #         size = 4, # 调整字体大小
  #         fontface = "bold") +
  #
  ### 添加 Padj 标签
  #annotate("text", 
#       x = 0.022, # 调整 Padj 的水平位置 (x轴)
#       y = max(as.numeric(factor(sig$metabolite))) + 0.5, # 调整 Padj 的垂直位置 (y轴) - 放在最上方
#       label = "Padj", 
#       size = 4, # 调整字体大小
#       fontface = "bold") +

# 移除原始的 subtitle 设置，使用 title
labs(title = "") +
  
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    
    # 调整 plot.margin 以适应新的顶部标签
    #plot.margin = ggplot2::margin(r = -80, t = 50), # 增加顶部边距 (t)
    
    # 移除 plot.subtitle 设置
    # plot.subtitle = element_text(hjust = 0.92, size = 12, face="bold"),
    
    plot.title.position = "plot",
    plot.title = element_text(size = 14, vjust = 2, face = "bold"),
    
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
  )
p1

p2=ggplot(data=sig, aes(x = ATE, y = reorder(metabolite,ATE))) +
  #annotate("rect",xmin=-Inf,xmax=Inf,ymin=0,ymax=14.5,fill="green4",alpha=0.05)+
  #annotate("rect",xmin=-Inf,xmax=Inf,ymin=14.5,ymax=25,fill="red4",alpha=0.05)+
  geom_point(data=sig, aes(x = ATE, y = metabolite, color = color),size=3) +
  scale_color_manual(values = c("blue4","red4"))+
  geom_errorbar(data=sig,mapping = aes(xmin = theta.lower, xmax = theta.upper,y=metabolite),color="black",width=0.3)  +
  geom_vline(xintercept = 0, color = "grey") +
  theme_test()+ ylab("") +
  xlab("ATE") +
  theme(legend.position = "none",
        text = element_text(size = 16),axis.title.x = element_text(size = 14,face = "bold"),
        axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        plot.margin = ggplot2::margin(10, 10, 10, -30),
        strip.text = element_text(face = "bold",size = 16),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA))
#geom_rect(data=sig[sig$ATE>0,],mapping = aes(xmin=x-0.01,xmax=x+0.04,y=metabolite))
#geom_text(data=sig[sig$ATE>0,],mapping=aes(x=x,y=metabolite,label=label),hjust=0)+
#geom_text(data=sig[sig$ATE<0,],mapping=aes(x=x,y=metabolite,label=label),hjust=0)+
#facet_grid(class~., scale='free',space = 'free')
p_double<-wrap_plots(p1,p2,nrow = 1,widths = c(0.5,0.5))
pdf('D:/microbiome/AD/1118/doubleml_top15.pdf',width = 6,height = 6)
p_double
dev.off()
########
# ==============================================================================
# 8. 绘制 FERM 子网络图 (FERM Subnetwork Visualization)
# ==============================================================================

# 提取 cluster 2  内部的差异边和节点
cluster1nodes<-cluster_1nodes
clu1net<-difnetwork[as.character(difnetwork$Var1)%in%cluster1nodes&as.character(difnetwork$Var2)%in%cluster1nodes,]
clus1no<-spedif_nodes[as.character(spedif_nodes$nodes)%in%as.character(cluster1nodes),]
c1_graph<-graph_from_data_frame(d=clu1net,vertices = clus1no)
V(c1_graph)$cluster<-factor(V(c1_graph)$cluster)
c1graph<-ggraph(c1_graph,
                layout = "manual",
                x = pos[names(V(c1_graph)), 1],
                y = pos[names(V(c1_graph)), 2]) +
  geom_edge_link0(aes(col = cor), width = 0.2) +
  geom_node_point(aes(fill = cluster), shape = 21, size = 3) +
  geom_mark_hull( # 注意这里哈
    aes(x, y, group = cluster, fill = cluster),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  scale_color_manual(values = c("orange3")) +
  scale_fill_manual(values = c("orange3")) +
  scale_edge_colour_gradientn(name = "difcor",
                              colors = c("grey90","skyblue1","mediumpurple"),
                              limits = c(0, 0.28),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "none")+
  guides(fill = guide_legend(ncol = 1))+
  #theme_graph(base_family = "sans")+theme(legend.position = "node")+
  scale_x_continuous(expand = expansion(add = 1))+
  scale_y_continuous(expand = expansion(add = 1))+
  theme(
    panel.background = element_rect(fill = "transparent"), # 设置绘图区背景为透明
    plot.background = element_rect(fill = "transparent", color = NA), # 设置整个图形背景为透明
    panel.grid.major = element_blank(), # 移除主要网格线
    panel.grid.minor = element_blank(), # 移除次要网格线
    legend.position = "none"
  )
c1graph
gsea<-wrap_plots(list(gseap[[1]],gseap[[2]],bar_species,gheat),ncol = 1,heights = c(0.38,0.05,0.25,1.2))
gsea<-ggdraw(gsea) +ggtitle("")+
  theme(plot.title = element_text(hjust = 0.01,vjust = 0),
        plot.margin = ggplot2::margin(t=15),
        title = element_text(size = 18,face="bold"))+
  draw_plot(c1graph,x = 0.3,y=0.8,width = 0.2,height = 0.12)
gsea
gseap[[1]]<-gseap[[1]]+theme(axis.title.y = element_text(size=12,face = "bold"),
                             axis.text.y = element_text(size=10))
gsea

pdf('D:/microbiome/AD/1118/gsea.pdf',width = 8,height = 17)
gsea
dev.off()

##############
manifest<-read.csv("D:/microbiome/AD/1118/manifest.csv",row.names = 1)
meta<-read.csv('D:/microbiome/AD/1118/meta.csv',row.names = 1)
str(manifest)
str(meta)

library(dplyr)
library(ggplot2)

# 1. 合并 manifest 与 meta，并按组别 + found_fraction 降序排列
plot_df <- manifest %>%
  mutate(Run.ID = sub("\\.pickle$", "", file)) %>%
  inner_join(meta, by = "Run.ID") %>%
  arrange(group, desc(found_fraction), Run.ID) %>%
  mutate(
    sample_idx = row_number(),
    # 若想固定组别顺序（例如先 Healthy 再 AD），可按需调整 levels
    group = factor(group, levels = c("Healthy", "AD"))
  )

# 2. 绘图
p <- ggplot(plot_df, aes(x = sample_idx, y = found_fraction, fill = group, color = group)) +
  geom_col(width = 1, linewidth = 0.1) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "black", linewidth = 0.7) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  scale_fill_manual(values = c("Healthy" = "#E8CE5A", "AD" = "#66CFD6")) +
  scale_color_manual(values = c("Healthy" = "#E8CE5A", "AD" = "#66CFD6")) +
  labs(
    x = "Samples",
    y = "AGORA Genus Percentage",
    fill = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal"
  )

print(p)
pdf('D:/microbiome/AD/1118/count.pdf',width = 6,height = 6)
p
dev.off()
# ==============================================================================
# 9. FERM 子网络详细分析和可视化 (FERM Detailed Analysis and Plot C)
# ==============================================================================

# 筛选 FERM 节点
subg_node<-names(membership(cluster))[membership(cluster)==2]
subg_node<-intersect(subg_node,colnames(agg_mat))
dif_subgraph<-difnetwork[difnetwork$Var1%in%subg_node&difnetwork$Var2%in%subg_node,]
subg_nodeattr<-nodeatr[nodeatr$nodes%in%subg_node,]
subg_nodeattr$Nsample<-myRankedlist[match(subg_nodeattr$nodes,names(myRankedlist))]
# 识别节点是否为已知肠道疾病相关物种（Documented）
known<-read.csv("D:/microbiome/AD/data/gutmdisorder.csv")
known<-unique(gsub(" .*","",known$Gut.Microbiota..ID.))
subg_nodeattr$nodetype<-ifelse(as.character(subg_nodeattr$nodes)%in%known,"Documented","Not Documented")
# 加载 TCG 分型数据（未使用）
# ... (TCG 数据加载和处理代码) ...

# 加载并处理差异丰度分析结果 (genus_dif.csv)
genus_dif<-read.csv("D:/microbiome/AD/diff_genus.csv",row.names = 1)
genus_dif<-na.omit(genus_dif)
genus_h<-genus_dif[genus_dif$Sign_group=="Healthy",]
genus_d<-genus_dif[genus_dif$Sign_group=="AD",]
rownames(genus_h)<-genus_h$genus
genus_h$sign<-case_when(genus_h$fdr<=0.001~"***",
                        genus_h$fdr<=0.01&genus_h$fdr>0.001~"**",
                        genus_h$fdr<=0.05&genus_h$fdr>0.01~"*",
                        .default = "")
rownames(genus_d)<-genus_d$genus
genus_d$sign<-case_when(genus_d$fdr<=0.001~"***",
                        genus_d$fdr<=0.01&genus_d$fdr>0.001~"**",
                        genus_d$fdr<=0.05&genus_d$fdr>0.01~"*",
                        .default = "")
sign_all<-rbind(genus_d,genus_h)
sign_all$nodes<-paste0(rownames(sign_all),sign_all$sign)
index<-which(subg_nodeattr$nodes%in%rownames(sign_all))
subg_nodeattr$nodes<-as.character(subg_nodeattr$nodes)
subg_nodeattr$nodes[index]<-sign_all$nodes[match(subg_nodeattr$nodes[index],rownames(sign_all))]

# 更新 FERM 节点的属性：丰度变化趋势 (Up/Down/Not)
subg_nodeattr$change<-ifelse(subg_nodeattr$nodes%in%paste0(rownames(genus_h),genus_h$sign),"Down",ifelse(subg_nodeattr$nodes%in%paste0(rownames(genus_d),genus_d$sign),"Up","Not"))
subg_nodeattr<-subg_nodeattr[order(subg_nodeattr$change),]
subg_nodeattr$change<-factor(subg_nodeattr$change,levels = c("Down","Up","Not"))

# 更新差异网络的节点名称（包含星号标记）
dif_subgraph$Var1<-as.character(dif_subgraph$Var1)
dif_subgraph$Var2<-as.character(dif_subgraph$Var2)
index<-which(dif_subgraph$Var1%in%rownames(sign_all))
dif_subgraph$Var1[index]<-sign_all$nodes[match(dif_subgraph$Var1[index],rownames(sign_all))]
index<-which(dif_subgraph$Var2%in%rownames(sign_all))
dif_subgraph$Var2[index]<-sign_all$nodes[match(dif_subgraph$Var2[index],rownames(sign_all))]

# 构建包含差异丰度信息的 FERM 子网络图对象
subgraph<-graph_from_data_frame(d=dif_subgraph,vertices = subg_nodeattr)
V(subgraph)$nodes<-subg_nodeattr$nodes
V(subgraph)$change<-factor(V(subgraph)$change,levels = c("Down","Up","Not"))

# 绘制圆形布局的 FERM 子网络图 (subgraph_p，图C)
legend_data <- data.frame(
  change = c("Down", "Up"),
  label = c("Down", "Up")
)
subgraph_p<-ggraph(subgraph, layout='linear', circular = TRUE) +
  geom_node_text(aes(x = 1.1 * x,
                     y = 1.1 * y,
                     label = nodes,color=change,
                     angle = -((-node_angle(x, y) + 90) %% 180) + 90),
                 hjust = "outward") +
  scale_color_manual(name="change",values=c("Up"="red","Down"="palegreen4","Not"="black"),
                     breaks = c("Up","Down"),
                     labels=c(expression(HT[Up]),expression(HT[Down])))+
  ggnewscale::new_scale_color()+
  geom_node_point(size=3,alpha=1,aes(color=nodetype))+
  geom_edge_arc(aes(colour = cor), width=0.1,
                alpha = 1.5) +
  scale_color_manual(name="node type",values=c("Documented"="red","Not Documented"="darkcyan"))+
  scale_edge_width_continuous(range = c(0,0.2)) +theme_void()+
  #theme_graph()+
  scale_edge_colour_gradientn(name = "Differential \ncorrelation",
                              colors = c("aquamarine3","grey90","mediumpurple"),
                              #limits = c(0, max(E(graph)$cor)),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "edge_colourbar")+
  expand_limits(x = c(-3.5, 2.5), y = c(-3, 3))+theme(plot.margin = ggplot2::margin(r=50))+
  ggtitle("C.")+labs(subtitle = "Cluster_2")+
  theme(plot.title = element_text(hjust = 0.06,vjust = -6),
        plot.subtitle = element_text(hjust = 0.5),
        plot.margin = ggplot2::margin(l=0,r=5,t=5,b=5),
        title = element_text(size = 16,face="bold"))
subgraph
data_prob<-data.frame(class=c("Cluster_2","other guilds"),n=c(sum(V(subgraph)$change!="Not"),nrow(sign_all)-sum(V(subgraph)$change!="Not")))
data_prob$prop<-data_prob$n/sum(data_prob$n)
data_prob <- data_prob %>%
  arrange(desc(class)) %>% # 重排序
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
data_prob$per<-paste0(round(data_prob$prop,3)*100,"%")
mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")

per_p<-ggplot(data_prob, aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = lab.ypos, label = per), color = "black",size=8)+
  scale_fill_manual(values = mycols) +
  theme_void()+labs(subtitle="Proportion of Differential Species")+
  ggtitle("D.")+
  theme(plot.title = element_text(hjust = -0.35,vjust = -2),
        plot.subtitle = element_text(hjust = 0.5,vjust=-5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 16,face="bold"))

spe_net_plot<-gridExtra::grid.arrange(difgraph,subgraph_p,ncol=1,heights=c(0.4,0.5))

