rm(list = ls())
library(philr)
library(data.table)
library(dplyr)
library(ape)
library(dplyr)
library(stringr)
library(tibble)
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(ggtree)
library(tidyverse)
library(ggClusterNet)
taxon<-fread('D:/microbiome/AD/data/bac120_taxonomy_r226.tsv',header = F)%>%as.data.frame()
#############
abundance_D<-read.csv('D:/microbiome/AD/data/AD/genus/abundance_genus.csv',row.names = 1)
colSums(abundance_D)%>%sort()%>%plot
colSums(abundance_D)%>%hist
temp<-colSums(abundance_D)
keep<-apply(abundance_D, 1, function(x){
  mean(x!=0)>0.05
})
abundance_D<-abundance_D[keep,]


#for(i in 1:ncol(abundance_D)){
#  abundance_D[,i]<-abundance_D[,i]/temp[i]
#}
#files<-list.files('D:/microbiome/micom/agora103_genus/890088c4-1a67-498b-8a54-3eb31ee06f63/data')
#for(i in 1:length(files)){
#  files[i]<-substr(files[i],1,nchar(files[i])-5)
#}
#abundance_D1<-abundance_D[intersect(rownames(abundance_D),files),]
#colSums(abundance_D1)%>%hist
#abundance_D1<-abundance_D1[,colSums(abundance_D1)>0.8]
#abundance_D<-abundance_D[,colnames(abundance_D1)]
#temp_D<-abundance_D1
#temp_D[temp_D!=0]<-1
#colSums(temp_D)%>%hist
abundance_H<-fread('D:/microbiome/AD/data/Healthy/genus/abundance_genus.csv')%>%as.data.frame()%>%column_to_rownames('V1')
colSums(abundance_H)%>%sort()%>%plot
colSums(abundance_H)%>%hist
abundance_H<-abundance_H[,colSums(abundance_H)<=100]
#abundance_H1<-abundance_H[intersect(rownames(abundance_H),files),]
#colSums(abundance_H1)%>%hist()
abundance_H<-abundance_H[,colSums(abundance_H)>90&colSums(abundance_H)<97]
abundance_H<-abundance_H[which(rowSums(abundance_H)!=0),]
keep<-apply(abundance_H, 1, function(x){
  mean(x!=0)>0.05
})
abundance_H<-abundance_H[keep,]

#temp_H<-abundance_H
#temp_H[temp_H!=0]<-1
#colSums(temp_H)%>%hist
#abundance_H<-abundance_H[,colnames(abundance_H1)]
data<-merge(abundance_D,abundance_H,by = "row.names",all = TRUE)
data[is.na(data)]<-0
rownames(data)<-data$Row.names
data$Row.names<-NULL
abundance<-data

write.csv(abundance,'D:/microbiome/AD/data/abundance_all_filtered.csv')
#rownames(abundance)[which(rownames(abundance)=='Erysipelatoclostridium')]<-'Clostridium'
#rownames(abundance)[which(rownames(abundance)=='Ruminiclostridium_E(f__Ruminococcaceae)')]<-'Ruminiclostridium'



#abundance<-abundance[-grep('unclassified',rownames(abundance)),,drop=F]
############
meta<-read.csv('D:/microbiome/AD/data/meta_958.csv',row.names = 1)
abundance<-abundance[,colnames(abundance)%in%rownames(meta)]
temp<-abundance
temp[temp!=0]<-1
pheatmap::pheatmap(temp,annotation_col = meta,
                   show_rownames = F,
                   show_colnames = F,
                   cluster_rows = F,
                   cluster_cols = F)

gtdb_tree<-read.tree('D:/microbiome/AD/data/bac120_r226.tree')
temp<-colSums(abundance)
temp<-temp[temp>0]
abundance<-abundance[,names(temp)]
for(i in 1:ncol(abundance)){
  abundance[,i]<-abundance[,i]/temp[i]
}
temp<-strsplit(taxon$V1,'\t')%>%unlist()
temp<-temp[-grep('d__',temp)]
rownames(taxon)<-temp
taxon$V1<-sub(
  pattern = ".*\t", 
  replacement = "", 
  x = taxon$V1
)
genus<-paste0('g__',rownames(abundance))
genus_no<-setdiff(genus,taxon$V6)
ncbi_gtdb<-read.csv('D:/microbiome/AD/ncbi_vs_gtdb_genus.csv')
genus_no<-genus_no[genus_no%in%ncbi_gtdb$NCBI_genus]
ncbi_gtdb$GTDB_genus<-sub("\\s.*", "", ncbi_gtdb$GTDB_genus)
ncbi_gtdb$GTDB_genus[match(genus_no,ncbi_gtdb$NCBI_genus)]
genus_no<-cbind(genus_no,ncbi_gtdb$GTDB_genus[match(genus_no,ncbi_gtdb$NCBI_genus)])%>%as.data.frame()
genus_no$genus_no<-gsub('g__','',genus_no$genus_no)
genus_no$V2<-gsub('g__','',genus_no$V2)
abundance$genus<-rownames(abundance)
abundance$genus[which(abundance$genus=='Erysipelatoclostridium')]<-'Clostridium'
abundance$genus[which(abundance$genus=='Ruminiclostridium')]<-'Ruminiclostridium_E(f__Ruminococcaceae)'
abundance$genus[which(abundance$genus=='Gemmiger')]<-'Gemmiger(f__Ruminococcaceae)'
abundance$genus[which(abundance$genus=='Anaerosinus')]<-'Anaerovibrio'




abundance$genus[match(genus_no$genus_no,abundance$genus)]<-genus_no$V2
setDT(abundance)
exp <- abundance[, lapply(.SD, function(x) sum(x, na.rm = TRUE)), 
             by = 'genus', 
             .SDcols = which(sapply(abundance, is.numeric))]
exp<-as.data.frame(exp)
exp<-exp%>%column_to_rownames('genus')
abundance<-exp
taxon<-taxon[taxon$V6%in%genus,]
genus_used<-unique(taxon$V6)
genus_used<-gsub('g__','',genus_used)
abundance_used<-abundance[genus_used,]
abundance_used<-na.omit(abundance_used)
rownames(abundance_used)<-paste0('g__',rownames(abundance_used))
setdiff(paste0('g__',rownames(abundance)),rownames(abundance_used))


(colSums(abundance_used)/colSums(abundance))%>%sort()%>%plot
(colSums(abundance_used)/colSums(abundance))%>%sort()%>%hist




tree_tip_ids <- gtdb_tree$tip.label
taxon_map <- taxon %>%
  as_tibble(rownames = "accession_id") %>%
  select(accession_id, genus = V6) %>%
  # 步骤 A: 过滤掉不在 GTDB 树中的 Accession ID
  filter(accession_id %in% tree_tip_ids) %>% # <--- 关键的新增筛选
  # 移除 g__ Genus 前缀
  mutate(genus = gsub("^g__", "", genus)) %>%
  # 移除空格/清理
  mutate(genus = gsub(" ", "_", genus))

representative_tips <- taxon_map %>%
  # 按属分组
  group_by(genus) %>%
  # 为每个属选择第一个 Accession ID 作为代表
  slice(1) %>%
  ungroup()
tips_to_keep <- representative_tips$accession_id
tips_to_drop <- gtdb_tree$tip.label[!gtdb_tree$tip.label %in% tips_to_keep]

genus_tree_pruned <- drop.tip(
  gtdb_tree, 
  tips_to_drop
)

rename_map <- representative_tips %>%
  filter(accession_id %in% genus_tree_pruned$tip.label) %>%
  select(accession_id, genus) %>%
  tibble::deframe()

genus_tree_pruned$tip.label <- rename_map[genus_tree_pruned$tip.label]
#######
genera_in_tree <- genus_tree_pruned$tip.label
genera_in_data <- gsub('g__','',rownames(abundance_used))
temp1<-length(intersect(genera_in_tree,genera_in_data))==length(genera_in_tree)&length(intersect(genera_in_tree,genera_in_data))==length(genera_in_data)
if(!temp1){
  common_genera <- intersect(genera_in_tree, genera_in_data)
  genus_tree_final <- drop.tip(
    genus_tree_pruned, 
    genera_in_tree[!genera_in_tree %in% common_genera]
  )
  rownames(abundance_used)<-gsub('g__','',rownames(abundance_used))
  abundance_used <- abundance_used[common_genera, ]
}
genus_tree_pruned<-genus_tree_final
#######
rownames(abundance_used)<-gsub('g__','',rownames(abundance_used))
abundance_matrix <- as.matrix(abundance_used)%>%t()
abundance_matrix[abundance_matrix == 0] <- min(abundance_matrix[abundance_matrix > 0]) / 100
if (!all(colnames(abundance_matrix) == genus_tree_pruned$tip.label)) {
  # 如果不一致，则需要重新排序
  abundance_matrix <- abundance_matrix[,genus_tree_pruned$tip.label]
  print("Warning: Abundance matrix was reordered to match tree tip labels.")
}
Ntip <- Ntip(genus_tree_pruned)
Nnode <- genus_tree_pruned$Nnode
unique_node_indices <- (Ntip + 1) : (Ntip + Nnode)
map_data<-cbind(unique_node_indices,genus_tree_pruned$node.label)%>%as.data.frame()
genus_tree_pruned$node.label <- as.character(unique_node_indices)
abundance_philr <- philr(
  abundance_matrix,
  genus_tree_pruned,
  part_weights = 'blw' # 使用分支长度加权，这是 PhILR 的标准做法
)
#saveRDS(genus_tree_pruned,'D:/microbiome/AD/genus_tree.rds')
#####ggClusterNet####



#########

id<-name.to.nn(genus_tree_pruned,colnames(abundance_philr))
subtree <- extract.clade(genus_tree_pruned, node = id[100])
plot(subtree)
treeplot<-ggtree(genus_tree_pruned, 
                 layout = "circular", # 更改为圆形布局
                 size = 0.5,         # 调整线条粗细
                 color = "darkgrey") +
  
  # 1. 添加叶子标签（属名）
  geom_tiplab(size = 3, 
              color = "black", 
              align = T, # 确保标签对齐
              linesize = 0.3)
treeplot
pdf('D:/microbiome/AD/tree.pdf')
treeplot
dev.off()


##########
library(caret)
abundance_philr<-as.data.frame(abundance_philr)
((rownames(meta)==rownames(abundance_philr))%>%sum())==nrow(meta)
abundance_philr$group<-meta$group
data<-abundance_philr
set.seed(123)
data_split <- createDataPartition(data$group, p = .70, list = FALSE)
training_data <- data[ data_split,]
testing_data  <- data[-data_split,]
train_index <- createFolds(training_data$group, k = 10)

randomForestFit <- training_data %>% train(group ~ .,
                                           method = "rf",
                                           data = .,
                                           tuneLength = 5,
                                           trControl = trainControl(method = "cv", indexOut = train_index))
randomForestFit
pr_rf<-predict(randomForestFit,testing_data)
confusionMatrix(pr_rf,reference = testing_data$group%>%as.factor())
pr_rf<-predict(randomForestFit,testing_data,type = 'prob')
library(pROC)
roc_curve<-roc(testing_data$group,pr_rf[,1])
plot(roc_curve)
sp.obj <- ci.sp(roc_curve, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj, type="shape", col="#5bd1d7")
text(0,0.5,paste0("AUC:",round(roc_curve[["auc"]],4)))
imp<-varImp(randomForestFit)[['importance']]
imp<-imp[order(imp$Overall,decreasing = T),,drop=F]
rownames(imp)<-gsub('`','',rownames(imp))
##########
meta_AD<-fread('D:/microbiome/AD/user_selected_run_list_W4as0en1k1.txt')%>%as.data.frame()%>%column_to_rownames('Run ID')
table(meta_AD$`Project ID`)
#meta_AD<-meta_AD[colnames(abundance_D),]
tr_id<-c(rownames(meta_AD)[which(meta_AD$`Project ID`=='PRJNA633959'|meta_AD$`Project ID`=='PRJNA554111')])
ts_id<-c(rownames(meta_AD)[which(meta_AD$`Project ID`=='PRJNA496408')])
va_id<-c(rownames(meta_AD)[which(meta_AD$`Project ID`=='PRJNA489760'|meta_AD$`Project ID`=='PRJEB11419')])
id_H<-rownames(data)[data$group=='Healthy']
set.seed(123)
id_H<-id_H[sample(1:length(id_H),length(id_H))]
tr_id_H<-id_H[1:503]
ts_id_H<-id_H[504:620]
va_id_H<-id_H[621:length(id_H)]

training_data<-data[match(c(tr_id,tr_id_H),rownames(data)),]
training_data<-na.omit(training_data)
testing_data<-data[c(ts_id,ts_id_H),]
testing_data<-na.omit(testing_data)
validation_data<-data[c(va_id,va_id_H),]
validation_data<-na.omit(validation_data)
train_index <- createFolds(training_data$group, k = 10)
randomForestFit <- training_data %>% train(group ~ .,
                                           method = "rf",
                                           data = .,
                                           tuneLength = 5,
                                           trControl = trainControl(method = "cv", indexOut = train_index))
randomForestFit










pr_rf_tr<-predict(randomForestFit,training_data,type = 'prob')
roc_curve<-roc(training_data$group,pr_rf_tr[,1])
auc(roc_curve)
roc_data_tr <- data.frame(
  TPR = rev(roc_curve$sensitivities), # 真阳性率
  FPR = rev(1 - roc_curve$specificities) # 假阳性率
)
roc_data_tr$label<-'Training Set, AUC = 1'


pr_rf_ts<-predict(randomForestFit,testing_data,type = 'prob')
pr_rf_ts$label<-'Testing Set, AUC = 1'
roc_curve<-roc(testing_data$group,pr_rf_ts[,1])
auc(roc_curve)
roc_data_ts <- data.frame(
  TPR = rev(roc_curve$sensitivities), # 真阳性率
  FPR = rev(1 - roc_curve$specificities) # 假阳性率
)
roc_data_ts$label<-'Testing Set, AUC = 1'
plot(roc_curve)
sp.obj <- ci.sp(roc_curve, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj, type="shape", col="#5bd1d7")
text(0,0.5,paste0("AUC:",round(roc_curve[["auc"]],4)))

pr_rf_va<-predict(randomForestFit,validation_data,type = 'prob')
pr_rf_va$label<-'Validation Set, AUC = 0.9996'
roc_curve<-roc(validation_data$group,pr_rf_va[,1])
auc(roc_curve)
roc_data_va <- data.frame(
  TPR = rev(roc_curve$sensitivities), # 真阳性率
  FPR = rev(1 - roc_curve$specificities) # 假阳性率
)
roc_data_va$label<-'Validation Set, AUC = 0.9996'
plot(roc_curve)
sp.obj <- ci.sp(roc_curve, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj, type="shape", col="#5bd1d7")
text(0,0.5,paste0("AUC:",round(roc_curve[["auc"]],4)))

df<-rbind(roc_data_tr,roc_data_ts,roc_data_va)
df$label<-factor(df$label,levels = c('Training Set, AUC = 1',
                                     'Testing Set, AUC = 1',
                                     'Validation Set, AUC = 0.9996'))
ggplot(df, aes(x = FPR, y = TPR, color = label)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(
    title = "ROC",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(color = guide_legend(title = NULL))

imp<-varImp(randomForestFit)[['importance']]
imp<-imp[order(imp$Overall,decreasing = T),,drop=F]
rownames(imp)<-gsub('`','',rownames(imp))



#########
library(ggtree)
library(ape)
library(phangorn)
n<-15
res<-list()
g_mean <- function(x) exp(mean(log(x)))
for(i in 1:n){#假设nrow(imp)>n
  subtree <- extract.clade(genus_tree_pruned, node = rownames(imp)[i]%>%as.numeric())
  #treeplot<-ggtree(subtree, 
  #                 #layout = "circular", # 更改为圆形布局
  #                 size = 0.5,         # 调整线条粗细
  #                 color = "darkgrey") +
  #  
  #  # 1. 添加叶子标签（属名）
  #  geom_tiplab(size = 4, 
  #              color = "black", 
  #              align = TRUE, # 确保标签对齐
  #              linesize = 0.2)
  #print(treeplot)
  current_node_name<-rownames(imp)[i]
  N <- length(subtree$tip.label)
  root_node <- N + 1
  child_edges <- which(subtree$edge[, 1] == root_node)
  first_layer_nodes <- subtree$edge[child_edges, 2]
  
  if (length(first_layer_nodes) != 2) {
    stop("错误: 树的根节点不是二叉的，或结构不清晰。")
  }
  
  node1 <- first_layer_nodes[1]
  node2 <- first_layer_nodes[2]
  
  # 3. 获取第一个分支下的所有叶子节点名称
  tips_in_clade1_indices <- Descendants(subtree, node1, type = "tips")
  species_in_clade1 <- subtree$tip.label[tips_in_clade1_indices[[1]]]
  
  # 4. 获取第二个分支下的所有叶子节点名称
  tips_in_clade2_indices <- Descendants(subtree, node2, type = "tips")
  species_in_clade2 <- subtree$tip.label[tips_in_clade2_indices[[1]]]
  species<-list(species_in_clade1,species_in_clade2)
  names(species)<-c('clade1','clade2')
  res[[current_node_name]]<-species
  test_val <- abundance_philr[1, current_node_name]
  test_ratio <- log(g_mean(abundance_matrix[1, species_in_clade1]) / 
                      g_mean(abundance_matrix[1, species_in_clade2]))
  
  if (sign(test_val) == sign(test_ratio)) {
    res[[current_node_name]]$formula <- "log(clade1/clade2)"
  } else {
    res[[current_node_name]]$formula <- "log(clade2/clade1)"
  }
}


diff_species <- read.csv('D:/microbiome/AD/diff_genus.csv',row.names = 1) 
diff_species <- diff_species$genus
# 创建一个数据框来存储结果
species_node_mapping <- data.frame(
  Species = diff_species,
  Node = NA,
  Position = NA, # 标记是在分母还是分子
  stringsAsFactors = FALSE
)

# 遍历每一个物种
for (i in 1:nrow(species_node_mapping)) {
  current_sp <- species_node_mapping$Species[i]
  
  found_nodes <- c()
  positions <- c()
  
  # 遍历 res 中的每一个 Balance 节点
  for (node_name in names(res)) {
    clade1 <- res[[node_name]]$clade1
    clade2 <- res[[node_name]]$clade2
    formula <- res[[node_name]]$formula
    
    # 检查物种是否在 Clade1 中
    if (current_sp %in% clade1) {
      found_nodes <- c(found_nodes, node_name)
      # 根据公式判断是分子还是分母
      pos <- if(formula == "log(clade1/clade2)") "Numerator (分子)" else "Denominator (分母)"
      positions <- c(positions, pos)
    } 
    # 检查物种是否在 Clade2 中
    else if (current_sp %in% clade2) {
      found_nodes <- c(found_nodes, node_name)
      pos <- if(formula == "log(clade2/clade1)") "Numerator (分子)" else "Denominator (分母)"
      positions <- c(positions, pos)
    }
  }
  
  # 如果找到了匹配的节点，填入信息（一个物种可能属于多个祖先节点）
  if (length(found_nodes) > 0) {
    species_node_mapping$Node[i] <- paste(found_nodes, collapse = "; ")
    species_node_mapping$Position[i] <- paste(positions, collapse = "; ")
  }
}

# 查看结果
print(species_node_mapping)






library(tidyverse)
library(ggpubr)
meta$group<-factor(meta$group,levels = c('Healthy','AD'))
abundance_philr_used<-abundance_philr[,c(rownames(imp)[1:n],'group')]
pheatmap::pheatmap(abundance_philr_used[,1:ncol(abundance_philr_used)-1]%>%t(),
                   show_rownames = T,
                   show_colnames = F,
                   cluster_rows = T,
                   cluster_cols = F,
                   annotation_col = meta)

genus_used<-unlist(res)%>%unique()
temp<-abundance_used
temp[temp!=0]<-1

pheatmap::pheatmap(temp[genus_used,],
                   show_rownames = T,
                   show_colnames = F,
                   cluster_rows = T,
                   cluster_cols = F,
                   annotation_col = meta)


df_long <- abundance_philr_used %>%
  # 确保 group 是 factor
  mutate(group = factor(group)) %>%
  # 将所有特征列（除了 group）聚集到两列：'Feature' 和 'Value'
  pivot_longer(
    cols = -group, # 不包含 group 列
    names_to = "Feature",
    values_to = "Value"
  )
COLOR_PALETTE <- c('#E7B800', '#00AFBB')
df_long$Feature<-factor(df_long$Feature,levels = colnames(abundance_philr_used))
theme_custom <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      # 移除主要的网格线，保留分面的背景
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # 设置轴线颜色
      axis.line = element_line(color = "grey50"),
      # 调整标题和轴标签的字体
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      # 调整图例位置和外观
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      # 设置分面标签 (strip) 的外观
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text = element_text(face = "bold",size = 14)
    )
}
df_long$Feature<-paste0('Node: ',df_long$Feature)
p <- ggplot(df_long, aes(x = group, y = Value, fill = group)) +
  
  # 1. 绘制小提琴图
  geom_violin(
    trim = TRUE,
    scale = "width",
    alpha = 0.8
  ) +
  
  # 2. 添加箱线图 (用于展示中位数和四分位数)
  geom_boxplot(
    width = 0.15, 
    color = "black", 
    alpha = 0.7, 
    outlier.shape = NA
  ) +
  
  # 3. 添加显著性检验标记 (核心修改)
  # method = "wilcox.test": 适用于两个独立样本的非参数检验
  # label = "p.signif": 用星号 (e.g., ****, ***) 来标记显著性水平
  # hide.ns = TRUE: 隐藏 p > 0.05 的不显著标记 (即 'ns')
  # tip.length: 调整连接显著性标记的线的长度
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif", # 显示星号 (*, **, ***, ****)
    hide.ns = TRUE, 
    label.y.npc = 0.9, # 将标记放在分面顶部 (90% Y轴高度)
    tip.length = 0.01,
    label.x = 1.35
  ) +
  
  # 4. 按 'Feature' 分面
  facet_wrap(~ Feature, scales = "free_y", ncol = 3) +
  
  # 5. 应用新的清新颜色和主题
  scale_fill_manual(values = COLOR_PALETTE) +
  theme_custom(base_size = 12) +
  
  # 6. 设置标签
  labs(
    title = "",
    x = "",
    y = "Value",
    fill = "Group"
  )
p
pdf('D:/microbiome/AD/violin.pdf',width = 10,height = 20)
p
dev.off()

p_tree<-list()

for(i in 1:length(res)){
  target_node_id <- as.numeric(names(res)[i])
  subtree <- extract.clade(genus_tree_pruned, node = target_node_id)
  
  clade1_species <- res[[i]]$clade1
  clade2_species <- res[[i]]$clade2
  
  # 使用 groupOTU 分组
  subtree_grouped <- groupOTU(subtree, 
                              list(Clade1 = clade1_species, 
                                   Clade2 = clade2_species))
  
  # --- 核心逻辑：根据 formula 分配颜色 ---
  # 目标：分母(Down) = 红色 (#E41A1C), 分子(Up) = 蓝色 (#377EB8)
  
  if (res[[i]]$formula == "log(clade1/clade2)") {
    # 此时 clade1 是分子(蓝)，clade2 是分母(红)
    color_map <- c("0" = "grey", "Clade1" = "#377EB8", "Clade2" = "#E41A1C")
  } else {
    # 此时 clade2 是分子(蓝)，clade1 是分母(红)
    color_map <- c("0" = "grey", "Clade1" = "#E41A1C", "Clade2" = "#377EB8")
  }
  # ---------------------------------------
  
  p_tree[[i]] <- ggtree(subtree_grouped, 
                        aes(color = group), 
                        size = 1.5) +
    geom_tiplab(size = 5, color = "black", offset = 0.01, align = FALSE) +
    # 应用动态颜色映射
    scale_color_manual(values = color_map) +
    theme_tree() +
    ggtitle(paste("Node:", target_node_id)) + # 建议保留标题方便识别
    xlim(0, max(subtree$edge.length) * 2.5) +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5, size = 14))
  print(p_tree[[i]])
}
library(patchwork)
pdf('D:/microbiome/AD/combined_tree_1.pdf',width = 10,height = 20)
wrap_plots(p_tree, ncol = 3)
dev.off()

c2_genus<-read.csv('D:/microbiome/AD/cluster_2.csv')
genus_used
intersect(genus_used,c2_genus$x)




#feature_stats <- df_long %>%
#  group_by(Feature) %>%
#  # 1. 计算每个特征在两组中的中位数
#  summarise(
#    Median_AD = median(Value[group == "AD"], na.rm = TRUE),
#    Median_Healthy = median(Value[group == "Healthy"], na.rm = TRUE),
#    # 2. 对每个特征进行 Wilcoxon 秩和检验
#    # 使用 ggpubr::compare_means 来获取 P 值
#    p_value = compare_means(Value ~ group, data = cur_data(), method = "wilcox.test")$p
#  ) %>%
#  ungroup()
#
#P_THRESHOLD <- 0.05
#feature_grouped <- feature_stats %>%
#  mutate(
#    # 确定特征的方向
#    Direction = case_when(
#      Median_AD > Median_Healthy ~ "up",
#      Median_AD < Median_Healthy ~ "down",
#      TRUE ~ "none"
#    ),
#    # 判断是否显著
#    Significance = if_else(p_value < P_THRESHOLD, "yes", "no")
#  )%>%as.data.frame()
##############
#res_long_df <- map_dfr(res, ~{
#  # 提取 clade1 数据
#  clade1_df <- tibble(
#    Clade = "clade1",
#    Species = as.character(.x$clade1)
#  )
#  
#  # 提取 clade2 数据
#  clade2_df <- tibble(
#    Clade = "clade2",
#    Species = as.character(.x$clade2)
#  )
#  
#  # 合并 clade1 和 clade2 的结果
#  bind_rows(clade1_df, clade2_df)
#  
#}, .id = "Feature") %>%
#  # 确保 Feature ID 是字符类型
#  mutate(Feature = as.character(Feature))
#
#
## -------------------- 步骤 2: 连接统计方向 --------------------
#
## 连接 species_df 和 feature_grouped，获取每个 Feature 的 Direction
#species_with_direction <- res_long_df %>%
#  left_join(
#    feature_grouped %>% select(Feature, Direction),
#    by = "Feature"
#  ) %>%
#  # 移除没有明确方向的特征（如果存在）
#  filter(Direction %in% c("up", "down"))
#
## -------------------- 步骤 3: 分组计数并重塑 --------------------
#
## 1. 按 Species 和 Direction 分组计数
#species_counts_long <- species_with_direction %>%
#  group_by(Species, Direction) %>%
#  summarise(
#    Count = n(), # 计算每个物种在每个方向（up/down）出现的次数
#    .groups = 'drop'
#  )
#
## 2. 将长格式重塑为宽格式 (pivot_wider)
#final_count_table <- species_counts_long %>%
#  pivot_wider(
#    names_from = Direction,  # 将 Direction 的值 ('up'/'down') 作为新的列名
#    values_from = Count,     # 将 Count 作为新的列值
#    values_fill = 0          # 如果某个物种在某个方向没有出现，填充为 0
#  ) %>%
#  # 对结果进行排序 (可选)
#  arrange(desc(up + down), Species)
#
#
#
#
#
########根据balance计算频繁项，得到物种共现规则#########
#library(arules)
#transactions_down <- species_with_direction %>%
#  filter(Direction == "down") %>%
#  # 将每个 Feature ID 视为一个事务，Species 为事务中的项
#  select(Feature, Species) %>%
#  # 将数据框转换为 arules 所需的 transactions 对象
#  # 每个 Feature ID 必须是唯一的 ID
#  split(.$Feature) %>% 
#  map(~.$Species) %>%
#  as(., "transactions")
#transactions_up <- species_with_direction %>%
#  filter(Direction == "up") %>%
#  select(Feature, Species) %>%
#  split(.$Feature) %>% 
#  map(~.$Species) %>%
#  as(., "transactions")
######down#####
#rules_down_full <- apriori(transactions_down, 
#                           parameter = list(support = 0.05, 
#                                            confidence = 0.5, 
#                                            maxlen = 3, 
#                                            target = "rules"))
#rules_with_pvalue_down <- interestMeasure(rules_down_full, 
#                                     measure = c("lift", "chi", "fishersExactTest", "phi"), 
#                                     transactions = transactions_down)
#
#rules_base_df_down <- as(rules_down_full, "data.frame") 
#
## 2. 从统计结果中选择需要的列并重命名
#stats_df_down <- rules_with_pvalue_down %>%
#  as.data.frame() %>%
#  # 仅选择我们需要的列，并立即重命名 P 值
#  select(
#    lift_stat = lift,             # 重命名 lift，避免与 rules_base_df 中的 lift 冲突
#    chiSquared, 
#    p_value = fishersExactTest    # 直接重命名 P 值列
#  )
#
## 3. 将基础规则数据和统计数据进行列绑定
#rules_final_df_down <- bind_cols(rules_base_df_down, stats_df_down) %>%
#  
#  # 筛选出显著的规则 (例如，P值 < 0.05)
#  filter(p_value < 0.05) %>%
#  
#  # 排序：使用基础数据中的 lift 进行排序 (也可以使用 lift_stat)
#  # 我们使用 rules_base_df 中原有的 lift
#  arrange(desc(lift), p_value)
#
#
######up#####
#rules_up_full <- apriori(transactions_up, 
#                           parameter = list(support = 0.05, 
#                                            confidence = 0.5, 
#                                            maxlen = 3, 
#                                            target = "rules"))
#rules_with_pvalue_up <- interestMeasure(rules_up_full, 
#                                          measure = c("lift", "chi", "fishersExactTest", "phi"), 
#                                          transactions = transactions_up)
#
#rules_base_df_up <- as(rules_up_full, "data.frame") 
#
## 2. 从统计结果中选择需要的列并重命名
#stats_df_up <- rules_with_pvalue_up %>%
#  as.data.frame() %>%
#  # 仅选择我们需要的列，并立即重命名 P 值
#  select(
#    lift_stat = lift,             # 重命名 lift，避免与 rules_base_df 中的 lift 冲突
#    chiSquared, 
#    p_value = fishersExactTest    # 直接重命名 P 值列
#  )
#
## 3. 将基础规则数据和统计数据进行列绑定
#rules_final_df_up <- bind_cols(rules_base_df_up, stats_df_up) %>%
#  
#  # 筛选出显著的规则 (例如，P值 < 0.05)
#  filter(p_value < 0.05) %>%
#  
#  # 排序：使用基础数据中的 lift 进行排序 (也可以使用 lift_stat)
#  # 我们使用 rules_base_df 中原有的 lift
#  arrange(desc(lift), p_value)

