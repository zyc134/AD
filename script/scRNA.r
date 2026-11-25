rm(list = ls())
library(Seurat)
library(SingleR)
library(celldex)
library(scRNAseq)
library(scater)
library(foreach)
library(stringr)
library(CellChat)
library(harmony)
base_dir<-"D:/microbiome/AD/data/scRNA/GSE157827_RAW"
setwd(base_dir)
sample_names <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
sample_names <- sample_names[sample_names != ""] 

seurat_list <- list()
for (sample in sample_names) {
  message("Processing: ", sample)
  
  # 关键修改：文件路径现在指向子目录下的标准化文件名
  # file.path() 确保路径分隔符在不同操作系统上正确（Windows/Linux）
  
  # 矩阵文件位于：base_dir/sample/matrix.mtx.gz
  mtx_file <- file.path(base_dir, sample, "matrix.mtx.gz")
  
  # 特征文件位于：base_dir/sample/features.tsv.gz
  feature_file <- file.path(base_dir, sample, "features.tsv.gz")
  
  # 条形码文件位于：base_dir/sample/barcodes.tsv.gz
  barcode_file <- file.path(base_dir, sample, "barcodes.tsv.gz")
  
  # 检查文件是否存在，防止因文件缺失导致脚本中断
  if (!all(file.exists(mtx_file, feature_file, barcode_file))) {
    warning(paste("警告: 样本", sample, "的文件缺失，跳过此样本。"))
    next
  }
  
  # 使用 ReadMtx 读取数据
  counts <- ReadMtx(mtx = mtx_file,
                    features = feature_file,
                    cells = barcode_file,
                    # 通常 scRNA-seq 数据中基因名在第二列
                    feature.column = 2) 
  
  # 创建 Seurat 对象
  seurat_obj <- CreateSeuratObject(counts = counts,
                                   project = "GSE157827_AD", # 建议将 project 名称更新为实际的 GSE 编号
                                   min.cells = 3,
                                   min.features = 200)
  
  # 添加样本元数据
  seurat_obj$sample <- sample
  
  # 将对象添加到列表中
  seurat_list[[sample]] <- seurat_obj
  
  message(paste("样本", sample, "创建成功。细胞数:", ncol(seurat_obj)))
}
combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_names)
combined$group <- ifelse(grepl("AD", combined$sample), "AD", "HC")
rm(seurat_list)
rm(seurat_obj)
rm(counts)
gc()
saveRDS(combined,'../combined.rds')
rm(list = ls())
gc()
combined<-readRDS('../combined.rds')
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
#qc_no<-VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 20)
#qc<-VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
gc()
symbol<-rownames(combined)
library(org.Hs.eg.db)
library(clusterProfiler)
res <- bitr(symbol,
            fromType = "SYMBOL",
            toType   = "ENSEMBL",
            OrgDb    = org.Hs.eg.db)
colnames(res)<-c('symbol','ensembl_id')
tmp<-table(res$ensembl_id)%>%as.data.frame()
tmp<-tmp[tmp$Freq==1,]
res<-res[res$ensembl_id%in%tmp$Var1,]
combined<-combined[rownames(combined)%in%res$symbol,]
symbol<-rownames(combined)
res <- bitr(symbol,
            fromType = "SYMBOL",
            toType   = "ENSEMBL",
            OrgDb    = org.Hs.eg.db)
colnames(res)<-c('symbol','ensembl_id')
tmp<-table(res$symbol)%>%as.data.frame()
tmp<-tmp[tmp$Freq==1,]
res<-res[res$symbol%in%tmp$Var1,]



combined<-combined[rownames(combined)%in%res$symbol,]
rownames(combined)<-res$ensembl_id
library(SeuratDisk)
SaveH5Seurat(
  combined, 
  filename = "../integrated_data.h5seurat", 
  overwrite = TRUE
)

# 步骤 B: 转换为 H5AD 格式
Convert(
  "../integrated_data.h5seurat", 
  dest = "h5ad", 
  overwrite = TRUE
)



write.csv(res,'../gene_map_df.csv')









combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(combined))
combined <- RunHarmony(combined, group.by.vars = "sample", reduction.use = "pca")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
DimPlot(combined, reduction = "umap",group.by = "seurat_clusters",label = TRUE)
DimPlot(combined, reduction = "umap",group.by = "group",label = TRUE)
ref_data <- celldex::HumanPrimaryCellAtlasData()
expr_data <- GetAssayData(combined, assay = "RNA", slot = "data")
combined.hesc <- SingleR(test = expr_data, ref = ref_data, labels = ref_data$label.main)
combined@meta.data$labels<-combined.hesc$labels
temp<-table(combined@meta.data$labels)%>%as.data.frame()
match(temp$Var1[temp$Freq<500],combined@meta.data$labels)
combined@meta.data$labels[combined@meta.data$labels%in%temp$Var1[temp$Freq<1000]]<-'others'
DimPlot(combined, reduction = "umap",group.by = "labels",label = TRUE)
saveRDS(combined,'../combined1.rds')
DotPlot(combined, features = 'PPARG', dot.scale = 5) + RotatedAxis()
FeaturePlot(combined, features = 'PPARG', split.by = "group", max.cutoff = 3,reduction = "umap")
FeaturePlot(combined, features = 'PPARG', split.by = "labels", max.cutoff = 3,reduction = "umap")
FeaturePlot(combined, features = 'IL1B', split.by = "cell_type", max.cutoff = 3,reduction = "umap")


combined<-readRDS('../combined1.rds')
dim(combined)
cell_label<-read.csv('D:/microbiome/AD/data/scRNA/output.csv')
unique(cell_label$supercluster_name)
unique(cell_label$cluster_name)
unique(cell_label$subcluster)
match(combined@active.ident%>%names(),cell_label$cell_id)
combined@meta.data[["cell_type"]]<-cell_label$supercluster_name
tmp<-table(combined@meta.data[["cell_type"]])%>%as.data.frame()
tmp$Freq[tmp$Var1=='Miscellaneous'|tmp$Var1=='Splatter']<-1
combined@meta.data[["cell_type"]][combined@meta.data[["cell_type"]]%in%tmp$Var1[tmp$Freq<500]]<-'Others'
gc()
combined<-subset(combined, subset = cell_type != 'Others')
dim(combined)

cellmeta<-combined@meta.data
prop_df <- cellmeta %>%
  group_by(group, cell_type) %>%
  summarise(n = n(), .groups = "drop_last") %>%  # 统计每组内各 cluster 数量
  mutate(proportion = n / sum(n)) %>%             # 计算在该组内的占比
  ungroup()
library(RColorBrewer)
cell_types <- levels(as.factor(combined$cell_type))
num_cell_types <- length(cell_types)
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_cell_types)
names(my_colors) <- cell_types
prop<-ggplot(prop_df, aes(x = group, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  labs(x = "Group", y = "Relative Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
prop

umap_cluster<-DimPlot(combined, reduction = "umap", 
                  group.by = 'seurat_clusters',
                  label = TRUE, 
                  #cols = my_colors,
                  repel = TRUE)

umap_all<-DimPlot(combined, reduction = "umap", 
        group.by = 'cell_type',
        label = F, 
        cols = my_colors,
        repel = TRUE)
pdf('../umap_cluster.pdf',height = 9.5,width = 9.5)
umap_cluster
dev.off()
pdf('../umap_celltype.pdf',height = 9.5,width = 9.5)
umap_all
dev.off()
p <- umap_cluster | umap_all | prop |
  plot_layout(widths = c(4, 4, 1))
p
pdf('../umap_prop.pdf',width = 16,height = 6)
p
dev.off()
Oligodendrocyte<-subset(combined, subset = cell_type == 'Oligodendrocyte')
FeaturePlot(Oligodendrocyte, features = 'PPARG', split.by = "group", max.cutoff = 3,reduction = "umap")



saveRDS(Oligodendrocyte,'../Oligodendrocyte.rds')

DimPlot(combined, reduction = "umap",group.by = "cell_type",label = TRUE)
pdf('../PPARG_group.pdf',width = 14,height = 9.5)
FeaturePlot(combined, features = 'PPARG', split.by = "group", max.cutoff = 3,reduction = "umap")
dev.off()
cell_types<-unique(combined@meta.data[["cell_type"]])
deg_list <- list()
for (celltype in cell_types) {
  print(celltype)
  subset_obj <- subset(combined, subset = cell_type == celltype)
  Idents(subset_obj) <- subset_obj$group
  deg <- FindMarkers(subset_obj, ident.1 = "AD", ident.2 = "HC",
                     min.pct = 0, logfc.threshold = 0)
  deg$cell_type <- celltype
  deg$gene <- rownames(deg)
  deg_list[[celltype]] <- deg
  gc()
}
deg_all <- do.call(rbind, deg_list)
deg_all<-deg_all[deg_all$p_val_adj<0.05&abs(deg_all$avg_log2FC)>1,]
write.csv(deg_all,'../deg_allcells.csv')

Idents(combined) <- "cell_type"
all_cell_markers <- FindAllMarkers(
  object = combined, 
  only.pos = TRUE,               # 仅报告上调的基因（这是我们通常想要的 Marker）
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  assay = "RNA"               # 推荐使用 'RNA' Assay 的数据进行差异分析
  # future::plan("multisession", workers = 4) # 可选：用于加速计算（如果数据量大）
)

markers_top3<-all_cell_markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

marker_all<-DotPlot(combined, features = unique(markers_top3$gene), dot.scale = 5) + RotatedAxis()

pdf('../markers_all.pdf',height = 6,width = 16)
marker_all
dev.off()







PPARG<-deg_all[grep('PPARG',deg_all$gene),]

Microglia<-subset(combined, subset = cell_type == 'Microglia')
saveRDS(Microglia,'../Microglia.rds')
Microglia <- FindClusters(Microglia, resolution = 0.3)
table(Microglia@meta.data[["RNA_snn_res.0.3"]])
#Microglia <- subset(Microglia, subset = RNA_snn_res.0.5 != '8')#cluster13细胞数量少于10
DimPlot(Microglia, reduction = "umap",group.by = "RNA_snn_res.0.3",label = TRUE)
cellmeta<-Microglia@meta.data
prop_df <- cellmeta %>%
  group_by(group, RNA_snn_res.0.3) %>%
  summarise(n = n(), .groups = "drop_last") %>%  # 统计每组内各 cluster 数量
  mutate(proportion = n / sum(n)) %>%             # 计算在该组内的占比
  ungroup()



Microglia_clustrs<-unique(Microglia@meta.data[["RNA_snn_res.0.5"]])
deg_Microglia <- list()
for (clustr in Microglia_clustrs) {
  print(clustr)
  subset_obj <- subset(Microglia, subset = RNA_snn_res.1 == clustr)
  Idents(subset_obj) <- subset_obj$group
  deg <- FindMarkers(subset_obj, ident.1 = "AD", ident.2 = "HC",
                     min.pct = 0, logfc.threshold = 0)
  deg$Microglia_clustrs <- clustr
  deg$gene <- rownames(deg)
  deg_Microglia[[clustr]] <- deg
}
deg_all_Microglia <- do.call(rbind, deg_Microglia)
deg_all_Microglia<-deg_all_Microglia[deg_all_Microglia$p_val_adj<0.05&abs(deg_all_Microglia$avg_log2FC)>1,]

#####cellchat#####

seurat.list <- SplitObject(combined, split.by = "group")
rm(combined)
rm(cell_label)
gc()
CellChatDB <- CellChatDB.human

data.input1 <- GetAssayData(seurat.list[["AD"]], slot = "data", assay = "RNA")
meta1 <- seurat.list[["AD"]]@meta.data
cellchat1 <- createCellChat(object = data.input1, meta = meta1, group.by = "cell_type")
cellchat1@DB <- CellChatDB


data.input2 <- GetAssayData(seurat.list[["HC"]], slot = "data", assay = "RNA")
meta2 <- seurat.list[["HC"]]@meta.data
cellchat2 <- createCellChat(object = data.input2, meta = meta2, group.by = "cell_type")
cellchat2@DB <- CellChatDB

cellchat1 <- subsetData(cellchat1)
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
cellchat1 <- computeCommunProb(cellchat1)
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
saveRDS(cellchat1,'../cellchat_AD.rds')

cellchat2 <- subsetData(cellchat2)
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- computeCommunProb(cellchat2)
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
saveRDS(cellchat2,'../cellchat_HC.rds')

cellchat1<-readRDS('../cellchat_AD.rds')
cellchat2<-readRDS('../cellchat_HC.rds')
cellchat1@options$ident <- "AD"
cellchat2@options$ident <- "HC"

object.list <- list(AD = cellchat1, HC = cellchat2)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
pdf('../cellchat_bar.pdf')
compareInteractions(cellchat_merged, show.legend = F, group = c(1, 2))
dev.off()
pdf('../cellchat_net_GSE269906.pdf',width = 11.5,height = 14)
netVisual_diffInteraction(cellchat, 
                          weight.scale = TRUE,
                          measure = 'count',
                          comparison = c(1, 2),
                          label.edge = F)
dev.off()

pdf('../cellchat_heatmap.pdf',width = 10,height = 6)
netVisual_heatmap(cellchat1)+netVisual_heatmap(cellchat2)
dev.off()

cellchat_AD<-netVisual_circle(cellchat1@net$count, weight.scale = T, label.edge = F, title.name = "AD Total Interactions")
cellchat_HC<-netVisual_circle(cellchat2@net$count, weight.scale = T, label.edge = F, title.name = "HC Total Interactions")
pdf('../cellchat_AD.pdf',width = 9,height = 9)
cellchat_AD
dev.off()
pdf('../cellchat_HC.pdf',width = 9,height = 9)
cellchat_HC
dev.off()



###############
rm(list = ls())
gc()
Microglia<-readRDS('../Microglia.rds')
library(SCENIC)
library(SCopeLoomR)

##normalized_data_matrix <- GetAssayData(Microglia, assay = "RNA", slot = "data")
#raw_count_matrix <- GetAssayData(Microglia, assay = "RNA", slot = "counts")
##gene_expression_sums <- colSums(raw_count_matrix%>%as.matrix())
##threshold_sum <- ncol(raw_count_matrix) * 0.03
##genes_over_threshold <- gene_expression_sums > threshold_sum
##genes_expressed_in_cells <- Matrix::rowSums(raw_count_matrix > 0) > (ncol(raw_count_matrix) * 0.01)
##filtered_genes <- genes_over_threshold & genes_expressed_in_cells
#exprMat_filtered <- raw_count_matrix
#loom <- SCopeLoomR::build_loom(
#  file.name = "./Mcell.loom",
#  dgem = exprMat_filtered,
#  default.embedding = NULL
#)
#loom$close()
#data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
#motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
#scenicOptions <- initializeScenic(org="hgnc", dbDir="../scenic_db", nCores=1)
#
#genesKept <- geneFiltering(exprMat_filtered%>%as.matrix(), scenicOptions)
#exprMat_filtered <- exprMat_filtered[genesKept, ]
#runCorrelation(exprMat_filtered%>%as.matrix(), scenicOptions)
#exprMat_filtered_log <- log2(exprMat_filtered+1) 
#runGenie3(exprMat_filtered_log%>%as.matrix(), scenicOptions)
#loom <- open_loom('../scenic_out/out_SCENIC1.loom') 
#regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
#regulons <- regulonsToGeneLists(regulons_incidMat)
#regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
#regulonAucThresholds <- get_regulon_thresholds(loom)
#embeddings <- get_embeddings(loom)
#close_loom(loom)
#identical(colnames(regulonAUC), colnames(Microglia))
#cellClusters <- data.frame(row.names = colnames(Microglia), 
#                           seurat_clusters = as.character(Microglia@meta.data[["cell_type"]]))
#cellTypes <- data.frame(row.names = colnames(Microglia), 
#                        celltype = Microglia@meta.data[["cell_type"]])

##########

##########
#loom <- open_loom("D:/microbiome/AD/data/scRNA/scenic_out/out_SCENIC1.loom")
#regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
#regulons <- regulonsToGeneLists(regulons_incidMat)
#regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
#regulonAucThresholds <- get_regulon_thresholds(loom)
#tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
#
#embeddings <- get_embeddings(loom)  
#close_loom(loom)
#rownames(regulonAUC)
#names(regulons)
#if ("PPARG" %in% rownames(Microglia)) {
#  message("正在可视化 PPARG 基因表达水平...")
#  
#  # UMAP 降维图上显示 PPARG 基因的表达
#  p1_gene <- FeaturePlot(Microglia, features = "PPARG", reduction = "umap") +
#    ggtitle("PPARG Gene Expression (UMAP)")
#  
#  # 小提琴图展示 PPARG 基因在不同细胞亚群中的表达
#  p2_gene <- VlnPlot(Microglia, features = "PPARG", group.by = "group", pt.size = 0) +
#    NoLegend() +
#    ggtitle("PPARG Gene Expression (VlnPlot)")
#  
#  # 组合展示
#  print(p1_gene | p2_gene)
#  
#} else {
#  warning("PPARG 基因不在 Seurat 对象的表达矩阵中。请检查基因名是否正确或数据准备过程。")
#}
#pparg_gene_index <- which(colnames(regulons_incidMat) == "PPARG")
#if (length(pparg_gene_index) > 0) {
#  message("PPARG 基因在 Regulon 矩阵中找到，正在查找调控它的转录因子...")
#  
#  # 提取 PPARG 基因在所有 Regulons 中的状态 (TRUE/FALSE)
#  pparg_targets_status <- regulons_incidMat[, pparg_gene_index, drop = FALSE]
#  
#  # 找出包含 PPARG 作为靶基因的 Regulons (即值为 TRUE 的行)
#  # 结果的行名就是调控 PPARG 的转录因子 (Regulon Leader)
#  tfs_regulating_pparg_incid <- rownames(pparg_targets_status)[as.vector(pparg_targets_status)]
#  
#  if (length(tfs_regulating_pparg_incid) > 0) {
#    # 结果可能包含扩展的 Regulon (例如 TF(+) 和 TF(E))，我们只关注 TF 名称
#    tfs_regulating_pparg <- unique(gsub("\\(\\+.*\\)", "", tfs_regulating_pparg_incid))
#    
#    cat("\n==================================================\n")
#    cat(" 调控 PPARG 的上游转录因子（Regulon Leaders）：\n")
#    cat("==================================================\n")
#    print(tfs_regulating_pparg)
#    
#    cat("\n完整包含 PPARG 的 Regulon 名称（包括特异性后缀）：\n")
#    print(tfs_regulating_pparg_incid)
#    cat("\n==================================================\n")
#    
#  } else {
#    message("在 pySCENIC 的结果中，没有发现任何 Regulon 将 PPARG 列为靶基因。")
#    message("这意味着 SCENIC 分析未推断出直接调控 PPARG 的上游转录因子。")
#  }
#} else {
#  warning("PPARG 基因不在 Regulon 靶基因列表中。请检查基因名是否正确，或 pySCENIC 的输入数据。")
#}
#regulon_name <- "AHRR(+)"
#if (regulon_name %in% rownames(regulonAUC)) {
#  message(paste0("正在可视化 ", regulon_name, " 调节因子的活性..."))
#  
#  # 1. 提取 AHRR(+) Regulon 的活性分数
#  ahrr_auc <- regulonAUC[regulon_name, ]
#  
#  # 2. 将 AUC 分数添加到 Seurat 对象的 metadata 中
#  # 确保细胞顺序与 Seurat 对象的列名一致
#  # 假设您的细胞分组信息存储在 seurat.data$seurat_annotations
#  #seurat.data@meta.data[[regulon_name]] <- ahrr_auc[colnames(seurat.data)]
#  
#  # ================================================================
#  # 3. 可视化
#  # ================================================================
#  
#  # A. UMAP 降维图上显示 AHRR(+) Regulon 的活性
#  p1_ahrr <- FeaturePlot(Microglia, features = gsub('\\(\\+\\)','',regulon_name), reduction = "umap") +
#    ggtitle(paste0(regulon_name, " Regulon Activity (UMAP)"))
#  
#  # B. 小提琴图展示 AHRR(+) Regulon 在不同细胞亚群中的活性分布
#  # 'group.by' 参数需要指定您的细胞分组列名，这里假设是 "seurat_annotations"
#  p2_ahrr <- VlnPlot(Microglia, features = gsub('\\(\\+\\)','',regulon_name), group.by = "group", pt.size = 0, assay = "RNA") +
#    NoLegend() +
#    ggtitle(paste0(regulon_name, " Regulon Activity (VlnPlot)"))
#  
#  # C. 点图展示不同亚群的平均活性和表达细胞百分比
#  p3_ahrr <- DotPlot(seurat.data, features = regulon_name, group.by = "seurat_annotations") +
#    RotatedAxis() +
#    ggtitle(paste0(regulon_name, " Activity (DotPlot)"))
#  
#  # 组合展示所有图
#  print(p1_ahrr / (p2_ahrr | p3_ahrr))
#  
#} else {
#  warning(paste0(regulon_name, " 调节因子不在 SCENIC 结果的 AUC 矩阵中。请检查 pySCENIC 运行结果。"))
#}
######Microglia通路#####
#library(WGCNA)
#datExpr <- as.data.frame(t(GetAssayData(Microglia, assay = "RNA", slot = "data")))
#gc()
#Microglia<-FindVariableFeatures(Microglia, 
#                                selection.method = "vst", 
#                                nfeatures = 2000)
#top10 <- head(VariableFeatures(Microglia), 10)
#
#plot1 <- VariableFeaturePlot(Microglia)
#plot2 <- LabelPoints(plot = plot1, 
#                     points = top10, 
#                     repel = TRUE, 
#                     xnudge = 0, 
#                     ynudge = 0)
#
## 组合并显示图
#plot1 + plot2
#hvg.genes <- VariableFeatures(Microglia)
#datExpr<-datExpr[,hvg.genes]
#PPARG_expression <- datExpr$PPARG
#AD_status <- as.numeric(Microglia@meta.data$group == "AD")
#datTraits <- data.frame(
#  PPARG_Expr = PPARG_expression,
#  AD_Status = AD_status,
#  row.names = rownames(datExpr) # 确保行名是细胞ID
#)
#powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#sft
#softPower = sft$powerEstimate
#softPower
#net = blockwiseModules(datExpr, 
#                       power = softPower,
#                       TOMType = "unsigned", # 无符号网络 (正负相关都视为连接)
#                       minModuleSize = 30,  # 模块最小基因数
#                       reassignThreshold = 0, 
#                       mergeCutHeight = 0.25, # 模块合并的阈值
#                       numericLabels = TRUE,
#                       pamStage = FALSE, 
#                       verbose = 3)
#moduleLabels = net$colors
#moduleColors = labels2colors(net$colors)
#MEs = net$MEs # Module Eigengenes (模块特征基因)
#geneTree = net$dendro[[1]]
#
#moduleTraitCor = cor(MEs, datTraits, use = "p")
#moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr))
#sizeGrWindow(10, 6)
#textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
#                   signif(moduleTraitPvalue, 1), ")", sep = "")
#dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3));
#labeledHeatmap(Matrix = moduleTraitCor,
#               xLabels = names(datTraits),
#               yLabels = names(MEs),
#               ySymbols = names(MEs),
#               colorLabels = FALSE,
#               colors = blueWhiteRed(50),
#               textMatrix = textMatrix,
#               setStdMargins = FALSE,
#               cex.text = 0.5,
#               zlim = c(-1,1),
#               main = paste("Module-trait relationships"))
##########
#deg_all<-read.csv('../deg_allcells.csv',row.names = 1)
#deg_Microglia<-deg_all[deg_all$cell_type=='Microglia',]
#pvalueCutoff <- 0.05
#library(clusterProfiler)
#library(org.Hs.eg.db)
#go_result <- enrichGO(gene = deg_Microglia$gene, 
#                      OrgDb = org.Hs.eg.db,
#                      keyType = "SYMBOL",
#                      ont = "all", 
#                      pvalueCutoff = pvalueCutoff)
#dotplot(go_result)
#gene_id<-bitr(deg_Microglia$gene, 
#              fromType = "SYMBOL", 
#              toType = "ENTREZID", 
#              OrgDb = org.Hs.eg.db)
#
#kegg_res <- enrichKEGG(gene = gene_id$ENTREZID, 
#                       keyType = 'kegg',
#                       organism = "hsa", 
#                       pvalueCutoff = pvalueCutoff)
#dotplot(kegg_res)
######
#diff_geneset<-deg_Microglia[order(deg_Microglia$avg_log2FC,decreasing = T),]
#ranks <- diff_geneset$avg_log2FC
#names(ranks)<-diff_geneset$gene
#Go_gseresult <- gseGO(ranks, 
#                      'org.Hs.eg.db', 
#                      keyType = "SYMBOL", 
#                      ont="all", 
#                      nPerm = 1000, 
#                      minGSSize = 10, 
#                      maxGSSize = 1000, 
#                      pAdjustMethod = "none",
#                      pvalueCutoff=0.05,
#)
#library(enrichplot)
#gseaplot2(Go_gseresult, 2, pvalue_table = T)
######hdWGCNA######
#rm(list = ls())
gc()
Microglia<-readRDS('../Microglia.rds')
Microglia <- FindVariableFeatures(Microglia, 
                                        selection.method = "vst", 
                                        nfeatures = 2000)
library(hdWGCNA)
obj.sub<-Microglia
obj.sub <- SetupForWGCNA(
  obj.sub,
  gene_select = "variable", 
  wgcna_name = "Microglia_Sub"
)
obj.sub <- MetacellsByGroups(
  obj.sub,
  group.by    = c('group',"cell_type"),
  ident.group = "group",
  reduction = 'pca',
  k = 25,
  max_shared = 10
)
obj.sub <- NormalizeMetacells(obj.sub)
obj.sub <- SetDatExpr(
  obj.sub,
  group_name = "Microglia", # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
obj.sub <- TestSoftPowers(
  obj.sub,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(obj.sub)
library(patchwork)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)
obj.sub <- ConstructNetwork(
  obj.sub,
  soft_power = 7,
  overwrite_tom = TRUE,
  tom_name = 'Microglia' # name of the topoligical overlap matrix written to disk
)
pdf('../Dendrogram_7.pdf')
PlotDendrogram(obj.sub, main='Microglia hdWGCNA Dendrogram')
dev.off()
obj.sub <- ModuleEigengenes(
  obj.sub,
  group.by.vars = "group" # 这里填入你 meta.data 中用于批次校正或区分的变量，如 group
)
obj.sub <- ModuleConnectivity(
  obj.sub,
  group.by = 'cell_type', 
  group_name = 'Microglia'
)
# rename the modules
obj.sub <- ResetModuleNames(
  obj.sub,
  new_name = "Microglia-M"
)

# plot genes ranked by kME for each module
PlotKMEs(obj.sub, ncol=5)



plot_list <- ModuleFeaturePlot(
  obj.sub,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
pdf('../umap_M2.pdf')
wrap_plots(plot_list[[2]], ncol=1)
dev.off()



hub_df <- GetHubGenes(obj.sub, n_hubs = 999999)
modules <- GetModules(obj.sub)
hub_df$color<-modules$color[match(hub_df$module,modules$module)]
###########
# ---------------------------------------------------------
# 1. 批量计算所有模块的活性差异 (Wilcoxon Test)
# ---------------------------------------------------------
library(ggpubr)
library(clusterProfiler)

# 提取所有模块的 Eigengenes (MEs) 矩阵
ME_df <- obj.sub@misc$Microglia_Sub$wgcna_net$multiMEs$Microglia_Sub$data
meta_mc <- obj.sub@misc$Microglia_Sub$wgcna_metacell_obj@meta.data

# 移除 grey 模块（未分配到任何模块的基因），通常名为 MEgrey 或以 grey 结尾
ME_df <- ME_df[, !grepl("grey", colnames(ME_df), ignore.case = TRUE)]

# 创建一个数据框来存储所有模块的P值
diff_results <- data.frame(Module_ME = character(), 
                           Module_Clean = character(),
                           Pval = numeric(), 
                           stringsAsFactors = FALSE)

# 循环遍历每个模块进行 Wilcoxon 检验
for (mod in colnames(ME_df)) {
  plot_df <- data.frame(
    ME_Value = ME_df[[mod]],
    group = meta_mc$group
  )
  
  # 进行两组间的非参数检验 (假设你的 group 只有两组)
  test_res <- wilcox.test(ME_Value ~ group, data = plot_df)
  
  # 清理模块名称 (例如把 "MEMicroglia-M1" 变成 "Microglia-M1"，以便与 modules 数据框匹配)
  clean_name <- gsub("^ME", "", mod) 
  
  diff_results <- rbind(diff_results, data.frame(
    Module_ME = mod, 
    Module_Clean = clean_name,
    Pval = test_res$p.value
  ))
}

# ---------------------------------------------------------
# 2. 筛选具有显著活性差异的模块并提取对应基因
# ---------------------------------------------------------
# 设定显著性阈值，例如 p < 0.05
pvalueCutoff <- 0.05
sig_modules_info <- diff_results[diff_results$Pval < pvalueCutoff, ]

print("具有显著活性差异的模块有：")
print(sig_modules_info$Module_Clean)

# 从总的 modules 数据框中提取这些显著模块的所有基因
modules <- GetModules(obj.sub)
# 过滤：只保留属于差异显著模块的基因
sig_modules_genes <- hub_df[hub_df$color%in%sig_modules_info$Module_Clean,c(1,2),drop=F]
sig_modules_genes
# ---------------------------------------------------------
# 3. 使用 gene_used 基因集对显著模块进行富集分析
# ---------------------------------------------------------
# 【注意】请确保 gene_used 是一个包含基因名的字符向量，例如：
# gene_used <- c("GeneA", "GeneB", "GeneC", ...)

# 构建 TERM2GENE 格式的数据框，clusterProfiler::enricher 需要这种两列格式：
# 第一列是Term(分类名/模块名)，第二列是Gene(基因名)
term2gene <- sig_modules_genes[, c("module", "gene_name")]

# 使用 enricher 进行超几何分布检验 (ORA)
head(gene_used)
head(term2gene)
unique(term2gene$module)
module_enrichment <- enricher(
  gene = gene_used,
  TERM2GENE = term2gene,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
module_enrichment@result$ID
# ---------------------------------------------------------
# 4. 可视化富集结果
# ---------------------------------------------------------
if (!is.null(module_enrichment) && nrow(module_enrichment@result[module_enrichment@result$p.adjust < 0.05, ]) > 0) {
  print("富集分析成功，正在生成结果图...")
  
  # 打印具体的富集表格查看
  print(head(module_enrichment@result))
  
  # 绘制气泡图
  pdf('../gene_used_module_enrichment.pdf', width = 7, height = 6)
  p <- dotplot(module_enrichment, showCategory = 15, title = "Enrichment of 'gene_used' in Active Modules")
  print(p)
  dev.off()
  
} else {
  print("基因集在差异模块中没有发现具有统计学显著性的富集结果。")
}


###########
modules<-modules[order(modules$`kME_Microglia-M2`,decreasing = T),]
modules$gene_name[modules$color=='black']


target_mod <- modules[modules$gene_name == "PPARG", ]
hub_genes <- GetHubGenes(obj.sub, n_hubs = 30)
M1_gene<-modules$gene_name[modules$module=='Microglia-M2']
pvalueCutoff <- 0.05
library(clusterProfiler)
library(org.Hs.eg.db)
go_result <- enrichGO(gene = M1_gene, 
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "all", 
                      pvalueCutoff = pvalueCutoff)
pdf('../go_Microglia.pdf')
dotplot(go_result,showCategory = 11)
dev.off()
go_Microglia<-dotplot(go_result,showCategory = 11)
ME_df <- obj.sub@misc$Microglia_Sub$wgcna_net$multiMEs$Microglia_Sub$data
meta_mc <- obj.sub@misc$Microglia_Sub$wgcna_metacell_obj@meta.data
plot_df <- data.frame(
  ME_M1 = ME_df$MEblack,
  group = meta_mc$group
)
library(ggpubr)
M2_act<-ggplot(plot_df, aes(x = group, y = ME_M1, fill = group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",   # 非参数检验，单细胞/宏基因组都推荐
    label = "p.signif",        # 显示 *, **, ***
    label.x = 1.5
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "M2 module activity",
    y = "M2 module eigengene",
    x = NULL
  )
pdf('../M2_activity.pdf')
ggplot(plot_df, aes(x = group, y = ME_M1, fill = group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",   # 非参数检验，单细胞/宏基因组都推荐
    label = "p.signif",        # 显示 *, **, ***
    label.x = 1.5
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "M2 module activity across groups",
    y = "M2 module eigengene",
    x = NULL
  )
dev.off()
cor(obj.sub@misc[["Microglia_Sub"]][["hMEs"]][["Microglia-M2"]], GetAssayData(obj.sub)["PPARG",])
#cor(obj.sub@misc[["Microglia_Sub"]][["hMEs"]][["Microglia-M1"]], GetAssayData(obj.sub)["MITF",])
#cor(obj.sub@misc[["Microglia_Sub"]][["hMEs"]][["Microglia-M1"]], GetAssayData(obj.sub)["MSR1",])
#cor(obj.sub@misc[["Microglia_Sub"]][["hMEs"]][["Microglia-M1"]], GetAssayData(obj.sub)["GPNMB",])
#cor(obj.sub@misc[["Microglia_Sub"]][["hMEs"]][["Microglia-M1"]], GetAssayData(obj.sub)["LPL",])

library(AUCell)
expr_matrix <- GetAssayData(Microglia, slot = "data")
cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 4, plotStats = TRUE)

library(stringr)
library(ggrastr)
id<-go_result@result[["ID"]]
p_aucell_go_Microglia<-list()
for(i in 1:length(id)){
  geneSets <- list(hdWGCNA_M1 = go_result@geneSets[[id[i]]])
  # 4. 计算 AUC 得分
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  # 提取 M1 的 AUC 分数并整合到 Seurat meta.data
  M1_auc_scores <- as.numeric(getAUC(cells_AUC)["hdWGCNA_M1", ])
  Microglia$M1_AUC <- M1_auc_scores
  hME_M1 <- obj.sub@misc[["Microglia_Sub"]][["hMEs"]][["Microglia-M2"]]
  #common_cells <- intersect(colnames(Oligodendrocyte), rownames(obj.sub@misc$Oligodendrocyte_Sub$hMEs))
  #hME_val <- hME_M1[common_cells]
  AUC_val <- Microglia@meta.data[, "M1_AUC"]
  cor_res <- cor.test(hME_M1, AUC_val, method = "spearman")
  library(ggplot2)
  df_cor <- data.frame(hME = hME_M1, AUC = AUC_val)
  p_aucell_go_Microglia[[i]]<-ggplot(df_cor, aes(x = hME, y = AUC)) +
    geom_point_rast(alpha = 0.3, size = 0.5, color = "steelblue") +
    geom_smooth(method = "lm", color = "red") +
    annotate("text", x = min(hME_M1), y = max(AUC_val), 
             label = paste0("rho = ", round(cor_res$estimate, 3), 
                            "\nP < 2.2e-16"), hjust = 0) +
    theme_classic() +
    ggtitle(id[i])
  #print(p_aucell_go_Oligodendrocyte[[i]])
  print(cor_res)
}
pdf('../aucell_Microglia.pdf',width = 10,height = 5)
wrap_plots(p_aucell_go_Microglia[[1]],
           p_aucell_go_Microglia[[2]],
           p_aucell_go_Microglia[[8]],ncol = 3)
dev.off()

pdf('../act_go_Microglia.pdf',width = 10,height = 7)
M2_act+go_Microglia
dev.off()
################
#gene_id<-bitr(M1_gene, 
#              fromType = "SYMBOL", 
#              toType = "ENTREZID", 
#              OrgDb = org.Hs.eg.db)
#library(ReactomePA)
#res_reactome <- enrichPathway(gene = gene_id$ENTREZID,
#                              organism = "human")
#pdf('../reactome_Microglia_M7.pdf')
#dotplot(res_reactome)
#dev.off()
ModuleNetworkPlot(
  obj.sub,
  mods = "Microglia-M2",
  n_inner = 30, # 标注前10个枢纽基因
  n_outer = 15,      # 展示前20个基因的连接
  plot_size = c(15,15)
)
##
#modules_M1<-modules[modules$module=='Microglia-M1',,drop=F]
#modules_M1<-modules_M1[order(modules_M1$`kME_Microglia-M1`,decreasing = T),,drop=F]
#ranks <- modules_M1$`kME_Microglia-M1`
#names(ranks)<-modules_M1$gene_name
#Go_gseresult <- gseGO(ranks, 
#                      'org.Hs.eg.db', 
#                      keyType = "SYMBOL", 
#                      ont="all", 
#                      nPerm = 1000, 
#                      minGSSize = 10, 
#                      maxGSSize = 1000, 
#                      pAdjustMethod = "none",
#                      pvalueCutoff=0.05)
#gseaplot2(Go_gseresult, 1, pvalue_table = T)

#####
#diff_geneset<-deg_Oligodendrocyte[order(deg_Oligodendrocyte$p_val_adj,decreasing = T),]
#ranks <- diff_geneset$p_val_adj
#names(ranks)<-diff_geneset$gene
#Go_gseresult <- gseGO(ranks, 
#                      'org.Hs.eg.db', 
#                      keyType = "SYMBOL", 
#                      ont="all", 
#                      nPerm = 1000, 
#                      minGSSize = 10, 
#                      maxGSSize = 1000, 
#                      pAdjustMethod = "none",
#                      pvalueCutoff=0.05,
#)
#library(enrichplot)
#gseaplot2(Go_gseresult, 1:12, pvalue_table = T)
#########
Microglia<-readRDS('../Microglia.rds')
Microglia <- FindVariableFeatures(Microglia, 
                                  selection.method = "vst", 
                                  nfeatures = 3000)
DefaultAssay(Microglia) <- "RNA"

obj.sub <- Microglia
obj.sub<-subset(obj.sub, subset = seurat_clusters == '8')
obj.sub@meta.data%>%colnames()
obj.sub <- FindVariableFeatures(obj.sub, nfeatures = 1000)
obj.sub <- ScaleData(obj.sub) 
obj.sub <- RunPCA(obj.sub, features = VariableFeatures(obj.sub))
Idents(obj.sub) <- "cell_type"
obj.sub <- FindClusters(obj.sub, resolution = 0.3)
all_markers <- FindAllMarkers(obj.sub,min.pct = 0.25, only.pos = TRUE,logfc.threshold = 0.25)
all_markers <- all_markers[-grep('^MT-',all_markers$gene),]
markers_top5<-all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

DoHeatmap(obj.sub, features = markers_top5$gene) + NoLegend()
DotPlot(obj.sub, features = unique(markers_top5$gene), dot.scale = 5) + RotatedAxis()


Microglia<-subset(obj.sub, subset = RNA_snn_res.0.3 != 4)
new_labels <- c(
  "0" = "Mature / State-specific Microglia",
  "1" = "Homeostatic Microglia (hMG)",
  "2" = "Surveilling / Complement-rich Microglia",
  "3" = "Ribosomal-high (Ribo-high) Microglia"
)
Microglia <- RenameIdents(Microglia, new_labels)
DimPlot(Microglia, reduction = "umap", label = TRUE, repel = TRUE)+
  ylim(c(3,8))+
  xlim(c(5,10))


FeaturePlot(Microglia, features = markers_top5$gene, 
            #split.by = "group", 
            max.cutoff = 3,
            reduction = "umap")






library(monocle)
library(ggplot2)
library(dplyr)
expression_matrix <- GetAssayData(object = Microglia, assay = "RNA", slot = "counts")
cell_metadata <- Microglia@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(Microglia[["RNA"]]), 
                              row.names = rownames(Microglia[["RNA"]]))
pd <- new('AnnotatedDataFrame', data = cell_metadata)
fd <- new('AnnotatedDataFrame', data = gene_annotation)
cds <- newCellDataSet(
  expression_matrix, 
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5, # 适用于 UMI 或单细胞 RNA-seq 数据
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 10)
ordering_genes <- VariableFeatures(Microglia)
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, 
                       max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
p1<-plot_cell_trajectory(cds, color_by = "Pseudotime")
p2<-plot_cell_trajectory(cds, 
                         color_by = "State") + 
  ggtitle("")
p1/p2
pData(cds)$`RNA_snn_res.0.3` <- Microglia@meta.data[rownames(pData(cds)), "RNA_snn_res.0.3"]
plot_cell_trajectory(cds, color_by = "RNA_snn_res.0.3")

df <- pData(cds)
ggplot(df, aes(x = Pseudotime, y = cell_type, fill = `RNA_snn_res.0.3`)) +
  ggridges::geom_density_ridges(alpha = 0.7) +
  theme_classic() +
  labs(title = "Cell Type Distribution along Pseudotime", x = "Pseudotime", y = "Annotated Clusters")



diff_test_res <- differentialGeneTest(cds[ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                      cores = 4)
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
top_genes <- head(sig_gene_names, 5)
plot_genes_in_pseudotime(cds[top_genes,], color_by = "Pseudotime")






######Oligodendrocyte hdWGCNA######
rm(list = ls())
library(hdWGCNA)
gc()
Oligodendrocyte<-readRDS('../Oligodendrocyte.rds')
Oligodendrocyte <- FindVariableFeatures(Oligodendrocyte, 
                                        selection.method = "vst", 
                                        nfeatures = 1000)
DefaultAssay(Oligodendrocyte) <- "RNA"

obj.sub <- Oligodendrocyte
#obj.sub<-subset(obj.sub, subset = seurat_clusters == '8')
obj.sub@meta.data%>%colnames()
obj.sub <- FindVariableFeatures(obj.sub, nfeatures = 1000)
obj.sub <- ScaleData(obj.sub) 
obj.sub <- RunPCA(obj.sub, features = VariableFeatures(obj.sub))
Idents(obj.sub) <- "cell_type"
obj.sub <- SetupForWGCNA(
  obj.sub,
  gene_select = "variable", 
  wgcna_name = "Oligodendrocyte_Sub"
)
obj.sub <- MetacellsByGroups(
  obj.sub,
  group.by    = c('group',"cell_type"),
  ident.group = "group",
  reduction = 'pca',
  k = 25,
  max_shared = 10
)
obj.sub <- NormalizeMetacells(obj.sub)
obj.sub <- SetDatExpr(
  obj.sub,
  group_name = "Oligodendrocyte", # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)
obj.sub <- TestSoftPowers(
  obj.sub,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(obj.sub)
library(patchwork)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)
obj.sub <- ConstructNetwork(
  obj.sub,
  soft_power = 6,
  overwrite_tom = TRUE,
  tom_name = 'Oligodendrocyte' # name of the topoligical overlap matrix written to disk
)
pdf('../Dendrogram_Oligodendrocyte.pdf')
PlotDendrogram(obj.sub, main='Oligodendrocyte hdWGCNA Dendrogram')
dev.off()
Oligodendrocyte_tree<-PlotDendrogram(obj.sub, main='Oligodendrocyte hdWGCNA Dendrogram')

obj.sub <- ModuleEigengenes(
  obj.sub,
  group.by.vars = "group" # 这里填入你 meta.data 中用于批次校正或区分的变量，如 group
)
obj.sub <- ModuleConnectivity(
  obj.sub,
  group.by = 'cell_type', group_name = 'Oligodendrocyte'
)
# rename the modules
obj.sub <- ResetModuleNames(
  obj.sub,
  new_name = "Oligodendrocyte-M"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(obj.sub, ncol=2)
pdf('../kMEs_Oligodendrocyte.pdf')
p
dev.off()
plot_list <- ModuleFeaturePlot(
  obj.sub,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
pdf('../umap_Oligodendrocyte_hdwgcna.pdf')
wrap_plots(plot_list, ncol=3)
dev.off()

hub_df <- GetHubGenes(obj.sub, n_hubs = 5)
modules <- GetModules(obj.sub)
M1_gene<-modules$gene_name[modules$module=='Oligodendrocyte-M1']
M2_gene<-modules$gene_name[modules$module=='Oligodendrocyte-M2']
M3_gene<-modules$gene_name[modules$module=='Oligodendrocyte-M3']
#target_mod <- modules[modules$gene_name == "PPARG", ]
#hub_genes <- GetHubGenes(obj.sub, n_hubs = 30)
#M1_gene<-modules$gene_name[modules$module=='Oligodendrocyte-M1']
pvalueCutoff <- 0.05
library(clusterProfiler)
library(org.Hs.eg.db)
go_result_1 <- enrichGO(gene = M1_gene, 
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "all", 
                        pvalueCutoff = pvalueCutoff)

go_result_2 <- enrichGO(gene = M2_gene, 
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "all", 
                        pvalueCutoff = pvalueCutoff)

go_result_3 <- enrichGO(gene = M3_gene, 
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "all", 
                        pvalueCutoff = pvalueCutoff)
dotplot(go_result_1)/dotplot(go_result_2)
GO_ME1<-dotplot(go_result_1)
ME_df <- obj.sub@misc$Oligodendrocyte_Sub$wgcna_net$multiMEs$Oligodendrocyte_Sub$data
meta_mc <- obj.sub@misc$Oligodendrocyte_Sub$wgcna_metacell_obj@meta.data
plot_df_1 <- data.frame(
  ME_M1 = ME_df$MEblue,
  group = meta_mc$group
)

plot_df_2 <- data.frame(
  ME_M2 = ME_df$MEturquoise,
  group = meta_mc$group
)

plot_df_3 <- data.frame(
  ME_M3 = ME_df$MEbrown,
  group = meta_mc$group
)

library(ggpubr)
p_M1<-ggplot(plot_df_1, aes(x = group, y = ME_M1, fill = group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",   # 非参数检验，单细胞/宏基因组都推荐
    label = "p.signif",        # 显示 *, **, ***
    label.x = 1.5
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Module activity across groups",
    y = "Module eigengene",
    x = NULL
  )+
  theme(legend.position = 'none')

p_M2<-ggplot(plot_df_2, aes(x = group, y = ME_M2, fill = group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",   # 非参数检验，单细胞/宏基因组都推荐
    label = "p.signif",        # 显示 *, **, ***
    label.x = 1.5
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Module activity across groups",
    y = "Module eigengene",
    x = NULL
  )

p_M3<-ggplot(plot_df_3, aes(x = group, y = ME_M3, fill = group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",   # 非参数检验，单细胞/宏基因组都推荐
    label = "p.signif",        # 显示 *, **, ***
    label.x = 1.5
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Module activity across groups",
    y = "Module eigengene",
    x = NULL
  )

pdf('../Oligodendrocyte_activity.pdf')
p_M1+p_M2+p_M3
dev.off()

pdf('../Oligodendrocyte_activity_m1.pdf',width = 6,height = 6)
p_M1
dev.off()

library(ReactomePA)
gene_id<-bitr(M1_gene, 
              fromType = "SYMBOL", 
              toType = "ENTREZID", 
              OrgDb = org.Hs.eg.db)
res_reactome <- enrichPathway(gene = gene_id$ENTREZID, organism = "human")
react_M1<-dotplot(res_reactome)
react_M1
pdf('reactome_Oligodendrocyte.pdf')
dotplot(res_reactome)
dev.off()

targetgene<-c('RASGRF1','ACTN2')
pdf('../2gene.pdf')
FeaturePlot(Oligodendrocyte, features = targetgene)
dev.off()
cor(obj.sub@misc[["Oligodendrocyte_Sub"]][["hMEs"]][["Oligodendrocyte-M1"]], GetAssayData(obj.sub)["RASGRF1",])
cor(obj.sub@misc[["Oligodendrocyte_Sub"]][["hMEs"]][["Oligodendrocyte-M1"]], GetAssayData(obj.sub)["ACTN2",])
cor(GetAssayData(obj.sub)["RASGRF1",], GetAssayData(obj.sub)["ACTN2",])


p_overlap <- FeaturePlot(
  obj.sub, 
  features = c("RASGRF1", "ACTN2"), 
  blend = TRUE, 
  blend.threshold = 0.1, # 调整阈值以过滤背景低表达噪音
  cols = c("lightgrey", "red", "green"), # 基础色、基因1色、基因2色
  order = TRUE # 让高表达的点显示在最上层
)
pdf('../overlap_Oligodendrocyte.pdf',width = 10,height = 5)
p_overlap
dev.off()

library(AUCell)
expr_matrix <- GetAssayData(Oligodendrocyte, slot = "data")
cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 4, plotStats = TRUE)


###########
library(ggrastr)
id<-go_result_1@result[["ID"]]
p_aucell_go_Oligodendrocyte<-list()
for(i in 1:length(id)){
  geneSets <- list(hdWGCNA_M1 = go_result_1@geneSets[[id[i]]])
  # 4. 计算 AUC 得分
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  # 提取 M1 的 AUC 分数并整合到 Seurat meta.data
  M1_auc_scores <- as.numeric(getAUC(cells_AUC)["hdWGCNA_M1", ])
  Oligodendrocyte$M1_AUC <- M1_auc_scores
  hME_M1 <- obj.sub@misc$Oligodendrocyte_Sub$hMEs$`Oligodendrocyte-M1`
  #common_cells <- intersect(colnames(Oligodendrocyte), rownames(obj.sub@misc$Oligodendrocyte_Sub$hMEs))
  #hME_val <- hME_M1[common_cells]
  AUC_val <- Oligodendrocyte@meta.data[, "M1_AUC"]
  cor_res <- cor.test(hME_M1, AUC_val, method = "spearman")
  library(ggplot2)
  df_cor <- data.frame(hME = hME_M1, AUC = AUC_val)
  p_aucell_go_Oligodendrocyte[[i]]<-ggplot(df_cor, aes(x = hME, y = AUC)) +
    geom_point_rast(alpha = 0.3, size = 0.5, color = "steelblue") +
    geom_smooth(method = "lm", color = "red") +
    annotate("text", x = min(hME_M1), y = max(AUC_val), 
             label = paste0("rho = ", round(cor_res$estimate, 3)), hjust = 0) +
    theme_classic() +
    ggtitle(id[i])
  #print(p_aucell_go_Oligodendrocyte[[i]])
}
wrap_plots(p_aucell_go_Oligodendrocyte, ncol=1)


tmp<-res_reactome@result
tmp<-tmp[tmp$p.adjust<0.05,]


id_reac<-tmp$ID[c(10,6,8,1,2,3,4,5,7,9)]
p_aucell_reac_Oligodendrocyte<-list()
for(i in 1:length(id_reac)){
  geneSets <- list(hdWGCNA_M1 = res_reactome@geneSets[[id_reac[i]]])
  gene_id<-bitr(geneSets$hdWGCNA_M1, 
                fromType = "ENTREZID", 
                toType = "SYMBOL", 
                OrgDb = org.Hs.eg.db)
  geneSets$hdWGCNA_M1<-gene_id$SYMBOL
  # 4. 计算 AUC 得分
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  # 提取 M1 的 AUC 分数并整合到 Seurat meta.data
  M1_auc_scores <- as.numeric(getAUC(cells_AUC)["hdWGCNA_M1", ])
  Oligodendrocyte$M1_AUC <- M1_auc_scores
  hME_M1 <- obj.sub@misc$Oligodendrocyte_Sub$hMEs$`Oligodendrocyte-M1`
  #common_cells <- intersect(colnames(Oligodendrocyte), rownames(obj.sub@misc$Oligodendrocyte_Sub$hMEs))
  #hME_val <- hME_M1[common_cells]
  AUC_val <- Oligodendrocyte@meta.data[, "M1_AUC"]
  cor_res <- cor.test(hME_M1, AUC_val, method = "spearman")
  library(ggplot2)
  df_cor <- data.frame(hME = hME_M1, AUC = AUC_val)
  p_aucell_reac_Oligodendrocyte[[i]]<-ggplot(df_cor, aes(x = hME, y = AUC)) +
    geom_point_rast(alpha = 0.3, size = 0.5, color = "steelblue") +
    geom_smooth(method = "lm", color = "red") +
    annotate("text", x = min(hME_M1), y = max(AUC_val), 
             label = paste0("rho = ", round(cor_res$estimate, 3)), hjust = 0) +
    theme_classic() +
    ggtitle(id_reac[i])
  #print(p_aucell_go_Oligodendrocyte[[i]])
  print(cor_res)
}
library(patchwork)
wrap_plots(p_aucell_reac_Oligodendrocyte[[1]], 
           p_aucell_reac_Oligodendrocyte[[8]],
           p_aucell_reac_Oligodendrocyte[[9]],ncol = 3)

go1<-p_aucell_go_Oligodendrocyte[[1]]
go2<-p_aucell_go_Oligodendrocyte[[2]]
go3<-p_aucell_go_Oligodendrocyte[[3]]

reac3<-p_aucell_reac_Oligodendrocyte[[3]]
reac4<-p_aucell_reac_Oligodendrocyte[[4]]

pdf('../aucell_cor.pdf',width = 6,height = 6)
(go1+go2+go3)/(reac3+reac4)
dev.off()

pdf('aucell_go.pdf',width = 6,height = 8)



pdf('../GO_react_M1.pdf',width = 6,height = 12)
grid.arrange(arrangeGrob(arrangeGrob(GO_ME1,react_M1,nrow=2,heights = c(1,2))))
dev.off()
pdf('../aucell_go_react.pdf',width = 6,height = 12)
wrap_plots(p_aucell_go_Oligodendrocyte, ncol = 3)/wrap_plots(p_aucell_reac_Oligodendrocyte, ncol = 3)+plot_layout(heights = c(1, 4))
dev.off()
ggsave('../aucell_go_react.pdf',width = 6,height=12)
pdf('../Dendrogram_Oligodendrocyte.pdf',width = 6,height = 6)
PlotDendrogram(obj.sub, main='Oligodendrocyte hdWGCNA Dendrogram')
dev.off()
pdf('../p_M1_Oligodendrocyte.pdf',width = 6,height = 6)
p_M1
dev.off()
#################






################
library(SCopeLoomR)
library(SCENIC)
loom_path <- "../scenic_out/out_SCENIC_Oligodendrocyte.loom"
loom <- open_loom(loom_path)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)
library(ComplexHeatmap)
cells_AUC <- getAUC(regulonAUC)
cells_AUC_zscore <- t(scale(t(cells_AUC)))


library(dplyr)
library(tidyr)
library(purrr)

# 1. 准备数据 (假设你的数据框名为 df_scenic)
# 行名为转录因子，列名为细胞名 (如 AD1_cell_1, NC1_cell_1)
df_scenic_t <- as.data.frame(t(cells_AUC)) # 转置后列为 TF，行为细胞

# 2. 根据列名前缀添加分组标签
df_scenic_t$group <- NA
df_scenic_t$group[grep('AD',rownames(df_scenic_t))]<-'AD'
df_scenic_t$group[grep('NC',rownames(df_scenic_t))]<-'NC'
table(df_scenic_t$group)

# 3. 过滤掉无法识别分组的细胞 (如有)
df_scenic_t <- df_scenic_t %>% filter(!is.na(group))

# 4. 批量运行 Wilcoxon 检验
# 遍历除了 group 以外的所有列 (即所有 TF)
de_results <- map_dfr(colnames(df_scenic_t)[-ncol(df_scenic_t)], function(tf) {
  test <- wilcox.test(df_scenic_t[[tf]] ~ df_scenic_t$group)
  
  # 计算均值以确定差异倍数 (Log2 Fold Change)
  avg_ad <- mean(df_scenic_t[[tf]][df_scenic_t$group == "AD"])
  avg_nc <- mean(df_scenic_t[[tf]][df_scenic_t$group == "NC"])
  log2FC <- log2((avg_ad + 1e-9) / (avg_nc + 1e-9)) # 添加微小值避免除以0
  
  data.frame(
    TF = tf,
    Mean_AD = avg_ad,
    Mean_NC = avg_nc,
    Log2FC = log2FC,
    p_val = test$p.value
  )
})

# 5. 校正 P 值并筛选显著 TF
de_results$p_adj <- p.adjust(de_results$p_val, method = "BH")
de_results <- de_results %>% arrange(p_adj)

# 查看显著性差异的 TF (如 p_adj < 0.05 且 Log2FC > 0 为 AD 高表达)
significant_tfs <- de_results %>% filter(p_adj < 0.05)
sig_tf<-cells_AUC[significant_tfs$TF,]
plot_matrix <- as.matrix(sig_tf)
plot_matrix <- t(scale(t(plot_matrix)))
group_info <- data.frame(
  Cell = colnames(plot_matrix),
  Group = ifelse(grepl("^AD1", colnames(plot_matrix)), "AD", "NC")
)
col_fun = list(Group = c("AD" = "#E41A1C", "NC" = "#377EB8"))
top_ann <- HeatmapAnnotation(
  Group = group_info$Group,
  col = col_fun,
  show_annotation_name = TRUE
)
library(circlize)
pdf('../heatmap_Oligodendrocyte.pdf')
Heatmap(plot_matrix, 
        name = "Relative Activity", 
        top_annotation = top_ann,
        column_split = group_info$Group, # 根据分组切开热图，增加可读性
        show_column_names = FALSE,       # 细胞太多，隐藏列名
        show_row_names = TRUE,
        cluster_columns = TRUE,          # 组内进行聚类
        cluster_rows = TRUE,
        column_title = "Differential TF Activity (AD vs NC)",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
dev.off()
##############
library(CIBERSORT)
rm(list = ls())
gc()
###CIBERSORTX数据准备
Epi <- readRDS("../combined1.rds")

library(Seurat)
library(dplyr)
library(readr)
library(data.table)
library(stringr)
library(stringr)

## 2.2 获得矩阵
count_matrix <- GetAssayData(Epi, assay = "RNA", slot = "counts")
#count_matrix <- as.data.frame(count_matrix)
count_matrix[1:3,1:12]

fwrite(count_matrix%>%as.data.frame(),file = "D:/microbiome/AD/data/scRNA/Step1.CibersortX_scRNAassay.txt",sep = "\t",row.names = F)

#load("F:\\work\\xiehe\\TCGA-STAD\\cancer_tpm.RData")
#fwrite(cancer,file = "Step1.CibersortX_TCGA_matrix_TPM.txt",sep = "\t",row.names = T)

#制作细胞类型的系数表
X <- AverageExpression(Epi)[[1]]    #这就是我们要的细胞类型系数表
#write.table(X,"sig.txt",sep = "\t",col.names = T,row.names = T)     #保存为sig.txt

results=CIBERSORT("sig.txt","test.txt", perm=1000, QN=TRUE)   #还是用那一个test.txt表达矩阵，同样运行

source("Cibersort.R")

library(tibble)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(e1071)             #必须加载，因为后面反卷积的算法是基于这个包里的数据

scRNA_harmony <- readRDS("scRNAsub.rds")  #导入单细胞数据
#先Idents
Idents(scRNA_harmony) <- "celltype"
#制作细胞类型的系数表
X <- AverageExpression(scRNA_harmony)[[1]]    #这就是我们要的细胞类型系数表
write.table(X,"sig.txt",sep = "\t",col.names = T,row.names = T)     #保存为sig.txt

results=CIBERSORT("sig.txt","Step1.CibersortX_TCGA_matrix_TPM.txt", perm=1000, QN=TRUE)   #还是用那一个test.txt表达矩阵，同样运行

##############
rm(list = ls())
gc()
library(Scissor)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(tibble)
base_dir<-"D:/microbiome/AD/data/scRNA/GSE157827_RAW"
setwd(base_dir)
combined<-readRDS('../combined1.rds')
names(combined)
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
phenotype <- as.numeric(meta$group == "AD")
names(phenotype) <- rownames(meta)
bulk_data <- log2(edgeR::cpm(exp) + 1)
set.seed(123)

infos <- Scissor(bulk_dataset = as.matrix(bulk_data), 
                 sc_dataset = combined, 
                 phenotype = phenotype, 
                 tag = "AD_Scissor",
                 alpha = 0.05, 
                 family = "binomial")
write.csv(bulk_data,'D:/microbiome/AD/data/scRNA/bulk_data.csv')
write.csv(meta,'D:/microbiome/AD/data/scRNA/meta.csv')

infos<-readRDS('D:/microbiome/AD/data/scRNA/scissor_out_Oligodendrocyte.rds')
Scissor_res <- rep("Background", ncol(combined))
names(Scissor_res) <- colnames(combined)
Scissor_res[infos$Scissor_pos] <- "Scissor+"
Scissor_res[infos$Scissor_neg] <- "Scissor-"
combined$Scissor_res <- Scissor_res
table(combined$Scissor_res)
my_cols1 <- c("Scissor+" = "#E41A1C", "Scissor-" = "#D3D3D3", "Background" = "#D3D3D3")
my_cols2 <- c("Scissor+" = "#D3D3D3", "Scissor-" = "#377EB8", "Background" = "#D3D3D3")

scissor_AD<-DimPlot(combined, group.by = "Scissor_res", cols = my_cols1, pt.size = 0.5, order = T) +
  #theme_dr() + # 或者使用 theme_classic()
  labs(title = "Scissor Identified Cells",
       subtitle = paste0("Scissor+: ", length(infos$Scissor_pos), " cells; ",
                         "Scissor-: ", length(infos$Scissor_neg), " cells"))
scissor_H<-DimPlot(combined, group.by = "Scissor_res", cols = my_cols2, pt.size = 0.5, order = T)+
  ##theme_dr() + # 或者使用 theme_classic()
  labs(title = "")
  #     subtitle = paste0("Scissor+: ", length(infos$Scissor_pos), " cells; ",
  #                       "Scissor-: ", length(infos$Scissor_neg), " cells"))
scissor_AD+scissor_H






infos<-readRDS('D:/microbiome/AD/data/scRNA/scissor_out_Oligodendrocyte.rds')
Scissor_res <- rep("Background", ncol(combined))
names(Scissor_res) <- colnames(combined)
Scissor_res[infos$Scissor_pos] <- "Scissor+"
#Scissor_res[infos$Scissor_neg] <- "Scissor-"
combined$Scissor_res <- Scissor_res
Microglia<-subset(combined, subset = cell_type == 'Microglia')
Microglia_scissor1<-DimPlot(Microglia, group.by = "Scissor_res", cols = my_cols1, pt.size = 0.5, order = T) +
  #theme_dr() + # 或者使用 theme_classic()
  labs(title = "Scissor Identified Cells")+
       #subtitle = paste0("Scissor+: ", length(infos$Scissor_pos), " cells; ",
        #                 "Scissor-: ", length(infos$Scissor_neg), " cells"))+
  xlim(c(3,12))+
  ylim(c(3,10))





infos<-readRDS('D:/microbiome/AD/data/scRNA/scissor_out_Oligodendrocyte.rds')
Scissor_res <- rep("Background", ncol(combined))
names(Scissor_res) <- colnames(combined)
#Scissor_res[infos$Scissor_pos] <- "Scissor+"
Scissor_res[infos$Scissor_neg] <- "Scissor-"
combined$Scissor_res <- Scissor_res
Microglia<-subset(combined, subset = cell_type == 'Microglia')
Microglia_scissor2<-DimPlot(Microglia, group.by = "Scissor_res", cols = my_cols2, pt.size = 0.5, order = T) +
  #theme_dr() + # 或者使用 theme_classic()
  labs(title = "Scissor Identified Cells")+
  #subtitle = paste0("Scissor+: ", length(infos$Scissor_pos), " cells; ",
  #                 "Scissor-: ", length(infos$Scissor_neg), " cells"))+
  xlim(c(3,12))+
  ylim(c(3,10))
pdf('../Microglia_scissor.pdf',width = 8,height = 4)
Microglia_scissor1+Microglia_scissor2
dev.off()
##########
infos<-readRDS('D:/microbiome/AD/data/scRNA/scissor_out_Oligodendrocyte.rds')
Scissor_res <- rep("Background", ncol(combined))
names(Scissor_res) <- colnames(combined)
Scissor_res[infos$Scissor_pos] <- "Scissor+"
Scissor_res[infos$Scissor_neg] <- "Scissor-"
combined$Scissor_res <- Scissor_res
Microglia<-subset(combined, subset = cell_type == 'Microglia')
seurat.list <- SplitObject(obj.sub, split.by = "group")
rm(combined)
rm(cell_label)
gc()
CellChatDB <- CellChatDB.human

data.input1 <- GetAssayData(seurat.list[["AD"]], slot = "data", assay = "RNA")
meta1 <- seurat.list[["AD"]]@meta.data
meta1$RNA_snn_res.0.2 <- paste0("C", meta1$RNA_snn_res.0.2)
cellchat1 <- createCellChat(object = data.input1, meta = meta1, group.by = "RNA_snn_res.0.2")
cellchat1@DB <- CellChatDB


data.input2 <- GetAssayData(seurat.list[["HC"]], slot = "data", assay = "RNA")
meta2 <- seurat.list[["HC"]]@meta.data
meta2$RNA_snn_res.0.2 <- paste0("C", meta2$RNA_snn_res.0.2)
cellchat2 <- createCellChat(object = data.input2, meta = meta2, group.by = "RNA_snn_res.0.2")
cellchat2@DB <- CellChatDB

cellchat1 <- subsetData(cellchat1)
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
cellchat1 <- computeCommunProb(cellchat1)
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
#saveRDS(cellchat1,'../cellchat_AD.rds')

cellchat2 <- subsetData(cellchat2)
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- computeCommunProb(cellchat2)
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
#saveRDS(cellchat2,'../cellchat_HC.rds')

cellchat1<-readRDS('../cellchat_AD.rds')
cellchat2<-readRDS('../cellchat_HC.rds')
cellchat1@options$ident <- "AD"
cellchat2@options$ident <- "HC"

object.list <- list(AD = cellchat1, HC = cellchat2)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
pdf('../cellchat_bar.pdf')
compareInteractions(cellchat_merged, show.legend = F, group = c(1, 2))
dev.off()
pdf('../cellchat_net.pdf',width = 11.5,height = 14)
netVisual_diffInteraction(cellchat_merged, 
                          weight.scale = TRUE,
                          measure = 'count',
                          comparison = c(1, 2),
                          label.edge = F)
dev.off()

pdf('../cellchat_heatmap.pdf',width = 10,height = 6)
netVisual_heatmap(cellchat1)+netVisual_heatmap(cellchat2)
dev.off()

cellchat_AD<-netVisual_circle(cellchat1@net$count, weight.scale = T, label.edge = F, title.name = "AD Total Interactions")
cellchat_HC<-netVisual_circle(cellchat2@net$count, weight.scale = T, label.edge = F, title.name = "HC Total Interactions")
pdf('../cellchat_AD.pdf',width = 9,height = 9)
cellchat_AD
dev.off()
pdf('../cellchat_HC.pdf',width = 9,height = 9)
cellchat_HC
dev.off()


##########

Oligodendrocyte<-subset(combined, subset = cell_type == 'Oligodendrocyte')
DimPlot(Oligodendrocyte, group.by = "Scissor_res", cols = my_cols, pt.size = 0.5, order = T) +
  #theme_dr() + # 或者使用 theme_classic()
  labs(title = "Scissor Identified Cells (AD vs Control)",
       subtitle = paste0("Scissor+: ", length(infos$Scissor_pos), " cells; ",
                         "Scissor-: ", length(infos$Scissor_neg), " cells"))
#########
Oligodendrocyte<-subset(combined, subset = cell_type == 'Oligodendrocyte')
scissor_Oligodendrocyte1<-DimPlot(Oligodendrocyte, group.by = "Scissor_res", cols = my_cols1, pt.size = 0.1, order = T,alpha = 0.5)
scissor_Oligodendrocyte2<-DimPlot(Oligodendrocyte, group.by = "Scissor_res", cols = my_cols2, pt.size = 0.1, order = T,alpha = 0.5)
scissor_Oligodendrocyte1+scissor_Oligodendrocyte2


Oligodendrocyte <- FindVariableFeatures(Oligodendrocyte, 
                                  selection.method = "vst", 
                                  nfeatures = 2000)
DefaultAssay(Oligodendrocyte) <- "RNA"
obj.sub <- Oligodendrocyte
obj.sub@meta.data%>%colnames()
obj.sub <- ScaleData(obj.sub) 
obj.sub <- RunPCA(obj.sub, features = VariableFeatures(obj.sub))
obj.sub <- FindClusters(obj.sub, resolution = 0.2)
obj.sub<-subset(obj.sub,RNA_snn_res.0.2 != 6)
Oligodendrocyte_02<-DimPlot(obj.sub, group.by = "RNA_snn_res.0.2", pt.size = 0.1, order = T,alpha = 0.5)
Idents(obj.sub) <- "RNA_snn_res.0.2"
all_markers <- FindAllMarkers(obj.sub,min.pct = 0.25, only.pos = TRUE,logfc.threshold = 0.25)
all_markers <- all_markers[-grep('^MT-',all_markers$gene),]
markers_top3<-all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)
Oligodendrocyte_02
dot_Oligodendrocyte<-DotPlot(obj.sub, features = unique(markers_top3$gene), dot.scale = 5) + RotatedAxis()
heatmap_Oligodendrocyte<-DoHeatmap(obj.sub, features = markers_top3$gene) + NoLegend()


pdf('../Oligodendrocyte_scissor.pdf',width = 16,height = 10)
grid.arrange(arrangeGrob(arrangeGrob(scissor_Oligodendrocyte1,scissor_Oligodendrocyte2, Oligodendrocyte_02,nrow=1),arrangeGrob(heatmap_Oligodendrocyte, dot_Oligodendrocyte,nrow=2,heights = c(1,1.5)),heights = c(1,2)))
dev.off()
library(patchwork)
pdf('../Oligodendrocyte_scissor_top3.pdf',width = 16,height = 3.3)
FeaturePlot(obj.sub, features = markers_top3$gene[markers_top3$cluster==5], max.cutoff = 3,reduction = "umap",ncol = 3)
dev.off()
pdf('../umap_scissor_Oligodendrocyte.pdf',width = 16,height = 3.3)
scissor_Oligodendrocyte1+scissor_Oligodendrocyte2+Oligodendrocyte_02
dev.off()

clustr5_gene<-all_markers[all_markers$cluster==5,]
diff_geneset<-clustr5_gene[order(clustr5_gene$avg_log2FC,decreasing = T),]
ranks <- diff_geneset$avg_log2FC
names(ranks)<-diff_geneset$gene
library(clusterProfiler)
Go_gseresult <- gseGO(ranks, 
                      'org.Hs.eg.db', 
                      keyType = "SYMBOL", 
                      ont="all", 
                      nPerm = 1000, 
                      minGSSize = 10, 
                      maxGSSize = 1000, 
                      pAdjustMethod = "none",
                      pvalueCutoff=0.05,
)
library(enrichplot)
gsea_1<-gseaplot2(Go_gseresult, 1:5, pvalue_table = F)
pdf('../GSEA_c5_Oligodendrocyte.pdf',width = 16,height = 6)
gsea_1
dev.off()



##########
Microglia <- FindVariableFeatures(Microglia, 
                                  selection.method = "vst", 
                                  nfeatures = 2000)
DefaultAssay(Microglia) <- "RNA"
obj.sub <- Microglia
obj.sub@meta.data%>%colnames()
obj.sub <- FindVariableFeatures(obj.sub, nfeatures = 2000)
obj.sub <- ScaleData(obj.sub) 
obj.sub <- RunPCA(obj.sub, features = VariableFeatures(obj.sub))
Idents(obj.sub) <- "cell_type"
obj.sub <- FindClusters(obj.sub, resolution = 0.2)
obj.sub<-subset(obj.sub,RNA_snn_res.0.2 != 4)







all_markers <- FindAllMarkers(obj.sub,min.pct = 0.25, only.pos = TRUE,logfc.threshold = 0.25)
all_markers <- all_markers[-grep('^MT-',all_markers$gene),]
markers_top5<-all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)
my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628")
Microglia_0.2<-DimPlot(obj.sub, group.by = 'RNA_snn_res.0.2', order = T, cols = my_colors) +
  xlim(c(3, 12)) +
  ylim(c(3, 10))
Microglia_0.2
Microglia_heatmap_0.2<-DoHeatmap(obj.sub, features = markers_top5$gene) + NoLegend()
Microglia_dot_0.2<-DotPlot(obj.sub, features = unique(markers_top5$gene), dot.scale = 5) + RotatedAxis()
library(gridExtra)
pdf('../Microglia_scissor_and_cluster.pdf',width = 12,height = 10)
grid.arrange(arrangeGrob(arrangeGrob(Microglia_scissor1,Microglia_scissor2, Microglia_0.2,nrow=1),arrangeGrob(Microglia_heatmap_0.2, Microglia_dot_0.2,nrow=2,heights = c(1,1.5)),heights = c(1,2)))
dev.off()
library(gridExtra)
clustr1_gene<-all_markers[all_markers$cluster==1,]
diff_geneset<-clustr1_gene[order(clustr1_gene$avg_log2FC,decreasing = T),]
ranks <- diff_geneset$avg_log2FC
names(ranks)<-diff_geneset$gene
library(clusterProfiler)
Go_gseresult <- gseGO(ranks, 
                      'org.Hs.eg.db', 
                      keyType = "SYMBOL", 
                      ont="all", 
                      nPerm = 1000, 
                      minGSSize = 10, 
                      maxGSSize = 1000, 
                      pAdjustMethod = "none",
                      pvalueCutoff=0.05,
)
library(enrichplot)
gsea_1<-gseaplot2(Go_gseresult, 1:5, pvalue_table = F,base_size = 8)
gsea_1

clustr2_gene<-all_markers[all_markers$cluster==2,]
diff_geneset<-clustr2_gene[order(clustr2_gene$avg_log2FC,decreasing = T),]
ranks <- diff_geneset$avg_log2FC
names(ranks)<-diff_geneset$gene
library(clusterProfiler)
Go_gseresult <- gseGO(ranks, 
                      'org.Hs.eg.db', 
                      keyType = "SYMBOL", 
                      ont="all", 
                      nPerm = 1000, 
                      minGSSize = 10, 
                      maxGSSize = 1000, 
                      pAdjustMethod = "none",
                      pvalueCutoff=0.05,
)
library(enrichplot)
gsea_2<-gseaplot2(Go_gseresult, 1:5, pvalue_table = F,base_size = 8)
library(aplot)
pdf('../gsea_mi_c1_c2.pdf',width = 12,height=6)
plot_list(gsea_1, gsea_2, ncol = 2)
dev.off()
gsea_plot<-plot_list(gsea_1, gsea_2, ncol = 1)

ggsave("../gsea_mi_c1_c2.pdf", plot = gsea_plot, width = 10, height = 12)
