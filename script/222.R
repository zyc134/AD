gutmgene<-read.csv('D:/microbiome/AD/data/gutmgene.csv')
gutmgene$Metabolite..ID.<-sub("\\([^()]*\\)$", "", gutmgene$Metabolite..ID.)
#gutmgene<-gutmgene[gutmgene$Metabolite..ID.%in%mets_data$V2,]

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