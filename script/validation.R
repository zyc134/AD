library(curatedMetagenomicData)
library(ggplot2)
library(ggpubr)

load("data/sampleMetadata.rda")
hypertension_meta<-sampleMetadata[which(sampleMetadata$disease=='hypertension'),]
healthy_meta<-hypertension_meta[which(hypertension_meta$age_category=='adult'),]
healthy_meta<-sampleMetadata[which(sampleMetadata$disease=='healthy'),]
healthy_meta<-healthy_meta[which(healthy_meta$age_category=='adult'),]
healthy_meta<-healthy_meta[which(healthy_meta$country=='CHN'|healthy_meta$country=='ITA'|healthy_meta$country=='AUT'),]
healthy_meta<-healthy_meta[which(healthy_meta$BMI>18.5&healthy_meta$BMI<23.9),]
healthy_meta<-healthy_meta[which(healthy_meta$age>55&healthy_meta$age<70),]
hypertension_meta<-rbind(hypertension_meta,healthy_meta)

ab_hypertension<-hypertension_meta |>
  dplyr::filter(body_site == 'stool') |>
  returnSamples("relative_abundance",rownames = "long")
abundance_hypertension<-ab_hypertension@assays@data@listData[["relative_abundance"]]

ab<-abundance_hypertension|>as.data.frame()
taxa<-strsplit(rownames(ab),"\\|")
genus<-unlist(lapply(taxa, function(x){x[6]}))
genus<-gsub("g__","",genus)
ab<-data.frame(genus=genus,ab)
ab_genus<-aggregate(.~genus,data=ab,sum)
rownames(ab_genus)<-ab_genus$genus
ab_genus<-ab_genus[,-1]
c3<-read.csv("output/cluster_3.csv")[,1]
c3<-intersect(c3,rownames(ab_genus)
pc<-WGCNA::moduleEigengenes(t(ab_genus[c3,]), colors = rep("FERM",length(c3)))
pc<-pc$eigengenes
pc$group<-hypertension_meta$disease
wilcox.test(MEFERM~group,data = pc)

ggplot(pc,aes(group,MEFERM,fill = group))+
  geom_boxplot()+
  scale_fill_manual(values = c("blue4","red4"))+
  coord_cartesian(ylim = c(min(pc$MEFERM),0.1))+
  stat_compare_means(comparisons = list(c("hypertension","healthy")),
                     label.y = 0.0,tip.length = 0)+
  theme_classic()+labs(x="",y="FERM(eigenvector)")+
  theme(
    text = element_text(size = 16),axis.title.x = element_text(size = 14,face = "bold"),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold",size = 16),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.text.x = element_text(angle = 45,hjust=1),
    plot.title.position = "plot",title = element_text(size = 16),
    plot.title = element_text(vjust = 5,hjust = 0.04,face = "bold"),
  )
ggsave("output/supplementary3.pdf",width = 6,height = 6)


