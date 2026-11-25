getdata <- function(runid_path = "./data/runid_test.csv", save_path = "./", session = NULL) {
  print(paste0('The data will be saved in ', save_path))
  runid <- read.csv(runid_path, header = FALSE)
  dir.create(file.path(save_path, 'species'), showWarnings = FALSE)
  dir.create(file.path(save_path, 'genus'), showWarnings = FALSE)
  
  total <- nrow(runid)
  
  # Initialize progress bar based on the context (Shiny app or console)
  if (!is.null(session)) {
    progress <- shiny::Progress$new(session, min = 0, max = total)
    progress$set(message = "Downloading data", value = 0)
    on.exit(progress$close())
  } else {
    progress <- utils::txtProgressBar(min = 0, max = total, style = 3)
    on.exit(close(progress))
  }
  
  genus_list <- list()
  species_list <- list()
  
  for (i in 1:total) {
    query <- POST("https://gmrepo.humangut.info/api/getFullTaxonomicProfileByRunID", 
                  body = list("run_id" = runid[i, 1]), 
                  encode = "json",
                  config = config(ssl_verifypeer = FALSE))
    retrieved_contents <- content(query)
    data <- fromJSON(xml_text(retrieved_contents))
    
    # Process genus
    if (!is.null(data$genus$relative_abundance)) {
      genus <- data$genus$relative_abundance %>% as.data.frame()
      if (ncol(genus) > 0) {
        colnames(genus) <- runid[i, 1]
        genus$name <- data$genus$scientific_name
        genus_name <- unique(data$genus$scientific_name)
        genus <- genus[match(genus_name, genus$name), ]
        rownames(genus) <- genus_name
        genus$name <- NULL
        genus_list[[i]] <- genus
      }
    }
    
    # Process species
    if (!is.null(data$species$relative_abundance)) {
      species <- data$species$relative_abundance %>% as.data.frame()
      if (ncol(species) > 0) {
        colnames(species) <- runid[i, 1]
        species$name <- data$species$scientific_name
        species_name <- unique(data$species$scientific_name)
        species <- species[match(species_name, species$name), ]
        rownames(species) <- species_name
        species$name <- NULL
        species_list[[i]] <- species
      }
    }
    
    # Update progress bar
    if (!is.null(session)) {
      progress$inc(1, detail = paste("Processing", i, "of", total))
    } else {
      utils::setTxtProgressBar(progress, i)
    }
  }
  
  # Combine genus and species data
  merge_fill <- function(data_list) {
    # Get all unique row names
    all_rows <- unique(unlist(lapply(data_list, rownames)))
    # Align all data frames to have the same row names
    aligned_list <- lapply(data_list, function(df) {
      if (is.null(df)) return(data.frame(matrix(0, nrow = length(all_rows), ncol = 0)))
      missing_rows <- setdiff(all_rows, rownames(df))
      if (length(missing_rows) > 0) {
        # Handle empty data frames correctly
        col_count <- ifelse(is.null(ncol(df)), 0, ncol(df))
        missing_df <- data.frame(matrix(0, nrow = length(missing_rows), ncol = col_count))
        if (col_count > 0) colnames(missing_df) <- colnames(df)
        rownames(missing_df) <- missing_rows
        df <- rbind(df, missing_df)
      }
      df[all_rows, , drop = FALSE]
    })
    # Combine using bind_cols
    result <- Reduce(function(x, y) dplyr::bind_cols(x, y), aligned_list)
    rownames(result) <- all_rows
    result
  }
  
  if (length(genus_list) > 0) {
    merged_df_genus <- merge_fill(genus_list)
    merged_df_genus <- merged_df_genus[!rownames(merged_df_genus) %in% 'Unknown', ]
    write.csv(merged_df_genus, paste0(save_path, '/genus/abundance_genus.csv'))
  }
  
  if (length(species_list) > 0) {
    merged_df_species <- merge_fill(species_list)
    merged_df_species <- merged_df_species[!rownames(merged_df_species) %in% 'Unknown', ]
    write.csv(merged_df_species, paste0(save_path, '/species/abundance_species.csv'))
  }
  
  print('Download completed.')
}
install_required_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

getdata('D:/microbiome/AD/user_selected_unique_run_ids_9sVBlwSCY0.csv','D:/microbiome/AD/data/')
getdata('D:/microbiome/AD/runid_H.csv','D:/microbiome/AD/data/Healthy/')
required_packages <- c("shiny", "MatchIt", 'caret',"ggplot2", "dplyr", "leaflet", 'tibble',"igraph", 'dplyr','tidyr',"httr", "jsonlite", "xml2", "shinyFiles", 'data.table',"progress")
install_required_packages(required_packages)
#######
#rm(list = ls())
abundance_D<-read.csv('D:/microbiome/AD/data/AD/genus/abundance_genus.csv',row.names = 1)
colSums(abundance_D)%>%sort()%>%plot
colSums(abundance_D)%>%hist
temp<-colSums(abundance_D)
for(i in 1:ncol(abundance_D)){
  abundance_D[,i]<-abundance_D[,i]/temp[i]
}
files<-list.files('D:/microbiome/micom/agora103_genus/890088c4-1a67-498b-8a54-3eb31ee06f63/data')
for(i in 1:length(files)){
  files[i]<-substr(files[i],1,nchar(files[i])-5)
}
abundance_D1<-abundance_D[intersect(rownames(abundance_D),files),]
colSums(abundance_D1)%>%hist
abundance_D1<-abundance_D1[,colSums(abundance_D1)>0.8]
abundance_D<-abundance_D[,colnames(abundance_D1)]
temp_D<-abundance_D1
temp_D[temp_D!=0]<-1
colSums(temp_D)%>%hist
abundance_H<-fread('D:/microbiome/AD/data/Healthy/genus/abundance_genus.csv')%>%as.data.frame()%>%column_to_rownames('V1')
colSums(abundance_H)%>%sort()%>%plot
colSums(abundance_H)%>%hist
abundance_H<-abundance_H[,colSums(abundance_H)<=100]
abundance_H1<-abundance_H[intersect(rownames(abundance_H),files),]
colSums(abundance_H1)%>%hist()
abundance_H1<-abundance_H1[,colSums(abundance_H1)>90&colSums(abundance_H1)<95]
temp_H<-abundance_H1
temp_H[temp_H!=0]<-1
colSums(temp_H)%>%hist
abundance_H<-abundance_H[,colnames(abundance_H1)]
data<-merge(abundance_D1,abundance_H1,by = "row.names",all = TRUE)
data[is.na(data)]<-0
rownames(data)<-data$Row.names
data$Row.names<-NULL
abundance<-data
data<-matrix(0,nrow(abundance)*ncol(abundance),4)%>%as.data.frame()
colnames(data)<-c('id','sample_id','abundance','genus')
id<-rep(rownames(abundance),ncol(abundance))
sample_id<-c()
for(i in 1:ncol(abundance)){
  sample_id<-c(sample_id,rep(colnames(abundance)[i],nrow(abundance)))
}
data$id<-id
data$sample_id<-sample_id
data$abundance<-abundance %>% as.matrix() %>% as.numeric()
data$genus<-id
write.csv(data,'D:/microbiome/AD/data/data.csv')
#####
meta<-c(colnames(abundance))%>%as.data.frame()
meta$group<-0
meta$group[match(colnames(abundance_D),meta$.)]<-'AD'
meta$group[match(colnames(abundance_H),meta$.)]<-'Healthy'
meta<-meta%>%column_to_rownames('.')
#exchange_mat<-fread('D:/microbiome/AD/data/res_grow/exchanges.csv')%>%as.data.frame()
exchange_mat<-fread('D:/microbiome/AD/data/res_grow_wd/exchanges.csv')%>%as.data.frame()
exchange_mat$metabolite<-paste0(exchange_mat$metabolite,'_',exchange_mat$direction)
exchange_mat<-exchange_mat %>%
  group_by(sample_id, metabolite) %>%
  summarise(total_flux = sum(flux, na.rm = TRUE), .groups = "drop")
exchange_mat1<-exchange_mat %>%
  pivot_wider(names_from = sample_id, 
              values_from = total_flux,
              values_fill = list(total_flux = 0))%>%as.data.frame()%>%column_to_rownames('metabolite')
rownames(exchange_mat1)<-gsub('\\[|\\]','_',rownames(exchange_mat1))
exchange_mat1<-t(exchange_mat1)%>%as.data.frame()
data<-cbind(
  exchange_mat1
)
#data<-cbind(shannon,simpson,exchange_mat,topo)

data$group<-meta$group[match(rownames(data),rownames(meta))]
#data$Disease.name<-NULL
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
varImp(randomForestFit)
ggplot(data,aes(x=group,y=`gthrd_e__import`,fill=group))+
  geom_boxplot()
#####
meta_AD<-fread('D:/microbiome/AD/user_selected_run_list_W4as0en1k1.txt')%>%as.data.frame()%>%column_to_rownames('Run ID')
table(meta_AD$`Project ID`)
meta_AD<-meta_AD[colnames(abundance_D),]
table(meta_AD$`Project ID`)
tr_id<-c(rownames(meta_AD)[which(meta_AD$`Project ID`=='PRJNA633959'|meta_AD$`Project ID`=='PRJNA554111')])
ts_id<-c(rownames(meta_AD)[which(meta_AD$`Project ID`=='PRJNA496408')])
va_id<-c(rownames(meta_AD)[which(meta_AD$`Project ID`=='PRJNA489760'|meta_AD$`Project ID`=='PRJEB11419')])
set.seed(123)
x <- colnames(abundance_H)
shuffled_x <- sample(x)
total_ratio <- 116 + 32 + 33
group_sizes <- floor(c(116, 32, 33) / total_ratio * length(x))
group_sizes[3] <- length(x) - sum(group_sizes[1:2])
tr_id_H <- shuffled_x[1:group_sizes[1]]
ts_id_H <- shuffled_x[(group_sizes[1] + 1):(group_sizes[1] + group_sizes[2])]
va_id_H <- shuffled_x[(group_sizes[1] + group_sizes[2] + 1):length(x)]
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
pr_rf_va$label<-'Validation Set, AUC = 0.9143'
roc_curve<-roc(validation_data$group,pr_rf_va[,1])
auc(roc_curve)
roc_data_va <- data.frame(
  TPR = rev(roc_curve$sensitivities), # 真阳性率
  FPR = rev(1 - roc_curve$specificities) # 假阳性率
)
roc_data_va$label<-'Validation Set, AUC = 0.9143'
plot(roc_curve)
sp.obj <- ci.sp(roc_curve, sensitivities=seq(0, 1, .01), boot.n=100) 
plot(sp.obj, type="shape", col="#5bd1d7")
text(0,0.5,paste0("AUC:",round(roc_curve[["auc"]],4)))

df<-rbind(roc_data_tr,roc_data_ts,roc_data_va)
df$label<-factor(df$label,levels = c('Training Set, AUC = 1',
                                     'Testing Set, AUC = 1',
                                     'Validation Set, AUC = 0.9143'))
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
