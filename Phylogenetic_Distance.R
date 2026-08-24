library(taxize)
manifest<-read.csv('D:/microbiome/APOLLO/manifest.csv',row.names = 1)
genus<-manifest$genus
out<-list()
for(i in 1:length(genus)){
  print(i)
  out[[i]]<-get_uid(genus[i])
}
out<-get_uid(manifest$genus)
uid<-matrix(0,length(out),2)%>%as.data.frame()
colnames(uid)<-c('genus','uid')
for(i in 1:length(out)){
  uid$genus[i]<-genus[i]
  uid$uid[i]<-out[[i]][1]
}
uid<-na.omit(uid)
write.csv(uid,'D:/microbiome/AD/data/uid_agora2.csv')


out<-classification(uid$genus, db='ncbi')
tr<-class2tree(out)
distance_mat<-tr[["distmat"]]%>%as.matrix()


