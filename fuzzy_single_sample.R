set.seed(412)
require(doSNOW)
require(foreach)
library(cluster)

scaling<-T
#####GROUND CORES####

no_out_data<-lapply(no_out_data,function(x) {
  spectra<-x
  colnames(spectra)<-as.numeric(seq(from=500,to=2450,by=10))
  spectra})

no_out_details<-lapply(no_out_ground_DATA,function(x) {
  spectra<-x[1:6]
  spectra})

pr_spectra_DATA <-lapply(no_out_data,function(x) prcomp(x, center=T,scale=scaling)) 
pr_scores_DATA <- mapply(function(y,z) data.frame(group=y$Sample,z$x[,1:5]),no_out_details,pr_spectra_DATA,SIMPLIFY=F,USE.NAMES=F) # principal component scores I used 3 to avoid overfitting the data
pr_scores_DATAFRAME <- do.call(rbind,pr_scores_DATA)


num_clusters <- c(2:7)

#   i="gku17"

cl <-makeCluster(8)
setMKLthreads(1)
registerDoSNOW(cl)

fanny_data_by_sample_pc_euc_no_out <- sapply(levels(pr_scores_DATAFRAME$group),simplify = F,USE.NAMES = T,FUN = function(i){
lapply(seq(1,2,0.1),function(phi){
  test_data<-pr_scores_DATAFRAME[,2:ncol(pr_scores_DATAFRAME)][pr_scores_DATAFRAME$group==i,]
    foreach(nc=num_clusters,.packages='cluster') %dopar% fanny(x=test_data,k=nc,metric='euclidean',memb.exp=phi)
  })
})
stopCluster(cl)
setMKLthreads(4)
for(i in 1:length(fanny_data_by_sample_pc_euc_no_out)){
  names(fanny_data_by_sample_pc_euc_no_out[[i]])<-paste0('phi',seq(1,2,0.1))
}

for(i in 1:length(fanny_data_by_sample_pc_euc_no_out)){
  for(z in 1:length(seq(1,2,0.1)))
  names(fanny_data_by_sample_pc_euc_no_out[[i]][[z]])<-paste0('cluster',seq(2,7))
}



