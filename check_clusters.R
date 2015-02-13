extract_FPI<-function(fuzz_data) {
  fpis<-sapply(fuzz_data,FUN = function(y){
    sapply(y,FUN = function(z){
    memberships<-z$membership
    part_coef<-sum(memberships^2)/nrow(memberships)
    FPI<-1-(((ncol(memberships)*part_coef)-1)/(ncol(memberships)-1))
  },USE.NAMES=T)
},USE.NAMES=T)
  return(fpis)
}
#####CORES####
if(data_type!='pit'){
pdf('Plots/check_clusters.pdf',height = 8,width = 12)
par(mfrow = c(2,2))
for (i in names(fanny_data_by_sample_pc_euc_no_out)){
tmp_fpi <-extract_FPI(fanny_data_by_sample_pc_euc_no_out[[i]])
plot(2:7,tmp_fpi[,1],type='l',ylim = range(tmp_fpi[,-1]),main=i,col='white')
for (z in 2:length(names(fanny_data_by_sample_pc_euc_no_out[[1]]))) {
  lines(2:7,tmp_fpi[,z],type='l')
  text(x=2,y=tmp_fpi[sample(1,1),z],label = colnames(tmp_fpi)[z],cex=.7)
}

test_levels<-unique(no_out_details[[i]]$b.master_hor)
limits<-sapply(test_levels,function(x){
  sum(no_out_details[[i]]$b.master_hor==x)*2
})

hor_lim <- c(limits[1],sum(limits[1:2])) 
if (length(unique(no_out_details[[i]]$b.master_hor))==3) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3])) 
if (length(unique(no_out_details[[i]]$b.master_hor))==4) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]))
if (length(unique(no_out_details[[i]]$b.master_hor))==5) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5])) 
if (length(unique(no_out_details[[i]]$b.master_hor))==6) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5]),sum(limits[1:6])) 
if (length(unique(no_out_details[[i]]$b.master_hor))==7) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5]),sum(limits[1:6]),sum(limits[1:7])) 

plot(1:hor_lim[length(hor_lim)],col='white',xlab = 'Depth',yaxt='n',ylab='' )
for (z in 1:length(hor_lim)) {
  abline(v = hor_lim[z])
}

text(hor_lim-limits/2,50,labels = unique(no_out_details[[i]]$b.master_hor))
}


dev.off()
setwd('Plots/')
shell.exec(file = "check_clusters.pdf")
setwd('../')

}
#####PITS#### 
if (data_type=='pit'){
  
  pdf('Plots/check_clusters.pdf',height = 8,width = 12)
  par(mfrow = c(2,2))
  for (i in names(fanny_data_by_sample_pc_euc_no_out)){
    tmp_fpi <-extract_FPI(fanny_data_by_sample_pc_euc_no_out[[i]])
    plot(2:4,tmp_fpi[,1],type='l',ylim = range(tmp_fpi[,-1]),main=i,col='white')
    for (z in 2:length(names(fanny_data_by_sample_pc_euc_no_out[[1]]))) {
      lines(2:4,tmp_fpi[,z],type='l')
      text(x=2,y=tmp_fpi[sample(1,1),z],label = colnames(tmp_fpi)[z],cex=.7)
    }
    
    test_levels<-unique(no_out_details[[i]]$b.master_hor)
    limits<-sapply(test_levels,function(x){
      sum(no_out_details[[i]]$b.master_hor[1:(length(no_out_details[[i]]$b.master_hor)/2)]==x)*5
    })
    
    hor_lim <- c(limits[1],sum(limits[1:2])) 
    if (length(unique(no_out_details[[i]]$b.master_hor))==3) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3])) 
    if (length(unique(no_out_details[[i]]$b.master_hor))==4) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]))
    
    plot(1:hor_lim[length(hor_lim)],col='white',xlab = 'Depth',yaxt='n',ylab='' )
    for (z in 1:length(hor_lim)) {
      abline(v = hor_lim[z])
    }
    
    text(hor_lim-limits/2,10,labels = unique(no_out_details[[i]]$b.master_hor))
  }
  
  
  dev.off()
  setwd('Plots/')
  shell.exec(file = "check_clusters.pdf")
  setwd('../')
}
