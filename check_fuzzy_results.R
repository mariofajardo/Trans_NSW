#####CORES####
if(data_type!='pit'){
# ##dunn's coefficient
# a<-c()
# for (i in names(fuzzy_data)) a<-c(a,fuzzy_data[[i]]$Clustering$coeff[1])
# hist(a,main= 'Dunn\'s coefficient for all the soil core clusters')
extract_CI<-function(fuzz_data) {
  aaply(fuzz_data$Clustering$membership,1,function(y){
    max_mem<-y[which.max(y)]
    sec_max<-y[order(y,decreasing = T)==2]
    CI <- sec_max/max_mem
  })
}
testing_CI<-list()
for (i in names(fuzzy_data)){
testing_CI[[i]]<-extract_CI(fuzzy_data[[i]])
}


extract_entropy<-function(fuzz_data) {
  aaply(fuzz_data$Clustering$membership,1,function(y){
    entropy<-abs(1/log(fuzz_data$Clustering$k.crisp)*sum(y*log(y)))
  })}

testing_entropy<-list()
for (i in names(fuzzy_data)){
testing_entropy[[i]]<-extract_entropy(fuzzy_data[[i]])
}

pdf(file = 'Plots/check_fuzzy_results.pdf',height = 8,width = 12)
par(mfrow=c(3,1))
for (i in names(fuzzy_data)){
  
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
  
plot(x=c(1:hor_lim[length(hor_lim)]),col='white',xlab = 'Depth',yaxt='n',ylab='',xlim=c(1,hor_lim[length(hor_lim)]),main = i)
for (z in 1:length(hor_lim)) {
  abline(v = hor_lim[z])
}
  
text(hor_lim-limits/2,50,labels = unique(no_out_details[[i]]$b.master_hor))
# #
# 
# test_levels<-unique(fuzzy_data[[i]]$Clustering$clustering)
# limits<-sapply(test_levels,function(x){
#   sum(fuzzy_data[[i]]$Clustering$clustering==x)*2
# })
# 
# hor_lim <- c(limits[1],sum(limits[1:2])) 
# if (length(unique(fuzzy_data[[i]]$Clustering$clustering))==3) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3])) 
# if (length(unique(fuzzy_data[[i]]$Clustering$clustering))==4) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]))
# if (length(unique(fuzzy_data[[i]]$Clustering$clustering))==5) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5])) 
# if (length(unique(fuzzy_data[[i]]$Clustering$clustering))==6) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5]),sum(limits[1:6])) 
# if (length(unique(fuzzy_data[[i]]$Clustering$clustering))==7) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5]),sum(limits[1:6]),sum(limits[1:7])) 
# 
# plot(x=c(1:hor_lim[length(hor_lim)]),col='white',xlab = 'Depth',yaxt='n',ylab='',xlim=c(1,hor_lim[length(hor_lim)]),main = i)
# for (z in 1:length(hor_lim)) {
#   abline(v = hor_lim[z])
# }
# text(hor_lim-limits/2,50,labels = unique(fuzzy_data[[i]]$Clustering$clustering))
# 
# #
plot(rep(fuzzy_data[[i]]$Clustering$membership[,1],each=2),type='l',main=paste0('Membership with phi = ',fuzzy_data[[i]]$Clustering$memb.exp),xlab='',xlim=c(1,hor_lim[length(hor_lim)]),xaxt = 'n')
lines(rep(fuzzy_data[[i]]$Clustering$membership[,2],each=2),type='l',col='blue')
if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 3){ 
lines(rep(fuzzy_data[[i]]$Clustering$membership[,3],each=2),type='l',col='red')
}
if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 4){
  lines(rep(fuzzy_data[[i]]$Clustering$membership[,4],each=2),type='l',col='violet')
}
if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 5){
  lines(rep(fuzzy_data[[i]]$Clustering$membership[,5],each=2),type='l',col='green')
}
if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 6){
  lines(rep(fuzzy_data[[i]]$Clustering$membership[,6],each=2),type='l',col='yellow')
}
if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 7){
  lines(rep(fuzzy_data[[i]]$Clustering$membership[,7],each=2),type='l',col='orange')
}
plot(rep(testing_entropy[[i]],each=2),main='Entropy',type='l',xlab='',xlim=c(1,hor_lim[length(hor_lim)]),xaxt = 'n')
abline(h=.5,col='red')

}
dev.off()

setwd('Plots/')
shell.exec(file = "check_fuzzy_results.pdf")
setwd('../')
}

#####PITS####
if (data_type=='pit'){
extract_CI<-function(fuzz_data) {
  aaply(fuzz_data$Clustering$membership,1,function(y){
    max_mem<-y[which.max(y)]
    sec_max<-y[order(y,decreasing = T)==2]
    CI <- sec_max/max_mem
  })
}
testing_CI<-list()
for (i in names(fuzzy_data)){
  testing_CI[[i]]<-extract_CI(fuzzy_data[[i]])
}


extract_entropy<-function(fuzz_data) {
  aaply(fuzz_data$Clustering$membership,1,function(y){
    entropy<-abs(1/log(fuzz_data$Clustering$k.crisp)*sum(y*log(y)))
  })}

testing_entropy<-list()
for (i in names(fuzzy_data)){
  testing_entropy[[i]]<-extract_entropy(fuzzy_data[[i]])
}

pdf(file = 'Plots/check_fuzzy_results.pdf',height = 8,width = 12)
par(mfrow=c(3,1))
for (i in names(fuzzy_data)){
  test_levels<-unique(no_out_details[[i]]$b.master_hor)
  limits<-sapply(test_levels,function(x){
    sum(no_out_details[[i]]$b.master_hor==x)*5
  })
  
  hor_lim <- c(limits[1],sum(limits[1:2])) 
  if (length(unique(no_out_details[[i]]$b.master_hor))==3) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3])) 
  if (length(unique(no_out_details[[i]]$b.master_hor))==4) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]))
  if (length(unique(no_out_details[[i]]$b.master_hor))==5) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5])) 
  if (length(unique(no_out_details[[i]]$b.master_hor))==6) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5]),sum(limits[1:6])) 
  if (length(unique(no_out_details[[i]]$b.master_hor))==7) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5]),sum(limits[1:6]),sum(limits[1:7])) 
  
  plot(x=c(1:hor_lim[length(hor_lim)]),col='white',xlab = 'Depth',yaxt='n',ylab='',xlim=c(1,hor_lim[length(hor_lim)]),main = i)
  for (z in 1:length(hor_lim)) {
    abline(v = hor_lim[z])
  }
  
  text(hor_lim-limits/2,25,labels = unique(no_out_details[[i]]$b.master_hor))
  #
  plot(rep(fuzzy_data[[i]]$Clustering$membership[,1],each=5),type='l',main=paste0('Membership with phi = ',fuzzy_data[[i]]$Clustering$memb.exp),xlab='',xaxt = 'n')
  lines(rep(fuzzy_data[[i]]$Clustering$membership[,2],each=5),type='l',col='blue')
  if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 3){ 
    lines(rep(fuzzy_data[[i]]$Clustering$membership[,3],each=5),type='l',col='red')
  }
  if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 4){
    lines(rep(fuzzy_data[[i]]$Clustering$membership[,4],each=5),type='l',col='violet')
  }
  if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 5){
    lines(rep(fuzzy_data[[i]]$Clustering$membership[,5],each=5),type='l',col='green')
  }
  if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 6){
    lines(rep(fuzzy_data[[i]]$Clustering$membership[,6],each=5),type='l',col='yellow')
  }
  if(ncol(fuzzy_data[[i]]$Clustering$membership) >= 7){
    lines(rep(fuzzy_data[[i]]$Clustering$membership[,7],each=5),type='l',col='orange')
  }
  plot(rep(testing_entropy[[i]],each=5),main='Entropy',type='l',xlab='',xlim=c(1,hor_lim[length(hor_lim)]),xaxt = 'n')
  
}
dev.off()

setwd('Plots/')
shell.exec(file = "check_fuzzy_results.pdf")
setwd('../')

}

#return the original functions
extract_FPI<-function(fuzz_data) {
  fpis<-sapply(fuzz_data,FUN = function(y){
    memberships<-y$membership
    part_coef<-sum(memberships^2)/nrow(memberships)
    FPI<-1-(((ncol(memberships)*part_coef)-1)/(ncol(memberships)-1))
  },USE.NAMES=T)
  return(fpis)
  names(fpis)
}
extract_CI<-function(fuzz_data) {
  aaply(fuzz_data$membership,1,function(y){
    max_mem<-y[which.max(y)]
    sec_max<-y[order(y,decreasing = T)==2]
    CI <- sec_max/max_mem
  })
}
extract_entropy<-function(fuzz_data) {
  aaply(fuzz_data$membership,1,function(y){
    entropy<-abs(1/log(which.min(tmp_fpi)+1)*sum(y*log(y))) #which min fpi is the number of cluster that has the lower fpi (starts from 2)
  })
}
