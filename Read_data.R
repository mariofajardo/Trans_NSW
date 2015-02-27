####This script will read the data from the asd instrument#####
#####################core_samples_3_rep_2cm_resolution###
library(pbapply)
library(reshape)
library(spectroscopy)
library(plyr)
setwd("..//Functions")#####    SETWD
load('correct_steps.RData')
load('check_plots.RData')
load('pr_varExp.RData')
prev_dir <- getwd()
setwd("../../Spectra/Transect/Core_samples/")#####    SETWD


hor <- sort_df(read.csv('NIR_hor.csv',header=T))     
list_df <- as.list.data.frame(split(hor,hor$id))
hor_names <- lapply(list_df,function(x){
  tmp <- data.frame(hor = rep(x$hor,times=x$bottom-x$top+1),master_hor = rep(x$master_horizon,times=x$bottom-x$top+1))
  tmp
})


###########read the spectra and assign horizon names###########
require(naturalsort)
pboptions(type='txt',style=3,title='Reading files')
files <- dir(pattern='txt',recursive=T)
files <- files[naturalorder(files)]
tmp_data <- pblapply(files,function(x){   
tmp <- read.csv(x,header=T)
data <- tmp[-1]
rep <- as.factor(rep(1:nrow(data), each=3, length.out=nrow(data)))
tmp2 <- as.data.frame(sapply(data, tapply, rep, mean))
tmp2$File.Name <- paste(regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)),paste0('Spectrum',sprintf('%05i',c(1:nrow(tmp2)))),sep='_')
sample <- regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T))
top <- seq(0,by=2,length.out=nrow(tmp2))
bottom <- seq(2,by=2,length.out=nrow(tmp2))
tmp3 <- cbind(File.Name=tmp2$File.Name,Sample=sample,top=top,bottom=bottom,tmp2[-length(tmp2)])
final <-tmp3[order(tmp3$File.Name,decreasing=F),]
})  


#####check the spectral differences####
pdf('../../../Codes/transect_n_s/Plots/check_spectra.pdf',width = 10,height = 6)
require(ggplot2)
require(reshape2)
require(spectroscopy)
require(signal)
load('../../../Codes/transect_n_s/RData/filter_sg2.RData')
sample<-1:length(files)
for (i in sample){
  spectra <-data.frame(tmp_data[[i]][,-c(1:4,2156)][rep(1:nrow(tmp_data[[i]][,-c(1:4,2156)]), each = 2), ])
  details <-c(seq(0,tmp_data[[i]]$top[length(tmp_data[[i]]$top)]),tmp_data[[i]]$top[length(tmp_data[[i]]$top)]+1)
  spectra<-spectra[naturalorder(details),]
  colnames(spectra)<-seq(350,2500)
  levels(spectra$sample) <- naturalsort(levels(spectra$sample))
  ####make some process####
  spectra <-1/log(spectra)
  spectra <-filter_sg(spectra,n=11,p=2,m=0)
  spectra <- trimSpec(spectra,wavlimits = c(500,2450),as.numeric(colnames(spectra)))
  spectra <-strip_spectra(spectra,as.numeric(colnames(spectra),c(500,2450),which=5))
  spectra <-snvBLC(spectra)
#   
  ###and some wavlength filtering ####
  spectra<-filter_sg2(spectra,n=5)
  spectra <- data.frame(sample=details,spectra)
  input_data<-melt(spectra,id.vars = 'sample')
  print(ggplot(input_data,aes(variable,sample, fill = value)) + 
    geom_tile() +
    labs(x='Wavelenght (nm)',
         y='Depth (cm)',
         title = paste0("Main absorbance differences in depth sample : ",tmp_data[[i]]$Sample[1])) +
    scale_x_discrete(breaks=c(500,1000,1500,2000,2450))+
      scale_y_reverse(breaks=c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_fill_gradient2(low = "green", high = "red",mid='white')+
    theme(axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          title=element_text(size=20)))
}
dev.off()
setwd('../../../Codes/transect_n_s/Plots/')
shell.exec('check_spectra.pdf')
setwd(prev_dir)



cores_NSW_DATA <- tmp_data


rawg_spectra<-lapply(cores_NSW_DATA,function(x) {
  spectra<-x[5:2155]
  spectra})

sample_details<-lapply(cores_NSW_DATA,function(x) {
  spectra<-x[1:4]
  spectra})

setwd(prev_dir)

View(rawg_spectra[[16]])
load('../../Codes/spectra_codes/RData/huntervalley_2cm_dry_cores.RDATA')
load('../../Codes/spectra_codes/RData/huntervalley_2cm_ground.RDATA')

plot(seq(350,2500),ground_DATA[[1]][1,7:2157],type='l')
lines(seq(350,2500),DATA[[1]][1,7:2157],col='red')
lines(seq(350,2500),rawg_spectra[[1]][1,],col='green')

pc_dried <- do.call(rbind,DATA)[,-c(1:6,2158)]
pc_ground <- do.call(rbind,ground_DATA[-37])[,7:2157] ###check sample 37
pc_sample <-do.call(rbind,rawg_spectra)

scores_dry <- prcomp(pc_dried)$x
scores_ground <- prcomp(pc_ground)$x
scores_sample <- prcomp(pc_sample)$x

data <-rbind(data.frame(dataset='ground',scores_ground[,1:10]),
             data.frame(dataset='sample',scores_sample[,1:10]),
             data.frame(dataset='dry',scores_dry[,1:10]))




find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hulls <- ddply(data,'dataset', find_hull)


  


require(ggplot2)
ggplot(data,aes(PC1,PC2,colour=dataset)) +
  geom_point() +
  geom_polygon(data = hulls, alpha = 0.2)


snv_data <-data.frame(snvBLC(rbind(pc_dried,pc_ground,pc_sample)))

scores_dry <- prcomp(snv_data[data$dataset=='dry',])$x
scores_ground <- prcomp(snv_data[data$dataset=='ground',])$x
scores_sample <- prcomp(snv_data[data$dataset=='sample',])$x

snv_scores <-rbind(data.frame(dataset='ground',scores_ground[,1:10]),
                          data.frame(dataset='sample',scores_sample[,1:10]),
                          data.frame(dataset='dry',scores_dry[,1:10]))


find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hulls <- ddply(snv_scores,'dataset', find_hull)

require(ggplot2)
ggplot(snv_scores,aes(PC1,PC2,colour=dataset)) +
  geom_point() +
  geom_polygon(data = hulls, alpha = 0.2)


#####plot3d####
# require(rgl)
# require(spectroscopy)
# input_data<-rawg_spectra
# open3d()
# for (i in 1:length(input_data)){
#   x<-i
#   y<-input_data
#   
#   sample3d<-strip_spectra(as.matrix(y[[x]][1:nrow(y[[x]]),]),datawavs=seq(350,2500),wavlimits=c(350,2500),which=2)
#   sample3d<-filter_sg(sample3d, n = 11, p = 2, m = 0)
#   sample3d_details<-sample_details[[x]][1:nrow(sample_details[[x]]),]
#   levels(sample3d_details$b.hor)<-terrain.colors(n=length(levels(sample3d_details$b.hor)),alpha=.75)
#   col <- matrix(sample3d_details$b.hor,nrow=nrow(sample3d),ncol=ncol(sample3d))
# 
#   
#   persp3d(x=1:nrow(sample3d)*2,y=as.numeric(colnames(sample3d)),z=sample3d,
#           aspect=c(2, 1, .5),
#           col=col,
#           xlab='Depth',
#           ylab='wavelength',
#           main=paste0('sample: ', regmatches(as.character(sample3d_details[1,1]),(regexpr(pattern="(?!=[_])[[:alnum:]]+",as.character(sample3d_details[1,1]),perl=T))))
#   )
#   
#   play3d(spin3d(axis=c(0,0,1), rpm=.001), duration=15)
# }
# 

#####fuzzy clustering####
cont_DATA<-pblapply(rawg_spectra,correct_step)    
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(500:2450))) 
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(500:2450),which=10))
snv_DATA <- pblapply(strip_DATA,snvBLC)


# check_plots(rawg_spectra,'Raw Spectra','Reflectance')
# check_plots(trim_DATA,'Trimmed Spectra','Reflectance')
# check_plots(abs_DATA,'Abs Spectra','Reflectance')
# check_plots(filt_DATA,'Filt Spectra','Reflectance')

pr_spectra<-prcomp(do.call(rbind,snv_DATA), center=T,scale=F) 
screeplot(pr_spectra)#visualize the PC
pr_varExp(do.call(rbind,snv_DATA))#check the acummulative variation on the data
pr_scores <- pr_spectra$x 

#####convex.hull analysis#####
require(tripack)
rand_tr <-tri.mesh(pr_scores[,1],pr_scores[,2])
rand.ch <- convex.hull(rand_tr,plot.it=F)
pr_poly <-cbind(x=c(rand.ch$y),y=c(rand.ch$y))

plot(pr_scores[,1],pr_scores[,2])
lines(c(rand.ch$x,rand.ch$x[1]),c(rand.ch$y,rand.ch$y[1]),col='blue')
#####outlier detection #####
mean_pcaA<- (colMeans(pr_scores[,1:5])) # mean of each of the components
mean_pcaA
cov_pcaA<- cov(pr_scores[,1:5])         # covariance matrix of the componets
cov_pcaA
chiMat<- matrix(NA,ncol=3,nrow=nrow(pr_scores)) 
chiMat[,1]<-mahalanobis(pr_scores[,1:5], mean_pcaA, cov_pcaA) # calculate mahalanobis distances

# Fit chi squared distribution and determine which rows can be removed based on the pcrit value (generally a pcrit of 0.975 is acceptable)
chiMat[,2]<- pchisq(c(chiMat[,1]), df =  5)       # chi squared distribution with 5 degrees of freedom (df = number of components)
pcrit<- 1-((0.24-0.003*5)/sqrt(nrow(pr_scores)))  # critical probability for outlier removal
pcrit
for (i in 1: nrow(chiMat)){
  if (chiMat[i,2] >= pcrit)
  {chiMat[i,3]<-0}else{chiMat[i,3]<-1}}           # which spectra are outliers
par(mfrow=c(1,2))
plot(chiMat[,1],chiMat[,2], xlab= "distance", ylab="cumlative prob")
points(chiMat[which(chiMat[,3]==0),1:2],pch='X',col='red')
abline(v=qchisq(pcrit,5),col="green")
plot(pr_scores[,1], pr_scores[,2], xlab="PCA 1", ylab="PCA 2", xlim=c(min(pr_scores[,1:2]), max(pr_scores[,1:2])),ylim=c(min(pr_scores[,1:2]), max(pr_scores[,1:2])))
points(pr_scores[which(chiMat[,3]==0),1:2],pch='X',col='red')

####Remove outliers from original dataset####
original_data <-do.call(rbind,snv_DATA)
original_details <- data.frame(do.call(rbind,sample_details))

new_data <- original_data[chiMat[,3] == 1,]
new_details<- original_details[chiMat[,3] == 1,]

####exclude the samples with too more than 10 outliers (20% of the the observations on a soil profile of 1m)####
count <- table(original_details[chiMat[,3] == 0,]$Sample)
exclude <- names(count)[count > 10]
no_out_details <- new_details[!(new_details$Sample %in% exclude),]
no_out_details$Sample <- droplevels(no_out_details$Sample)

######Dataset A filtered and transformed without outliers####
no_out_ground_DATA <- split.data.frame(cbind(no_out_details,new_data[!(new_details$Sample %in% exclude),]),no_out_details$Sample,drop=T)

no_out_data<-lapply(no_out_ground_DATA,function(x) {
  spectra<-x[5:ncol(x)]
  colnames(spectra)<-as.numeric(seq(from=500,to=2450,by=10))
  spectra})

new_pr_scores<-prcomp(do.call(rbind,no_out_data),center = T)$x

length(snv_DATA)              # original samples
length(no_out_ground_DATA)    # no outliers samples

#Data preparation
# no_out_data<-lapply(no_out_ground_DATA,function(x) {
#   spectra<-x[7:ncol(x)]
#   colnames(spectra)<-as.numeric(seq(from=500,to=2450,by=10))
#   spectra})
# 
# no_out_details<-lapply(no_out_ground_DATA,function(x) {
#   spectra<-x[1:6]
#   spectra})
# 
# no_out_pr_spectra<-prcomp(do.call(rbind,no_out_data), center=T,scale=F) 
# no_out_pr_Exp <- pr_varExp(do.call(rbind,no_out_data))
# no_out_pr_scores <- no_out_pr_spectra$x 
# 
# ####convex hull ploygon####
# rand.tr<-tri.mesh(no_out_pr_scores[,1],no_out_pr_scores[,2])
# rand.ch<-convex.hull(rand.tr, plot.it=F) 
# pr_poly = cbind(x=c(rand.ch$x),y=c(rand.ch$y))
# 
# 
# ####Visualize Dataset A without the outliers#### 
# 
# 
# sample_details_pc<-data.frame(do.call(rbind,no_out_details))$b.master_hor
# #####visualize the different horizons on pca space####
# 
# data <- data.frame(Soil_Types=sample_details_pc, no_out_pr_scores)
# ggplot(data, aes_string(x='PC1', y='PC2', col='Soil_Types')) + geom_point(alpha=.9,size=5) +
#   labs(x=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'), y=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'),title=paste0('Cumulative explanation :',round(no_out_pr_Exp[2,2],2),'%'))
# 
# ####create confidence ellipses####
# conf_ellipse <- data.frame()
# for(i in levels(data$Soil_Types)){
#   conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(data[data$Soil_Types==i,], 
#                                                                ellipse(cor(PC1,PC2), 
#                                                                        scale=c(sd(PC1),sd(PC2)), 
#                                                                        centre=c(mean(PC1),mean(PC2)),
#                                                                        level=.95))),Soil_Types=i))
# }
# #####plot confidence ellipses####
# ggplot(data, aes_string(x='PC1', y='PC2', colour='Soil_Types')) + geom_point(alpha=.6,size=6) +
#   labs(x=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'), 
#        y=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'),
#        title=paste0('Cumulative explanation :',
#                     round(no_out_pr_Exp[2,2],2),'%')) +
#   geom_path(data=conf_ellipse, aes_string(x='x', y='y',colour='Soil_Types'), size=1, linetype=1)
# 

#####ploting convex hulls####
# 
# find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
# hulls <- ddply(data,'Soil_Types', find_hull)
# 
# ggplot(data, aes_string(x='PC1', y='PC2', colour='Soil_Types')) + geom_point(alpha=.6,size=6) +
#   labs(x=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'), 
#        y=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'),
#        title=paste0('Cumulative explanation :',
#                     round(no_out_pr_Exp[2,2],2),'%')) +
#   geom_polygon(data = hulls, alpha = 0.2)
# 

# 
# ####and in 3d pca space###
# sample3d_details <- do.call(rbind,no_out_details)
# open3d(cex=.001)
# pca3d(no_out_pr_scores[,1:3],col=as.numeric(sample3d_details$b.master_hor),
#       radius=.15,
#       group=sample3d_details$b.master_hor,
#       show.plane=F,
#       show.group.labels=T
# )
# 
# #confidence spheres 70% conf#
# number_of_classes<-length(levels(sample3d_details$b.master_hor))
# palette(terrain.colors(number_of_classes))
# for (i in 1:number_of_classes){
#   c1<- mean(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1])
#   c2<- mean(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],2])
#   c3<- mean(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],3])
#   centremean <- c(c1,c2,c3)
#   plot3d(ellipse3d(cov(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.7),add=T,alpha=.12)
#   #segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col=palette()[i])
#   segments3d(c(0,c1),c(0,c2),c(0,c3),lwd=2,col='black')
#   #wire3d(ellipse3d(cov(no_out_pr_scores[sample3d_details$b.master_hor==levels(sample3d_details$b.master_hor)[i],1:3]),col=palette()[i],centre=centremean,level=.5,alpha=.7))
# }
# 
# 
# print(paste0('variance explained :', round(no_out_pr_Exp[3,2],2),'%'))

#####Bring in the fuzzy cluster of Dataset A and do some processing#####
setwd(prev_dir)
data_type<-'dry'
no_out_DATA<-no_out_data
source('fuzzy_single_sample.R',local = T)
source('check_clusters.R',local = T)


extract_FPI<-function(fuzz_data) {
  fpis<-sapply(fuzz_data,FUN = function(y){
    memberships<-y$membership
    part_coef<-sum(memberships^2)/nrow(memberships)
    FPI<-1-(((ncol(memberships)*part_coef)-1)/(ncol(memberships)-1))
  },USE.NAMES=T)
  return(fpis)
  names(fpis)
} #new functions (need to run before creating fuzzy data)
extract_CI<-function(fuzz_data) {
  aaply(fuzz_data$membership,1,function(y){
    max_mem<-y[which.max(y)]
    sec_max<-y[order(y,decreasing = T)==2]
    CI <- sec_max/max_mem
  })
}#new functions (need to run before creating fuzzy data)
extract_entropy<-function(fuzz_data) {
  aaply(fuzz_data$membership,1,function(y){
    entropy<-abs(1/log(which.min(tmp_fpi)+1)*sum(y*log(y))) #which min fpi is the number of cluster that has the lower fpi (starts from 2)
  })
}#new functions (need to run before creating fuzzy data)

fuzzy_data <-list()
phi <-c('phi1.5')
for (i in names(fanny_data_by_sample_pc_euc_no_out)){
  tmp_fpi <-extract_FPI(fanny_data_by_sample_pc_euc_no_out[[i]][[phi]])
  x<-names(tmp_fpi)[which.min(tmp_fpi)]
  fuzzy_data[[i]][['Clustering']]<- fanny_data_by_sample_pc_euc_no_out[[i]][[phi]][[x]]
  fuzzy_data[[i]][['Fpi']]<- tmp_fpi[which.min(tmp_fpi)]
  fuzzy_data[[i]][['Confusion_Index']]<- extract_CI(fanny_data_by_sample_pc_euc_no_out[[i]][[phi]][[x]])
  fuzzy_data[[i]][['Entropy']]<- extract_entropy(fanny_data_by_sample_pc_euc_no_out[[i]][[phi]][[x]])
}

source('check_fuzzy_results.R')

#####Digital Gradient####

pdf('Plots/Digital_gradient.pdf',height = 6,width = 10)
horizon_class<-list()
for (i in names(fuzzy_data)){
  sample<-i
  lag <-4 #in cm
  threshold<-.4
  source('Digital_Gradient.R') # inputs are the sample name and the desired lag (in cm)
}
dev.off()
setwd('Plots/')
shell.exec('Digital_gradient.pdf')
setwd(prev_dir)

