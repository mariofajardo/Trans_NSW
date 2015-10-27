CORES<-readRDS('RData/NSW_CORES_with_lab_and_predictionsThu_Oct_22_16_23_20_2015.rds')

###############################################
#..script to use NIR spectroscopy for soil....# 
#.......profiles horizons recognition.........#
#..SOME FUNCTIONS ARE MULTICORE (8cores)......#
#..AND USE LOTS OF RESOURCES AVOID USE OTHER .#
#..........PROGRAMS WHEN RUNNING..............#
###############################################

#####Read spectra files and outlier removal####
require(pbapply)
require(spectroscopy)
require(tripack)
require(plyr)
require(signal)
smoothing<-T
outlier <- T

#####Read data and outlier removal####
###############################################
#....script for reading txt spectra.......#####
#....and remove outliers..................#####
###############################################


load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/clustering_horizons/correct_steps.RData')
load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/clustering_horizons/pr_varExp.RData')
load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/clustering_horizons/check_plots.RData')
load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/filter_sg2.RData')

raw_DATA <- by(CORES$spectra,INDICES = CORES$responses$Sample,function(x) x)

cont_DATA<-pblapply(raw_DATA,correct_step)    
print('Trimming spectra')
trim_DATA<- pblapply(cont_DATA,function(x) trimSpec(x, wavlimits=range(500:2450))) 
print('Converting to absorbance')
abs_DATA<-pblapply(trim_DATA,function(x) absorbance<-log(1/x))  
print('S-G filter')
filt_DATA<- pblapply(abs_DATA,function(x) filter_sg(x, n = 11, p = 2, m = 0)) 
#####filter by sample#### 
if (smoothing){
  print('S-G filter by sample')
  # filt_DATA <- pblapply(filt_DATA,function(x) filter_sg2(x, n = 3, p = 2, m = 0))
  filt_DATA<- pblapply(filt_DATA,function(x) 
    t(aaply(x,2,function(y) c(runmed(y,k = 5)))))
}
# check_plots(filt_DATA,'Filt Spectra','Reflectance')
print('Stripping spectra')
strip_DATA<- pblapply(filt_DATA,function(x) strip_spectra(x,c(500:2450),which=10))
print('SNV transformation')
snv_DATA <- pblapply(strip_DATA,snvBLC)


#use this to check different types of spectral pre-processing effect#

# check_plots(raw_spectra,'Raw Spectra','Reflectance')
# check_plots(cont_DATA,'step corrected Spectra','Reflectance')
# check_plots(trim_DATA,'Trimmed Spectra','Reflectance')
# check_plots(abs_DATA,'Abs Spectra','Reflectance')
# check_plots(filt_DATA,'Filt Spectra','Reflectance')
# check_plots(snv_DATA,'Normalized Spectra','Reflectance')

if(outlier){
  pr_spectra<-prcomp(do.call(rbind,snv_DATA), center=T,scale=F) 
  screeplot(pr_spectra)#visualize the PC
  pr_varExp(do.call(rbind,snv_DATA))#check the acummulative variation on the data
  pr_scores <- pr_spectra$x 
  
  #####convex.hull analysis#####
  rand_tr <-tri.mesh(pr_scores[,1],pr_scores[,2],duplicate = 'remove')
  rand.ch <- convex.hull(rand_tr,plot.it=F)
  pr_poly <-cbind(x=c(rand.ch$y),y=c(rand.ch$y))
  
  plot(pr_scores[,1],
       pr_scores[,2],
       xlab="PCA 1", ylab="PCA 2",
       xlim=c(min(pr_scores[,1:2]),
              max(pr_scores[,1:2])),
       ylim=c(min(pr_scores[,1:2]),
              max(pr_scores[,1:2])))
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
  
  plot(pr_scores[,1],
       pr_scores[,2],
       xlab="PCA 1", ylab="PCA 2",
       xlim=c(min(pr_scores[,1:2]),
              max(pr_scores[,1:2])),
       ylim=c(min(pr_scores[,1:2]),
              max(pr_scores[,1:2])))
  
  points(pr_scores[which(chiMat[,3]==0),1:2],pch='X',col='red')
  
  ####Remove outliers from original dataset####
  original_data <-do.call(rbind,snv_DATA)
  original_details <- CORES$responses
  
  new_data <- original_data[chiMat[,3] == 1,]
  new_details<- original_details[chiMat[,3] == 1,]
  
  ####exclude the samples with too more than 10 outliers (20% of the the observations on a soil profile of 1m)####
  count <- table(original_details[chiMat[,3] == 0,]$Sample)
  exclude <- names(count)[count > 10]
  no_out_details <- new_details[!(new_details$Sample %in% exclude),]
  no_out_details$Sample <- droplevels(no_out_details$Sample)
  
  ######Dataset filtered and transformed without outliers####
  no_out_DATA <- split.data.frame(cbind(no_out_details,new_data[!(new_details$Sample %in% exclude),]),no_out_details$Sample,drop=T)
  no_out_data<-lapply(no_out_DATA,function(x) {
    spectra<-x[7:ncol(x)]
    colnames(spectra)<-as.numeric(seq(from=500,to=2450,by=10))
    spectra})
  
  new_pr_scores<-prcomp(do.call(rbind,no_out_data),center = T)$x
  
  length(snv_DATA)              # original samples
  length(no_out_DATA)           # outliers free dataset 
  
  rm(list=as.character(ls()[!ls()%in%c('no_out_ground_DATA','no_out_DATA','prev_dir')]))
  gc(reset=T)
} else {
  no_out_DATA <- snv_DATA
  names(no_out_DATA)<-sapply(sample_details,function(x) unique(x$Sample))
  rm(list=as.character(ls()[!ls()%in%c('no_out_ground_DATA','no_out_DATA','prev_dir')]))
}



#####Figure Main differences in depth####
# sample 0C

pdf(file = 'Plots/Main_differences_in_depth_1.pdf',width = 8,height = 8)
require(reshape2)
require(ggplot2)
sample<-'0C'
spectra <-data.frame(sample=seq(1,(2*nrow(no_out_DATA[[sample]][,-c(1:17)]))),no_out_DATA[[sample]][,-c(1:17)][rep(1:nrow(no_out_DATA[[sample]][,-c(1:17)]), each = 2), ])
names(spectra)<-c('sample',seq(500,2450,10))
input_data<-melt(spectra,measure.vars=c(names(no_out_DATA[[sample]][,-c(1:17)])),value.name='abs')
ggplot(input_data, aes(variable,-sample, fill = abs)) + 
  geom_tile() +
  labs(x='Wavelength (nm)',
       y='Depth (cm)',
       title = "Absorbance in depth") +
  scale_x_discrete(breaks=c(500,1000,1500,2000,2450)) +
  scale_fill_gradient2(low = "green", high = "red",mid='white')+
  theme(axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        legend.text=element_text(size=20),
        legend.title=element_text(size=30),
        title=element_text(size=30),
        axis.text=element_text(size=20))

dev.off()


#####
no_out_details <- lapply(no_out_DATA,function(x) x<-x[,c(1:17)])
new_pr_scores<-prcomp(do.call(rbind,no_out_DATA)[,-c(1:17)],center = T)$x



###create fuzzy classes for soil materials with whole dataset#### 
# it takes a while (about five to six hours with makeCluster(8)) to run this... 
# the previuosly calculated result can be loaded instead#   

# require(foreach)
# require(doSNOW)
# 
# scores<-prcomp(do.call(rbind,no_out_DATA)[,-c(1:17)])$x[,1:5]
# 
# 
# cl <-makeCluster(8)
# setMKLthreads(1) #if using Revolutions R
# registerDoSNOW(cl)
# 
# num_clusters <-2:20
# fanny_data_by_sample_pc_euc_no_out <- lapply(seq(1.3,1.5,0.05),function(phi){
#   test_data<-scores
#   foreach(nc=num_clusters,.packages='cluster') %dopar% fanny(x=test_data,k=nc,metric='euclidean',memb.exp=phi)
# })
# 
# stopCluster(cl)
# setMKLthreads(4)
# 
# JUST_IN_CASE <-fanny_data_by_sample_pc_euc_no_out
# names(fanny_data_by_sample_pc_euc_no_out)<-paste0('phi',seq(1.3,1.5,0.05))
# 
# 
# for(z in 1:length(seq(1.3,1.5,0.05))){
#   names(fanny_data_by_sample_pc_euc_no_out[[z]])<-paste0('cluster',num_clusters)
# }
# 
# 
# extract_FPI<-function(fuzz_data) {
#   fpis<-sapply(fuzz_data,FUN = function(y){
#     memberships<-y$membership
#     part_coef<-sum(memberships^2)/nrow(memberships)
#     FPI<- 1-((ncol(memberships)/(ncol(memberships)-1))*(1-part_coef))
#   },USE.NAMES=T)
#   return(fpis)
# }
# 
# tmp_fpi <-do.call(rbind,lapply(names(fanny_data_by_sample_pc_euc_no_out),function(x) extract_FPI(fanny_data_by_sample_pc_euc_no_out[[x]])))
# plot(2:20,tmp_fpi[5,],type='l',ylim = c(0,1),main='FPI values',col='white')
# 
# 
# for (z in 1:5) {
#   lines(2:20,tmp_fpi[z,],type='l')
#   text(x=2,y=tmp_fpi[z,1],label = names(fanny_data_by_sample_pc_euc_no_out)[z],cex=.7)
# }
# 
# 
# plot(rep(fanny_data_by_sample_pc_euc_no_out$phi1.5$cluster7$membership[48:107,1],each=2),type='l',ylim = c(0,1))
# 
# for (i in 2:13){
#   lines(rep(fanny_data_by_sample_pc_euc_no_out$phi1.4$cluster13$membership[48:107,i],each=2),type='l',col=i)
# }
# 
# fuzzy_2_20_whole <- fanny_data_by_sample_pc_euc_no_out
# # 
# save(fuzzy_2_20_whole,file='RData/fuzzy_2_20_whole.RData')
load('RData/fuzzy_2_20_whole.RData')
fanny_data_by_sample_pc_euc_no_out <-fuzzy_2_20_whole

#####Comparison of confuzion index#####

extract_FPI<-function(fuzz_data) {
  fpis<-sapply(fuzz_data,FUN = function(y){
    memberships<-y$membership
    part_coef<-sum(memberships^2)/nrow(memberships)
    FPI<- 1-((ncol(memberships)/(ncol(memberships)-1))*(1-part_coef))
  },USE.NAMES=T)
  return(fpis)
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
require(zoo)

pdf('Plots/Confusion_index_whole.pdf',height = 10,width = 16)
fuzzy_data <-list()
phi <-c('phi1.3')
tmp_fpi <-extract_FPI(fanny_data_by_sample_pc_euc_no_out[[phi]])
#   tmp_fpi<- as.zoo(tmp_fpi)
#   x<-names(rollapply(tmp_fpi, 3, function(x) which.min(x)==2)[rollapply(tmp_fpi, 3, function(x) which.min(x)==2)==TRUE][1])
x <- 'cluster13'
fuzzy_data[['Clustering']]<- fanny_data_by_sample_pc_euc_no_out[[phi]][[x]]
fuzzy_data[['Fpi']]<- tmp_fpi[names(tmp_fpi)==x]
fuzzy_data[['Confusion_Index']]<- extract_CI(fanny_data_by_sample_pc_euc_no_out[[phi]][[x]])
fuzzy_data[['Entropy']]<- extract_entropy(fanny_data_by_sample_pc_euc_no_out[[phi]][[x]])

par(mfrow = c(2,1))

hist(fuzzy_data$Confusion_Index,
     main=paste0('Confusion Index values for 13 classes with a membership exponent of : ','1.3'),
     xlab='Confusion index',
     pch=20,
     cex.lab=1.5,
     cex.main=2,
     ylim=c(0,1800)
)
rug(fuzzy_data$Confusion_Index)

fuzzy_data <-list()
phi <-c('phi1.4')
tmp_fpi <-extract_FPI(fanny_data_by_sample_pc_euc_no_out[[phi]])
#   tmp_fpi<- as.zoo(tmp_fpi)
#   x<-names(rollapply(tmp_fpi, 3, function(x) which.min(x)==2)[rollapply(tmp_fpi, 3, function(x) which.min(x)==2)==TRUE][1])
x <- 'cluster13'
fuzzy_data[['Clustering']]<- fanny_data_by_sample_pc_euc_no_out[[phi]][[x]]
fuzzy_data[['Fpi']]<- tmp_fpi[names(tmp_fpi)==x]
fuzzy_data[['Confusion_Index']]<- extract_CI(fanny_data_by_sample_pc_euc_no_out[[phi]][[x]])
fuzzy_data[['Entropy']]<- extract_entropy(fanny_data_by_sample_pc_euc_no_out[[phi]][[x]])

hist(fuzzy_data$Confusion_Index,
     main=paste0('Confusion Index values for 13 classes with a membership exponent of : ','1.4'),
     xlab='Confusion index',
     pch=20,
     cex.lab=1.5,
     cex.main=2,
     ylim=c(0,1800)
)
rug(fuzzy_data$Confusion_Index)

dev.off()
setwd('Plots/')
shell.exec(file = 'Confusion_index_whole.pdf')
setwd(prev_dir)


#####check clusters####
pdf('Plots/FPI_whole.pdf',height = 7,width = 16)

extract_FPI<-function(fuzz_data) {
  fpis<-sapply(fuzz_data,FUN = function(y){
    memberships<-y$membership
    part_coef<-sum(memberships^2)/nrow(memberships)
    FPI<- 1-((ncol(memberships)/((ncol(memberships)-1))*(1-part_coef)))
  },USE.NAMES=T)
  return(fpis)
}

tmp_fpi <-do.call(rbind,lapply(names(fanny_data_by_sample_pc_euc_no_out),function(x) extract_FPI(fanny_data_by_sample_pc_euc_no_out[[x]])))
plot(2:20,
     tmp_fpi[1,],
     type='l',
     ylim = c(0,1),
     main='Fuzziness performance index (FPI) values',
     col='white',
     ylab='FPI',
     xlab='Classes',
     pch=20,
     cex.lab=1.5,
     cex.main=2,
     cex.axis=1.5)

for (z in 1:5) {
  lines(2:20,tmp_fpi[z,],type='l')
  text(x=2,y=tmp_fpi[z,1],label = paste0('exp ',substring(names(fanny_data_by_sample_pc_euc_no_out)[z],first = 4)),cex=1.2)
}


abline(v=13,lty=2,h=tmp_fpi[2,12])

dev.off()
setwd('Plots/')
shell.exec(file = 'FPI_whole.pdf')

#####
require(plyr)
phi <- 'phi1.4'
classes <-'cluster13'
class_num<-13

extract_CI<-function(fuzz_data) {
  aaply(fuzz_data,1,function(y){
    max_mem<-y[which.max(y)]
    sec_max<-y[order(y,decreasing = T)==2]
    CI <- sec_max/max_mem
  })
}
extract_entropy<-function(fuzz_data) {
  aaply(fuzz_data,1,function(y){
    entropy<-abs(1/log(class_num)*sum(y*log(y))) #which min fpi is the number of cluster that has the lower fpi (starts from 2)
  })
}

fuzzy_data <- list()
data <- data.frame(Sample=do.call(rbind,no_out_DATA)$Sample)

for (i in names(no_out_DATA)){
  fuzzy_data[[i]][['Clustering']][['membership']]<- fanny_data_by_sample_pc_euc_no_out[[phi]][[classes]]$membership[data$Sample==i,]
  fuzzy_data[[i]][['Clustering']][['clustering']]<- fanny_data_by_sample_pc_euc_no_out[[phi]][[classes]]$clustering[data$Sample==i]
  fuzzy_data[[i]][['Confusion_Index']]<- extract_CI(fanny_data_by_sample_pc_euc_no_out[[phi]][[classes]]$membership[data$Sample==i,])
  fuzzy_data[[i]][['Entropy']]<- extract_entropy(fanny_data_by_sample_pc_euc_no_out[[phi]][[classes]]$membership[data$Sample==i,])
  fuzzy_data[[i]][['Clustering']][['k.crisp ']] <-as.numeric(fanny_data_by_sample_pc_euc_no_out[[phi]][[classes]]$k.crisp)
  fuzzy_data[[i]][['Clustering']][['memb.exp']] <- as.numeric(fanny_data_by_sample_pc_euc_no_out[[phi]][[classes]]$memb.exp)
  
}

data_type<-'dry'

no_out_details <- split(do.call(rbind,no_out_DATA)[,1:17],do.call(rbind,no_out_DATA)[,2])

check_classes <- data.frame(do.call(rbind,no_out_details),class=do.call(c,sapply(fuzzy_data,function(x) x$Clustering$clustering,USE.NAMES = F)),do.call(rbind,no_out_DATA)[,-c(1:17)],row.names = NULL)

#####Check centroids#####
my_col <-palette(terrain.colors(14))
pdf(file = 'Centroids.pdf',width = 16,height = 7)
require(naturalsort)
centroids <- ddply(check_classes,'class',function(x) {
  colMeans(x[,-c(1:18)])
},.drop=T)

centroids <- melt(centroids,id.vars = 'class')
centroids$variable<- as.numeric(substring(centroids$variable,first = 2))
centroids$class <- paste0('class',centroids$class) 


paper_theme <- theme_update(title = element_text(size=16,face='bold'), 
                            axis.text = element_text(colour='black'), 
                            axis.text.x = element_text(size=14), 
                            axis.title.x = element_text(vjust=-0.2, size=14, face='bold'), 
                            axis.text.y = element_text(size=14),
                            axis.title.y = element_text(size=14, face='bold', angle=90, vjust=0.5), legend.text = element_text(size=12,face='bold'), 
                            legend.key = element_blank() ,
                            legend.title = element_text(size=12), 
                            legend.background=element_blank(),
                            panel.grid.minor = element_blank(), 
                            panel.grid.major = element_blank(), 
                            panel.background = element_rect(fill = 'transparent',colour = 'black'))

centroids$class <-factor(centroids$class,levels = unique(centroids$class))  #correct order for ggplot
ggplot(centroids,aes(x = variable,y = value,group=class,colour=class))+
  geom_line(size=1.2)+
  scale_colour_manual(values = my_col)+
  xlab(label = 'Wavelength')+
  ylab(label = 'Absorbance')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.text=(element_text(size=20)),
        legend.title=element_text(size=20))
dev.off()
shell.exec(file = 'Centroids.pdf')


#####Check centroids relative distribution####
centroids <- ddply(check_classes,'class',function(x) {
  colMeans(x[,-c(1:18)])
},.drop=T)


dis_mat<-as.matrix(dist(prcomp(centroids[,-1])$x[,1:5]))
rownames(dis_mat) <- rownames(centroids)


a<-as.dist(dis_mat)
str(a)
# 
a<-hclust(a)
# 
par(mfrow=c(1,1))
# 
library(ape)
b<-as.phylo(a)
b$tip.label <- as.character(b$tip.label)

tip_col <-as.character(b$tip.label)
levels(tip_col) <-my_col

pdf(file = 'Centroids_dendro_numbers.pdf',width = 18,height = 7)
plot(b, 
     type = "unrooted", 
     cex = 0.1,
     label.offset = 5,
     show.tip.label=F,
     rotate.tree=90
)

tiplabels(text=rep(' ',length(tip_col)), 
          frame = "circle",
          col = as.character(tip_col),
          bg=as.character(tip_col),
          cex = 2)

legend(17,8.3,
       legend = b$tip.label,
       fill = unique(as.character(tip_col)),
       y.intersp=1,
       cex = 1.7)

dev.off()

shell.exec(file = 'Centroids_dendro_numbers.pdf')
setwd(prev_dir)


# #####Reorder centroids and their respective colors based on its distance####
# check_classes <- data.frame(do.call(rbind,no_out_details),class=do.call(c,sapply(fuzzy_data,function(x) x$Clustering$clustering,USE.NAMES = F)),do.call(rbind,no_out_DATA)[,-c(1:6)],row.names = NULL)
# check_classes$class <-factor(check_classes$class,levels = c(1,2,11,3,4,5,10,7,12,6,9,8))
# levels(check_classes$class)<-letters[1:12]
# 
# #and re order my colors#
# my_col <-palette(value = RColorBrewer::brewer.pal(12,name = 'Paired'))
# my_col <-my_col[c(1,2,11,3,4,5,10,7,12,6,9,8)]
# 
# 
# #####Plot again####
# ##### Dendrogram####
# centroids <- ddply(check_classes,'class',function(x) {
#   colMeans(x[,-c(1:7)])
# },.drop=T)
# 
# 
# dis_mat<-as.matrix(dist(prcomp(centroids[,-1])$x[,1:5]))
# rownames(dis_mat) <- centroids$class
# 
# 
# a<-as.dist(dis_mat)
# str(a)
# # 
# a<-hclust(a)
# # 
# par(mfrow=c(1,1))
# # 
# library(ape)
# b<-as.phylo(a)
# 
# pdf(file = '../Papers/Fuzzy_horizons/Reviews/Centroids_dendro.pdf',width = 19,height = 8)
# plot(b, 
#      type = "unrooted", 
#      cex = 0.1,
#      label.offset = 5,
#      show.tip.label=F,
#      rotate.tree=90
# )
# 
# tiplabels(text=rep(' ',length(tip_col)), 
#           frame = "circle",
#           col = my_col,
#           bg=my_col,
#           cex = 2)
# 
# legend(15,6.5,
#        legend = b$tip.label,
#        fill = my_col,
#        y.intersp=1,
#        cex = 2.5,
#        box.col='white')
# 
# dev.off()
# 
# setwd('../Papers/Fuzzy_horizons/Reviews/')
# shell.exec(file = 'Centroids_dendro.pdf')
# setwd(prev_dir)


#####Spectra centroids####
#####Centroids spectra####
# pdf(file = '../Papers/Fuzzy_horizons/Reviews/Centroids.pdf',width = 19,height = 8)
# 
# centroids <- ddply(check_classes,'class',function(x) {
#   colMeans(x[,-c(1:7)])
# },.drop=T)
# 
# centroids <- melt(centroids,id.vars = c('class'))
# centroids$variable<- as.numeric(substring(centroids$variable,first = 2))
# 
# 
# paper_theme <- theme_update(title = element_text(size=16,face='bold'), 
#                             axis.text = element_text(colour='black'), 
#                             axis.text.x = element_text(size=14), 
#                             axis.title.x = element_text(vjust=-0.2, size=14, face='bold'), 
#                             axis.text.y = element_text(size=14),
#                             axis.title.y = element_text(size=14, face='bold', angle=90, vjust=0.5), 
#                             legend.text = element_text(size=12), 
#                             legend.key = element_blank() ,
#                             legend.title = element_text(size=12), 
#                             legend.background=element_blank(),
#                             panel.grid.minor = element_blank(), 
#                             panel.grid.major = element_blank(), 
#                             panel.background = element_rect(fill = 'transparent',colour = 'black'))
# 
# centroids$class <-factor(centroids$class,levels = unique(centroids$class))  #correct order for ggplot
# ggplot(centroids,aes(x = variable,y = value,group=class,colour=class))+
#   geom_line(size=1.5)+
#   scale_colour_manual(values = my_col)+
#   xlab(label = 'Wavelength')+
#   ylab(label = 'Absorbance')+
#   theme(axis.text.x=element_text(size=20),
#         axis.text.y=element_text(size=20),
#         legend.text=(element_text(size=35)),
#         legend.title=element_text(size=40),
#         axis.title.x=element_text(size=30),
#         axis.title.y=element_text(size=30))
# dev.off()
# setwd('../Papers/Fuzzy_horizons/Reviews/')
# shell.exec(file = 'Centroids.pdf')
# setwd(prev_dir)

# #####Relative distribution by horizons reordered####
# 
# pdf(file = 'Centroids_distribution.pdf',width = 16,height = 10)
# ggplot(check_classes,aes(x = b.master_hor,group=class))+
#   geom_histogram(fill='gray')+
#   geom_density()+
#   theme(axis.text.x=element_text(size=20),
#         axis.text.y=element_text(size=20),
#         legend.text=(element_text(size=35)),
#         legend.title=element_text(size=40),
#         axis.title.x=element_text(size=30),
#         axis.title.y=element_text(size=30),
#         strip.text=element_text(size=22))+
#   xlab(label = 'Horizon type')+
#   ylab(label = 'Counts by class')+
#   ggtitle(label = 'Relative distribution of classes by horizons type')+
#   facet_wrap(~class,nrow=6)
# 
# dev.off()
# setwd('../Papers/Fuzzy_horizons/Reviews/')
# shell.exec(file = 'Centroids_distribution.pdf')
# setwd(prev_dir)

#####Figure for PhD_presentation_not included in paper#####

#####


# 
# 
# ####take branches ####
# 
# groups<-cutree(a,k=7)
# 
# plot(b, 
#      type = "unrooted", 
#      cex = 0.1,
#      label.offset = 5,
#      show.tip.label=F
# )
# tip_col <-as.factor(groups)
# 
# tiplabels(text=rep(' ',length(tip_col)), 
#           frame = "circle",
#           col = as.character(tip_col),
#           bg=as.character(tip_col),
#           cex = 0.3)
# 
# 
# 
# 
# legend(15,18,legend = unique(as.character(tip_col)),fill = unique(as.character(tip_col)))
# 
# View(data.frame(hor=no_out_DATA$ch12$b.hor,class=groups[1:47]))
# 
# scores<-prcomp(do.call(rbind,no_out_DATA)[,-c(1:6)])$x[,1:5]
# plot(scores[,1],scores[,2],col=groups)
# 
# 
# data <- data.frame(clus=groups,do.call(rbind,no_out_DATA))
# 
# 
# #theoretical centroid#
# theor_centroids <- ddply(.data = data,.variables = 'clus',.fun = function(x) colMeans(x[,-c(1:7)]))
# 
# plot(seq(500,2450,10),theor_centroids[7,-1],type='l')
# 
# for (i in 1:6){
#   lines(seq(500,2450,10),theor_centroids[i,-1],type='l',col=i)
#    }
# 
# data_for_plot<-rbind(data.frame(type='data',do.call(rbind,no_out_DATA)[,-c(1:6)]),
#                      data.frame(type='centroid',theor_centroids[,-1]))
# 
# scores<-prcomp(data_for_plot[,-1])$x[,1:5]
# 
# plot(scores[data_for_plot$type=='data',1],scores[data_for_plot$type=='data',2],col=groups)
# points(scores[data_for_plot$type=='centroid',1],scores[data_for_plot$type=='centroid',2],pch=10,lwd = 10)
# text(scores[data_for_plot$type=='centroid',1],scores[data_for_plot$type=='centroid',2],labels = rownames(theor_centroids),pch=10,lwd = 10,adj = c(2,2))
# 
# 

#####Check relative distribution of classes by horizons####


# source('check_fuzzy_results.R')

#####Horizons recognition with Digital gradients####
whole_dataset <-T
pdf(file = 'Digital_gradient.pdf',height = 10,width = 16)
horizon_class<-list()
for (i in names(fuzzy_data)){
  sample<-i
  lag <-4 #in cm
  threshold<-.2 #may be different for each profile
  source('../Digital_Gradient.R') # inputs are the sample name and the desired lag (in cm)
}
dev.off()

shell.exec('Digital_gradient.pdf')

#creation of an horizon classification#
str(horizon_class)

horizon_class <- sapply(horizon_class, function(x) x[seq(1,length(x),2)])

no_out_details_final <- no_out_details[sapply(horizon_class,function(x) unique(x!='homogeneus'),USE.NAMES = F)]
no_out_data_final <- no_out_DATA[sapply(horizon_class,function(x) unique(x!='homogeneus'),USE.NAMES = F)]

data <- data.frame(hor_class=do.call(c,horizon_class[horizon_class!='homogeneus']),do.call(rbind,no_out_data_final[names(horizon_class)]),row.names = NULL) 

thickness_of_first_spd_hor <- unlist(by(data,data$Sample,function(x) regmatches(x$hor_class[1],regexpr('[[:digit:]]{2}$',x$hor_class[1],perl=T))))

first_spd_hor_thickness <-data.frame(Sample=names(horizon_class),thickness_of_first_spd_hor=NA)
first_spd_hor_thickness[!names(horizon_class)%in%names(thickness_of_first_spd_hor),'thickness_of_first_spd_hor'] <- 100
first_spd_hor_thickness[first_spd_hor_thickness$Sample%in%names(thickness_of_first_spd_hor),'thickness_of_first_spd_hor'] <- thickness_of_first_spd_hor


number_of_spd_hor <- unlist(by(data,data$Sample,function(x) length(unique(x$hor_class)),simplify = F))

spd_hor_number <-data.frame(Sample=names(horizon_class),number_of_spd_hor=NA)
spd_hor_number[!names(horizon_class)%in%names(number_of_spd_hor),'number_of_spd_hor'] <- 1
spd_hor_number[names(horizon_class)%in%names(number_of_spd_hor),'number_of_spd_hor'] <- number_of_spd_hor

spd_hor_number

horizons <- data.frame(Sample=spd_hor_number[,1],num_spd_hor=spd_hor_number[,2],first_spd_hor_thickness=first_spd_hor_thickness[,2]) 


# date<-gsub(' |:','_',date())
# saveRDS(horizons,paste0('RData/','horizons',date,'.rds'))
