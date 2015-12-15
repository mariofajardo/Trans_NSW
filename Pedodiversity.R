# ####This script will read the data from the asd instrument#####
# #####################core_samples_3_rep_2cm_resolution###
# library(pbapply)
# library(reshape)
# library(spectroscopy)
# library(plyr)
# 
# load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/clustering_horizons/correct_steps.RData')
# load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/clustering_horizons/pr_varExp.RData')
# load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/clustering_horizons/check_plots.RData')
# load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/filter_sg2.RData')
# 
# 
# setwd("../../Spectra/Transect/Core_samples/")#####    SETWD
# 
# 
# 
# ###########read the spectra and assign horizon names###########
# require(naturalsort)
# pboptions(type='txt',style=3,title='Reading files')
# files <- dir(pattern='txt',recursive=T)
# files <- files[naturalorder(files)]
# tmp_data <- pblapply(files,function(x){   
#   tmp <- read.csv(x,header=T)
#   data <- tmp[-1]
#   rep <- rep(1:nrow(data), each=3, length.out=nrow(data))
#   File.Name <- paste(regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)),paste0('Spectrum',sprintf('%05i',c(1:nrow(data)))),sep='_')
#   sample <- rep(regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)),length.out=nrow(data))
#   top <- rep(seq(0,length.out = floor(length(sample)/3),by = 2),each=3)
#   bottom <- rep(seq(2,length.out = floor(length(sample)/3),by = 2),each=3)
#   tmp3 <- data.frame(File.Name=File.Name,Sample=sample,top=top,bottom=bottom,Rep=rep,data[-length(data)])
#   })  
# 
# tmp_data<-lapply(tmp_data,function(x){
#   x$system<-gsub('[0-9]','',x$Sample)
#   x}
# )
# 
# tmp_data<-lapply(tmp_data,function(x){
#   x$system[x$system=='C']<-'CROP'
#   x}
# )
# 
# tmp_data<-lapply(tmp_data,function(x){
#   x$system[x$system=='N']<-'NAT'
#   x}
# )
# 
# tmp_data<-lapply(tmp_data,function(x){
#   x$Sample<-regmatches(x$Sample,regexpr('[[:digit:]]+',x$Sample,perl = T))
#   x}
# )
# pedodiversity_dataset <-do.call(rbind,tmp_data)
# pedodiversity_dataset$index<-with(pedodiversity_dataset,paste(Sample,system,top,bottom,sep='_'))
# pedodiversity_dataset <- pedodiversity_dataset[pedodiversity_dataset$top<10,]
# saveRDS(pedodiversity_dataset,file = '../../../Codes/transect_n_s/RData/pedodiversity_dataset.RDS')
# #####Predict properties and pedodiversity by selected depth increments####
require(Cubist)
require(spectroscopy)
require(plyr)
require(foreach)
require(doSNOW)
require(clhs)

# # 
# #####CEC####
setwd("~/University of Sydney/PhD/Data/Lab_work/Codes/transect_n_s")
CORES<-readRDS('RData/NSW_CORES_with_lab_values_18_8_2105.RDS')
pedodiversity_dataset <- readRDS('RData/pedodiversity_dataset.RDS') 


#####CEC####
cl <-makeCluster(5)
registerDoSNOW(cl)
input <- CORES$calibration[,-c(1:20)]
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_CEC_surf_cores.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_CEC_surf_cores.RData')
source('Katoomba_Properties_pedodiversity.R')

goof(CORES$calibration$ECEC,predictions$Mean,main='CEC')


input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_CEC_surf_cores.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_CEC_surf_cores.RData')
source('Katoomba_Properties_pedodiversity.R')



pedodiversity_dataset$approx_top<-rep(c(rep(0,6),rep('out',3),rep(5,6)),length.out=length(pedodiversity_dataset$top))
pedodiversity_dataset$approx_bottom<-rep(c(rep(5,6),rep('out',3),rep(10,6)),length.out=length(pedodiversity_dataset$top))

pedodiversity_dataset$index<-with(pedodiversity_dataset,toupper(paste(Sample,system,approx_top,approx_bottom,sep='_')))

CEC_diversity<-by(predictions$Mean,pedodiversity_dataset$index,function(x) sd(x))
# CEC_diversity1<-by(predictions$Mean,pedodiversity_dataset$index,function(x) diversity(x,'simpson'))

pedodiversity <- data.frame(Sample=names(CEC_diversity),CEC_diversity=unlist(as.list(CEC_diversity),use.names = F))

Chemical<-readRDS('RData/NSW_CHEMICAL_with_lab_and_predictionTue_Nov_10_16_03_58_2015.rds')

Chemical$predictions$CEC_diversity <- pedodiversity[match(with(Chemical$responses,toupper(paste(site,system,top,bottom,sep='_'))),pedodiversity$Sample),'CEC_diversity']

# #####pH####

input <- CORES$calibration[,-c(1:20)]
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_pH_surf_cores.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_pH_surf_cores.RData')
source('Katoomba_Properties_pedodiversity.R')
goof(CORES$calibration$pH.Level..H2O.,predictions$Mean,main='pH')


input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_pH_surf_cores.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_pH_surf_cores.RData')
source('Katoomba_Properties_pedodiversity.R')


pedodiversity_dataset$approx_top<-rep(c(rep(0,6),rep('out',3),rep(5,6)),length.out=length(pedodiversity_dataset$top))
pedodiversity_dataset$approx_bottom<-rep(c(rep(5,6),rep('out',3),rep(10,6)),length.out=length(pedodiversity_dataset$top))

pedodiversity_dataset$index<-with(pedodiversity_dataset,toupper(paste(Sample,system,approx_top,approx_bottom,sep='_')))

pH_diversity<-by(predictions$Mean,pedodiversity_dataset$index,function(x) sd(x))

pedodiversity <- data.frame(Sample=names(pH_diversity),pH_diversity=unlist(as.list(pH_diversity),use.names = F))

Chemical$predictions$pH_diversity <- pedodiversity[match(with(Chemical$responses,toupper(paste(site,system,top,bottom,sep='_'))),pedodiversity$Sample),'pH_diversity']

#####Clay####
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

cl <-makeCluster(5)
registerDoSNOW(cl)
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Clay_surf_transect.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Clay_surf_transect.RData')
source('Katoomba_Properties_pedodiversity.R')
stopCluster(cl)


pedodiversity_dataset$approx_top<-rep(c(rep(0,6),rep('out',3),rep(5,6)),length.out=length(pedodiversity_dataset$top))
pedodiversity_dataset$approx_bottom<-rep(c(rep(5,6),rep('out',3),rep(10,6)),length.out=length(pedodiversity_dataset$top))

pedodiversity_dataset$index<-with(pedodiversity_dataset,toupper(paste(Sample,system,approx_top,approx_bottom,sep='_')))

Clay_diversity<-by(predictions$Mean,pedodiversity_dataset$index,function(x) sd(x))

pedodiversity <- data.frame(Sample=names(Clay_diversity),Clay_diversity=unlist(as.list(Clay_diversity),use.names = F))

Chemical$predictions$Clay_diversity <- pedodiversity[match(with(Chemical$responses,toupper(paste(site,system,top,bottom,sep='_'))),pedodiversity$Sample),'Clay_diversity']

#####EC####
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

cl <-makeCluster(5)
registerDoSNOW(cl)
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_EC_surf_transect.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_EC_surf_transect.RData')
source('Katoomba_Properties_pedodiversity.R')
stopCluster(cl)


pedodiversity_dataset$approx_top<-rep(c(rep(0,6),rep('out',3),rep(5,6)),length.out=length(pedodiversity_dataset$top))
pedodiversity_dataset$approx_bottom<-rep(c(rep(5,6),rep('out',3),rep(10,6)),length.out=length(pedodiversity_dataset$top))

pedodiversity_dataset$index<-with(pedodiversity_dataset,toupper(paste(Sample,system,approx_top,approx_bottom,sep='_')))

EC_diversity<-by(predictions$Mean,pedodiversity_dataset$index,function(x) sd(x))

pedodiversity <- data.frame(Sample=names(EC_diversity),EC_diversity=unlist(as.list(EC_diversity),use.names = F))

Chemical$predictions$EC_diversity <- pedodiversity[match(with(Chemical$responses,toupper(paste(site,system,top,bottom,sep='_'))),pedodiversity$Sample),'EC_diversity']

#####TC####
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

cl <-makeCluster(5)
registerDoSNOW(cl)
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_TC_surf_cores.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_TC_surf_cores.RData')
source('Katoomba_Properties_pedodiversity.R')
stopCluster(cl)


pedodiversity_dataset$approx_top<-rep(c(rep(0,6),rep('out',3),rep(5,6)),length.out=length(pedodiversity_dataset$top))
pedodiversity_dataset$approx_bottom<-rep(c(rep(5,6),rep('out',3),rep(10,6)),length.out=length(pedodiversity_dataset$top))

pedodiversity_dataset$index<-with(pedodiversity_dataset,toupper(paste(Sample,system,approx_top,approx_bottom,sep='_')))

TC_diversity<-by(predictions$Mean,pedodiversity_dataset$index,function(x) sd(x))

pedodiversity <- data.frame(Sample=names(TC_diversity),TC_diversity=unlist(as.list(TC_diversity),use.names = F))

Chemical$predictions$TC_diversity <- pedodiversity[match(with(Chemical$responses,toupper(paste(site,system,top,bottom,sep='_'))),pedodiversity$Sample),'TC_diversity']

#####Slaking_a####
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

cl <-makeCluster(5)
registerDoSNOW(cl)
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Slaking_a_coef.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Slaking_a_coef.RData')
source('Katoomba_Properties_pedodiversity.R')
stopCluster(cl)


pedodiversity_dataset$approx_top<-rep(c(rep(0,6),rep('out',3),rep(5,6)),length.out=length(pedodiversity_dataset$top))
pedodiversity_dataset$approx_bottom<-rep(c(rep(5,6),rep('out',3),rep(10,6)),length.out=length(pedodiversity_dataset$top))

pedodiversity_dataset$index<-with(pedodiversity_dataset,toupper(paste(Sample,system,approx_top,approx_bottom,sep='_')))

Slaking_a_diversity<-by(predictions$Mean,pedodiversity_dataset$index,function(x) sd(x))

pedodiversity <- data.frame(Sample=names(Slaking_a_diversity),Slaking_a_diversity=unlist(as.list(Slaking_a_diversity),use.names = F))

Chemical$predictions$Slaking_a_diversity <- pedodiversity[match(with(Chemical$responses,toupper(paste(site,system,top,bottom,sep='_'))),pedodiversity$Sample),'Slaking_a_diversity']


#####Update with spectral-hull diversity####
Chemical$predictions <- as.data.frame(Chemical$predictions)

pedodiversity_dataset <- readRDS('RData/pedodiversity_dataset_with_katoomba_preds.RDS') 


pedodiversity_dataset$approx_top<-rep(c(rep(0,6),rep('out',3),rep(5,6)),length.out=length(pedodiversity_dataset$top))
pedodiversity_dataset$approx_bottom<-rep(c(rep(5,6),rep('out',3),rep(10,6)),length.out=length(pedodiversity_dataset$top))
pedodiversity_dataset$index<-with(pedodiversity_dataset,toupper(paste(Sample,system,approx_top,approx_bottom,sep='_')))

require(tripack)
require(naturalsort)

#####In_hull diversity function####
inhull_diversity_spectra <- function(variable_matrix=NULL,index=NULL){
  require(splancs)
  cha<-function(x,y){
    chull(x,y)->i
    return(areapl(cbind(x[i],y[i])))
  }
  components<-prcomp(variable_matrix,center = T)
  scores <- data.frame(PC1=components$x[,'PC1'],PC2=components$x[,'PC2'])
  test<-split(scores,index)
  hull_areas<-data.frame(index=names(test),In_hull_diversity_spectra=sapply(test,function(x) cha(x[,1],x[,2])))
}

inhull_diversity_prop <- function(variable_matrix=NULL,index=NULL){
  require(splancs)
  cha<-function(x,y){
    chull(x,y)->i
    return(areapl(cbind(x[i],y[i])))
  }
  components<-prcomp(variable_matrix,scale. = T,center = T)
  scores <- data.frame(PC1=components$x[,'PC1'],PC2=components$x[,'PC2'])
  test<-split(scores,index)
  hull_areas<-data.frame(index=names(test),In_hull_diversity_prop=sapply(test,function(x) cha(x[,1],x[,2])))
}


diversity_spectra<- inhull_diversity_spectra(pedodiversity_dataset[!grepl('OUT',pedodiversity_dataset$index),grepl('X',colnames(pedodiversity_dataset))]
                                    ,pedodiversity_dataset$index[!grepl('OUT',pedodiversity_dataset$index)])

diversity_prop<- inhull_diversity_prop(pedodiversity_dataset[!grepl('OUT',pedodiversity_dataset$index),grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))]
                                    ,pedodiversity_dataset$index[!grepl('OUT',pedodiversity_dataset$index)])

View(diversity_spectra[match(Chemical$responses$index,paste(diversity_spectra$index,'R1',sep='_')),])

Chemical$predictions$diversity_spectra <- diversity_spectra[match(Chemical$responses$index,paste(diversity_spectra$index,'R1',sep='_')),'In_hull_diversity_spectra']
Chemical$predictions$diversity_prop <- diversity_prop[match(Chemical$responses$index,paste(diversity_prop$index,'R1',sep='_')),'In_hull_diversity_prop']

#####Normalized centroid distance#####


Normalized_centroid_diversity_prop <- function(variable_matrix=pedodiversity_dataset[!grepl('OUT',pedodiversity_dataset$index),grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))]
                                          ,index=pedodiversity_dataset$index[!grepl('OUT',pedodiversity_dataset$index)]){
  
  scores<-as.data.frame(prcomp(variable_matrix,scale. = T,center = T)$x[,1:6])
  test<-split(scores,index)
  Normalized_centroid_diversity<-data.frame(index=names(test),Normalized_centroid_diversity_prop=sapply(test,function(x){
  sqrt(sum(apply(x,2,function(y) (y-mean(y))^2))/nrow(x))
  }))
}


Normalized_centroid_diversity_spectra <- function(variable_matrix=pedodiversity_dataset[!grepl('OUT',pedodiversity_dataset$index),grepl('X',colnames(pedodiversity_dataset))]
                                               ,index=pedodiversity_dataset$index[!grepl('OUT',pedodiversity_dataset$index)]){
  
  scores<-as.data.frame(prcomp(variable_matrix,center = T)$x[,1:10])
  test<-split(scores,index)
  Normalized_centroid_diversity<-data.frame(index=names(test),Normalized_centroid_diversity_spectra=sapply(test,function(x){
    sqrt(sum(apply(x,2,function(y) (y-mean(y))^2))/nrow(x))
  }))
}

Norm_centroid_diversity_spectra<- Normalized_centroid_diversity_spectra(pedodiversity_dataset[!grepl('OUT',pedodiversity_dataset$index),grepl('X',colnames(pedodiversity_dataset))]
                                    ,pedodiversity_dataset$index[!grepl('OUT',pedodiversity_dataset$index)])

require(car)
# pedodiversity_dataset$predicted_Clay <- logit(pedodiversity_dataset$predicted_Clay)
# pedodiversity_dataset$predicted_Clay <- null

Norm_centroid_diversity_prop<- Normalized_centroid_diversity_prop(pedodiversity_dataset[!grepl('OUT',pedodiversity_dataset$index),grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))]
                                                              ,pedodiversity_dataset$index[!grepl('OUT',pedodiversity_dataset$index)])


Chemical$predictions$Norm_centroid_diversity_spectra <- Norm_centroid_diversity_spectra[match(Chemical$responses$index,paste(Norm_centroid_diversity_spectra$index,'R1',sep='_')),'Normalized_centroid_diversity_spectra']
Chemical$predictions$Norm_centroid_diversity_prop <- Norm_centroid_diversity_prop[match(Chemical$responses$index,paste(Norm_centroid_diversity_prop$index,'R1',sep='_')),'Normalized_centroid_diversity_prop']


######Check how in -hull distances work

#with properties
# pedodiversity_dataset <- pedodiversity_dataset[!grepl('out',pedodiversity_dataset$approx_top),]
# components <- prcomp(pedodiversity_dataset[,grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))
# ],scale. = T,center = T)

#with spectra
pedodiversity_dataset <- pedodiversity_dataset[!grepl('out',pedodiversity_dataset$approx_top),]
components <- prcomp(pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))],center = T)


sample<-c(5,6,7)
system <- c('CROP','NAT')
depth <- 0  
scores <- data.frame(components$x)

logical_query <-pedodiversity_dataset$Sample%in%sample&
  pedodiversity_dataset$system%in%system&
  pedodiversity_dataset$approx_top%in%depth
scores_dataframe <- data.frame(pedodiversity_dataset[,c('index','top','system','Sample')],scores)
scores_plot <- data.frame(pedodiversity_dataset[logical_query,c('index','top','system','Sample')],scores[logical_query,])

scores_plot$index <- factor(scores_plot$index)
scores_plot$Sample <- reorder(scores_plot$Sample,naturalorder(scores_plot$Sample))
# ggplot(scores,aes(PC1,PC2,colour=index.index))+
#   geom_point()

require(plyr)
require(ggplot2)
find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hulls <- ddply(scores_plot,'index', find_hull)

whole <- ggplot(scores_plot, aes_string(x='PC1', y='PC2',colour='Sample',shape='system',group='index')) + 
  geom_point(alpha=.6,size=6) +
  geom_polygon(data = hulls, alpha = 0.2)+
  ggtitle('Convex hulls in spectra principal component space')+
  coord_equal()

whole
####Pick the three closest points####
selection <- as.list(by(scores[,1:2],scores_dataframe$index,function(x) dist(x)))

solutions<-unlist(lapply(selection,function(x) {
dist_mat<-as.matrix(x) #distance matrix

condition_1<-colSums(dist_mat==dist_mat[dist_mat!=0][which.min(dist_mat[dist_mat!=0])])

two_closest_points<-colnames(dist_mat)[as.logical(condition_1)]
third_closest_point_distance <- sort(dist_mat[,two_closest_points])[5]

condition_2 <- rowSums(dist_mat[,two_closest_points]==third_closest_point_distance)

third_closest_point <- rownames(dist_mat[,two_closest_points])[as.logical(condition_2)]

solution <- c(two_closest_points,third_closest_point)
}),use.names = F)


new_scores <- scores[solutions,]


new_scores_dataframe<- data.frame(pedodiversity_dataset[solutions,c('index','top','system','Sample')],scores[solutions,])


logical_query1 <-pedodiversity_dataset[solutions,'Sample']%in%sample&
  pedodiversity_dataset[solutions,'system']%in%system&
  pedodiversity_dataset[solutions,'approx_top']%in%depth

new_scores_plot <- data.frame(pedodiversity_dataset[solutions,][logical_query1,c('index','top','system','Sample')],new_scores[logical_query1,])

new_scores_plot$index <- factor(new_scores_plot$index)
new_scores_plot$Sample <- reorder(new_scores_plot$Sample,naturalorder(new_scores_plot$Sample))


require(plyr)
require(ggplot2)
find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hulls <- ddply(new_scores_plot,'index', find_hull)

whole + 
  geom_point(data=new_scores_plot,alpha=.6,size=6,colour='black') +
  geom_polygon(data = hulls, alpha = 0.2)


#####

diversity_spectra_closest<- inhull_diversity_spectra(pedodiversity_dataset[solutions,grepl('X',colnames(pedodiversity_dataset))],pedodiversity_dataset[solutions,'index'])
diversity_prop_closest<- inhull_diversity_prop(pedodiversity_dataset[solutions,grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))]
                                                                  ,pedodiversity_dataset[solutions,'index'])

diversity_centroid_spectra_closest<- Normalized_centroid_diversity_spectra(pedodiversity_dataset[solutions,grepl('X',colnames(pedodiversity_dataset))],pedodiversity_dataset[solutions,'index'])
diversity_centroid_prop_closest<- Normalized_centroid_diversity_prop(pedodiversity_dataset[solutions,grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))]
                                               ,pedodiversity_dataset[solutions,'index'])



Chemical$predictions$diversity_spectra_closest <- diversity_spectra_closest[match(Chemical$responses$index,paste(diversity_spectra_closest$index,'R1',sep='_')),'In_hull_diversity_spectra']
Chemical$predictions$diversity_prop_closest <- diversity_prop_closest[match(Chemical$responses$index,paste(diversity_prop_closest$index,'R1',sep='_')),'In_hull_diversity_prop']

Chemical$predictions$diversity_centroid_spectra_closest <- diversity_centroid_spectra_closest[match(Chemical$responses$index,paste(diversity_centroid_spectra_closest$index,'R1',sep='_')),'Normalized_centroid_diversity_spectra']
Chemical$predictions$diversity_centroid_prop_closest <- diversity_centroid_prop_closest[match(Chemical$responses$index,paste(diversity_centroid_prop_closest$index,'R1',sep='_')),'Normalized_centroid_diversity_prop']


# 
# ####Remove outliers first (by sample)####
# sample<-0:48
# system <- c('CROP','NAT')
# depth <- c(0,5)  
# scores <- data.frame(components$x)
# 
# logical_query <-pedodiversity_dataset$Sample%in%sample&
#   pedodiversity_dataset$system%in%system&
#   pedodiversity_dataset$approx_top%in%depth
# scores <- data.frame(pedodiversity_dataset[logical_query,c('index','top','system','Sample')],scores[logical_query,])
# 
# scores$index <- factor(scores$index)
# scores$Sample <- reorder(scores$Sample,naturalorder(scores$Sample))
# 
# new_scores_no_outliers <- do.call(rbind,as.list(by(scores,scores$index,function(x) {
# pr_scores <- x[,5:9] #scores[,5:9][scores$index=='1_CROP_0_5',]
# 
# #####convex.hull analysis#####
# rand_tr <-tri.mesh(pr_scores[,1],pr_scores[,2],duplicate = 'remove')
# rand.ch <- convex.hull(rand_tr,plot.it=F)
# pr_poly <-cbind(x=c(rand.ch$y),y=c(rand.ch$y))
# 
# plot(pr_scores[,1],
#      pr_scores[,2],
#      xlab="PCA 1", ylab="PCA 2",
#      xlim=c(min(pr_scores[,1:2]),
#             max(pr_scores[,1:2])),
#      ylim=c(min(pr_scores[,1:2]),
#             max(pr_scores[,1:2])))
# lines(c(rand.ch$x,rand.ch$x[1]),c(rand.ch$y,rand.ch$y[1]),col='blue')
# #####outlier detection #####
# mean_pcaA<- (colMeans(pr_scores[,1:3])) # mean of each of the components
# mean_pcaA
# cov_pcaA<- cov(pr_scores[,1:3])         # covariance matrix of the componets
# cov_pcaA
# chiMat<- matrix(NA,ncol=3,nrow=nrow(pr_scores)) 
# chiMat[,1]<-mahalanobis(pr_scores[,1:3], mean_pcaA, cov_pcaA) # calculate mahalanobis distances
# 
# # Fit chi squared distribution and determine which rows can be removed based on the pcrit value (generally a pcrit of 0.975 is acceptable)
# chiMat[,2]<- pchisq(c(chiMat[,1]), df =  5)       # chi squared distribution with 5 degrees of freedom (df = number of components)
# pcrit<- 1-((0.24-0.003*5)/sqrt(nrow(pr_scores)))  # critical probability for outlier removal
# pcrit
# for (i in 1: nrow(chiMat)){
#   if (chiMat[i,2] >= pcrit)
#   {chiMat[i,3]<-0}else{chiMat[i,3]<-1}}           # which spectra are outliers
# par(mfrow=c(1,2))
# plot(chiMat[,1],chiMat[,2], xlab= "distance", ylab="cumlative prob")
# points(chiMat[which(chiMat[,3]==0),1:2],pch='X',col='red')
# abline(v=qchisq(pcrit,5),col="green")
# 
# plot(pr_scores[,1],
#      pr_scores[,2],
#      xlab="PCA 1", ylab="PCA 2",
#      xlim=c(min(pr_scores[,1:2]),
#             max(pr_scores[,1:2])),
#      ylim=c(min(pr_scores[,1:2]),
#             max(pr_scores[,1:2])))
# 
# points(pr_scores[which(chiMat[,3]==0),1:2],pch='X',col='red')
# 
# ####Remove outliers from original dataset####
# 
# 
# solution <- x[chiMat[,3] == 1,]
# })))
# 
# 
# require(plyr)
# require(ggplot2)
# find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
# hulls <- ddply(new_scores_no_outliers,'index', find_hull)
# 
# whole + 
#   geom_point(data=new_scores_no_outliers,alpha=.6,size=6,colour='black') +
#   geom_polygon(data = hulls, alpha = 0.2)
# 
# 

####And by dataset (better)####
####Remove outliers first (by sample)####
sample<-0:48
system <- c('CROP','NAT')
depth <- c(0,5)
scores <- data.frame(components$x)

logical_query <-pedodiversity_dataset$Sample%in%sample&
  pedodiversity_dataset$system%in%system&
  pedodiversity_dataset$approx_top%in%depth
scores <- data.frame(pedodiversity_dataset[logical_query,c('index','top','system','Sample')],scores[logical_query,])

scores$index <- factor(scores$index)
scores$Sample <- reorder(scores$Sample,naturalorder(scores$Sample))

remove_outliers <- function(x) {
  pr_scores <- x[,5:9] #scores[,5:9][scores$index=='1_CROP_0_5',]
  
  #####convex.hull analysis#####
  rand_tr <-tri.mesh(pr_scores[,1],pr_scores[,2],duplicate = 'remove')
  rand.ch <- convex.hull(rand_tr,plot.it=F)
  pr_poly <-cbind(x=c(rand.ch$y),y=c(rand.ch$y))
  
#   plot(pr_scores[,1],
#        pr_scores[,2],
#        xlab="PCA 1", ylab="PCA 2",
#        xlim=c(min(pr_scores[,1:2]),
#               max(pr_scores[,1:2])),
#        ylim=c(min(pr_scores[,1:2]),
#               max(pr_scores[,1:2])))
#   lines(c(rand.ch$x,rand.ch$x[1]),c(rand.ch$y,rand.ch$y[1]),col='blue')
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
#   par(mfrow=c(1,2))
#   plot(chiMat[,1],chiMat[,2], xlab= "distance", ylab="cumlative prob")
#   points(chiMat[which(chiMat[,3]==0),1:2],pch='X',col='red')
#   abline(v=qchisq(pcrit,5),col="green")
  
#   plot(pr_scores[,1],
#        pr_scores[,2],
#        xlab="PCA 1", ylab="PCA 2",
#        xlim=c(min(pr_scores[,1:2]),
#               max(pr_scores[,1:2])),
#        ylim=c(min(pr_scores[,1:2]),
#               max(pr_scores[,1:2])))
#   
#   points(pr_scores[which(chiMat[,3]==0),1:2],pch='X',col='red')
  
  ####Remove outliers from original dataset####
  
  
  solution <- x[chiMat[,3] == 1,]
}

no_out_scores <- remove_outliers(scores)


require(plyr)
require(ggplot2)
find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hulls <- ddply(no_out_scores,'index', find_hull)


ggplot(no_out_scores, aes_string(x='PC1', y='PC2',colour='Sample',shape='system',group='index')) + 
  geom_point(alpha=.6,size=6) +
  geom_polygon(data = hulls, alpha = 0.2)+
  ggtitle('Convex hulls in spectra principal component space')+
  coord_equal()

sample<-seq(0,48,3)
system <- c('CROP','NAT')
depth <- c(0,5)

plot_scores <- data.frame(no_out_scores)

logical_query <-plot_scores$Sample%in%sample&
  plot_scores$system%in%system&
  plot_scores$index%in%no_out_scores$index
plot_scores <- data.frame(plot_scores[logical_query,])

plot_scores$index <- factor(plot_scores$index)
plot_scores$Sample<-as.numeric(plot_scores$Sample)
plot_scores$Sample <- reorder(plot_scores$Sample,naturalorder(plot_scores$Sample))


find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hulls <- ddply(plot_scores,'index', find_hull)

ggplot(plot_scores, aes_string(x='PC1', y='PC2',colour='Sample',shape='system',group='index')) + 
  geom_point(alpha=.6,size=6) +
  geom_polygon(data = hulls, alpha = 0.2)


diversity_spectra<- inhull_diversity_spectra(pedodiversity_dataset[!grepl('OUT',pedodiversity_dataset$index),grepl('X',colnames(pedodiversity_dataset))]
                                             ,pedodiversity_dataset$index[!grepl('OUT',pedodiversity_dataset$index)])


condition <- rownames(pedodiversity_dataset)%in%rownames(no_out_scores)&!grepl('OUT',pedodiversity_dataset$index)

no_out_dataset_hull_spectra <- inhull_diversity_spectra(pedodiversity_dataset[condition,grepl('X',colnames(pedodiversity_dataset))],pedodiversity_dataset[condition,'index'])
no_out_dataset_hull_prop <- inhull_diversity_prop(pedodiversity_dataset[condition,grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))],pedodiversity_dataset[condition,'index'])


no_out_dataset_centroid_spectra <- Normalized_centroid_diversity_spectra(pedodiversity_dataset[condition,grepl('X',colnames(pedodiversity_dataset))],pedodiversity_dataset[condition,'index'])
no_out_dataset_centroid_prop <- Normalized_centroid_diversity_prop(pedodiversity_dataset[condition,grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))],pedodiversity_dataset[condition,'index'])

 

#####add the data ####

Chemical$predictions$diversity_inhull_spectra_dataset_without_outliers <- no_out_dataset_hull_spectra[match(Chemical$responses$index,paste(no_out_dataset_hull_spectra$index,'R1',sep='_')),'In_hull_diversity_spectra']
Chemical$predictions$diversity_inhull_prop_dataset_without_outliers <- no_out_dataset_hull_prop[match(Chemical$responses$index,paste(no_out_dataset_hull_prop$index,'R1',sep='_')),'In_hull_diversity_prop']


Chemical$predictions$diversity_centroid_spectra_dataset_without_outliers <- no_out_dataset_centroid_spectra[match(Chemical$responses$index,paste(no_out_dataset_centroid_spectra$index,'R1',sep='_')),'Normalized_centroid_diversity_spectra']
Chemical$predictions$diversity_centroid_prop_dataset_without_outliers <- no_out_dataset_centroid_prop[match(Chemical$responses$index,paste(no_out_dataset_centroid_prop$index,'R1',sep='_')),'Normalized_centroid_diversity_prop']


colnames(Chemical$predictions)
pedodiversity<-data.frame(Site=as.numeric(Chemical$responses$site),
                            system=Chemical$responses$system,
                            top=Chemical$responses$top,
                            index=Chemical$responses$index,
                            Chemical$predictions[,13:24])

pedodiversity<-pedodiversity[pedodiversity$Site%in%0:26&pedodiversity$top%in%0,]
require(reshape2)
pedodiversity <- melt(pedodiversity,id.vars = c('Site','system','top','index'))

ggplot(pedodiversity,aes(as.numeric(Site),value,group=variable,colour=variable))+
  geom_line(alpha=.5)+
  geom_smooth(size=2)+
  ggtitle('Diversity comparison')+
  facet_wrap(~variable+system)


pedodiversity<-data.frame(Site=as.numeric(Chemical$responses$site),
                          system=Chemical$responses$system,
                          top=Chemical$responses$top,
                          index=Chemical$responses$index,
                          Chemical$predictions[,13:24])

pedodiversity<-pedodiversity[pedodiversity$Site%in%0:26&pedodiversity$top%in%0,]

pdf('Plots/pedodiversities.pdf',width = 20,height = 15)

for (i in colnames(pedodiversity)[-c(1:4)]){
print(ggplot(pedodiversity,aes_string('Site',i,colour='system'))+
  geom_line(alpha=.5)+
  geom_smooth(size=2)+
  ggtitle('Diversity comparison')+
  ylim(-1,5))
    
}
dev.off()

setwd('Plots/')

shell.exec('pedodiversities.pdf')
# date<-gsub(' |:','_',date())
# saveRDS(Chemical,paste0('RData/','NSW_CHEMICAL_with_predictions_surface_and_pedodiversity',date,'.rds'))
# # 

####Explore the newly derived pedodiversity###
##at different scales###
depth<-5
ncomp<-4
 



Scalewise_Normalized_centroid_diversity_pro<-function(variable_matrix=NULL){
  variable_matrix<-variable_matrix[,grepl('predicted',colnames(variable_matrix))&!grepl('sd',colnames(variable_matrix))]
  scores<-as.data.frame(prcomp(variable_matrix,scale. = T,center = T)$x[,1:ncomp])
  sqrt(sum(apply(scores,2,function(y) (y-mean(y))^2))/nrow(scores))
  }


#First all of them#

lag<-seq(1,25,by = 1)

pedodivs<-list()
for(x in lag){
distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
sites<-sites[sites$Var2%in%0:26&sites$Var1%in%0:26,]

conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$Sample%in%as.numeric(i))

mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
  Scalewise_Normalized_centroid_diversity_pro(pedodiversity_dataset[x,])
}))
pedodivs[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs)

####Now for crop####
lag<-seq(1,25,by = 1)

pedodivs_crop<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var2%in%0:26&sites$Var1%in%0:26,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='CROP'&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_Normalized_centroid_diversity_pro(pedodiversity_dataset[x,])
  }))
  pedodivs_crop[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs_crop)


####Now for nat####
lag<-seq(1,25,by = 1)

pedodivs_nat<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var2%in%0:26&sites$Var1%in%0:26,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='NAT'&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_Normalized_centroid_diversity_pro(pedodiversity_dataset[x,])
  }))
  pedodivs_nat[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs_nat)



# pedodivs_NS<-data.frame(all=unlist(pedodivs),
#                         nat=unlist(pedodivs_nat),
#                         crop=unlist(pedodivs_crop))
pedodivs_NS<-data.frame(nat=unlist(pedodivs_nat),
                        crop=unlist(pedodivs_crop))

pedodivs_NS$lag<-1:25*50


require(reshape2)
pedodivs_NS<-melt(pedodivs_NS,variable_name = 'Type',id.vars = 'lag')

NS<-ggplot(pedodivs_NS,aes(lag,value,colour=Type))+
  geom_point(size=5)+
  geom_line()+
  geom_smooth(alpha=.2)+
xlab('Lag in Km')

NS



#####And for EW#####
depth<-5
ncomp<-4


#First all of them#

lag<-seq(1,21,by = 1)

pedodivs<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%27:48&sites$Var2%in%27:48,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_Normalized_centroid_diversity_pro(pedodiversity_dataset[x,])
  }))
  pedodivs[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs)

####Now for crop####
lag<-seq(1,21,by = 1)

pedodivs_crop<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%27:48&sites$Var2%in%27:48,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='CROP'&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_Normalized_centroid_diversity_pro(pedodiversity_dataset[x,])
  }))
  pedodivs_crop[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs_crop)


####Now for nat####
lag<-seq(1,21,by = 1)

pedodivs_nat<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%27:48&sites$Var2%in%27:48,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='NAT'&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_Normalized_centroid_diversity_pro(pedodiversity_dataset[x,])
  }))
  pedodivs_nat[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs_nat)



# pedodivs_EW<-data.frame(all=unlist(pedodivs),
#                         nat=unlist(pedodivs_nat),
#                         crop=unlist(pedodivs_crop))
pedodivs_EW<-data.frame(nat=unlist(pedodivs_nat),
                        crop=unlist(pedodivs_crop))

pedodivs_EW$lag<-1:21*50


require(reshape)
require(plyr)
pedodivs_EW<-melt(pedodivs_EW,variable_name = 'Type',id.vars = 'lag')

EW<-ggplot(pedodivs_EW,aes(lag,value,colour=Type))+
  geom_point(size=5)+
  geom_line()+
  geom_smooth(alpha=.2)+
  xlab('Lag in Km')

EW
NS



#####Now for HULL DISTANCE####
##at different scales###
depth<-5
ncomp<-4




# Scalewise_Normalized_centroid_diversity_pro<-function(variable_matrix=NULL){
#   variable_matrix<-variable_matrix[,grepl('predicted',colnames(variable_matrix))&!grepl('sd',colnames(variable_matrix))]
#   scores<-as.data.frame(prcomp(variable_matrix,scale. = T,center = T)$x[,1:ncomp])
#   sqrt(sum(apply(scores,2,function(y) (y-mean(y))^2))/nrow(scores))
# }

Scalewise_inhull_diversity_prop <- function(variable_matrix=NULL){
  require(splancs)
  cha<-function(x,y){
    chull(x,y)->i
    return(areapl(cbind(x[i],y[i])))
  }
  variable_matrix<-variable_matrix[,grepl('predicted',colnames(variable_matrix))&!grepl('sd',colnames(variable_matrix))]
  components<-prcomp(variable_matrix,scale. = T,center = T)
  scores <- data.frame(PC1=components$x[,'PC1'],PC2=components$x[,'PC2'])
  cha(scores[,1],scores[,2])
}



# # 
# # WEIRD BEHAVIOUR!..INTERESTING....#
# #First all of them#
# 
# 
# 
# 
# lag<-seq(1,25,by = 1)
# 
# pedodivs<-list()
# for(x in lag){
#   distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
#   sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
#   sites<-sites[sites$Var1%in%0:26,]
#   
#   conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$Sample%in%as.numeric(i))
#   
#   mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
#     Scalewise_inhull_diversity_prop(pedodiversity_dataset[x,])
#   }))
#   pedodivs[[x]]<-mean_pedodiv_in_lag
# }
# 
# 
# plot(lag*50,pedodivs)
# 
# ####Now for crop####
# lag<-seq(1,25,by = 1)
# 
# pedodivs_crop<-list()
# for(x in lag){
#   distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
#   sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
#   sites<-sites[sites$Var1%in%0:26,]
#   
#   conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='CROP'&pedodiversity_dataset$Sample%in%as.numeric(i))
#   
#   mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
#     Scalewise_inhull_diversity_prop(pedodiversity_dataset[x,])
#   }))
#   pedodivs_crop[[x]]<-mean_pedodiv_in_lag
# }
# 
# 
# plot(lag*50,pedodivs_crop)
# 
# ####Now for nat####
# lag<-seq(1,25,by = 1)
# 
# pedodivs_nat<-list()
# for(x in lag){
#   distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
#   sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
#   sites<-sites[sites$Var1%in%0:26,]
#   
#   conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='NAT'&pedodiversity_dataset$Sample%in%as.numeric(i))
#   
#   mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
#     Scalewise_inhull_diversity_prop(pedodiversity_dataset[x,])
#   }))
#   pedodivs_nat[[x]]<-mean_pedodiv_in_lag
# }
# 
# 
# plot(lag*50,pedodivs_nat)
# 
# 
# 
# # pedodivs_NS<-data.frame(all=unlist(pedodivs),
# #                         nat=unlist(pedodivs_nat),
# #                         crop=unlist(pedodivs_crop))
# pedodivs_NS<-data.frame(nat=unlist(pedodivs_nat),
#                         crop=unlist(pedodivs_crop))
# 
# pedodivs_NS$lag<-1:25*50
# 
# 
# require(reshape2)
# pedodivs_NS<-melt(pedodivs_NS,variable_name = 'Type',id.vars = 'lag')
# 
# NS<-ggplot(pedodivs_NS,aes(lag,value,colour=Type))+
#   geom_point(size=5)+
#   geom_line()+
#   geom_smooth(alpha=.2)+
#   xlab('Lag in Km')
# 
# NS
# 
# 
# 
# #####And for EW#####
# depth<-5
# ncomp<-4
# 
# 
# # 
# # 
# # Scalewise_Normalized_centroid_diversity_pro<-function(variable_matrix=NULL){
# #   variable_matrix<-variable_matrix[,grepl('predicted',colnames(variable_matrix))&!grepl('sd',colnames(variable_matrix))]
# #   scores<-as.data.frame(prcomp(variable_matrix,scale. = T,center = T)$x[,1:ncomp])
# #   sqrt(sum(apply(scores,2,function(y) (y-mean(y))^2))/nrow(scores))
# # }
# # 
# 
# #First all of them#
# 
# lag<-seq(1,21,by = 1)
# 
# pedodivs<-list()
# for(x in lag){
#   distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
#   sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
#   sites<-sites[sites$Var1%in%27:48,]
#   
#   conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$Sample%in%as.numeric(i))
#   
#   mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
#     Scalewise_inhull_diversity_prop(pedodiversity_dataset[x,])
#   }))
#   pedodivs[[x]]<-mean_pedodiv_in_lag
# }
# 
# 
# plot(lag*50,pedodivs)
# 
# ####Now for crop####
# lag<-seq(1,21,by = 1)
# 
# pedodivs_crop<-list()
# for(x in lag){
#   distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
#   sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
#   sites<-sites[sites$Var1%in%27:48,]
#   
#   conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='CROP'&pedodiversity_dataset$Sample%in%as.numeric(i))
#   
#   mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
#     Scalewise_inhull_diversity_prop(pedodiversity_dataset[x,])
#   }))
#   pedodivs_crop[[x]]<-mean_pedodiv_in_lag
# }
# 
# 
# plot(lag*50,pedodivs_crop)
# 
# 
# ####Now for nat####
# lag<-seq(1,21,by = 1)
# 
# pedodivs_nat<-list()
# for(x in lag){
#   distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
#   sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
#   sites<-sites[sites$Var1%in%27:48,]
#   
#   conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='NAT'&pedodiversity_dataset$Sample%in%as.numeric(i))
#   
#   mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
#     Scalewise_inhull_diversity_prop(pedodiversity_dataset[x,])
#   }))
#   pedodivs_nat[[x]]<-mean_pedodiv_in_lag
# }
# 
# 
# plot(lag*50,pedodivs_nat)
# 
# 
# 
# # pedodivs_EW<-data.frame(all=unlist(pedodivs),
# #                         nat=unlist(pedodivs_nat),
# #                         crop=unlist(pedodivs_crop))
# pedodivs_EW<-data.frame(nat=unlist(pedodivs_nat),
#                         crop=unlist(pedodivs_crop))
# 
# pedodivs_EW$lag<-1:21*50
# 
# 
# require(reshape)
# require(plyr)
# pedodivs_EW<-melt(pedodivs_EW,variable_name = 'Type',id.vars = 'lag')
# 
# EW<-ggplot(pedodivs_EW,aes(lag,value,colour=Type))+
#   geom_point(size=5)+
#   geom_line()+
#   geom_smooth(alpha=.2)+
#   xlab('Lag in Km')
# 
# EW
# NS

#####Possible reason #
#This function depends of the initial principal componen space

depth<-5
ncomp<-4

variable_matrix_for_hull<-pedodiversity_dataset[,grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))]
components<-prcomp(variable_matrix_for_hull,scale. = T,center = T)
scores <- data.frame(PC1=components$x[,'PC1'],PC2=components$x[,'PC2'])



Scalewise_inhull_diversity_prop <- function(scores=NULL){
  require(splancs)
  cha<-function(x,y){
    chull(x,y)->i
    return(areapl(cbind(x[i],y[i])))
  }
  cha(scores[,1],scores[,2])
}


#First all of them#




lag<-seq(1,25,by = 1)

pedodivs<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%0:26&sites$Var2%in%0:26,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs)

####Now for crop####
lag<-seq(1,25,by = 1)

pedodivs_crop<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%0:26&sites$Var2%in%0:26,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='CROP'&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs_crop[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs_crop)

####Now for nat####
lag<-seq(1,25,by = 1)

pedodivs_nat<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%0:26&sites$Var2%in%0:26,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='NAT'&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs_nat[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs_nat)



# pedodivs_NS<-data.frame(all=unlist(pedodivs),
#                         nat=unlist(pedodivs_nat),
#                         crop=unlist(pedodivs_crop))
pedodivs_NS<-data.frame(nat=unlist(pedodivs_nat),
                        crop=unlist(pedodivs_crop))

pedodivs_NS$lag<-1:25*50


require(reshape2)
pedodivs_NS<-melt(pedodivs_NS,variable_name = 'Type',id.vars = 'lag')

NS<-ggplot(pedodivs_NS,aes(lag,value,colour=Type))+
  geom_point(size=5)+
  geom_line()+
  geom_smooth(alpha=.2)+
  xlab('Lag in Km')

NS



#####And for EW#####
#First all of them#

lag<-seq(1,21,by = 1)

pedodivs<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%27:48&sites$Var2%in%27:48,]

  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs)

####Now for crop####
lag<-seq(1,21,by = 1)

pedodivs_crop<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%27:48&sites$Var2%in%27:488,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='CROP'&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs_crop[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs_crop)


####Now for nat####
lag<-seq(1,21,by = 1)

pedodivs_nat<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%27:48&sites$Var2%in%27:48,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='NAT'&pedodiversity_dataset$Sample%in%as.numeric(i))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs_nat[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs_nat)



# pedodivs_EW<-data.frame(all=unlist(pedodivs),
#                         nat=unlist(pedodivs_nat),
#                         crop=unlist(pedodivs_crop))
pedodivs_EW<-data.frame(nat=unlist(pedodivs_nat),
                        crop=unlist(pedodivs_crop))

pedodivs_EW$lag<-1:21*50


require(reshape)
require(plyr)
pedodivs_EW<-melt(pedodivs_EW,variable_name = 'Type',id.vars = 'lag')

EW<-ggplot(pedodivs_EW,aes(lag,value,colour=Type))+
  geom_point(size=5)+
  geom_line()+
  geom_smooth(alpha=.2)+
  xlab('Lag in Km')

EW
NS


#####################################################################
  ###And with richness vs area (ruled by biodiversity Power law)###
depth<-5
ncomp<-3
 
variable_matrix_for_hull<-pedodiversity_dataset[,grepl('predicted',colnames(pedodiversity_dataset))&!grepl('sd',colnames(pedodiversity_dataset))]
components<-prcomp(variable_matrix_for_hull,scale. = T,center = T)
scores <- data.frame(PC1=components$x[,'PC1'],PC2=components$x[,'PC2'])



# Scalewise_inhull_diversity_prop <- function(scores=NULL){
#   require(splancs)
#   cha<-function(x,y){
#     chull(x,y)->i
#     return(areapl(cbind(x[i],y[i])))
#   }
#   cha(scores[,1],scores[,2])
# }

Scalewise_Normalized_centroid_diversity_pro<-function(variable_matrix=NULL){
  variable_matrix<-variable_matrix[,grepl('predicted',colnames(variable_matrix))&!grepl('sd',colnames(variable_matrix))]
  scores<-as.data.frame(prcomp(variable_matrix,scale. = T,center = T)$x[,1:ncomp])
  sqrt(sum(apply(scores,2,function(y) (y-mean(y))^2))/nrow(scores))
}


#First all of them#




lag<-seq(1,26,by = 1)

pedodivs<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%0:26&sites$Var2%in%0:26,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$Sample%in%as.numeric(i[1]):as.numeric(i[2]))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs)

####Now for crop####
lag<-seq(1,26,by = 1)

pedodivs_crop<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%0:26&sites$Var2%in%0:26,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='CROP'&pedodiversity_dataset$Sample%in%as.numeric(i[1]):as.numeric(i[2]))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs_crop[[x]]<-mean_pedodiv_in_lag
}

plot(lag*50,pedodivs_crop)

####Now for nat####
lag<-seq(1,26,by = 1)

pedodivs_nat<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%0:26&sites$Var2%in%0:26,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='NAT'&pedodiversity_dataset$Sample%in%as.numeric(i[1]):as.numeric(i[2]))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs_nat[[x]]<-mean_pedodiv_in_lag
}

plot(lag*50,pedodivs_nat)


# pedodivs_NS<-data.frame(all=unlist(pedodivs),
#                         nat=unlist(pedodivs_nat),
#                         crop=unlist(pedodivs_crop))
pedodivs_NS<-data.frame(nat=unlist(pedodivs_nat),
                        crop=unlist(pedodivs_crop))

pedodivs_NS$lag<-1:26*50


require(reshape2)
pedodivs_NS<-melt(pedodivs_NS,variable_name = 'Type',id.vars = 'lag')

NS<-ggplot(pedodivs_NS,aes(lag,value,colour=Type))+
  geom_point(size=5)+
  geom_line()+
  geom_smooth(alpha=.2)+
  xlab('Area measured')+
  ylab('Norm_centroid_dist_prop')

NS



#####And for EW#####
#First all of them#
lag<-seq(1,21,by = 1)

pedodivs<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%27:48&sites$Var2%in%27:48,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$Sample%in%as.numeric(i[1]):as.numeric(i[2]))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs[[x]]<-mean_pedodiv_in_lag
}


plot(lag*50,pedodivs)

####Now for crop####
lag<-seq(1,21,by = 1)

pedodivs_crop<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%27:48&sites$Var2%in%27:48,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='CROP'&pedodiversity_dataset$Sample%in%as.numeric(i[1]):as.numeric(i[2]))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs_crop[[x]]<-mean_pedodiv_in_lag
}

plot(lag*50,pedodivs_crop)

####Now for nat####
lag<-seq(1,21,by = 1)

pedodivs_nat<-list()
for(x in lag){
  distances<-expand.grid(pedodiversity_dataset$Sample,pedodiversity_dataset$Sample)
  sites<-unique(distances[as.numeric(distances$Var1)-as.numeric(distances$Var2)==x,])
  sites<-sites[sites$Var1%in%27:48&sites$Var2%in%27:48,]
  
  conditions<-apply(sites,1,function(i) pedodiversity_dataset$top<depth&pedodiversity_dataset$system=='NAT'&pedodiversity_dataset$Sample%in%as.numeric(i[1]):as.numeric(i[2]))
  
  mean_pedodiv_in_lag<-mean(apply(conditions,2,function(x){
    Scalewise_inhull_diversity_prop(scores[x,])
  }))
  pedodivs_nat[[x]]<-mean_pedodiv_in_lag
}

plot(lag*50,pedodivs_nat)


# pedodivs_NS<-data.frame(all=unlist(pedodivs),
#                         nat=unlist(pedodivs_nat),
#                         crop=unlist(pedodivs_crop))
pedodivs_EW<-data.frame(nat=unlist(pedodivs_nat),
                        crop=unlist(pedodivs_crop))

pedodivs_EW$lag<-1:21*50


require(reshape2)
pedodivs_EW<-melt(pedodivs_EW,variable_name = 'Type',id.vars = 'lag')

EW<-ggplot(pedodivs_EW,aes(lag,value,colour=Type))+
  geom_point(size=5)+
  geom_line()+
  geom_smooth(alpha=.2)+
  xlab('Area measured')+
  ylab('Norm_centroid_dist_prop')

EW
NS


