####This script wil read the data from the asd instrument#####
#####samples_N_S####
#####################chemical_5_rep  - 2_depths###
library(pbapply)
library(reshape)
library(spectroscopy)
library(naturalsort)
prev_dir <- getwd()
setwd("..//Functions")#####    SETWD
load('correct_steps.RData')
load('check_plots.RData')
setwd("../../Spectra/Transect/Chemical_samples_n_s/")#####    SETWD
###########read the spectra and assign horizon names###########
pboptions(type='txt',style=3,title='Reading files')
files <- dir(pattern='txt',recursive=T)
files <- naturalsort(files)
tmp_data <- pblapply(files,function(x){   
  tmp <- read.csv(x,header=T)[,-2153] #for some reason the last variables are just NA (ASD instrument)
  tmp<-tmp[naturalorder(tmp$File.Name),]
  data <- tmp[-1]
  rep <- as.factor(rep(1:nrow(data), each=5, length.out=nrow(data)))
  tmp2 <- data.frame(sapply(data, tapply, rep, mean,USE.NAMES = F))
  File_Name <- paste('sample',regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)),
                     paste0('Spectrum',
                            sprintf('%05i',c(1:nrow(tmp2)))),
                     sep='_')
  
  sample <- as.character(regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)))
  top <- c(0,0,5,5,0,0,5,5)
  bottom <- c(5,5,10,10,5,5,10,10)
  system <- c('CROP','CROP','CROP','CROP','NAT','NAT','NAT','NAT')
  tmp3 <- data.frame(File.Name=File_Name,system=system,Sample=sample,top=top,bottom=bottom,tmp2,row.names = NULL)
  final <-tmp3[order(tmp3$File.Name,decreasing=F),]
}) 


rawg_spectra<-lapply(tmp_data,function(x) {
  spectra<-x[6:2156]
  spectra})

sample_details<-lapply(tmp_data,function(x) {
  spectra<-x[1:5]
  spectra})

names(sample_details)<-paste0('sample_',mapply(function(x,y) regmatches(y,regexpr('[0-9]+(?=[.])',x,perl = T)),files,files))
names(rawg_spectra)<-paste0('sample_',mapply(function(x,y) regmatches(y,regexpr('[0-9]+(?=[.])',x,perl = T)),files,files))

rawg_spectra<-lapply(rawg_spectra,function(x) {
  names(x)<-seq(350,2500)
  x})

setwd(prev_dir)

require(spectroscopy)
spectra <- lapply(rawg_spectra,function(x) x<-correct_step(x))
spectra <- lapply(rawg_spectra,function(x) x<-filter_sg(x,n = 11,p = 2,m = 0))
spectra <- lapply(rawg_spectra,function(x) x<-trimSpec(x,c(450,2450),350:2500))
spectra <- lapply(rawg_spectra,function(x) x<-strip_spectra(x,450:2450,c(450,2450),which = 5))
# check_plots(spectra,Title = 'CHECK',Ylab = 'CHECK')

NSW_NS <- list()


NSW_NS$spectra <- do.call(rbind,spectra)
NSW_NS$responses <-do.call(rbind,sample_details)

require(ggplot2)

###take spectrum from 2, 9, 14 and 23 (i.e. site 1,8,13,22####

samples_poster <- c(2,9,14,23)

data_poster<-data.frame(Site=rep(paste0('Site',c(1,8,13,22)),each=2),System=rep(c('Crop','Nat'),4),do.call(rbind,lapply(spectra[samples_poster],function(x) rbind(colMeans(x[1:2,]),colMeans(x[5:6,])))))
require(reshape2)

data_poster_melt <- melt(data_poster,id.vars = colnames(data_poster)[1:2])

data_poster_melt$Site <-factor(data_poster_melt$Site,levels(data_poster_melt$Site)[c(1,4,2,3)])

spectra_poster <- ggplot(data_poster_melt,aes(as.numeric(variable),y=value))+
  geom_line(aes(linetype=System),size=1,colour='red')+
  theme_minimal()+
  ylab('Reflectance')+
  xlab('Wavelength (nm)')+
  facet_wrap(~Site,nrow=4)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))
  
spectra_poster

# ggsave(filename = '../../../Codes/transect_n_s/Plots/spectra_poster.pdf',spectra_poster,width = 28,height = 42,units='cm')

#####PCA for poster samples ####

scores <-as.data.frame(prcomp(NSW_NS$spectra,center = T,scale. = T)$x)

scores_poster <- ggplot(scores,aes(PC1,PC2))+
  geom_point(aes(colour=NSW_NS$responses$system))+
  theme_bw()+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))+
  geom_text(aes(PC1,PC2,label=NSW_NS$responses$Sample),vjust=2)
  

scores_poster

require(ellipse)
conf_ellipse <-data.frame()
conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(scores[NSW_NS$responses$Sample%in%0:5,], 
                                                               ellipse(cor(PC1,PC2), 
                                                                       scale=c(sd(PC1),sd(PC2)), 
                                                                       centre=c(mean(PC1),mean(PC2)),
                                                                       level=.90))),Sample='Sites 0 to 5'))
conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(scores[NSW_NS$responses$Sample%in%6:26,], 
                                                             ellipse(cor(PC1,PC2), 
                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                     centre=c(mean(PC1),mean(PC2)),
                                                                     level=.90))),Sample='Sites 6 to 26'))


#####plot confidence ellipses####
setwd(prev_dir)
load('pr_varExp.RData')
no_out_pr_Exp<-pr_varExp(NSW_NS$spectra)

scores$System <-NSW_NS$responses$system

scores_poster<-ggplot(scores, aes(x=PC1, y=PC2,colour=System)) + 
  geom_point(alpha=.6,size=7) +
  labs(x=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'), 
       y=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'),
       title=paste0('Cumulative explanation :',
                    round(no_out_pr_Exp[2,2],2),'%')) +
  geom_path(data=conf_ellipse, aes_string(x='x', y='y',colour='Sample'), size=1, linetype=1)+
  geom_text(aes(x=PC1, y=PC2,label=NSW_NS$responses$Sample),vjust=2,size=5)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        title=element_text(size=30),
        legend.text=element_text(size = 20))


scores_poster

# ggsave(filename = 'Plots/spectra_scores_poster.pdf',scores_poster,width = 42,height = 35,units='cm')



#####samples_E_W####
library(pbapply)
library(reshape)
library(spectroscopy)
library(naturalsort)
prev_dir <- getwd()
setwd("..//Functions")#####    SETWD
load('correct_steps.RData')
load('check_plots.RData')
setwd("../../Spectra/Transect/chemical_samples_e_w/")#####    SETWD
###########read the spectra and assign horizon names###########
pboptions(type='txt',style=3,title='Reading files')
files <- dir(pattern='txt',recursive=T)
files <- naturalsort(files)
tmp_data <- pblapply(files,function(x){   
  tmp <- read.csv(x,header=T)[,-2153] #for some reason the last variables are just NA (ASD instrument)
  tmp<-tmp[naturalorder(tmp$Wavelength),]
  data <- tmp[-1]
  rep <- as.factor(rep(1:nrow(data), each=5, length.out=nrow(data)))
  tmp2 <- data.frame(sapply(data, tapply, rep, mean,USE.NAMES = F))
  File_Name <- paste('sample',regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)),
                     paste0('Spectrum',
                            sprintf('%05i',c(1:nrow(tmp2)))),
                     sep='_')
  
  sample <- as.character(regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)))
  top <- c(0,0,5,5,0,0,5,5)
  bottom <- c(5,5,10,10,5,5,10,10)
  system <- c('CROP','CROP','CROP','CROP','NAT','NAT','NAT','NAT')
  tmp3 <- data.frame(File.Name=File_Name,system=system,Sample=sample,top=top,bottom=bottom,tmp2,row.names = NULL)
  final <-tmp3[order(tmp3$File.Name,decreasing=F),]
}) 


rawg_spectra<-lapply(tmp_data,function(x) {
  spectra<-x[6:2156]
  spectra})

sample_details<-lapply(tmp_data,function(x) {
  spectra<-x[1:5]
  spectra})

names(sample_details)<-paste0('sample_',mapply(function(x,y) regmatches(y,regexpr('[0-9]+(?=[.])',x,perl = T)),files,files))
names(rawg_spectra)<-paste0('sample_',mapply(function(x,y) regmatches(y,regexpr('[0-9]+(?=[.])',x,perl = T)),files,files))

rawg_spectra<-lapply(rawg_spectra,function(x) {
  names(x)<-seq(350,2500)
  x})

setwd(prev_dir)

require(spectroscopy)
spectra <- lapply(rawg_spectra,function(x) x<-correct_step(x))
spectra <- lapply(rawg_spectra,function(x) x<-filter_sg(x,n = 11,p = 2,m = 0))
spectra <- lapply(rawg_spectra,function(x) x<-trimSpec(x,c(450,2450),350:2500))
spectra <- lapply(rawg_spectra,function(x) x<-strip_spectra(x,450:2450,c(450,2450),which = 5))
# check_plots(spectra,Title = 'CHECK',Ylab = 'CHECK')

NSW_EW <- list()


NSW_EW$spectra <- do.call(rbind,spectra)
NSW_EW$responses <-do.call(rbind,sample_details)



scores <-rbind(data.frame(Sample=NSW_NS$responses$Sample,
                          Transect='N_S',
                          system=NSW_NS$responses$system,
                          top=NSW_NS$responses$top,
                          bottom=NSW_NS$responses$bottom,
                          prcomp(NSW_NS$spectra,center = T,scale. = T)$x)[,1:10],
               data.frame(Sample=NSW_EW$responses$Sample,
                          Transect='E_W',
                          system=NSW_EW$responses$system,
                          top=NSW_EW$responses$top,
                          bottom=NSW_EW$responses$bottom,
                          prcomp(NSW_EW$spectra,center = T,scale. = T)$x)[,1:10]
               )


scores_poster <- ggplot(scores,aes(PC1,PC2))+
  geom_point(aes(colour=Transect),size=5)+
  theme_bw()+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))+
  geom_text(aes(PC1,PC2,label=Sample),vjust=2)+
  facet_wrap(~system)


scores_poster

require(ellipse)
conf_ellipse <-data.frame()
conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(scores[NSW_NS$responses$Sample%in%0:5,], 
                                                             ellipse(cor(PC1,PC2), 
                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                     centre=c(mean(PC1),mean(PC2)),
                                                                     level=.90))),Sample='Sites 0 to 5'))
conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(scores[NSW_NS$responses$Sample%in%6:26,], 
                                                             ellipse(cor(PC1,PC2), 
                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                     centre=c(mean(PC1),mean(PC2)),
                                                                     level=.90))),Sample='Sites 6 to 26'))


#####plot confidence ellipses####
setwd(prev_dir)
load('pr_varExp.RData')
no_out_pr_Exp<-pr_varExp(NSW_NS$spectra)

scores$System <-NSW_NS$responses$system

scores_poster<-ggplot(scores, aes(x=PC1, y=PC2,colour=System)) + 
  geom_point(alpha=.6,size=7) +
  labs(x=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'), 
       y=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'),
       title=paste0('Cumulative explanation :',
                    round(no_out_pr_Exp[2,2],2),'%')) +
  geom_path(data=conf_ellipse, aes_string(x='x', y='y',colour='Sample'), size=1, linetype=1)+
  geom_text(aes(x=PC1, y=PC2,label=NSW_NS$responses$Sample),vjust=2,size=5)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        title=element_text(size=30),
        legend.text=element_text(size = 20))


scores_poster


# 
# require(clhs)
# set.seed(134)
# calibration_set<-clhs(as.data.frame(scores[,1:10]),size=round(nrow(scores)*.3,0),iter=50000,progress=T)
# 
# 
# write.csv(NSW_NS$responses[calibration_set,],file='../../../Codes/transect_n_s/DATA/calibration_samples.csv')
# 
