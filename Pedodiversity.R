####This script will read the data from the asd instrument#####
#####################core_samples_3_rep_2cm_resolution###
library(pbapply)
library(reshape)
library(spectroscopy)
library(plyr)

load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/clustering_horizons/correct_steps.RData')
load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/clustering_horizons/pr_varExp.RData')
load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/clustering_horizons/check_plots.RData')
load('../../../../../PhD/Project_Docs/Thesis/thesis_scripts/RData/filter_sg2.RData')


setwd("../../Spectra/Transect/Core_samples/")#####    SETWD



###########read the spectra and assign horizon names###########
require(naturalsort)
pboptions(type='txt',style=3,title='Reading files')
files <- dir(pattern='txt',recursive=T)
files <- files[naturalorder(files)]
tmp_data <- pblapply(files,function(x){   
  tmp <- read.csv(x,header=T)
  data <- tmp[-1]
  rep <- rep(1:nrow(data), each=3, length.out=nrow(data))
  File.Name <- paste(regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)),paste0('Spectrum',sprintf('%05i',c(1:nrow(data)))),sep='_')
  sample <- rep(regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)),length.out=nrow(data))
  top <- rep(seq(0,length.out = floor(length(sample)/3),by = 2),each=3)
  bottom <- rep(seq(2,length.out = floor(length(sample)/3),by = 2),each=3)
  tmp3 <- data.frame(File.Name=File.Name,Sample=sample,top=top,bottom=bottom,Rep=rep,data[-length(data)])
  })  

tmp_data<-lapply(tmp_data,function(x){
  x$system<-gsub('[0-9]','',x$Sample)
  x}
)

tmp_data<-lapply(tmp_data,function(x){
  x$system[x$system=='C']<-'CROP'
  x}
)

tmp_data<-lapply(tmp_data,function(x){
  x$system[x$system=='N']<-'NAT'
  x}
)

tmp_data<-lapply(tmp_data,function(x){
  x$Sample<-regmatches(x$Sample,regexpr('[[:digit:]]+',x$Sample,perl = T))
  x}
)
pedodiversity_dataset <-do.call(rbind,tmp_data)
pedodiversity_dataset$index<-with(pedodiversity_dataset,paste(Sample,system,top,bottom,sep='_'))
pedodiversity_dataset <- pedodiversity_dataset[pedodiversity_dataset$top<10,]
# saveRDS(pedodiversity_dataset,file = '../../../Codes/transect_n_s/RData/pedodiversity_dataset.RDS')
#####Predict properties and pedodiversity by selected depth increments####
require(Cubist)
require(spectroscopy)
require(plyr)
require(foreach)
require(doSNOW)
require(clhs)


#####CEC####
setwd("~/University of Sydney/PhD/Data/Lab_work/Codes/transect_n_s")
CORES<-readRDS('RData/NSW_CORES_with_lab_values_18_8_2105.RDS')
pedodiversity_dataset <- readRDS('RData/pedodiversity_dataset.RDS') 

props<-CORES$calibration[,6:20]
spectra<-CORES$calibration[,21:2171]

#####CEC####
input <- CORES$calibration[,-c(6:20)]
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_CEC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_CEC.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')
goof(CORES$calibration$ECEC,predictions$Mean,main='CEC')


input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

cl <-makeCluster(5)
registerDoSNOW(cl)
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_CEC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_CEC.RData')
source('Katoomba_Properties_pedodiversity.R')
stopCluster(cl)


pedodiversity_dataset$approx_top<-rep(c(rep(0,6),rep('out',3),rep(5,6)),length.out=length(pedodiversity_dataset$top))
pedodiversity_dataset$approx_bottom<-rep(c(rep(5,6),rep('out',3),rep(10,6)),length.out=length(pedodiversity_dataset$top))

pedodiversity_dataset$index<-with(pedodiversity_dataset,toupper(paste(Sample,system,approx_top,approx_bottom,sep='_')))

CEC_diversity<-by(predictions$Mean,pedodiversity_dataset$index,function(x) sd(x))
# CEC_diversity1<-by(predictions$Mean,pedodiversity_dataset$index,function(x) diversity(x,'simpson'))

pedodiversity <- data.frame(Sample=names(CEC_diversity),CEC_diversity=unlist(as.list(CEC_diversity),use.names = F))

Chemical<-readRDS('RData/NSW_CHEMICAL_with_lab_and_predictionThu_Oct_22_17_12_51_2015.rds')

Chemical$predictions$CEC_diversity <- pedodiversity[match(with(Chemical$responses,toupper(paste(site,system,top,bottom,sep='_'))),pedodiversity$Sample),'CEC_diversity']

#####pH####

input <- CORES$calibration[,-c(6:20)]
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_pH.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_pH.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')
goof(CORES$calibration$pH.Level..CaCl2.,predictions$Mean,main='pH')


input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

cl <-makeCluster(5)
registerDoSNOW(cl)
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_pH.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_pH.RData')
source('Katoomba_Properties_pedodiversity.R')
stopCluster(cl)


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
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Clay.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Clay.RData')
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
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_EC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_EC.RData')
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
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_TC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_TC.RData')
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

# date<-gsub(' |:','_',date())
# saveRDS(Chemical,paste0('RData/','NSW_CHEMICAL_with_lab_and_prediction',date,'.rds'))
