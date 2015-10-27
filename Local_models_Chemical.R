Chemical<-readRDS('RData/NSW_CHEMICAL_with_lab_values.RDS')

Chemical$responses[,9:23] <- apply(Chemical$responses[,9:23],2,as.numeric)
View(Chemical$responses)


# #####Create models#####
# ## For the Chemical dataset    ##
# require(Cubist)
# require(spectroscopy)
# require(plyr)
# require(foreach)
# require(doSNOW)
# require(clhs)
# 
# props<-Chemical$responses
# 
# spectra<-Chemical$spectra
# components<-prcomp(spectra)
# cal_samples <- clhs(as.data.frame(components$x[,1:10]),size=80,iter=15000)
# set.seed(1)
# cal_iter <- lapply(1:50,function(x) sample(cal_samples,size=79))
# 
# #####spectra pretreatments####
# spectra <- filter_sg(spectra,n = 11,p = 2,m = 0)
# spectra <- trimSpec(spectra,c(450,2450),350:2500)
# spectra <- strip_spectra(spectra,450:2450,c(450,2450),which = 10)
# spectra <- snvBLC(spectra)
# 
# #provide parallel backend first#
# cl <-makeCluster(4)
# setMKLthreads(2)
# registerDoSNOW(cl)
# 
# models<-foreach(properties=colnames(props),                
#                 .packages = 'Cubist') %:% 
#   foreach(iter=cal_iter) %dopar% {
#     cubist(spectra[iter,],props[iter,properties],committees = 5)
#   }
# 
# names(models)<-colnames(props)
# 
# stopCluster(cl)
# 
# cl <-makeCluster(8)
# setMKLthreads(1)
# registerDoSNOW(cl)
# 
# preds<-foreach(model=models,                
#                .packages = 'Cubist') %dopar% {
#                  
#                  sapply(model,FUN = function(model_iter,iter){
#                    predict(model_iter,spectra[-cal_samples,])
#                  })
#                }
# names(preds)<-colnames(props)
# 
# 
# stopCluster(cl)
# 
# View(preds[[1]])
# for (i in colnames(props)){
#   print(i)
#   print(goof(props[-cal_samples,i],rowMeans(preds[[i]])))
# }
# #####Update dataset ####
# 
# details<-rbind(data.frame(dataset='Chemical',rbind(CHEMICAL$EW$responses,
#                                                    CHEMICAL$NS$responses)),
#                data.frame(dataset='Cores',rbind(CORES$EW$details,
#                                                 CORES$NS$details)))
# 
# details<-cbind(details,predictions)
# 
# details<-split(details,details$dataset)
# 
# CHEMICAL$NS$responses<-details$Chemical[gsub('[A-a]','',details$Chemical$Sample)%in%c(0:26),] 
# CHEMICAL$EW$responses<-details$Chemical[!gsub('[A-a]','',details$Chemical$Sample)%in%c(0:26),] 
# 
# CORES$NS$details<-details$Cores[gsub('[A-a]','',details$Cores$Sample)%in%c(0:26),] 
# CORES$EW$details<-details$Cores[!gsub('[A-a]','',details$Cores$Sample)%in%c(0:26),] 
# 
# CORES<-saveRDS(CORES,file = 'RData/NSW_CORES_with_predictions_Katoomba.rds',compress = 'xz')
# CHEMICAL <-saveRDS(CHEMICAL,file = 'RData/NSW_CHEMICAL_with_predictions_Katoomba.rds',compress = 'xz')


#####Predictions using Katoomba####
# CORES<-readRDS('RData/NSW_CORES_with_lab_values_18_8_2105.RDS')

require(Cubist)
require(spectroscopy)
require(plyr)
require(foreach)
require(doSNOW)
require(clhs)


#####CEC####
input <- Chemical$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_CEC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_CEC.RData')

source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')
goof(Chemical$responses$ECEC,predictions$Mean,main='CEC')


input <- Chemical$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_CEC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_CEC.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

Chemical$predictions$CEC_pred<-predictions$Mean
Chemical$predictions$CEC_pred_sd<-predictions$Standard_deviation


#####Clay####
input <- Chemical$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Clay.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Clay.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')


input <- Chemical$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Clay.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Clay.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

Chemical$predictions$Clay_pred<-predictions$Mean
Chemical$predictions$Clay_pred_sd<-predictions$Standard_deviation

#####Slaking_coef_a####

input <- Chemical$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Slaking_a_coef.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Slaking_a_coef.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')


Structchem<-readRDS('../usyd_spectral_lib/RData/Struct_chem_dataset.rds')
Chemical$predictions$Slaking_coef_a_pred<-predictions$Mean
Chemical$predictions$Slaking_coef_a_pred_sd<-predictions$Standard_deviation

goof(Structchem$responses$SI_a,predictions$Mean[Chemical$responses$index%in%substring(Structchem$responses$index,first=8)])

date<-gsub(' |:','_',date())
saveRDS(Chemical,paste0('RData/','NSW_CHEMICAL_with_lab_and_prediction',date,'.rds'))

