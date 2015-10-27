CORES<-readRDS('RData/NSW_CORES_with_lab_values_18_8_2105.RDS')


#####Create models#####
## For the Cores dataset    ##
require(Cubist)
require(spectroscopy)
require(plyr)
require(foreach)
require(doSNOW)
require(clhs)

props<-CORES$calibration[,6:20]

spectra<-CORES$calibration[,21:2171]
components<-prcomp(spectra)
cal_samples <- clhs(as.data.frame(components$x[,1:10]),size=80,iter=15000)
set.seed(1)
cal_iter <- lapply(1:50,function(x) sample(cal_samples,size=79))

#####spectra pretreatments####
spectra <- filter_sg(spectra,n = 11,p = 2,m = 0)
spectra <- trimSpec(spectra,c(450,2450),350:2500)
spectra <- strip_spectra(spectra,450:2450,c(450,2450),which = 10)
spectra <- snvBLC(spectra)

#provide parallel backend first#
cl <-makeCluster(4)
setMKLthreads(2)
registerDoSNOW(cl)

models<-foreach(properties=colnames(props),                
                .packages = 'Cubist') %:% 
  foreach(iter=cal_iter) %dopar% {
    cubist(spectra[iter,],props[iter,properties],committees = 5)
  }

names(models)<-colnames(props)

stopCluster(cl)

cl <-makeCluster(8)
setMKLthreads(1)
registerDoSNOW(cl)

preds<-foreach(model=models,                
                .packages = 'Cubist') %dopar% {
  
  sapply(model,FUN = function(model_iter,iter){
    predict(model_iter,spectra[-cal_samples,])
  })
}
names(preds)<-colnames(props)


stopCluster(cl)

View(preds[[1]])
for (i in colnames(props)){
  print(i)
print(goof(props[-cal_samples,i],rowMeans(preds[[i]])))
}
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
CORES<-readRDS('RData/NSW_CORES_with_lab_values_18_8_2105.RDS')

require(Cubist)
require(spectroscopy)
require(plyr)
require(foreach)
require(doSNOW)
require(clhs)

props<-CORES$calibration[,6:20]
spectra<-CORES$calibration[,21:2171]

#####CEC####
input <- spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_CEC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_CEC.RData')

source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')
goof(CORES$calibration$ECEC,predictions$Mean,main='CEC')


input <- CORES$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_CEC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_CEC.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

CORES$responses$CEC_pred<-predictions$Mean
CORES$responses$CEC_pred_sd<-predictions$Standard_deviation


#####pH####
props<-CORES$calibration[,6:20]
spectra<-CORES$calibration[,21:2171]


input <- spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_pH.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_pH.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

goof(CORES$calibration$pH.Level..H2O.,predictions$Mean,main='pH')


input <- CORES$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_pH.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_pH.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

CORES$responses$pH_pred<-predictions$Mean
CORES$responses$pH_pred_sd<-predictions$Standard_deviation



#####EC####
#** predictions for values lower than 3 ds/m**#
props<-CORES$calibration[,6:20]
spectra<-CORES$calibration[,21:2171]


input <- spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_EC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_EC.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

goof(CORES$calibration$Conductivity,predictions$Mean,main='EC')

input <- CORES$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_EC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_EC.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

CORES$responses$EC_pred<-predictions$Mean
CORES$responses$EC_pred_sd<-predictions$Standard_deviation

#####TC####

props<-CORES$calibration[,6:20]
spectra<-CORES$calibration[,21:2171]


input <- spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_TC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_TC.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')


input <- CORES$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_TC.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_TC.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

CORES$responses$TC_pred<-predictions$Mean
CORES$responses$TC_pred_sd<-predictions$Standard_deviation


#####Clay####
props<-CORES$calibration[,6:20]
spectra<-CORES$calibration[,21:2171]


input <- spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Clay.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Clay.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')


input <- CORES$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Clay.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Clay.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

CORES$responses$Clay_pred<-predictions$Mean
CORES$responses$Clay_pred_sd<-predictions$Standard_deviation

# date<-gsub(' |:','_',date())
# saveRDS(CORES,paste0('RData/','NSW_CORES_with_lab_and_predictions',date,'.rds'))

#####Slaking_coef_a####

CORES<-readRDS('RData/NSW_CORES_with_lab_and_predictionsFri_Aug_28_13_34_47_2015.rds')
input <- CORES$spectra
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Slaking_a_coef.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Slaking_a_coef.RData')
source('../usyd_spectral_lib/Katoomba_spectral_models_JAVA.R')

CORES$responses$Slaking_coef_a_pred<-predictions$Mean
CORES$responses$Slaking_coef_a_pred_sd<-predictions$Standard_deviation

# date<-gsub(' |:','_',date())
# saveRDS(CORES,paste0('RData/','NSW_CORES_with_lab_and_predictions',date,'.rds'))
# 


View(CORES$responses)


require(ggplot2)
require(reshape2)

pedons<-by(CORES$responses,CORES$responses$Sample,function(x) x)

pdf('Plots/check_Pedons_Cores_predictions.pdf')

for (i in pedons){

  preds<-i[,colnames(i)[c(1:5,6,8,10,12,14,16)]]
  
  preds<-melt(preds,id.vars = colnames(preds)[1:5])
  
  print(ggplot(preds,aes(top,value,colour=variable))+
  geom_line(size=2)+
  coord_flip()+
  scale_x_continuous('Depth',trans='reverse'))
}
dev.off()

