CORES<-readRDS('RData/NSW_CORES.rds')
CHEMICAL <-readRDS('RData/NSW_CHEMICAL.rds')

#####Predicting with Katoomba models#####
## For the Chemical dataset    ##
require(Cubist)
require(spectroscopy)
require(plyr)
require(foreach)
require(doSNOW)
####CEC####

prop<-lapply(dir('../usyd_spectral_lib/Katoomba_spectral_models/',pattern = 'model',
                 full.names = T),function(x) readRDS(x))
val <- lapply(dir('../usyd_spectral_lib/Katoomba_spectral_models/',pattern = 'val',
                  full.names = T),function(x) readRDS(x))
input<-rbind(CHEMICAL$EW$spectra,
             CHEMICAL$NS$spectra,
             CORES$EW$spectra,
             CORES$NS$spectra)
names(prop)<-c('CEC','Clay','EC','pH','TC')
#provide parallel backend first#
cl <-makeCluster(5)
setMKLthreads(1)
registerDoSNOW(cl)

source('../usyd_spectral_lib/Katoomba_predict.R')
predictions<-do.call(cbind,mapply(Katoomba_predict,prop,val,SIMPLIFY = F,USE.NAMES = T))

stopCluster(cl)


#####Update dataset ####

details<-rbind(data.frame(dataset='Chemical',rbind(CHEMICAL$EW$responses,
             CHEMICAL$NS$responses)),
             data.frame(dataset='Cores',rbind(CORES$EW$details,
             CORES$NS$details)))

details<-cbind(details,predictions)

details<-split(details,details$dataset)

CHEMICAL$NS$responses<-details$Chemical[gsub('[A-a]','',details$Chemical$Sample)%in%c(0:26),] 
CHEMICAL$EW$responses<-details$Chemical[!gsub('[A-a]','',details$Chemical$Sample)%in%c(0:26),] 

CORES$NS$details<-details$Cores[gsub('[A-a]','',details$Cores$Sample)%in%c(0:26),] 
CORES$EW$details<-details$Cores[!gsub('[A-a]','',details$Cores$Sample)%in%c(0:26),] 

CORES<-saveRDS(CORES,file = 'RData/NSW_CORES_with_predictions_Katoomba.rds',compress = 'xz')
CHEMICAL <-saveRDS(CHEMICAL,file = 'RData/NSW_CHEMICAL_with_predictions_Katoomba.rds',compress = 'xz')

