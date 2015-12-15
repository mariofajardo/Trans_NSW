Transect_NSW_surface_samples<-readRDS('RData/NSW_CHEMICAL_with_lab_values.RDS')

texture_lab<-read.csv('DATA/psa.csv')
Transect_NSW_surface_samples$responses$index <-paste0('SAMPLE_',Transect_NSW_surface_samples$responses$index) 
texture_lab$index <-paste0(toupper(texture_lab$.id),'_R1')

test<-merge(Transect_NSW_surface_samples$responses,texture_lab,by.x ='index',by.y ='index')

test[,c(10:24,26:29)]<-sapply(colnames(test)[c(10:24,26:29)],function(x) as.numeric(test[,x]))

test<- test[,c(1:28,42)]
colnames(test)

test <- test[match(Transect_NSW_surface_samples$responses$index,test$index),]
test <- test[complete.cases(test),]

samples_in<-Transect_NSW_surface_samples$responses$index%in%test$index
Transect_NSW_surface_samples$spectra<-Transect_NSW_surface_samples$spectra[samples_in,]
Transect_NSW_surface_samples$responses<-test
saveRDS(Transect_NSW_surface_samples,'RData/NSW_CHEMICAL_with_lab_values_13_8_2105_revised.RDS')
