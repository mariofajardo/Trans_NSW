####Data overview####

CORES<-readRDS('RData/NSW_CORES.rds')
CHEMICAL <-readRDS('RData/NSW_CHEMICAL.rds')

CORES_comp<-prcomp(rbind(CORES$EW$spectra,CORES$NS$spectra))

CORES_scores <-cbind(rbind(data.frame(dataset='CORES_EW',CORES$EW$details),
                     data.frame(dataset='CORES_NS',CORES$NS$details)),
                     CORES_comp$x[,1:20])

CORES_scores$system<-gsub('[0-9]','',CORES_scores$Sample)

CHEMICAL_comp<-prcomp(rbind(CHEMICAL$EW$spectra,CHEMICAL$NS$spectra))
CHEMICAL_scores <-cbind(rbind(data.frame(dataset='CHEMICAL_EW',CHEMICAL$EW$responses),
                           data.frame(dataset='CHEMICAL_NS',CHEMICAL$NS$responses)),
                        CHEMICAL_comp$x[,1:20])


require(ggplot2)

##Check cores dataset##

ggplot(CORES_scores[CORES_scores$top<2,],aes(PC1,PC2,colour=dataset))+
  geom_point(size=2)+
  geom_text(aes(label=Sample),vjust=3)+
  facet_wrap(~dataset*system)



ggplot(CHEMICAL_scores,aes(PC1,PC2,colour=dataset))+
  geom_point(size=4)+
  geom_text(aes(label=Sample),vjust=3)+
  facet_wrap(~dataset*system)

PROJECTED_CORES<-data.frame(CORES_scores[,1:5],predict(CHEMICAL_comp,rbind(CORES$EW$spectra,CORES$NS$spectra)))
PROJECTED_CORES$system<-gsub('[0-9]','',PROJECTED_CORES$Sample)


base_chemical<-ggplot(CHEMICAL_scores,aes(PC1,PC2,colour=dataset))+
  geom_point(size=4)+
  geom_text(aes(label=Sample),vjust=3)+
  facet_wrap(~dataset*system)

base_chemical+
  geom_point(data = PROJECTED_CORES[PROJECTED_CORES$top<4,])+
  geom_text(data = PROJECTED_CORES[PROJECTED_CORES$top<4,],
            aes(x=PC1,
                y=PC2,
                label=Sample),vjust=3)+
  facet_wrap(~dataset*system)









require(clhs)
set.seed(123)
train_data_cores<-clhs(CORES_scores[,6:25],
                          size = 100, 
                          iter = 15000, 
                          progress = T, 
                          simple = T)
# write.csv(CORES_scores[train_data_cores,1:5],file='DATA/calibration_samples_core.csv')


ggplot(data = CORES_scores,aes(PC1,PC2))+
geom_point(alpha=.2)+
geom_point(data=CORES_scores[train_data_cores,],aes(PC1,PC2),colour='red',size=3)

###93 out of the 100 samples did not have enough material ...so I mixed some samples
## now ...based in the spectral distribution I will pick 7 samples more 
##from the sub-surface observations in order to have higher representativity in B-C horizons

require(clhs)
set.seed(123)
new_scores<-CORES_scores[-train_data_cores,]

set.seed(111)
train_data_cores_2<-clhs(new_scores[new_scores$top>20,6:25],
                       size = 7, 
                       iter = 15000, 
                       progress = T, 
                       simple = T)


ggplot(data = CORES_scores,aes(PC1,PC2))+
  geom_point(alpha=.2)+
  geom_point(data=CORES_scores[train_data_cores,],aes(PC1,PC2),colour='red',size=3)+
  geom_point(data=new_scores[train_data_cores_2,],aes(PC1,PC2),colour='blue',size=3)

write.csv(file='DATA/calibration_samples_core2.csv',new_scores[train_data_cores_2,1:5])
