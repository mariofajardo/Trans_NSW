#Read lab analyses#

LAB_NSW_NS<-read.csv('DATA/NS_Calibration.csv',
                     header = T,
                     dec = '.',
                     stringsAsFactors=F,skip = 3)[-1,-c(3,6)]

LAB_NSW_EW<-read.csv('DATA/EW_Calibration.csv',
                     header = T,
                     dec = '.',
                     stringsAsFactors=F,skip = 3)[-1,-c(3,6)]

#fix "date error" ...probably an excel thing#

LAB_NSW_EW$Depth <- gsub('10-May','5-10',LAB_NSW_EW$Depth)
LAB_NSW_NS$Depth <- gsub('10-May','5-10',LAB_NSW_NS$Depth)

LAB_NSW_CORES<-read.csv('DATA/Cores_calibration.csv',
                     header = T,
                     dec = '.',
                     stringsAsFactors=F,skip = 3)[-1,-c(3,6)]

##Bring the spectra in

Spectra_NSW_chemical <- readRDS('RData/NSW_CHEMICAL.rds')


##Add Join samples NS and EW with chemical samples and predict the one that is missing plus the reps#
##and validate the models##

Spectra_NSW_chemical_responses<-rbind(Spectra_NSW_chemical$NS$responses,
                                      Spectra_NSW_chemical$EW$responses)
Spectra_NSW_chemical_responses$Sample<-as.character(Spectra_NSW_chemical_responses$Sample)
Spectra_NSW_chemical_responses$system<-as.character(Spectra_NSW_chemical_responses$system)

INDEX_SPECTRA_NS_EW<-with(Spectra_NSW_chemical_responses,paste(Sample,system,top,bottom,sep='_'))

INDEX_SPECTRA_NS_EW <- paste(INDEX_SPECTRA_NS_EW,c('R1','R2'),sep='_')

INDEX_LAB_NS_EW <- c(paste(regmatches(LAB_NSW_NS$Code,regexpr('[0-9]+',LAB_NSW_NS$Code,perl=T)),
                      toupper(regmatches(LAB_NSW_NS$Code,regexpr('[a-zA-Z]+',LAB_NSW_NS$Code,perl=T))),
                      gsub('-','_',LAB_NSW_NS$Depth),'R1',sep='_'),
                     paste(regmatches(LAB_NSW_EW$Code,regexpr('[0-9]+',LAB_NSW_EW$Code,perl=T)),
                           toupper(regmatches(LAB_NSW_EW$Code,regexpr('[a-zA-Z]+',LAB_NSW_EW$Code,perl=T))),
                           gsub('-','_',LAB_NSW_EW$Depth),'R1',sep='_'))

Spectra_NSW_chemical_responses$index<-INDEX_SPECTRA_NS_EW


Spectra_NSW_chemical_spectra <-rbind(Spectra_NSW_chemical$NS$spectra,
                                     Spectra_NSW_chemical$EW$spectra)
  
Spectra_NSW_chemical_spectra$index<-INDEX_SPECTRA_NS_EW


LAB_NSW<-rbind(LAB_NSW_NS,LAB_NSW_EW)


LAB_NSW$index<-INDEX_LAB_NS_EW
require(plyr)

Spectra_NSW_chemical_responses<-join(Spectra_NSW_chemical_responses,LAB_NSW,
                                        by = 'index')

View(Spectra_NSW_chemical_responses)

# ###make models to predict rep 2 chemical values ###
# 
# require(Cubist)
# require(clhs)
# require(spectroscopy)
# 
# scores <-as.data.frame(prcomp(Spectra_NSW_chemical_spectra[Spectra_NSW_chemical_spectra$index%in%INDEX_LAB_NS_EW,as.character(seq(350,2500))])$x[,1:20])
# set.seed(1)
# train <- clhs(scores,size=nrow(scores)*.9,iter=15000)
# 
# require(doSNOW)
# require(foreach)
# 
# 
# spectra <-Spectra_NSW_chemical_spectra[Spectra_NSW_chemical_spectra$index%in%INDEX_LAB_NS_EW,as.character(seq(350,2500))]
# spectra_predicted <- Spectra_NSW_chemical_spectra[!Spectra_NSW_chemical_spectra$index%in%INDEX_LAB_NS_EW,as.character(seq(350,2500))] 
# 
# pdf('Plots/check_models.pdf')
# models<-lapply(colnames(LAB_NSW)[6:20],function(x){
#   
#   cl<-makeCluster(8)
#   registerDoSNOW(cl)
#   
#   tmp_preds<-foreach(iters=1:16,
#                      .packages = c('Cubist','spectroscopy','plyr'),
#                      .export = c('spectra','train','LAB_NSW','spectra_predicted')) %dopar% {  
#   
#   set.seed(iters+1)  
#   iteration<-sample(train,size = floor(length(train)*.95))
# #   if(x%in%colnames(LAB_NSW_NS)[6]){
# #   model<-cubist(spectra[iteration,],as.numeric(LAB_NSW_NS[iteration,x]))
# #   }else{
#   model<-cubist(spectra[iteration,],as.numeric(LAB_NSW[iteration,x]),committees = 5,control = cubistControl(rules = 2))
# #   }
#   tmp_result<-list()
#   tmp_result$val_preds<-predict(model,spectra[-train,])
#   tmp_result$preds<-predict(model,spectra_predicted)
#   tmp_result
#   }
#   
#   acc <- goof(as.numeric(LAB_NSW[-train,x]),colMeans(ldply(tmp_preds,function(x) x$val_preds)))
#   text(as.numeric(LAB_NSW[-train,x]),colMeans(ldply(tmp_preds,function(x) x$val_preds)),labels = LAB_NSW[-train,'index'])
#   title(main=paste0(x,' R2:', round(acc$R2,2)))
#   preds<-colMeans(ldply(tmp_preds,function(x) x$preds))
#   stand_dev<-apply(ldply(tmp_preds,function(x) x$preds),2,sd)  
#   results<-list(preds=preds,sd=stand_dev,acc=acc)
#   stopCluster(cl)
#   results
#   })
# 
# dev.off()
# 
# ###not very good models.... :(   #

#Export the database of Chemical samples#


Transect_NSW_surface_samples <- list()
Transect_NSW_surface_samples$spectra<-Spectra_NSW_chemical_spectra[Spectra_NSW_chemical_spectra$index%in%LAB_NSW$index,as.character(350:2500)]
Transect_NSW_surface_samples$responses <-LAB_NSW

require(reshape2)
Transect_NSW_surface_samples$responses$system <- unlist(regmatches(Transect_NSW_surface_samples$responses$Code,regexec('[A-Za-z]+',Transect_NSW_surface_samples$responses$Code)))
Transect_NSW_surface_samples$responses$site <- unlist(regmatches(Transect_NSW_surface_samples$responses$Code,regexec('[0-9]+',Transect_NSW_surface_samples$responses$Code)))
Transect_NSW_surface_samples$responses$top <- substring(Transect_NSW_surface_samples$responses$Depth,1,1)
Transect_NSW_surface_samples$responses$bottom <- as.character(as.numeric(Transect_NSW_surface_samples$responses$top)+5)



require(rgdal)

coordinates<-readOGR('DATA/transects_NSEW.kml','WAYPOINT')
str(coordinates)

coordinates@data$Name<-toupper(gsub(' ','',coordinates@data$Name))

coordinates<-data.frame(index_coord=coordinates$Name,coordinates@coords[,-3])

coordinates$index_coord<-as.character(coordinates$index_coord)
coordinates$index_coord[c(90)]<-c('31CROP')

Transect_NSW_surface_samples$responses$index_coord <-toupper(with(Transect_NSW_surface_samples$responses,paste(site,system,sep = '')))

Transect_NSW_surface_samples$responses<-join(coordinates,Transect_NSW_surface_samples$responses,'index_coord')


saveRDS(Transect_NSW_surface_samples,'RData/NSW_CHEMICAL_with_lab_values.RDS')

