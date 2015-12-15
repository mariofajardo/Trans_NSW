# #Read lab analyses#
# 
# LAB_NSW_NS<-read.csv('DATA/NS_Calibration.csv',
#                      header = T,
#                      dec = '.',
#                      stringsAsFactors=F,skip = 3)[-1,-c(3,6)]
# 
# LAB_NSW_EW<-read.csv('DATA/EW_Calibration.csv',
#                      header = T,
#                      dec = '.',
#                      stringsAsFactors=F,skip = 3)[-1,-c(3,6)]
# 
# 
# LAB_NSW_CORES<-read.csv('DATA/cores_lab_data.csv',
#                         header = T,
#                         dec = '.',
#                         stringsAsFactors=F)[,-c(3,6)]
# 
# ##Bring the spectra in
# 
# Spectra_NSW_chemical <- readRDS('../usyd_spectral_lib/RData/Struct_chem_dataset.rds')
# Spectra_NSW_cores_kattomba <- readRDS('RData/NSW_CORES_with_lab_and_predictionsThu_Oct_22_16_23_20_2015.rds')
# 
# 
# ##Add Join samples NS and EW with chemical samples#
# 
# Spectra_NSW_chemical_responses<-rbind(Spectra_NSW_chemical$NS$responses,
#                                       Spectra_NSW_chemical$EW$responses)
# Spectra_NSW_chemical_responses$Sample<-as.character(Spectra_NSW_chemical_responses$Sample)
# Spectra_NSW_chemical_responses$system<-as.character(Spectra_NSW_chemical_responses$system)
# 
# INDEX_SPECTRA_NS_EW<-with(Spectra_NSW_chemical_responses,paste(Sample,system,top,bottom,sep='_'))
# 
# INDEX_SPECTRA_NS_EW <- paste(INDEX_SPECTRA_NS_EW,c('R1','R2'),sep='_')
# 
# INDEX_LAB_NS_EW <- c(paste(regmatches(LAB_NSW_NS$Code,regexpr('[0-9]+',LAB_NSW_NS$Code,perl=T)),
#                            toupper(regmatches(LAB_NSW_NS$Code,regexpr('[a-zA-Z]+',LAB_NSW_NS$Code,perl=T))),
#                            gsub('-','_',LAB_NSW_NS$Depth),'R1',sep='_'),
#                      paste(regmatches(LAB_NSW_EW$Code,regexpr('[0-9]+',LAB_NSW_EW$Code,perl=T)),
#                            toupper(regmatches(LAB_NSW_EW$Code,regexpr('[a-zA-Z]+',LAB_NSW_EW$Code,perl=T))),
#                            gsub('-','_',LAB_NSW_EW$Depth),'R1',sep='_'))
# 
# Spectra_NSW_chemical_responses$index<-INDEX_SPECTRA_NS_EW
# 
# 
# Spectra_NSW_chemical_spectra <-rbind(Spectra_NSW_chemical$NS$spectra,
#                                      Spectra_NSW_chemical$EW$spectra)
# 
# Spectra_NSW_chemical_spectra$index<-INDEX_SPECTRA_NS_EW
# 
# 
# LAB_NSW<-rbind(LAB_NSW_NS,LAB_NSW_EW)
# 
# 
# LAB_NSW$index<-INDEX_LAB_NS_EW
# require(plyr)
# 
# Spectra_NSW_chemical_responses<-join(Spectra_NSW_chemical_responses,LAB_NSW,
#                                      by = 'index')
# 
# View(Spectra_NSW_chemical_responses)
# 
# 
# 
# #Export the database of Chemical samples#
# 
# 
# Transect_NSW_surface_samples <- list()
# Transect_NSW_surface_samples$spectra<-Spectra_NSW_chemical_spectra[Spectra_NSW_chemical_spectra$index%in%LAB_NSW$index,as.character(350:2500)]
# Transect_NSW_surface_samples$responses <-LAB_NSW
# 
# require(reshape2)
# Transect_NSW_surface_samples$responses$system <- unlist(regmatches(Transect_NSW_surface_samples$responses$Code,regexec('[A-Za-z]+',Transect_NSW_surface_samples$responses$Code)))
# Transect_NSW_surface_samples$responses$site <- unlist(regmatches(Transect_NSW_surface_samples$responses$Code,regexec('[0-9]+',Transect_NSW_surface_samples$responses$Code)))
# Transect_NSW_surface_samples$responses$top <- substring(Transect_NSW_surface_samples$responses$Depth,1,1)
# Transect_NSW_surface_samples$responses$bottom <- as.character(as.numeric(Transect_NSW_surface_samples$responses$top)+5)
# 
# 
# 
# require(rgdal)
# 
# coordinates<-readOGR('DATA/transects_NSEW.kml','WAYPOINT')
# str(coordinates)
# 
# coordinates@data$Name<-toupper(gsub(' ','',coordinates@data$Name))
# 
# coordinates<-data.frame(index_coord=coordinates$Name,coordinates@coords[,-3])
# 
# coordinates$index_coord<-as.character(coordinates$index_coord)
# coordinates$index_coord[c(90)]<-c('31CROP')
# 
# Transect_NSW_surface_samples$responses$index_coord <-toupper(with(Transect_NSW_surface_samples$responses,paste(site,system,sep = '')))
# 
# Transect_NSW_surface_samples$responses<-join(coordinates,Transect_NSW_surface_samples$responses,'index_coord')
# 
# 
# # saveRDS(Transect_NSW_surface_samples,'RData/NSW_CHEMICAL_with_lab_values.RDS')

#####After this predict with kattomba... add the Structchem dataset and... predict the pedodiversity with soil cores
###$i.e. lots of parallel work which ended up like this###

#####From models predictions and laboratory analyses####
Chemical <-readRDS('RData/NSW_CHEMICAL_with_predictions_surface_and_pedodiversityTue_Nov_10_17_25_41_2015.rds')

##### From the slaking index investigations#####
StruTchem <-readRDS('../usyd_spectral_lib/RData/Struct_chem_dataset.rds')

#####From the fuzzy horizons investigation#####
Horizons <-readRDS('RData/horizonsWed_Oct_28_13_25_17_2015.rds')

#####From the MIR (Bruker instrument)#####
MIR_DATA <- readRDS('RData/MIR_DATA_Thu_Oct_29_12_31_07_2015.rds')

#####From Biom tables of microbial analysis#####



#####Build final dataset####
colnames(Chemical$responses)
colnames(Chemical$responses)[9:23]

Chemical$responses[,paste(colnames(Chemical$responses)[9:23],'method',sep='_')] <- 'Laboratory'
# View(Chemical$responses)

####Add slaking coefficient from laboratory (Structchem dataset) and from predictions (Katoomba)
Chemical$responses$Slaking_coef_a <- NA 
# View(Chemical$responses[match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),])
StruTchem$responses$index

Chemical$responses[match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Slaking_coef_a'] <- StruTchem$responses[,2]

# View(Chemical$responses)

Chemical$responses[match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Slaking_coef_a_method'] <- 'Laboratory'
# View(Chemical$responses)


Chemical$responses[-match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Slaking_coef_a']<-as.data.frame(Chemical$predictions)[-match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Slaking_coef_a_pred']

Chemical$responses[-match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Slaking_coef_a_method'] <- 'Predicted'
View(Chemical$responses)

####Add Clay percentage from laboratory (Structchem dataset) and from predictions (Katoomba)
Chemical$responses$Clay <- NA 
# View(Chemical$responses[match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),])
StruTchem$responses$index

Chemical$responses[match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Clay'] <- StruTchem$responses$Clay_perc

# View(Chemical$responses)

Chemical$responses[match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Clay_method'] <- 'Laboratory'
# View(Chemical$responses)


Chemical$responses[-match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Clay']<-as.data.frame(Chemical$predictions)[-match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Clay_pred']

Chemical$responses[-match(StruTchem$responses$index,paste('SAMPLE',Chemical$responses$index,sep='_')),'Clay_method'] <- 'Predicted'
View(Chemical$responses)



####Add horizons information####

Horizons[match(with(Chemical$responses,paste0(site,substring(system,first=1,last=1))),as.character(Horizons$Sample)),]

Chemical$predictions$num_spd_hor<-Horizons[match(with(Chemical$responses,paste0(site,substring(system,first=1,last=1))),as.character(Horizons$Sample)),'num_spd_hor']
Chemical$predictions$first_spd_hor_thickness<-Horizons[match(with(Chemical$responses,paste0(site,substring(system,first=1,last=1))),as.character(Horizons$Sample)),'first_spd_hor_thickness']


#####Add MIR data####

Chemical$MIR_spectra <- MIR_DATA[match(Chemical$responses$index,MIR_DATA$index),]


# #####Add Microbes data ####
# Chemical <- readRDS('RData/NSW_datasetFri_Oct_30_11_38_43_2015.rds')
# Bacterias <-readRDS('RData/OTUS_16S_Fri_Oct_30_12_20_25_2015.rds')
# 
# # 
# Chemical$Bacterias <- Bacterias

# 
# date<-gsub(' |:','_',date())
# saveRDS(Chemical,paste0('RData/','NSW_dataset',date,'.rds'))
# 
# 





#####further modifications####





