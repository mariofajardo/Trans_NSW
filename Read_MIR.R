####First check if everything is OK####
# ###olivier plate order (5 to 10 march)####
# Plate1<-data.frame(Letter=rep(toupper(letters[1:8]),each=12),
#                    Number=1:12,
#                    Sample=c(rep(0:47,each=2)[-96],
#                             'blank'),
#                    System=rep(c('Nat','Crop')),
#                    stringsAsFactors=F)
# Plate1[96,4] <- 'blank'  
# 
# Plate1$Position <- paste0(Plate1$Letter,Plate1$Number)
# Plate1$top <-5
# Plate1$bottom <-10
# 
# ###olivier plate order(5 to 10 march)####
# 
# Plate2 <- data.frame(Letter=toupper(c('a','a','a','h')),
#                      Number=c(1:3,12),
#                      Sample=c(47,48,48,
#                               'blank'),
#                      System=rep(c('Crop','Nat','Crop','blank')),
#                      stringsAsFactors=F)
# Plate2$Position <- paste0(Plate2$Letter,Plate2$Number)
# Plate2$top <-5
# Plate2$bottom <-10
# 
# ####samples from 0 to 5 october ####
# Plate3 <- read.csv('../../MIR/Third_run/Plate_3.csv')[,-7]
# Plate3$Position <- toupper(paste0(Plate3$Letter.psition,Plate3$Number.Position))
# colnames(Plate3) <- c("Letter","Number","top","bottom","Sample","System","Position") 
# 
# #####samples from 5 to 10 october####
# Plate4 <- Plate1
# 
# #####samples for testing assumptions####
# Plate5 <- read.csv('../../MIR/Fourth_run/Plate_5.csv')[,-7]
# Plate5$Position <- paste0(Plate5$Letter.psition,Plate5$Number.Position)
# colnames(Plate5) <- c("Letter","Number","top","bottom","Sample","System","Position") 
# 
# #####samples for moisture test####
# Plate_wet <- Plate1
# 
# 
# 
# #Read spectra#
# 
# files_MIR_bruker_Plate1 <- dir('../../MIR/Second_run/Plate1/',full.names = T) #5 to 10 march
# files_MIR_bruker_Plate2 <- dir('../../MIR/Second_run/Plate2/',full.names = T)
# files_MIR_bruker_Plate3 <- dir(path = '../../MIR/Third_run/',pattern = '.1')
# files_MIR_bruker_Plate4 <- dir(path = '../../MIR/Second_run_sec_blank/',pattern = '.2')
# files_MIR_bruker_Plate5 <- dir(path = '../../MIR/Fourth_run/',pattern = '.0')
# test_with_moisture <- dir('../../MIR/First_run/',full.names = T)
# 
# MIR_bruker1 <- lapply(files_MIR_bruker_Plate1,function(x){
#   tmp <- read.csv(x,header=F)
#   tmp2<-data.frame(Position=substring(regmatches(x,regexpr('([[:word:]])+(?=[.2])',x,perl = T)),first = 6),t(tmp[,2]))
#   colnames(tmp2)[-1] <- ceiling(tmp[,1])
#   tmp2
# })
# 
# MIR_bruker1<-do.call(rbind,MIR_bruker1)
# 
# MIR_bruker2 <- lapply(files_MIR_bruker_Plate2,function(x){
#   tmp <- read.csv(x,header=F)
#   tmp2<-data.frame(Position=substring(regmatches(x,regexpr('([[:word:]])+(?=[.0])',x,perl = T)),first = 7),t(tmp[,2]))
#   colnames(tmp2)[-1] <- ceiling(tmp[,1])
#   tmp2
# })
# 
# MIR_bruker2<-do.call(rbind,MIR_bruker2)
# 
# 
# setwd('../../MIR/Third_run/')
# MIR_bruker3 <- lapply(files_MIR_bruker_Plate3,function(x){
#   tmp <- read.csv(x,header=F)
#   tmp2<-data.frame(Position=substring(regmatches(x,regexpr('([[:word:]])+(?=[.0])',x,perl = T)),first = 7),t(tmp[,2]))
#   colnames(tmp2)[-1] <- ceiling(tmp[,1])
#   tmp2
# })
# 
# MIR_bruker3<-do.call(rbind,MIR_bruker3)
# 
# setwd('../../MIR/Second_run_sec_blank/')
# 
# MIR_bruker4 <- lapply(files_MIR_bruker_Plate4,function(x){
#   tmp <- read.csv(x,header=F)
#   tmp2<-data.frame(Position=substring(regmatches(x,regexpr('([[:word:]])+(?=[.0])',x,perl = T)),first = 7),t(tmp[,2]))
#   colnames(tmp2)[-1] <- ceiling(tmp[,1])
#   tmp2
# })
# 
# MIR_bruker4<-do.call(rbind,MIR_bruker4)
# 
# setwd('../../MIR/Fourth_run/')
# 
# MIR_bruker5 <- lapply(files_MIR_bruker_Plate5,function(x){
#   tmp <- read.csv(x,header=F)[1:1763,]
#   tmp2<-data.frame(Position=substring(regmatches(x,regexpr('([[:word:]])+(?=[.0])',x,perl = T)),first = 8),t(tmp[,2]))
#   colnames(tmp2)[-1] <- ceiling(tmp[,1])
#   tmp2
# })
# 
# MIR_bruker5<-do.call(rbind,MIR_bruker5)
# 
# setwd('../../MIR/First_run/')
# 
# MIR_bruker_wet <- lapply(test_with_moisture,function(x){
#   tmp <- read.csv(x,header=F)
#   tmp2<-data.frame(Position=substring(regmatches(x,regexpr('([[:word:]])+(?=[.2])',x,perl = T)),first = 6),t(tmp[,2]))
#   colnames(tmp2)[-1] <- ceiling(tmp[,1])
#   tmp2
# })
# 
# MIR_bruker_wet<-do.call(rbind,MIR_bruker_wet)
# 
# ####Now check####
# 
# 
# require(spectroscopy)
# 
# 
# #### mir1 = 5 -10 march in blue
# #### mir3 = 5 -10 october in red
# #### mir4 = 5 -10 march with new blank in green
# #### mir5 = 5 10 october new measurement in purple
# 
# #sample 45 crop 5 to 10
# plot(as.numeric(colnames(MIR_bruker1)[-1]),snvBLC(MIR_bruker1[MIR_bruker1$Position=='H8',-1]),type='l',col='blue') 
# lines(as.numeric(colnames(MIR_bruker3)[-1]),snvBLC(MIR_bruker3[MIR_bruker3$Position=='H10',-1]),type='l',col='red')
# lines(as.numeric(colnames(MIR_bruker4)[-1]),snvBLC(MIR_bruker4[MIR_bruker4$Position=='H8',-1]),type='l',col='green')
# lines(as.numeric(colnames(MIR_bruker5)[-1]),snvBLC(MIR_bruker5[MIR_bruker5$Position=='A8',-1]),type='l',col='purple')
# 
# plot(as.numeric(colnames(MIR_bruker3)[-1]),snvBLC(MIR_bruker3[MIR_bruker3$Position=='H10',-1]),type='l',col='red') 
# lines(as.numeric(colnames(MIR_bruker5)[-1]),snvBLC(MIR_bruker5[MIR_bruker5$Position=='A8',-1]),type='l',col='purple')
# 
# 
# 
# #sample 47 crop 5 to 10
# #### mir2 = 5 - 10 march in blue
# #### mir5 = 5 - 10 october in red
# #### mirwet = 5 - 10 march first testing violet
# 
# plot(as.numeric(colnames(MIR_bruker2)[-1]),snvBLC(MIR_bruker2[MIR_bruker2$Position=='A1',-1]),type='l',col='blue') 
# lines(as.numeric(colnames(MIR_bruker5)[-1]),snvBLC(MIR_bruker5[MIR_bruker5$Position=='A10',-1]),type='l',col='red')
# lines(as.numeric(colnames(MIR_bruker_wet)[-1]),snvBLC(MIR_bruker_wet[MIR_bruker_wet$Position=='A1',-1]),type='l',col='violet')
# 
# 
# #sample 48 nat 5 to 10
# #### mir2 = 5 - 10 march in blue
# #### mir5 = 5 - 10 october in red
# 
# plot(as.numeric(colnames(MIR_bruker2)[-1]),snvBLC(MIR_bruker2[MIR_bruker2$Position=='A2',-1]),type='l',col='blue') 
# lines(as.numeric(colnames(MIR_bruker5)[-1]),snvBLC(MIR_bruker5[MIR_bruker5$Position=='A12',-1]),type='l',col='red')
# 
# 
# ####Some samples from 0 - 5 
# ####plus some ground in march and scanned in october from 0 to 5 in red
# plot(as.numeric(colnames(MIR_bruker5)[-1]),snvBLC(MIR_bruker5[MIR_bruker5$Position=='A1',-1]),type='l',col='purple') 
# lines(as.numeric(colnames(MIR_bruker3)[-1]),snvBLC(MIR_bruker3[MIR_bruker3$Position=='H8',-1]),type='l',col='red',lwd=3) 
# lines(as.numeric(colnames(MIR_bruker2)[-1]),snvBLC(MIR_bruker2[MIR_bruker2$Position=='A1',-1]),type='l',col='blue',lwd=3) 
# 

#####Compile the dataset####

Plate1 <- read.csv('../../MIR/Third_run/Plate_3.csv')[,-7]
Plate1$Position <- toupper(paste0(Plate1$Letter.psition,Plate1$Number.Position))
colnames(Plate1) <- c("Letter","Number","top","bottom","Sample","System","Position") 


Plate2<-read.csv('../../MIR/Fourth_run/Plate_5.csv')[,-7]
Plate2$Position <- toupper(paste0(Plate2$Letter.psition,Plate2$Number.Position))
colnames(Plate2) <- c("Letter","Number","top","bottom","Sample","System","Position") 


Plate3<-data.frame(Letter=rep(toupper(letters[1:8]),each=12),
                   Number=1:12,
                   Sample=c(rep(0:47,each=2)[-96],
                            'blank'),
                   System=rep(c('Nat','Crop')),
                   stringsAsFactors=F)
Plate3[96,4] <- 'blank' 
Plate3$Position <- paste0(Plate3$Letter,Plate3$Number)
Plate3$top <-5
Plate3$bottom <-10

files_MIR_bruker_Plate1 <- dir('../../MIR/Third_run/',full.names = T,pattern = '.1') 
files_MIR_bruker_Plate2 <- dir('../../MIR/Fourth_run/',full.names = T,pattern = '.0')
files_MIR_bruker_Plate3 <- dir('../../MIR/Fifth_run/',full.names = T,pattern = '.0')



MIR_bruker1 <- lapply(files_MIR_bruker_Plate1,function(x){
    tmp <- read.csv(x,header=F)
    tmp2<-data.frame(Position=substring(regmatches(x,regexpr('([[:word:]])+(?=[.1])',x,perl = T)),first = 7),t(tmp[,2]))
    colnames(tmp2)[-1] <- ceiling(tmp[,1])
    tmp2
  })
  
MIR_bruker1<-do.call(rbind,MIR_bruker1)


MIR_bruker2 <- lapply(files_MIR_bruker_Plate2,function(x){
   tmp <- read.csv(x,header=F)[1:1763,]
   tmp2<-data.frame(Position=substring(regmatches(x,regexpr('([[:word:]])+(?=[.0])',x,perl = T)),first = 8),t(tmp[,2]))
   colnames(tmp2)[-1] <- ceiling(tmp[,1])
   tmp2
 })
  
MIR_bruker2<-do.call(rbind,MIR_bruker2)


MIR_bruker3 <- lapply(files_MIR_bruker_Plate3,function(x){
  tmp <- read.csv(x,header=F)[1:1763,]
  tmp2<-data.frame(Position=substring(regmatches(x,regexpr('([[:word:]])+(?=[.0])',x,perl = T)),first = 7),t(tmp[,2]))
  colnames(tmp2)[-1] <- ceiling(tmp[,1])
  tmp2
})

MIR_bruker3<-do.call(rbind,MIR_bruker3)


MIR_DATA_0_5 <- rbind(merge(Plate1[c(1:92),],MIR_bruker1,by = 'Position'),
                      merge(Plate2[c(1:6),],MIR_bruker2,by = 'Position'))

MIR_DATA_5_10 <- rbind(merge(Plate2[c(10:12),],MIR_bruker2,by = 'Position'),
                      merge(Plate3,MIR_bruker3,by = 'Position'))

require(naturalsort)
MIR_DATA_0_5 <- MIR_DATA_0_5[naturalorder(MIR_DATA_0_5$System),]
MIR_DATA_0_5 <- MIR_DATA_0_5[naturalorder(MIR_DATA_0_5$Sample),]

MIR_DATA_5_10 <- MIR_DATA_5_10[naturalorder(MIR_DATA_5_10$System),]
MIR_DATA_5_10 <- MIR_DATA_5_10[naturalorder(MIR_DATA_5_10$Sample),]


MIR_DATA <- rbind(MIR_DATA_0_5,MIR_DATA_5_10)

MIR_DATA$index <- with(MIR_DATA,toupper(paste(Sample,System,top,bottom,'R1',sep='_'))) 


# date<-gsub(' |:','_',date())
# saveRDS(MIR_DATA,paste0('RData/','MIR_DATA_',date,'.rds'))
