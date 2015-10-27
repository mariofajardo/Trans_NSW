####This script will read the data from the asd instrument#####
#####################core_samples_3_rep_2cm_resolution###
library(pbapply)
library(reshape)
library(spectroscopy)
library(plyr)
setwd("..//Functions")#####    SETWD
load('correct_steps.RData')
load('check_plots.RData')
load('pr_varExp.RData')
prev_dir <- getwd()
setwd("../../Spectra/Transect/Core_samples/")#####    SETWD


hor <- sort_df(read.csv('NIR_hor.csv',header=T))     
list_df <- as.list.data.frame(split(hor,hor$id))
hor_names <- lapply(list_df,function(x){
  tmp <- data.frame(hor = rep(x$hor,times=x$bottom-x$top+1),master_hor = rep(x$master_horizon,times=x$bottom-x$top+1))
  tmp
})


###########read the spectra and assign horizon names###########
require(naturalsort)
pboptions(type='txt',style=3,title='Reading files')
files <- dir(pattern='txt',recursive=T)
files <- files[naturalorder(files)]
tmp_data <- pblapply(files,function(x){   
tmp <- read.csv(x,header=T)
data <- tmp[-1]
rep <- as.factor(rep(1:nrow(data), each=3, length.out=nrow(data)))
tmp2 <- as.data.frame(sapply(data, tapply, rep, mean))
tmp2$File.Name <- paste(regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T)),paste0('Spectrum',sprintf('%05i',c(1:nrow(tmp2)))),sep='_')
sample <- regmatches(x,regexpr("[[:alnum:]]+(?=[.])",x,perl=T))
top <- seq(0,by=2,length.out=nrow(tmp2))
bottom <- seq(2,by=2,length.out=nrow(tmp2))
tmp3 <- cbind(File.Name=tmp2$File.Name,Sample=sample,top=top,bottom=bottom,tmp2[-length(tmp2)])
final <-tmp3[order(tmp3$File.Name,decreasing=F),]
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


#####check the spectral differences####
pdf('../../../Codes/transect_n_s/Plots/check_spectra.pdf',width = 10,height = 6)
require(ggplot2)
require(reshape2)
require(spectroscopy)
require(signal)
load('../../../Codes/transect_n_s/RData/filter_sg2.RData')
sample<-1:length(files)
for (i in sample){
  spectra <-data.frame(tmp_data[[i]][,-c(1:4,2156)][rep(1:nrow(tmp_data[[i]][,-c(1:4,2156)]), each = 2), ])
  details <-c(seq(0,tmp_data[[i]]$top[length(tmp_data[[i]]$top)]),tmp_data[[i]]$top[length(tmp_data[[i]]$top)]+1)
  spectra<-spectra[naturalorder(details),]
  colnames(spectra)<-seq(350,2500)
  levels(spectra$sample) <- naturalsort(levels(spectra$sample))
  ####make some process####
  spectra <- trimSpec(spectra,wavlimits = c(500,2450),as.numeric(colnames(spectra)))
  spectra <-1/log(spectra)
  spectra <-filter_sg(spectra,n=11,p=2,m=0)
  spectra <-strip_spectra(spectra,as.numeric(colnames(spectra),c(500,2450),which=5))
  spectra <-snvBLC(spectra)
#   
  ###and some wavlength filtering ####
  spectra<-filter_sg2(spectra,n=5)
  spectra <- data.frame(sample=details,spectra)
  input_data<-melt(spectra,id.vars = 'sample')
  print(ggplot(input_data,aes(variable,sample, fill = value)) + 
    geom_tile() +
    labs(x='Wavelength (nm)',
         y='Depth (cm)',
         title = paste0("Main absorbance differences in depth sample : ",tmp_data[[i]]$Sample[1])) +
    scale_x_discrete(breaks=c(500,1000,1500,2000,2450))+
      scale_y_reverse(breaks=c(0,10,20,30,40,50,60,70,80,90,100))+
    scale_fill_gradient2(low = "green", high = "red",mid='white')+
    theme(axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          title=element_text(size=20)))
}
dev.off()
setwd('../../../Codes/transect_n_s/Plots/')
shell.exec('check_spectra.pdf')
setwd(prev_dir)



cores_NSW_DATA <- tmp_data

rawg_spectra<-lapply(cores_NSW_DATA,function(x) {
  spectra<-x[5:2155]
  spectra})

sample_details<-lapply(cores_NSW_DATA,function(x) {
  spectra<-x[c(1:4,2157)]
  spectra})

# ####Save for later####
spectra <- lapply(rawg_spectra,function(x) x<-correct_step(x))
# spectra <- lapply(rawg_spectra,function(x) x<-filter_sg(x,n = 11,p = 2,m = 0))
# spectra <- lapply(rawg_spectra,function(x) x<-trimSpec(x,c(450,2450),350:2500))
# spectra <- lapply(rawg_spectra,function(x) x<-strip_spectra(x,450:2450,c(450,2450),which = 5))

details <-do.call(rbind,sample_details)
spectra<-do.call(rbind,spectra)

NSW_EW <- list()
NSW_NS <- list()

NSW_NS$spectra <- spectra[gsub('[A-a]','',details$Sample)%in%c(0:26),]
NSW_EW$spectra <- spectra[!gsub('[A-a]','',details$Sample)%in%c(0:26),]

NSW_NS$details <- details[gsub('[A-a]','',details$Sample)%in%c(0:26),]
NSW_EW$details <- details[!gsub('[A-a]','',details$Sample)%in%c(0:26),]

NSW_CORES<- list()
NSW_CORES$EW <- NSW_EW
NSW_CORES$NS <- NSW_NS

saveRDS(object =NSW_CORES,file = '../../../Codes/transect_n_s/RData/NSW_CORES.Rds')
#####
####Add chemical data and predictions####

NSW_CORES_import <-readRDS('RData/NSW_CORES.Rds')
lab_data <-read.csv('DATA/cores_lab_data.csv',stringsAsFactors=F)[,-c(1,3,6)]

require(naturalsort)
NSW_CORES<-list()
NSW_CORES$spectra <- rbind(NSW_CORES_import$NS$spectra,NSW_CORES_import$EW$spectra)
NSW_CORES$responses <- rbind(NSW_CORES_import$NS$details,NSW_CORES_import$EW$details)

NSW_CORES$calibration$responses <-lab_data[naturalorder(lab_data$Code),] 
NSW_CORES$calibration$responses$Depth<-gsub('-','_',NSW_CORES$calibration$responses$Depth)

NSW_CORES$calibration$responses$top<-regmatches(NSW_CORES$calibration$responses$Depth,regexpr('[[:print:]]+(?=[[:punct:]])',NSW_CORES$calibration$responses$Depth,perl = T))
NSW_CORES$calibration$responses$bottom<-regmatches(NSW_CORES$calibration$responses$Depth,regexpr('(?<=[[:punct:]])[[:print:]]+',NSW_CORES$calibration$responses$Depth,perl = T))

NSW_CORES$calibration$responses$site<-regmatches(NSW_CORES$calibration$responses$Code,regexpr('[0-9]+',NSW_CORES$calibration$responses$Code,perl = T)) 
NSW_CORES$calibration$responses$system<-regmatches(NSW_CORES$calibration$responses$Code,regexpr('[[:alpha:]]+',NSW_CORES$calibration$responses$Code,perl = T)) 
NSW_CORES$calibration$responses$sample <-toupper(with(NSW_CORES$calibration$responses,paste0(site,system)))
NSW_CORES$calibration$responses$system<-gsub('C','Crop',NSW_CORES$calibration$responses$system)
NSW_CORES$calibration$responses$system<-gsub('N','Nat',NSW_CORES$calibration$responses$system)

NSW_CORES$whole<-data.frame(NSW_CORES$responses[NSW_CORES$responses$Sample%in%NSW_CORES$calibration$responses$sample,],
                            NSW_CORES$spectra[NSW_CORES$responses$Sample%in%NSW_CORES$calibration$responses$sample,])

NSW_CORES$whole<-droplevels(NSW_CORES$whole)

NSW_CORES$calibration<-do.call(rbind,by(NSW_CORES$whole,NSW_CORES$whole$Sample,function(x) {
  selection_top<-as.numeric(naturalsort(NSW_CORES$calibration$responses$top[NSW_CORES$calibration$responses$sample%in%as.character(x$Sample)]))
  selection_bottom<-as.numeric(naturalsort(NSW_CORES$calibration$responses$bottom[NSW_CORES$calibration$responses$sample%in%as.character(x$Sample)]))
  
  properties <-NSW_CORES$calibration$response[NSW_CORES$calibration$response$sample%in%x$Sample,c(4:19)]
  properties<- properties[naturalorder(properties$Depth),]
  test<-data.frame(Sample=unique(x$Sample),
                   top=selection_top,
                   bottom=selection_bottom,
                   system=unique(x$system),
                   properties,
                   do.call(rbind,mapply(function(up,down) colMeans(x[x$top>=up&x$bottom<=down,grepl('X',colnames(x))]),selection_top,selection_bottom,SIMPLIFY = F)))
  })) 

NSW_CORES<-list(calibration=NSW_CORES$calibration,
                spectra=NSW_CORES$spectra,
                responses=NSW_CORES$responses)
# saveRDS(NSW_CORES,'RData/NSW_CORES_with_lab_values_18_8_2105.RDS')

