CORES<-readRDS('RData/NSW_CORES.rds')
require(spectroscopy)
require(ggplot2)
spectra<-CORES$NS$spectra[CORES$NS$details$Sample=='8C',]
spectra <- trimSpec(spectra,wavlimits = c(500,2450),as.numeric(colnames(spectra)))
spectra <-1/log(spectra)
spectra <-filter_sg(spectra,n=11,p=2,m=0)
spectra <-strip_spectra(spectra,as.numeric(colnames(spectra),c(500,2450),which=5))
spectra <-snvBLC(spectra)

spectra <- data.frame(sample=CORES$NS$details[CORES$NS$details$Sample=='8C','top'],spectra)
require(reshape2)
input_data<-melt(spectra,id.vars = 'sample')


ggplot(input_data,aes(variable,sample, fill = value)) + 
        geom_tile() +
        labs(x='Wavelength (nm)',
             y='Depth (cm)') +
        scale_x_discrete(breaks=c(500,1000,1500,2000,2450))+
        scale_y_reverse(breaks=c(0,10,20,30,40,50,60,70,80,90,100))+
        scale_fill_gradient2(low = "green", high = "red",mid='black')+
        theme(axis.title.x = element_text(size=20),
              axis.title.y = element_text(size=20),
              legend.text=element_text(size=15),
              legend.title=element_text(size=15),
              title=element_text(size=20))

to_plot <-rbind(data.frame(Hor_ID='A',Wave=seq(500,2450),value=as.numeric(spectra[1,-1])),
                data.frame(Hor_ID='Bt',Wave=seq(500,2450),value=as.numeric(spectra[25,-1])),
                data.frame(Hor_ID='BC',Wave=seq(500,2450),value=as.numeric(spectra[50,-1])))

ggplot(to_plot,aes(x=Wave,y=value,group=Hor_ID,colour=Hor_ID))+
  geom_line(size=1.3)+
  labs(x='Wavelength (nm)',
       y='Absorbance',
       title='Spectra of different horizons')+
  theme_bw()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        title=element_text(size=20))
  

