CORES<-readRDS('RData/NSW_CORES_with_predictions.rds')
CHEMICAL <-readRDS('RData/NSW_CHEMICAL_with_predictions.rds')

require(reshape2)
require(ggplot2)
pdf(file = 'Plots/check_Pedons_Cores_NS.pdf',width = 10,height = 16)
for(i in unique(CORES$NS$details$Sample)){
 pedon<-CORES$NS$details[CORES$NS$details$Sample==i,] 
 properties<-colnames(pedon)[grepl('Mean',colnames(pedon),perl = T)]
 property_name<-gsub('.Mean','',properties,perl = T)
 pedon[,gsub('Mean','Max',properties)]<-pedon[,properties]+pedon[,colnames(pedon)[grepl('Standard',colnames(pedon))]]
 pedon[,gsub('Mean','Min',properties)]<-pedon[,properties]-pedon[,colnames(pedon)[grepl('Standard',colnames(pedon))]]
 
 for (j in property_name){
 pedon_tmp<-melt(pedon,id.vars = c('top',colnames(pedon)[grepl(j,colnames(pedon))]))
 print(ggplot(pedon_tmp,aes_string('top',paste0(j,'.Mean')))+
  geom_line(size=1,color='green')+
  geom_line(aes_string('top',paste0(j,'.Min'),size=1),color='red')+
  geom_line(aes_string('top',paste0(j,'.Max'),size=1),,color='blue')+
  coord_flip()+
  scale_x_continuous('Depth',trans='reverse')+
  geom_ribbon(ymax=pedon_tmp[,paste0(j,'.Max')],
              ymin=pedon_tmp[,paste0(j,'.Min')],
              fill='lightgrey',
              alpha=.7)+
    ggtitle(label=paste0('Pedon: ',i,'; Predicted property: ',j))
 )
  }
}
dev.off()


####Cores EW ###

pdf(file = 'Plots/check_Pedons_Cores_EW.pdf',width = 10,height = 16)
for(i in unique(CORES$EW$details$Sample)){
  pedon<-CORES$EW$details[CORES$EW$details$Sample==i,] 
  properties<-colnames(pedon)[grepl('Mean',colnames(pedon),perl = T)]
  property_name<-gsub('.Mean','',properties,perl = T)
  pedon[,gsub('Mean','Max',properties)]<-pedon[,properties]+pedon[,colnames(pedon)[grepl('Standard',colnames(pedon))]]
  pedon[,gsub('Mean','Min',properties)]<-pedon[,properties]-pedon[,colnames(pedon)[grepl('Standard',colnames(pedon))]]
  
  for (j in property_name){
    pedon_tmp<-melt(pedon,id.vars = c('top',colnames(pedon)[grepl(j,colnames(pedon))]))
    print(ggplot(pedon_tmp,aes_string('top',paste0(j,'.Mean')))+
            geom_line(size=2,color='green')+
            geom_line(aes_string('top',paste0(j,'.Min'),size=1),color='red')+
            geom_line(aes_string('top',paste0(j,'.Max'),size=1),,color='blue')+
            coord_flip()+
            scale_x_continuous('Depth',trans='reverse')+
            geom_ribbon(ymax=pedon_tmp[,paste0(j,'.Max')],
                        ymin=pedon_tmp[,paste0(j,'.Min')],
                        fill='lightgrey',
                        alpha=.7)+
            ggtitle(label=paste0('Pedon: ',i,'; Predicted property: ',j))+
            scale_y_continuous(j)
    )
  }
}
dev.off()

####CHEMICAL NS ###

pdf(file = 'Plots/check_Pedons_CHEMICAL_NS.pdf',width = 10,height = 16)
for(i in unique(CHEMICAL$NS$responses$Sample)){
  for(y in unique(CHEMICAL$NS$responses$system)){
  pedon<-CHEMICAL$NS$responses[CHEMICAL$NS$responses$Sample==i&CHEMICAL$NS$responses$system==y,] 
  properties<-colnames(pedon)[grepl('Mean',colnames(pedon),perl = T)]
  property_name<-gsub('.Mean','',properties,perl = T)
  pedon[,gsub('Mean','Max',properties)]<-pedon[,properties]+pedon[,colnames(pedon)[grepl('Standard',colnames(pedon))]]
  pedon[,gsub('Mean','Min',properties)]<-pedon[,properties]-pedon[,colnames(pedon)[grepl('Standard',colnames(pedon))]]
  
  for (j in property_name){
    pedon_tmp<-melt(pedon,id.vars = c('top','system',colnames(pedon)[grepl(j,colnames(pedon))]))
    print(ggplot(pedon_tmp,aes_string('top',paste0(j,'.Mean')),color='green')+
            geom_line(size=1)+
            coord_flip()+
            scale_x_continuous('Depth',trans='reverse')+
            geom_ribbon(ymax=pedon_tmp[,paste0(j,'.Max')],
                        ymin=pedon_tmp[,paste0(j,'.Min')],
                        fill='lightgrey',
                        alpha=.7)+
            ggtitle(label=paste0('Pedon: ',i,'; Predicted property: ',j))+
            scale_y_continuous(j)

    )
  }
}
}
dev.off()


####Check predicted variables in space####

