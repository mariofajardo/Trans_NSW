
# input_table_tmp<-read.table('DATA/pear_joinedNW1_corediv/bdiv_even678/weighted_unifrac_pc.txt',blank.lines.skip = T,fill = T)
input_table_tmp<-read.table('Y:/VPino_hpc/ubuntu_qiime/ITS/Pipit_ITSvpNS/core_output/bdiv_even640/bray_curtis_pc.txt',blank.lines.skip = T,fill = T)


Scores <-data.frame(PC1=input_table_tmp[seq(7,165,by = 2),2],
                    PC2=input_table_tmp[seq(7,165,by = 2),3],
                    PC3=input_table_tmp[seq(7,165,by = 2),4])
require(naturalsort)

Scores$ID <-as.character(input_table_tmp[seq(7,165,by = 2),1])
Scores<-Scores[naturalorder(Scores$ID),]

Scores <- as.matrix(Scores[,-4])
Scores <-apply(Scores,2,as.numeric)

Scores<-data.frame(Scores,Sample_data_tmp[naturalorder(Sample_data_tmp$SampleID.x),])

Scores <- as.data.frame(Scores)

ggplot(Scores,aes(PC1,PC2,colour=Site))+
  geom_point()+
  geom_text(data = Scores,aes(x=PC1,y=PC2,label=Site,vjust=1))


conf_ellipse <-data.frame()
conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(Scores[Scores$Site%in%0:5,], 
                                                             ellipse(-cor(PC1,PC2), 
                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                     centre=c(mean(-PC1),mean(PC2)),
                                                                     level=.90))),Site='Sites 0 to 5'))
conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(Scores[Scores$Site%in%6:13,], 
                                                             ellipse(-cor(PC1,PC2), 
                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                     centre=c(mean(-PC1),mean(PC2)),
                                                                     level=.90))),Site='Sites 6 to 13'))

#####ploting confidence ellipses####
scores_bugs<-ggplot(Scores, aes(x=-PC1, y=-PC2))+
  geom_point(alpha=.6,size=6,aes(colour=System)) +
  labs(x='Var. Exp = 36%',y='Var. Exp = 9%',title='Cumulative explanation : 45%')+
  geom_path(data=conf_ellipse,aes(x=x, y=y,colour=Site),size=1,linetype=1)+
  geom_text(aes(x=-PC1, y=-PC2,label=Site,colour=System),vjust=2,size=5)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        title=element_text(size=30),
        legend.text=element_text(size = 20))
scores_bugs

# ggsave('../../../../../../../Dropbox/Mario-Vane/PCa_bichos.png',scores_bugs)

###Gradients of colors####

####pH####
Scores$pH.Level..H2O.<-as.numeric(Scores$pH.Level..H2O.)
colnames(Scores)[19] <- 'pH'
colnames(Scores)[19] <- 'pH'

pH_bugs<-ggplot(Scores, aes(x=-PC1, y=PC2,colour=pH))+
  scale_color_gradient(low = 'grey',high = 'black')+
  geom_point(size=6)+
  labs(x='Var. Exp = 36%',y='Var. Exp = 9%',title='pH gradient')+
  geom_text(aes(x=-PC1, y=PC2,label=Site),vjust=2,size=5)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        title=element_text(size=30),
        legend.text=element_text(size = 20))
pH_bugs

# ggsave('../../../../../../../Dropbox/Mario-Vane/pH_bugs.png',pH_bugs)

#####Magnessium####
Scores$Exc..Magnesium<-as.numeric(Scores$Exc..Magnesium)
colnames(Scores)[22] <- 'Exch.Mg'

Exch.Mg_bugs<-ggplot(Scores, aes(x=-PC1, y=PC2,colour=Exch.Mg))+
  scale_color_gradient(low = 'grey',high = 'black')+
  geom_point(size=6)+
  labs(x='Var. Exp = 36%',y='Var. Exp = 9%',title='Exchangeable Mg gradient (meq/100g)')+
  geom_text(aes(x=-PC1, y=PC2,label=Site),vjust=2,size=5)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        title=element_text(size=30),
        legend.text=element_text(size = 20))
Exch.Mg_bugs

# ggsave('../../../../../../../Dropbox/Mario-Vane/Exch.Mg.png',Exch.Mg_bugs)

#####Electric conductivity####
Scores$Conductivity<-as.numeric(Scores$Conductivity)

Conductivity_bugs<-ggplot(Scores, aes(x=-PC1, y=PC2,colour=Conductivity))+
  scale_color_gradient(low = 'grey',high = 'black')+
  geom_point(size=6)+
  labs(x='Var. Exp = 36%',y='Var. Exp = 9%',title='Electric conductivity gradient (dS/m)')+
  geom_text(aes(x=-PC1, y=PC2,label=Site),vjust=2,size=5)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        title=element_text(size=30),
        legend.text=element_text(size = 20))
Conductivity_bugs

# ggsave('../../../../../../../Dropbox/Mario-Vane/Conductivity.png',Conductivity_bugs)

#####ECEC####
Scores$ECEC<-as.numeric(Scores$ECEC)

ECEC_bugs<-ggplot(Scores, aes(x=-PC1, y=PC2,colour=ECEC))+
  scale_color_gradient(low = 'grey',high = 'black')+
  geom_point(size=6)+
  labs(x='Var. Exp = 36%',y='Var. Exp = 9%',title='ECEC gradient (meq/100g)')+
  geom_text(aes(x=-PC1, y=PC2,label=Site),vjust=2,size=5)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        title=element_text(size=30),
        legend.text=element_text(size = 20))
ECEC_bugs

# ggsave('../../../../../../../Dropbox/Mario-Vane/ECEC.png',ECEC_bugs)


#####Sites####
Scores$Site<-as.numeric(Scores$Site)

Site_bugs<-ggplot(Scores, aes(x=-PC1, y=PC2,colour=Site))+
  scale_color_gradient(low = 'grey',high = 'black',limits=c(0, 12))+
  geom_point(size=6)+
  labs(x='Var. Exp = 36%',y='Var. Exp = 9%',title='Sites')+
  geom_text(aes(x=-PC1, y=PC2,label=Site),vjust=2,size=5)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        title=element_text(size=30),
        legend.text=element_text(size = 20))
Site_bugs

# ggsave('../../../../../../../Dropbox/Mario-Vane/Site_bugs.png',Site_bugs)



PC1_gradient<-ggplot(Scores,aes(Site*50,PC1,group=System,colour=System))+
  geom_line(alpha=0.0)+
  geom_smooth(size=1.1)+
  labs(x='Kilometers')

PC1_gradient


# ggsave('../../../../../../../Dropbox/Mario-Vane/PC1_bichos.png',PC1_gradient)

#####nOW WITH SOIL####
NSW_NS <- readRDS('RData/NSW_CHEMICAL.rds')$NS
load('pr_varExp.RData')
no_out_pr_Exp<-pr_varExp(NSW_NS$spectra[NSW_NS$responses$Sample%in%0:13,])

scores <-data.frame(Sample=NSW_NS$responses$Sample[NSW_NS$responses$Sample%in%0:13],
                    Transect='N_S',
                    system=NSW_NS$responses$system[NSW_NS$responses$Sample%in%0:13],
                    top=NSW_NS$responses$top[NSW_NS$responses$Sample%in%0:13],
                    bottom=NSW_NS$responses$bottom[NSW_NS$responses$Sample%in%0:13],
                    prcomp(NSW_NS$spectra[NSW_NS$responses$Sample%in%0:13,],center = T,scale. = T)$x[,1:10])



scores_NSW <- ggplot(scores,aes(PC1,PC2))+
  geom_point(aes(colour=system),size=5)+
  theme_bw()+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))+
  geom_text(aes(PC1,PC2,label=Sample),vjust=2)
scores_NSW

require(ellipse)
conf_ellipse_soil <-data.frame()
conf_ellipse_soil <- rbind(conf_ellipse_soil, cbind(as.data.frame(with(scores[scores$Sample%in%0:5,], 
                                                                       ellipse(cor(PC1,PC2), 
                                                                               scale=c(sd(PC1),sd(PC2)), 
                                                                               centre=c(mean(PC1),mean(PC2)),
                                                                               level=.90))),Sample='Sites 0 to 5'))
conf_ellipse_soil <- rbind(conf_ellipse_soil, cbind(as.data.frame(with(scores[scores$Sample%in%6:13,], 
                                                                       ellipse(cor(PC1,PC2), 
                                                                               scale=c(sd(PC1),sd(PC2)), 
                                                                               centre=c(mean(PC1),mean(PC2)),
                                                                               level=.90))),Sample='Sites 6 to 13'))


#####plot confidence ellipses####
scores_poster<-ggplot(scores, aes(x=PC1, y=PC2,colour=system)) + 
  geom_point(alpha=.6,size=7) +
  labs(x=paste0('PC1 :',round(no_out_pr_Exp[1,1],2),'%'), 
       y=paste0('PC2 :',round(no_out_pr_Exp[2,1],2),'%'),
       title=paste0('Cumulative explanation :',
                    round(no_out_pr_Exp[2,2],2),'%')) +
  geom_path(data=conf_ellipse_soil, aes_string(x='x', y='y',colour='Sample'), size=1, linetype=1)+
  geom_text(aes(x=PC1, y=PC2,label=scores$Sample),vjust=2,size=5)+
  theme(strip.text=element_text(size=20,colour = 'white'),
        axis.text.x=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        title=element_text(size=30),
        legend.text=element_text(size = 20))


scores_poster

ggsave('../../../../../../../Dropbox/Mario-Vane/PCa_soil.png',scores_poster)


PC1_gradient_soil<-ggplot(scores,aes(as.numeric(Sample)*50,-PC1,group=system,colour=system))+
  geom_line(alpha=0.0)+
  geom_smooth(size=1.1)+
  labs(x='Kilometers')

ggsave('../../../../../../../Dropbox/Mario-Vane/PC1_soil.png',PC1_gradient_soil)


####take that sample 13###
Scores1<-Scores
Scores1 <- Scores1[-(79:80),]

index_spectra <-paste(NSW_NS$responses$Sample,NSW_NS$responses$system,NSW_NS$responses$top,NSW_NS$responses$bottom,'R1',sep='_') 
scores <-data.frame(Sample=NSW_NS$responses$Sample[index_spectra%in%Scores$index],
                    Transect='N_S',
                    system=NSW_NS$responses$system[index_spectra%in%Scores$index][NSW_NS$responses$Sample[index_spectra%in%Scores$index]%in%0:13],
                    top=NSW_NS$responses$top[index_spectra%in%Scores$index][NSW_NS$responses$Sample[index_spectra%in%Scores$index]%in%0:13],
                    bottom=NSW_NS$responses$bottom[index_spectra%in%Scores$index][NSW_NS$responses$Sample[index_spectra%in%Scores$index]%in%0:13],
                    prcomp(NSW_NS$spectra[index_spectra%in%Scores$index,][NSW_NS$responses$Sample[index_spectra%in%Scores$index]%in%0:13,],center = T,scale. = T)$x[,1:10])




scores<-scores[-c(53:54),]
scores<-scores[seq(1,52,2),]
require(plyr)
scores <-adply(scores,1,function(x) rbind(x,x,x))
require(naturalsort)

testing <- data.frame(FirstPC_soil=scores$PC1,FirstPC_bugs=Scores1$PC1)

testing$FirstPC_soil <- testing$FirstPC_soil[naturalorder(paste(scores$Sample,scores$system,sep='_'))]
testing$FirstPC_bugs <- testing$FirstPC_bugs[naturalorder(paste(Scores1$Site,Scores1$system,sep='_'))]



testing<-scale(testing)
testing <- as.data.frame(testing)



require(ggplot2)
ggplot(testing,aes(FirstPC_soil,FirstPC_bugs))+
  geom_point(size=2)+
  geom_abline(slope=1)+
  stat_abline()