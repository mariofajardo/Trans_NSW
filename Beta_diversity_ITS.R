####iNPUT otus###

require(phyloseq)
###Import BIOM table data###
OTUS<-import_biom(BIOMfilename = 'Y:/VPino_hpc/ubuntu_qiime/ITS-NSW/VP2_VP3/VP2_VP3/core_output_e1111/table_even1111_JSON.biom')

colnames(data.frame(OTUS@otu_table@.Data))

####Correct sample data from mapfiles####
Sample_data <- read.table('Y:/VPino_hpc/ubuntu_qiime/ITS-NSW/NSW_ITS_metadatfile_for_R.txt',sep='\t',header=T,stringsAsFactors = F)

Sample_data_tmp<-Sample_data[match(colnames(OTUS@otu_table@.Data),Sample_data$SampleID),]
rownames(Sample_data_tmp)<-Sample_data_tmp$SampleID

####How is the structure of sample data ????####
# data(GlobalPatterns) #example
# str(GlobalPatterns@sam_data)
#this command is from phyloseq package

Sample_data <-sample_data(Sample_data_tmp)
OTUS<-merge_phyloseq(OTUS,Sample_data)

colnames(tax_table(OTUS)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                               o = "Order", f = "Family", g = "Genus", s = "Species")
rank_names(OTUS)

sum(taxa_sums(OTUS))

####

input_table_tmp<-read.table('Y:/VPino_hpc/ubuntu_qiime/ITS-NSW/VP2_VP3/VP2_VP3/core_output_e1111/bdiv_even1111/bray_curtis_pc.txt',blank.lines.skip = T,fill = T)

Scores <-data.frame(PC1=input_table_tmp[seq(7,(nrow(input_table_tmp)-2),by = 2),2],
                    PC2=input_table_tmp[seq(7,(nrow(input_table_tmp)-2),by = 2),3],
                    PC3=input_table_tmp[seq(7,(nrow(input_table_tmp)-2),by = 2),4])

####And with ordination function####
####Ordination####
# 
# components <- ordinate(OTUS,'PCoA')
# Scores <- data.frame(PC1=components$vectors[,1],PC2=components$vectors[,2])
# 

####


require(naturalsort)

Scores$ID <-as.character(input_table_tmp[seq(7,(nrow(input_table_tmp)-2),by = 2),1])
Scores<-Scores[match(colnames(OTUS@otu_table@.Data),Scores$ID),]

Scores <- as.matrix(Scores[,-4])
Scores <-apply(Scores,2,as.numeric)

Scores<-data.frame(Scores,Sample_data_tmp)

Scores <- as.data.frame(Scores)

colnames(Scores)

require(ggplot2)
ggplot(Scores,aes(PC1,PC2,colour=System))+
  geom_point(size=3)+
  geom_text(data = Scores,aes(x=PC1,y=PC2,label=Site,vjust=1))

require(ellipse)
conf_ellipse <-data.frame()
conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(Scores[Scores$Site%in%0:5,], 
                                                             ellipse(cor(PC1,PC2), 
                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                     centre=c(mean(PC1),mean(PC2)),
                                                                     level=.90))),Site='Sites 0 to 5'))

conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(Scores[Scores$Site%in%6:13,], 
                                                             ellipse(cor(PC1,PC2), 
                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                     centre=c(mean(PC1),mean(PC2)),
                                                                     level=.90))),Site='Sites 6 to 13'))


conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(Scores[Scores$System%in%'CROP',], 
                                                             ellipse(cor(PC1,PC2), 
                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                     centre=c(mean(PC1),mean(PC2)),
                                                                     level=.90))),Site='System == Crop'))

conf_ellipse <- rbind(conf_ellipse, cbind(as.data.frame(with(Scores[Scores$System%in%'NAT',], 
                                                             ellipse(cor(PC1,PC2), 
                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                     centre=c(mean(PC1),mean(PC2)),
                                                                     level=.90))),Site='System == NAT'))


#####ploting confidence ellipses####
scores_bugs<-ggplot(Scores, aes(x=PC1, y=PC2))+
  geom_point(alpha=.6,size=6,aes(colour=System)) +
  labs(x='Var. Exp = 5%',y='Var. Exp = 4.2%',title='Cumulative explanation : 9.2%')+
  geom_path(data=conf_ellipse,aes(x=x, y=y,colour=Site),size=1,linetype=1)+
  geom_text(aes(x=PC1, y=PC2,label=Site,colour=System),vjust=2,size=5)+
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
colnames(Scores)[which(colnames(Scores)=='pH.Level..H2O.')] <- 'pH'

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
colnames(Scores)[which(colnames(Scores)=='Exc..Magnesium')] <- 'Exch.Mg'

Exch.Mg_bugs<-ggplot(Scores, aes(x=-PC1, y=PC2,colour=Exch.Mg))+
  scale_color_gradient(low = 'grey',high = 'black')+
  geom_point(size=6)+
  labs(x='Var. Exp = 5%',y='Var. Exp = 4.2%',title='Exchangeable Mg gradient (meq/100g)')+
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
  labs(x='Var. Exp = 5%',y='Var. Exp = 4.2%',title='Electric conductivity gradient (dS/m)')+
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
  labs(x='Var. Exp = 5%',y='Var. Exp = 4.2%',title='ECEC gradient (meq/100g)')+
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
  labs(x='Var. Exp = 5%',y='Var. Exp = 4.2%',title='Sites')+
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



PC1_gradient<-ggplot(Scores,aes(Km/50,PC1,group=System,colour=System))+
  geom_line(alpha=0.0)+
  geom_smooth(size=1.1)+
  labs(x='Site')+
  facet_wrap(~Transect)

PC1_gradient


# ggsave('../../../../../../../Dropbox/Mario-Vane/PC1_bichos.png',PC1_gradient)


