####iNPUT otus###
require(ggplot2)
require(phyloseq)
###Import BIOM table data###
OTUS<-import_biom(BIOMfilename = 'Y:/VPino_hpc/ubuntu_qiime/16S-NSW/8DIV1/table_mc5000_JSON.biom',
                  treefilename = 'Y:/VPino_hpc/ubuntu_qiime/16S-NSW/8DIV1/rep_set.tre')

# colnames(data.frame(OTUS@otu_table@.Data))

####Correct sample data from mapfiles####
Sample_data <- read.table('Y:/VPino_hpc/ubuntu_qiime/16S-NSW/8DIV1/merged_mapping_NSW_WITH_SOIL_INFO.txt',sep='\t',header=T,stringsAsFactors = F)

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

# saveRDS(OTUS,file = 'RData/BactNSW_OTU.RDS')
####Unifrac distance mae in the Cluster####
OTUS<-readRDS(file = 'RData/BactNSW_OTU.RDS')
####

input_table_tmp<-read.table('Y:/VPino_hpc/ubuntu_qiime/16S-NSW/8DIV1/bdiv_even5000/bray_curtis_pc.txt',blank.lines.skip = T,fill = T)

Scores <-data.frame(PC1=input_table_tmp[seq(7,(nrow(input_table_tmp)-2),by = 2),2],
                    PC2=input_table_tmp[seq(7,(nrow(input_table_tmp)-2),by = 2),3],
                    PC3=input_table_tmp[seq(7,(nrow(input_table_tmp)-2),by = 2),4])
require(naturalsort)

###
Sample_data <- read.table('Y:/VPino_hpc/ubuntu_qiime/16S-NSW/8DIV1/merged_mapping_NSW_WITH_SOIL_INFO.txt',sep='\t',header=T,stringsAsFactors = F)

Sample_data_tmp<-Sample_data[match(colnames(OTUS@otu_table@.Data),Sample_data$SampleID),]
rownames(Sample_data_tmp)<-Sample_data_tmp$SampleID

###


Scores$ID <-as.character(input_table_tmp[seq(7,(nrow(input_table_tmp)-2),by = 2),1])
Scores<-Scores[match(colnames(OTUS@otu_table@.Data),Scores$ID),]

Scores <- as.matrix(Scores[,-4])
Scores <-apply(Scores,2,as.numeric)

Scores<-data.frame(Scores,Sample_data_tmp)

Scores <- as.data.frame(Scores)


# ####With  my scores ####
# dist_mat <- as.dist(read.csv('Y:/VPino_hpc/ubuntu_qiime/16S/divbase10/bdiv_even10000/weighted_unifrac_dm.txt',sep='\t',header=T)[,-1])
# 
# pr_comp <- prcomp(dist_mat)
# 
# # View(pr_comp$x)
# 
# 
# # plot(pr_comp$x[,'PC1'],pr_comp$x[,'PC2'])
# 
# Scores <- pr_comp$x[,1:2]
# 
# Scores<-Scores[match(paste0('X',colnames(OTUS@otu_table@.Data)),rownames(Scores)),]
# Scores<-data.frame(Scores,Sample_data_tmp)
colnames(Scores)

require(ggplot2)
ggplot(Scores,aes(PC1,PC2,colour=Total.Nitrogen))+
  geom_point()+
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

#####ploting confidence ellipses####
scores_bugs<-ggplot(Scores, aes(x=PC1, y=PC2))+
  geom_point(alpha=.6,size=6,aes(colour=System)) +
  labs(x='Var. Exp = 24%',y='Var. Exp = 13%',title='Cumulative explanation : 37%')+
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



PC1_gradient<-ggplot(Scores,aes(Site,PC1,group=System,colour=System))+
  geom_line(alpha=0.0)+
  geom_smooth(size=1.1)+
  labs(x='Site')
PC1_gradient

#####Canonical Correlation Analysis####

OTUS <- readRDS(file = 'RData/BactNSW_OTU.RDS')


# load('RData/TRFLP_Height_bact_fung_s0_s13_5cm.RData')

# names(TRFLP_Height_bact_fung_s0_s13_5cm) <- c('Bact_F1.2','Fung_F1.3')
# 
# require(vegan)
# require(ggplot2)
# 
# test<-list()
# test$details<-TRFLP_Height_bact_fung_s0_s13_5cm[['Bact_F1.2']]$details
# test$data <- TRFLP_Height_bact_fung_s0_s13_5cm[['Bact_F1.2']]$data
# 
# 
# ####average replicates ####
# sum_data <- do.call(rbind,by(test$data,test$details$Site,function(x) ceiling(rbind(colSums(x[1:3,],na.rm = T),colSums(x[4:6,],na.rm = T)))) )
# sum_details <- data.frame(sites =rep(seq(0,13),each=2), system= rep(c('Natural','Crop'),14))
# 
# sum_details$system <- gsub('Natural','Nat',sum_details$system)
# sum_details$index<-toupper(with(sum_details,paste(sites,system,'0_5_R1',sep='_')))
# 
# sum_data <- sum_data[!sum_details$index%in%c('12_CROP_0_5_R1','12_NAT_0_5_R1'),]
# sum_details <- sum_details[!sum_details$index%in%c('12_CROP_0_5_R1','12_NAT_0_5_R1'),]
# 

# soil <-readRDS('../transect_n_s/RData/NSW_CHEMICAL_with_lab_values.RDS')
# soil <- soil$responses[soil$responses$index%in%sum_details$index,]
# soil <- soil[!soil$index_coord=='12CROP',]
# 
# 
# soil<-apply(soil[,9:23],2,as.numeric)
# 
# 
# require(vegan)
# 
# 
# sum_data<-sum_data[,!apply(sum_data==0,2,any)]
require(phyloseq)
require(vegan)
soil <- as(sample_data(OTUS), "data.frame")

bugs <- prune_taxa(taxa_sums(OTUS) > 2500, OTUS)

bugs <- as.data.frame(t(as(otu_table(bugs), "matrix")))

colnames(soil)
soil[,c(17:31,33,51,53,70,68,79,80)] <- apply(soil[,c(17:31,33,51,53,70,68,79,80)],2,function(x) as.numeric(x)) 
soil_info <- soil
soil <- soil[,c(17:31,33,51,53,70,68,79,80)]

#CCA calculated in the HPC#
cca_16s<-CCorA(bugs,soil,stand.X = T)
cca_16s
# cca_16s$p.Pillai
# 
# ####Check validity of canonical analysis####
# 
library(CCP)
## Define number of observations, number of dependent variables, number of independent variables.
N <-  dim(soil)[1]
p <-  dim(soil)[2]
q <-  dim(bugs)[2]

rho <- cancor(soil,bugs,xcenter = T)$cor

p.asym(rho, N, p, q, tstat = "Wilks")
p.asym(rho, N, p, q, tstat = "Hotelling")
p.asym(rho, N, p, q, tstat = "Pillai")
p.asym(rho, N, p, q, tstat = "Roy")

res1 <- p.asym(rho, N, p, q)
plt.asym(res1,rhostart=7)
plt.asym(res1,rhostart=2)
plt.asym(res1,rhostart=18)

# 
# 
# summary(cca_16s)
# # readRDS('RData\CCA_16s.RDS')
# 
biplot(cca_16s,xlabs=rownames(bugs))


cca_manual_loading_soil <-data.frame(soil_loadx=cca_16s$corr.X.Cx[,1],soil_loady=cca_16s$corr.X.Cx[,2])
cca_manual_loading_bacteria <-data.frame(bact_loadx=cca_16s$corr.Y.Cy[,1],bact_loady=cca_16s$corr.Y.Cy[,2])


cca_manual_scores_soil <- data.frame(Can_Axis_1=cca_16s$Cx[,1],Can_Axis_2=cca_16s$Cx[,2])
cca_manual_scores_bact <- data.frame(Can_Axis_1=cca_16s$Cy[,1],Can_Axis_2=cca_16s$Cy[,2])

cca_manual_scores_soil<-data.frame(cca_manual_scores_soil,soil_info)

rownames(cca_16s$corr.X.Cx)<-c('Am_N','Nit_N','P','K','EC','pH_CaCl','pH_H2O','Ex_Al','Ex_Ca',
                                   'Ex_Mg','Ex_K','Ex_Na','N','TC','ECEC','site','SI_a','Clay',
                               'Norm_pedodiv','sd_pedodiv','sd_hor','first_sd_hor')


rownames(cca_16s$corr.Y.Cx) <- paste0('OTU_',regmatches(rownames(cca_16s$corr.Y.Cx),regexpr('([0-9])+$',perl = T,rownames(cca_16s$corr.Y.Cx))))

require(grid)
require(ggplot2)
soil_bact_CCA <-ggplot(cca_manual_scores_soil)+
  geom_point(aes(Can_Axis_1,Can_Axis_2,colour=system),size=8)+
  ggtitle(label = paste0('First two Canonical Axis for for Soil-bacteria matrix ',' Pillai trace signif. = ',cca_16s$p.Pillai))+  
  geom_segment(data = cca_manual_loading_soil,aes(x=0,y=0,
                                                  xend=soil_loadx*3,
                                                  yend=soil_loady*3),
               arrow=arrow(length=unit(0.2,"cm")))+
  geom_text(data = cca_manual_loading_soil,aes(x=soil_loadx*3,
                                               y=soil_loady*3,
                                               label=rownames(cca_16s$corr.X.Cx)),
            size=8)+
  geom_text(aes(x=Can_Axis_1,y=Can_Axis_2,label=site,vjust=-1),size=8)+
  theme_bw()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20,vjust=2),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.text=element_text(size=25),
        title=element_text(size=25))
soil_bact_CCA



cca_manual_scores_bact<-data.frame(cca_manual_scores_bact,soil_info)

bact_soil_CCA <-ggplot(cca_manual_scores_bact)+
  geom_point(aes(Can_Axis_1,Can_Axis_2,colour=system),size=8)+
  ggtitle(label = paste0('First two Canonical Axis for Bacteria-soil matrix ',' Pillai trace signif. = ',cca_16s$p.Pillai))+  
  geom_segment(data = cca_manual_loading_bacteria,aes(x=0,y=0,
                                                      xend=bact_loadx*3,
                                                      yend=bact_loady*3),
               arrow=arrow(length=unit(0.2,"cm")))+
  geom_text(data = cca_manual_loading_bacteria,aes(x=bact_loadx*3,
                                                   y=bact_loady*3,
                                                   label=rownames(cca_16s$corr.Y.Cx)),
            size=8)+
  geom_text(aes(x=Can_Axis_1,y=Can_Axis_2,label=site,vjust=-1),size=8)+
  theme_bw()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20,vjust=2),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.text=element_text(size=25),
        title=element_text(size=25))
bact_soil_CCA


mix_bact_soil_CCA <-ggplot(cca_manual_scores_bact)+
  geom_point(aes(Can_Axis_1,Can_Axis_2,colour=system),size=8,alpha=0)+
  ggtitle(label = paste0('First two Canonical Axis for Bact-soil matrix ',' Pillai trace signif. = ',round(cca_16s$p.Pillai,5)))+  
  geom_segment(data = cca_manual_loading_bacteria,aes(x=0,y=0,
                                                      xend=bact_loadx*4,
                                                      yend=bact_loady*4),
               ,colour='blue',
               size=1,
               arrow=arrow(length=unit(0.2,"cm")))+
  geom_segment(data = cca_manual_loading_soil,aes(x=0,y=0,
                                                  xend=soil_loadx*6,
                                                  yend=soil_loady*6)
               ,colour='red',
               size=1,
               arrow=arrow(length=unit(0.2,"cm")))+
  geom_text(data = cca_manual_loading_bacteria,aes(x=bact_loadx*4,
                                                   y=bact_loady*4,
                                                   label=rownames(cca_16s$corr.Y.Cx)),
            size=7)+
  geom_text(data = cca_manual_loading_soil,aes(x=soil_loadx*6,
                                               y=soil_loady*6,
                                               label=rownames(cca_16s$corr.X.Cx),
                                               vjust=-1),
            size=7)+
  geom_text(aes(x=Can_Axis_1,y=Can_Axis_2,label=site),size=5,alpha=0)+
  theme_bw()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20,vjust=2),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.text=element_text(size=25),
        title=element_text(size=25))

mix_bact_soil_CCA
#New.0.ReferenceOTU1488
#####Now for the rare OTUS#####
require(phyloseq)
require(vegan)
soil <- as(sample_data(OTUS), "data.frame")

bugs <- prune_taxa(taxa_sums(OTUS) < 5, OTUS)

bugs <- as.data.frame(t(as(otu_table(bugs), "matrix")))

colnames(soil)
soil[,c(17:31,33,51,53,70,68,79,80)] <- apply(soil[,c(17:31,33,51,53,70,68,79,80)],2,function(x) as.numeric(x)) 
soil_info <- soil
soil <- soil[,c(17:31,33,51,53,70,68,79,80)]

#CCA calculated in the HPC#
cca_16s<-CCorA(bugs,soil,stand.X = T)
cca_16s
# cca_16s$p.Pillai
# 
# ####Check validity of canonical analysis####
# 
library(CCP)
## Define number of observations, number of dependent variables, number of independent variables.
N <-  dim(soil)[1]
p <-  dim(soil)[2]
q <-  dim(bugs)[2]

rho <- cancor(soil,bugs,xcenter = T)$cor

p.asym(rho, N, p, q, tstat = "Wilks")
p.asym(rho, N, p, q, tstat = "Hotelling")
p.asym(rho, N, p, q, tstat = "Pillai")
p.asym(rho, N, p, q, tstat = "Roy")

res1 <- p.asym(rho, N, p, q)
plt.asym(res1,rhostart=7)
plt.asym(res1,rhostart=2)
plt.asym(res1,rhostart=18)

# 
# 
# summary(cca_16s)
# # readRDS('RData\CCA_16s.RDS')
# 
biplot(cca_16s,xlabs=rownames(bugs))


cca_manual_loading_soil <-data.frame(soil_loadx=cca_16s$corr.X.Cx[,1],soil_loady=cca_16s$corr.X.Cx[,2])
cca_manual_loading_bacteria <-data.frame(bact_loadx=cca_16s$corr.Y.Cy[,1],bact_loady=cca_16s$corr.Y.Cy[,2])


cca_manual_scores_soil <- data.frame(Can_Axis_1=cca_16s$Cx[,1],Can_Axis_2=cca_16s$Cx[,2])
cca_manual_scores_bact <- data.frame(Can_Axis_1=cca_16s$Cy[,1],Can_Axis_2=cca_16s$Cy[,2])

cca_manual_scores_soil<-data.frame(cca_manual_scores_soil,soil_info)

rownames(cca_16s$corr.X.Cx)<-c('Am_N','Nit_N','P','K','EC','pH_CaCl','pH_H2O','Ex_Al','Ex_Ca',
                               'Ex_Mg','Ex_K','Ex_Na','N','TC','ECEC','site','SI_a','Clay',
                               'Norm_pedodiv','sd_pedodiv','sd_hor','first_sd_hor')


rownames(cca_16s$corr.Y.Cx) <- paste0('OTU_',regmatches(rownames(cca_16s$corr.Y.Cx),regexpr('([0-9])+$',perl = T,rownames(cca_16s$corr.Y.Cx))))

require(grid)
require(ggplot2)
soil_bact_CCA <-ggplot(cca_manual_scores_soil)+
  geom_point(aes(Can_Axis_1,Can_Axis_2,colour=system),size=8)+
  ggtitle(label = paste0('First two Canonical Axis for for Soil-bacteria matrix ',' Pillai trace signif. = ',cca_16s$p.Pillai))+  
  geom_segment(data = cca_manual_loading_soil,aes(x=0,y=0,
                                                  xend=soil_loadx*3,
                                                  yend=soil_loady*3),
               arrow=arrow(length=unit(0.2,"cm")))+
  geom_text(data = cca_manual_loading_soil,aes(x=soil_loadx*3,
                                               y=soil_loady*3,
                                               label=rownames(cca_16s$corr.X.Cx)),
            size=8)+
  geom_text(aes(x=Can_Axis_1,y=Can_Axis_2,label=site,vjust=-1),size=8)+
  theme_bw()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20,vjust=2),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.text=element_text(size=25),
        title=element_text(size=25))
soil_bact_CCA



cca_manual_scores_bact<-data.frame(cca_manual_scores_bact,soil_info)

bact_soil_CCA <-ggplot(cca_manual_scores_bact)+
  geom_point(aes(Can_Axis_1,Can_Axis_2,colour=system),size=8)+
  ggtitle(label = paste0('First two Canonical Axis for Bacteria-soil matrix ',' Pillai trace signif. = ',cca_16s$p.Pillai))+  
  geom_segment(data = cca_manual_loading_bacteria,aes(x=0,y=0,
                                                      xend=bact_loadx*3,
                                                      yend=bact_loady*3),
               arrow=arrow(length=unit(0.2,"cm")))+
  geom_text(data = cca_manual_loading_bacteria,aes(x=bact_loadx*3,
                                                   y=bact_loady*3,
                                                   label=rownames(cca_16s$corr.Y.Cx)),
            size=8)+
  geom_text(aes(x=Can_Axis_1,y=Can_Axis_2,label=site,vjust=-1),size=8)+
  theme_bw()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20,vjust=2),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.text=element_text(size=25),
        title=element_text(size=25))
bact_soil_CCA


mix_bact_soil_CCA <-ggplot(cca_manual_scores_bact)+
  geom_point(aes(Can_Axis_1,Can_Axis_2,colour=system),size=8,alpha=0)+
  ggtitle(label = paste0('First two Canonical Axis for Bact-soil matrix ',' Pillai trace signif. = ',round(cca_16s$p.Pillai,5)))+  
  geom_segment(data = cca_manual_loading_bacteria,aes(x=0,y=0,
                                                      xend=bact_loadx*4,
                                                      yend=bact_loady*4),
               ,colour='blue',
               size=1,
               arrow=arrow(length=unit(0.2,"cm")))+
  geom_segment(data = cca_manual_loading_soil,aes(x=0,y=0,
                                                  xend=soil_loadx*6,
                                                  yend=soil_loady*6)
               ,colour='red',
               size=1,
               arrow=arrow(length=unit(0.2,"cm")))+
  geom_text(data = cca_manual_loading_bacteria,aes(x=bact_loadx*4,
                                                   y=bact_loady*4,
                                                   label=rownames(cca_16s$corr.Y.Cx)),
            size=7)+
  geom_text(data = cca_manual_loading_soil,aes(x=soil_loadx*6,
                                               y=soil_loady*6,
                                               label=rownames(cca_16s$corr.X.Cx),
                                               vjust=-1),
            size=7)+
  geom_text(aes(x=Can_Axis_1,y=Can_Axis_2,label=site),size=5,alpha=0)+
  theme_bw()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20,vjust=2),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.text=element_text(size=25),
        title=element_text(size=25))

mix_bact_soil_CCA



OTUS_check <-prune_taxa(taxa_sums(OTUS) < 25, OTUS)

plot(phy_tree(OTUS_check))

