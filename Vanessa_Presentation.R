####OTUS Bacteria####
require(vegan)
require(phyloseq)
require(ggplot2)

OTUS_16S<-import_biom(BIOMfilename = 'Y:/VPino_hpc/ubuntu_qiime/16S/divbase10/table_mc10000_JSON.biom')


Sample_data <- read.table('Y:/VPino_hpc/ubuntu_qiime/16S/16S_ecometadata_NSW_WITH_SOIL_INFO_for_R.txt',sep='\t',header=T,stringsAsFactors = F)

Sample_data_tmp<-Sample_data[match(colnames(OTUS_16S@otu_table@.Data),Sample_data$SampleID),]
rownames(Sample_data_tmp)<-Sample_data_tmp$SampleID

####How is the structure of sample data ????####
# data(GlobalPatterns) #example
# str(GlobalPatterns@sam_data)
#this command is from phyloseq package

Sample_data <-sample_data(Sample_data_tmp)
OTUS_16S<-merge_phyloseq(OTUS_16S,Sample_data)

colnames(tax_table(OTUS_16S)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                               o = "Order", f = "Family", g = "Genus", s = "Species")
rank_names(OTUS_16S)



####OTUS fungi####


###Import BIOM table data###
OTUS_ITS<-import_biom(BIOMfilename = 'Y:/VPino_hpc/ubuntu_qiime/ITS-NSW/VP2_VP3/VP2_VP3/core_output_e1111/table_mc1111_JSON.biom')

colnames(data.frame(OTUS_ITS@otu_table@.Data))

####Correct sample data from mapfiles####
Sample_data <- read.table('Y:/VPino_hpc/ubuntu_qiime/ITS-NSW/NSW_ITS_metadatfile_for_R.txt',sep='\t',header=T,stringsAsFactors = F)

Sample_data_tmp<-Sample_data[match(colnames(OTUS_ITS@otu_table@.Data),Sample_data$SampleID),]
rownames(Sample_data_tmp)<-Sample_data_tmp$SampleID

####How is the structure of sample data ????####
# data(GlobalPatterns) #example
# str(GlobalPatterns@sam_data)
#this command is from phyloseq package

Sample_data <-sample_data(Sample_data_tmp)
OTUS_ITS<-merge_phyloseq(OTUS_ITS,Sample_data)

colnames(tax_table(OTUS_ITS)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                               o = "Order", f = "Family", g = "Genus", s = "Species")
rank_names(OTUS_ITS)

#####



# OTUS <- subset_samples(OTUS,grepl('NS',sample_names(OTUS)))


soil16_S <- as(sample_data(OTUS_16S), "data.frame")
bugs16_S <- as.data.frame(t(as(otu_table(OTUS_16S), "matrix")))

soil16_S_info <- soil16_S
colnames(soil16_S)
soil16_S <- soil16_S[,c(17:31,51,53,70)]

soil16_S$C_N <- soil16_S$Total.Carbon/soil16_S$Total.Nitrogen
soil16_S$Ca_Mg <- soil16_S$Exc..Calcium/soil16_S$Exc..Magnesium
colnames(soil16_S)[18] <- 'Pedodiv'

soilITS <- as(sample_data(OTUS_ITS), "data.frame")
bugsITS <- as.data.frame(t(as(otu_table(OTUS_ITS), "matrix")))

soilITS_info <- soilITS
soilITS <- soilITS[,c(17:31,51,53,70)]


soilITS$C_N <- soilITS$Total.Carbon/soilITS$Total.Nitrogen
soilITS$Ca_Mg <- soilITS$Exc..Calcium/soilITS$Exc..Magnesium
colnames(soilITS)[18] <- 'Pedodiv'

#diversities 16S#
soil16_S$Obs_OTUS_16S <- rowSums(bugs16_S)
soil16_S$Shannon_div_16S <- diversity(bugs16_S,1,index = 'shannon')

# soil16_S$chao1_16S <- chao1(bugs16_S,taxa.row = F)
# soil16_S$inv_sim_div_16S <- diversity(bugs16_S,1,index = 'invsimpson')


#diversities ITS#
soilITS$Obs_OTUS_ITS <- rowSums(bugsITS)
soilITS$Shannon_div_ITS <- diversity(bugsITS,1,index = 'shannon')

# soilITS$chao1_ITS <- chao1(bugsITS,taxa.row = F)
# soilITS$inv_sim_div_ITS <- diversity(bugsITS,1,index = 'invsimpson')



soil <- soil16_S
###Check###
View(data.frame(rownames(soil),rownames(soilITS[match(soil16_S_info$index,soilITS_info$index),])))
    
soil[,c('Shannon_div_ITS',
        'Obs_OTUS_ITS')] <- soilITS[match(soil16_S_info$index,soilITS_info$index),c('Shannon_div_ITS',
                                                                                       'Obs_OTUS_ITS')]
soil_info <- soil16_S_info[complete.cases(soil),]
soil <- soil[complete.cases(soil),]

#####fOLLOWING ANALYSES WILL BE JUST FOR n-s TRANSECT####
soil <- soil[soil_info$Transect=='NS',]
soil_info <- soil_info[soil_info$Transect=='NS',]

components<-prcomp(soil,center = T,scale. = T)


check_PCA<-data.frame(Ecosystem=soil_info$system,
                      site=soil_info$site,
                      top=soil_info$top,
                      bottom=soil_info$bottom,
                      components$x[,1:2])
Loadings<-data.frame(Property=colnames(soil),
                     PC1_rot=components$rotation[,1],
                     PC2_rot=components$rotation[,2])


require(grid)
require(ggplot2)

Vanessas_theme <-  theme_bw()+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30,vjust=2),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        strip.text=element_text(size=25),
        title=element_text(size=25))


check_PCA <- cbind(check_PCA,soil_info)

colnames(check_PCA)[23:37]<-c('Am_N','Nit_N','P','K','EC','pH_CaCl','pH_H2O','Ex_Al','Ex_Ca',
                               'Ex_Mg','Ex_K','Ex_Na','N','TC','ECEC')


rownames(Loadings)[1:15] <- c('Am_N','Nit_N','P','K','EC','pH_CaCl','pH_H2O','Ex_Al','Ex_Ca',
                        'Ex_Mg','Ex_K','Ex_Na','N','TC','ECEC')

levels(Loadings$Property)[7:11] <-c('Ex_Al','Ex_Ca','Ex_Mg','Ex_K','Ex_Na')

PCA_Obs_OTUS<-ggplot(check_PCA)+
  geom_point(aes(PC1,PC2,colour=Ecosystem),size=4)+
  ggtitle(label = paste0('Alpha diversity (N-S) + Soil properties Var. exp = ',round(summary(components)$importance[6],2)))+
  geom_text(aes(PC1,PC2,label=site,vjust=2),size=12,alpha=0.5)+  
  geom_segment(data = Loadings,aes(x=0,y=0,xend=PC1_rot*components$sdev[1]*sqrt(nrow(check_PCA)),
                                   yend=PC2_rot*components$sdev[2]*sqrt(nrow(check_PCA))),
               arrow=arrow(length=unit(0.2,"cm")),
               linetype=2)+
  geom_text(data = Loadings,aes(x=PC1_rot*components$sdev[1]*sqrt(nrow(check_PCA)),
                                y=PC2_rot*components$sdev[2]*sqrt(nrow(check_PCA)),
                                label=Property),size=12)
 

PCA_Obs_OTUS <- PCA_Obs_OTUS+Vanessas_theme
PCA_Obs_OTUS




ggsave(filename = 'Y:/VPino_hpc/Plots/PCA_Obs_OTUS.png',width = 30,height = 16,plot = PCA_Obs_OTUS)
# 
# 
# 




soil <- cbind(soil,check_PCA)

Pedodiv_OTUS_16S <- ggplot(soil,aes(Pedodiv,Obs_OTUS_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Pedodiv_OTUS_16S <- Pedodiv_OTUS_16S + Vanessas_theme 

Pedodiv_OTUS_16S

ggsave(filename = 'Y:/VPino_hpc/Plots/Pedodiv_OTUS_16S_mix.png',plot = Pedodiv_OTUS_16S)


Pedodiv_OTUS_ITS <- ggplot(soil,aes(Pedodiv,Obs_OTUS_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Pedodiv_OTUS_ITS <- Pedodiv_OTUS_ITS + Vanessas_theme 

Pedodiv_OTUS_ITS

ggsave(filename = 'Y:/VPino_hpc/Plots/Pedodiv_OTUS_ITS_mix.png',plot = Pedodiv_OTUS_ITS)



Pedodiv_div_16S <- ggplot(soil,aes(Pedodiv,Shannon_div_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Pedodiv_div_16S <- Pedodiv_div_16S + Vanessas_theme 

Pedodiv_div_16S

ggsave(filename = 'Y:/VPino_hpc/Plots/Pedodiv_div_16S_mix.png',plot = Pedodiv_div_16S)


Pedodiv_div_ITS <- ggplot(soil,aes(Pedodiv,Shannon_div_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Pedodiv_div_ITS <- Pedodiv_div_ITS + Vanessas_theme 

Pedodiv_div_ITS

ggsave(filename = 'Y:/VPino_hpc/Plots/Pedodiv_div_ITS_mix.png',plot = Pedodiv_div_ITS)






pH_OBS_OTUS_16S <- ggplot(soil,aes(pH.Level..H2O.,Obs_OTUS_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

pH_OBS_OTUS_16S <- pH_OBS_OTUS_16S + Vanessas_theme 

pH_OBS_OTUS_16S

ggsave(filename = 'Y:/VPino_hpc/Plots/pH_OBS_OTUS_16S_mix.png',plot = pH_OBS_OTUS_16S)


pH_OBS_OTUS_ITS <- ggplot(soil,aes(pH.Level..H2O.,Obs_OTUS_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

pH_OBS_OTUS_ITS <- pH_OBS_OTUS_ITS + Vanessas_theme 

pH_OBS_OTUS_ITS

ggsave(filename = 'Y:/VPino_hpc/Plots/pH_OBS_OTUS_ITS_mix.png',plot = pH_OBS_OTUS_ITS)


# PC1_gradient<-ggplot(check_PCA,aes(Km,PC1,group=Ecosystem,colour=Ecosystem))+
#   geom_line(alpha=0.0)+
#   geom_smooth(size=1.1)+
#   labs(x='Kilometers')+
#   ggtitle('First PC variation')+
#   facet_wrap(~Transect)
# 
# PC1_gradient <- PC1_gradient+ Vanessas_theme
# PC1_gradient

# ggsave(filename = 'Y:/VPino_hpc/Plots/PC1_gradient.pdf',plot = PC1_gradient)

colnames(soil)


Slacking_Obs_OTUS16S <- ggplot(soil,aes(Slaking_coef_a,Obs_OTUS_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Slacking_Obs_OTUS16S <- Slacking_Obs_OTUS16S + Vanessas_theme 
Slacking_Obs_OTUS16S


ggsave(filename = 'Y:/VPino_hpc/Plots/Slacking_Obs_OTUS_16S_mix.png',plot = Slacking_Obs_OTUS16S)


Slacking_Obs_OTUS_ITS <- ggplot(soil,aes(Slaking_coef_a,Obs_OTUS_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Slacking_Obs_OTUS_ITS <- Slacking_Obs_OTUS_ITS + Vanessas_theme 
Slacking_Obs_OTUS_ITS

# t.test(soil$pH.Level..H2O.[soil$system=='Crop'],soil$pH.Level..H2O.[soil$system=='Nat'],'greater')

ggsave(filename = 'Y:/VPino_hpc/Plots/Slacking_Obs_OTUS_ITS_mix.png',plot = Slacking_Obs_OTUS_ITS)



Clay_Obs_OTUS16S <- ggplot(soil,aes(Clay,Obs_OTUS_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')
Clay_Obs_OTUS16S <- Clay_Obs_OTUS16S + Vanessas_theme
Clay_Obs_OTUS16S

ggsave(filename = 'Y:/VPino_hpc/Plots/Clay_Obs_OTUS16S_mix.png',plot = Clay_Obs_OTUS16S)


Clay_Obs_OTUSITS <- ggplot(soil,aes(Clay,Obs_OTUS_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')
Clay_Obs_OTUSITS <- Clay_Obs_OTUSITS + Vanessas_theme
Clay_Obs_OTUSITS

ggsave(filename = 'Y:/VPino_hpc/Plots/Clay_Obs_OTUSITS_mix.png',plot = Clay_Obs_OTUSITS)


Shannon16S <- ggplot(soil,aes(site*50,Shannon_div_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')+
  xlab('Km')
Shannon16S <- Shannon16S + Vanessas_theme
Shannon16S

ggsave(filename = 'Y:/VPino_hpc/Plots/Shannon16S_mix.png',plot = Shannon16S)


Obs_OTUS16S <- ggplot(soil,aes(site*50,Obs_OTUS_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')+
  xlab('Km')
Obs_OTUS16S <- Obs_OTUS16S + Vanessas_theme
Obs_OTUS16S

ggsave(filename = 'Y:/VPino_hpc/Plots/Obs_OTUS16S_mix.png',plot = Obs_OTUS16S)



ShannonITS <- ggplot(soil,aes(site*50,Shannon_div_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')+
  xlab('Km')
ShannonITS <- ShannonITS + Vanessas_theme
ShannonITS

ggsave(filename = 'Y:/VPino_hpc/Plots/ShannonITS_mix.png',plot = ShannonITS)


Obs_OTUSITS <- ggplot(soil,aes(site*50,Obs_OTUS_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')+
  xlab('Km')
Obs_OTUSITS <- Obs_OTUSITS + Vanessas_theme
Obs_OTUSITS

ggsave(filename = 'Y:/VPino_hpc/Plots/Obs_OTUSITS_mix.png',plot = Obs_OTUSITS)

Pedodiv <- ggplot(soil,aes(site*50,Pedodiv,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')+
  xlab('Km')
Pedodiv <- Pedodiv + Vanessas_theme
Pedodiv

ggsave(filename = 'Y:/VPino_hpc/Plots/Pedodiv_mix.png',plot = Pedodiv)



# ggplot(soil,aes(inv_sim_div_ITS,pH.Level..H2O.))+
#   geom_point()+
#   geom_smooth(size=1.1,method='lm')+
#   facet_wrap(~system)

# ggplot(soil,aes(inv_sim_div_ITS,C_N))+
#   geom_point()+
#   geom_smooth(size=1.1,method='lm')+
#   facet_wrap(~system)

# ggplot(soil,aes(inv_sim_div_ITS,Exc..Calcium))+
#   geom_point()+
#   geom_smooth(size=1.1,method='lm')+
#   facet_wrap(~system)

# ggplot(soil,aes(abundance_ITS,Exc..Aluminium))+
#   geom_point()+
#   geom_smooth(size=1.1,method='lm')+
#   facet_wrap(~system)
# 
# ggplot(soil,aes((Exc..Sodium/sqrt((Exc..Calcium+Exc..Magnesium)/2)),abundance_ITS))+
#   geom_point()+
#   geom_smooth(size=1.1,method='lm')+
#   facet_wrap(~system)+
#   xlab('sodium adsorbption ratio')

ggplot(soil,aes(Obs_OTUS_16S,Conductivity))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')+
  facet_wrap(~system)

ggplot(soil,aes(Obs_OTUS_16S,ECEC))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')+
  facet_wrap(~system)

# ggplot(soil,aes(inv_sim_div_16S,ECEC))+
#   geom_point()+
#   geom_smooth(size=1.1,method='lm')

Total_Nitrogen16S <- ggplot(soil,aes(Total.Nitrogen,Obs_OTUS_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Total_Nitrogen16S <- Total_Nitrogen16S + Vanessas_theme
Total_Nitrogen16S

ggsave(filename = 'Y:/VPino_hpc/Plots/Total_Nitrogen16S_mix.png',plot = Total_Nitrogen16S)


Total_NitrogenITS <- ggplot(soil,aes(Total.Nitrogen,Obs_OTUS_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Total_NitrogenITS <- Total_NitrogenITS + Vanessas_theme
Total_NitrogenITS

ggsave(filename = 'Y:/VPino_hpc/Plots/Total_NitrogenITS_mix.png',plot = Total_NitrogenITS)

Total_CarbonITS <- ggplot(soil,aes(Total.Carbon,Obs_OTUS_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Total_CarbonITS <- Total_CarbonITS + Vanessas_theme
Total_CarbonITS

ggsave(filename = 'Y:/VPino_hpc/Plots/Total_CarbonITS_mix.png',plot = Total_CarbonITS)


Total_Carbon16S <- ggplot(soil,aes(Total.Carbon,Obs_OTUS_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Total_Carbon16S <- Total_Carbon16S + Vanessas_theme
Total_Carbon16S

ggsave(filename = 'Y:/VPino_hpc/Plots/Total_Carbon16S_mix.png',plot = Total_Carbon16S)



colnames(soil)

Total_Carbon_divITS <- ggplot(soil,aes(Total.Carbon,Shannon_div_ITS,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Total_Carbon_divITS <- Total_Carbon_divITS + Vanessas_theme
Total_Carbon_divITS

ggsave(filename = 'Y:/VPino_hpc/Plots/Total_CarbondivITS_mix.png',plot = Total_Carbon_divITS)


Total_Carbon16S <- ggplot(soil,aes(Total.Carbon,Shannon_div_16S,colour=Ecosystem))+
  geom_point()+
  geom_smooth(size=1.1,method='lm')

Total_Carbon16S <- Total_Carbon16S + Vanessas_theme
Total_Carbon16S

ggsave(filename = 'Y:/VPino_hpc/Plots/Total_Carbondiv16S_mix.png',plot = Total_Carbon16S)








summary((with(soil[soil_info$Ecosystem=='CROP',],lm(Total.Nitrogen~Obs_OTUS_16S))))

linear <- (with(soil[soil_info$Ecosystem=='CROP',],lm(Total.Nitrogen~Obs_OTUS_16S)))

predicted <- predict(with(soil,linear,Total.Nitrogen[soil_info$Ecosystem=='CROP']))

require(spectroscopy)

goof(predicted,soil$Total.Nitrogen[soil_info$Ecosystem=='CROP'])
     