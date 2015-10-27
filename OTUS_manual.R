require(phyloseq)
###Import BIOM table data###
OTUS<-import_biom(BIOMfilename = 'DATA/otus_join_nf_def/otu_table_mc2_w_tax_no_pynast_failures.biom',
                  treefilename ='DATA/otus_join_nf_def/rep_set.tre',
                  refseqfilename ='DATA/otus_join_nf_def/new_refseqs.fna',
                  version=1.9)

colnames(data.frame(OTUS@otu_table@.Data))

####Correct sample data from mapfiles####
mapfile_in_use <- read.table('DATA/Mapfiles/mapping_vane80.txt',sep='\t',header=T,stringsAsFactors = F)
mapfile_corrected <- read.table('DATA/Mapfiles/mapfile_16Svp1.csv_corrected.txt',sep='\t',header=T,stringsAsFactors = F)

as.character(mapfile_corrected$BarcodeSequence) %in% as.character(mapfile_in_use$BarcodeSequence)

require(plyr)
new_mapfile <-merge(mapfile_in_use,mapfile_corrected,'BarcodeSequence')
new_mapfile <- new_mapfile[order(new_mapfile$SampleID.y),]


new_mapfile$Site <- rep(0:13,each = 6)[-c(81:84)]
new_mapfile$Rep <- rep(1:3,length.out = 80)
new_mapfile$top <- 0 
new_mapfile$bottom <-5

new_mapfile$index <-paste(new_mapfile$Site,new_mapfile$System,new_mapfile$top,new_mapfile$bottom,'R1',sep='_')

#####Import Chemical Data ####

Lab_data <-readRDS('RData/NSW_CHEMICAL_with_lab_values.RDS')

####JUST FOR NOW####
#copy the data of Site 12 nat 5 to 10 to Site 12 0 to 5 ... because that sample was not analysed#

Lab_data_corrected_tmp <-rbind(Lab_data$responses,Lab_data$responses[Lab_data$responses$index=='12_NAT_5_10_R1',]) 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$Depth <- '0-5' 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$index <- '12_NAT_0_5_R1' 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$top <- '0'
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$bottom <- '5'

Sample_data <- merge(Lab_data_corrected_tmp,new_mapfile,by.x = 'index',by.y = 'index')

# write.csv(Sample_data,'DATA/Sample_data_export.csv')
####Now add that to the OTUS ####

#first order in the same order than the OTU table#

Sample_data<-Sample_data[match(colnames(OTUS@otu_table@.Data),Sample_data$SampleID.x),]
rownames(Sample_data)<-Sample_data$SampleID.x

####How is the structure of sample data ????####
# data(GlobalPatterns) #example
# str(GlobalPatterns@sam_data)
#this command is from phyloseq package

Sample_data <-sample_data(Sample_data)
OTUS<-merge_phyloseq(OTUS,Sample_data)

colnames(tax_table(OTUS)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                               o = "Order", f = "Family", g = "Genus", s = "Species")
rank_names(OTUS)

#####NOW create our own database####

require(vegan)
OTUs_table <-data.frame(OTUS@otu_table@.Data)
OTUs_table <-as.data.frame(t(OTUs_table[,naturalsort(colnames(OTUs_table))]))
OTUs_table$System<-rep(c('NAT','NAT','NAT','CROP','CROP','CROP'),length.out = 80)
OTUs_table$Site <- rep(0:13,each = 6)[-c(81:84)]
OTUs_table$Sample <-paste0(OTUs_table$Site,OTUs_table$System)

OTU_sum<-data.frame(do.call(rbind,by(OTUs_table[,c(ncol(OTUs_table)-5:ncol(OTUs_table))],OTUs_table$Sample,colSums,simplify = T)))
OTU_sum<-OTU_sum[naturalorder(row.names(OTU_sum)),]
tail(rownames(OTU_sum))
OTU_sum$System<-c(rep(c('CROP','NAT'),length.out = 26),'NAT')
tail(OTU_sum$System)
OTU_sum$Site <- rep(0:13,each=2)[-28]
OTU_sum$index <-paste(OTU_sum$Site,OTU_sum$System,'0_5_R1',sep='_')

lab_data<-readRDS('RData/NSW_CHEMICAL_with_lab_values.RDS')$responses

Lab_data_corrected_tmp <-rbind(lab_data,lab_data[lab_data$index=='12_NAT_5_10_R1',]) 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$Depth <- '0-5' 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$index <- '12_NAT_0_5_R1' 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$top <- '0'
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$bottom <- '5'

Sample_data <- merge(Lab_data_corrected_tmp,OTU_sum,by.x = 'index',by.y = 'index')

Sample_data$Shannon_D<-apply(Sample_data[,grepl('^X',colnames(Sample_data))],1,diversity)

ggplot(Sample_data,aes(Site,Shannon_D,colour=System,group=System))+
         geom_line()


ggplot(Sample_data,aes(System,Shannon_D,colour=System,group=System))+
  geom_boxplot()

Sample_data$Simpson_D<-apply(Sample_data[,grepl('^X',colnames(Sample_data))],1,diversity,'simpson')

ggplot(Sample_data,aes(Site,Simpson_D,colour=System,group=System))+
  geom_line()


ggplot(Sample_data,aes(System,Simpson_D,colour=System,group=System))+
  geom_boxplot()

Sample_data<-Sample_data[naturalorder(Sample_data$index),]

input<-Sample_data[,!grepl('^X',colnames(Sample_data))]
input<-data.frame(ID=input[,1],apply(input[,c(10:24,31,32)],2,as.numeric))


PCA <- prcomp(input[,-1],scale = T,center = T)


biplot(PCA,scale = F)

###check###
sum(PCA$rotation[,3]*as.data.frame(scale(input[,-1]))[17,])
PCA$x[17,3]

# and nicer#

check_PCA<-data.frame(system=Sample_data$System,
                      site=Sample_data$Site,
                      INDEX=Sample_data[,2],
                      PCA$x[,1:2])
Loadings<-data.frame(Property=colnames(input[,-1]),
                     PC1_rot=PCA$rotation[,1],
                     PC2_rot=PCA$rotation[,2])

sum(Loadings$PC1_rot*as.data.frame(scale(input[,-1]))[17,])
check_PCA$PC1[17]

require(grid)
PCA_analysis<-ggplot(check_PCA)+
  geom_point(aes(PC1,PC2,colour=system),size=3)+
  ggtitle(label = paste0('First two eigenvectors for all measured properties ',' Var. exp = ',round(summary(PCA)$importance[6],2)))+
  geom_text(aes(PC1,PC2,label=INDEX,vjust=2),size=5)+  
  geom_segment(data = Loadings,aes(x=0,y=0,xend=PC1_rot*PCA$sdev[1]*sqrt(nrow(check_PCA)),
                                   yend=PC2_rot*PCA$sdev[2]*sqrt(nrow(check_PCA))),
               arrow=arrow(length=unit(0.2,"cm")),
               linetype=2)+
  geom_text(data = Loadings,aes(x=PC1_rot*PCA$sdev[1]*sqrt(nrow(check_PCA)),
                                y=PC2_rot*PCA$sdev[2]*sqrt(nrow(check_PCA)),
                                label=Property),size=5)+
  theme_bw()+
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=5),
        axis.title.y = element_text(size=5,vjust=2),
        legend.text=element_text(size=5),
        legend.title=element_text(size=5),
        strip.text=element_text(size=15),
        title=element_text(size=15))

PCA_analysis
