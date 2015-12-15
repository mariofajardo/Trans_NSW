require(phyloseq)
###Import BIOM table data###
OTUS<-import_biom(BIOMfilename = 'Y:/VPino_hpc/ubuntu_qiime/ITSvpNS_output/7otus_hpc/',version = 1.9)

colnames(data.frame(OTUS@otu_table@.Data))

####Correct sample data from mapfiles####
mapfile_in_use <- read.table('DATA/mapfile_ITS1.csv_corrected.txt',sep='\t',header=T,stringsAsFactors = F)
# mapfile_corrected <- read.table('DATA/Mapfiles/mapfile_16Svp1.csv_corrected.txt',sep='\t',header=T,stringsAsFactors = F)
# 
# as.character(mapfile_corrected$BarcodeSequence) %in% as.character(mapfile_in_use$BarcodeSequence)
# 
# require(plyr)
# new_mapfile <-merge(mapfile_in_use,mapfile_corrected,'BarcodeSequence')
# new_mapfile <- new_mapfile[order(new_mapfile$SampleID.y),]


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