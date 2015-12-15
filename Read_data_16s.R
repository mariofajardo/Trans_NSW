require(phyloseq)
###Import BIOM table data###
OTUS_16S<-import_biom(BIOMfilename = 'Y:/VPino_hpc/ubuntu_qiime/16SvpNS_output/Diversity/7otus_16SvpNS/json_biom.biom',
                      treefilename ='Y:/VPino_hpc/ubuntu_qiime/16SvpNS_output/Diversity/7otus_16SvpNS/rep_set.tre',
                      #                   refseqfilename ='DATA/otus_join_nf_def/new_refseqs.fna',
                      version=1.9)

colnames(data.frame(OTUS_16S@otu_table@.Data)) ###check what I have

####join sample data from mapfiles####
mapfile_in_use <- read.table('Y:/VPino_hpc/ubuntu_qiime/16SvpNS_output/mapfile_16SvpNS_for_R.txt',sep='\t',stringsAsFactors = F,header=T)
mapfile_in_use <- mapfile_in_use[order(mapfile_in_use[,1]),]#### need to order first

mapfile_in_use$Site <- c(rep(0:26,each = 6),rep('Blank',16),rep('Empty',10),rep('MC',4))
mapfile_in_use$Rep <- c(rep(1:3,length.out = 162),c(1,10,11:16,2:9),c(1,10,2:9),1:4)
mapfile_in_use$top <- 0 
mapfile_in_use$bottom <-5
mapfile_in_use$system <- c(rep(c(rep('Nat',3),rep('Crop',3)),27),rep('Blank',16),rep('Empty',10),rep('MC',4))
mapfile_in_use$index <-toupper(with(mapfile_in_use,paste(Site,system,top,bottom,'R1',sep='_'))) #0_NAT_0_5_R1
#####Import Chemical Data ####

Lab_data <-readRDS('RData/NSW_datasetFri_Oct_30_11_38_43_2015.rds')[c('responses','predictions')]

Lab_data <- cbind(Lab_data$responses,Lab_data$predictions)

colnames(Lab_data)

Lab_data <- Lab_data[,-c(1,4:7)]

####JUST FOR NOW####
#copy the data of Site 12 nat 5 to 10 to Site 12 0 to 5 ... because that sample was not analysed#

Lab_data_corrected_tmp <-rbind(Lab_data,Lab_data[Lab_data$index=='12_NAT_5_10_R1',]) 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$Depth <- '0-5' 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$index <- '12_NAT_0_5_R1' 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$top <- '0'
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$bottom <- '5'

Sample_data <- merge(Lab_data_corrected_tmp,mapfile_in_use,by.x = 'index',by.y = 'index')

# write.csv(Sample_data,'DATA/Sample_data_export.csv')
####Now add that to the OTUS ####

#first order in the same order than the OTU table#

Sample_data<-Sample_data[match(colnames(OTUS_16S@otu_table@.Data),Sample_data$SampleID),]
rownames(Sample_data)<-colnames(OTUS_16S@otu_table@.Data)
Sample_data[is.na(Sample_data)]<-0

####How is the structure of sample data ????####
# data(GlobalPatterns) #example
# str(GlobalPatterns@sam_data)
#this command is from phyloseq package

Sample_data <-sample_data(Sample_data)


OTUS_16S<-merge_phyloseq(OTUS_16S,Sample_data)

colnames(tax_table(OTUS_16S)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                   o = "Order", f = "Family", g = "Genus", s = "Species")
rank_names(OTUS_16S)



OTUS_16S <- prune_samples(grepl('R1',as.data.frame(OTUS_16S@sam_data@.Data)[,1]),OTUS_16S)

# date<-gsub(' |:','_',date())
# saveRDS(OTUS_16S,paste0('RData/','OTUS_16S_',date,'.rds'))

