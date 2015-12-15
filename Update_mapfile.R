metadata_ITS_old <- read.table('Y:/VPino_hpc/ubuntu_qiime/ITS/ecometdata_ITSvpNS.txt',sep='\t',header=T)

NSW_dataset <- readRDS('RData/NSW_datasetTue_Nov_10_17_32_34_2015.rds')



NSW_dataset$responses$index #27_NAT_0_5_R1 (example)


NSW_dataset$responses <-rbind(NSW_dataset$responses,NSW_dataset$responses[NSW_dataset$responses$index=='12_NAT_5_10_R1',]) 
NSW_dataset$responses[nrow(NSW_dataset$responses),]$Depth <- '0-5' 
NSW_dataset$responses[nrow(NSW_dataset$responses),]$index <- '12_NAT_0_5_R1' 
NSW_dataset$responses[nrow(NSW_dataset$responses),]$top <- '0'
NSW_dataset$responses[nrow(NSW_dataset$responses),]$bottom <- '5'

NSW_dataset$predictions <-rbind(NSW_dataset$predictions,NSW_dataset$predictions[NSW_dataset$responses$index=='12_NAT_5_10_R1',]) 
View(NSW_dataset$predictions)

metadata_ITS_old$index <- with(metadata_ITS_old,paste(Site,Ecosystem,'0_5_R1',sep='_'))

order_for_join <- match(metadata_ITS_old$index,NSW_dataset$responses$index)

metadata_ITS_new <- data.frame(metadata_ITS_old,NSW_dataset$responses[order_for_join,],NSW_dataset$predictions[order_for_join,])

colnames(metadata_ITS_new)

col_to_keep <- c(1:14,16,17,22:37,39:87)
write.table(metadata_ITS_new[,col_to_keep],file = 'Y:/VPino_hpc/ubuntu_qiime/ITS/ecometdata_ITSvpNS_WITH_SOIL_INFO.txt',sep='\t',row.names = F,quote = F)



#####And for 16S####
metadata_16S_old <- read.table('Y:/VPino_hpc/ubuntu_qiime/16S-NSW/NSW_16S_metadatfile_R.txt',sep='\t',header=T)

NSW_dataset <- readRDS('RData/NSWdataset_Thu_Nov_12_11_21_41_2015.rds')



NSW_dataset$responses$index #27_NAT_0_5_R1 (example)


NSW_dataset$responses <-rbind(NSW_dataset$responses,NSW_dataset$responses[NSW_dataset$responses$index=='12_NAT_5_10_R1',]) 
NSW_dataset$responses[nrow(NSW_dataset$responses),]$Depth <- '0-5' 
NSW_dataset$responses[nrow(NSW_dataset$responses),]$index <- '12_NAT_0_5_R1' 
NSW_dataset$responses[nrow(NSW_dataset$responses),]$top <- '0'
NSW_dataset$responses[nrow(NSW_dataset$responses),]$bottom <- '5'

NSW_dataset$predictions <-rbind(NSW_dataset$predictions,NSW_dataset$predictions[NSW_dataset$responses$index=='12_NAT_5_10_R1',]) 
View(NSW_dataset$predictions)

metadata_16S_old$index <- with(metadata_16S_old,paste(Site,Ecosystem,'0_5_R1',sep='_'))

order_for_join <- match(metadata_16S_old$index,NSW_dataset$responses$index)

metadata_16S_new <- data.frame(metadata_16S_old,NSW_dataset$responses[order_for_join,],NSW_dataset$predictions[order_for_join,])

colnames(metadata_16S_new)

identical(metadata_16S_new$Ammonium.Nitrogen,metadata_16S_new$Ammonium.Nitrogen.1)
col_to_keep <- c(1:80)

View(metadata_16S_new[,col_to_keep])

write.table(metadata_16S_new[,col_to_keep],file = 'Y:/VPino_hpc/ubuntu_qiime/16S-NSW/8DIV1/merged_mapping_NSW_WITH_SOIL_INFO.txt',sep='\t',row.names = F,quote = F)

