require(phyloseq)
###Import BIOM table data###
OTUS_16S<-import_biom(BIOMfilename = 'Y:/VPino_hpc/ubuntu_qiime/16SvpNS_output/Diversity/7otus_16SvpNS/json_biom.biom',
                  treefilename ='Y:/VPino_hpc/ubuntu_qiime/16SvpNS_output/Diversity/7otus_16SvpNS/rep_set.tre',
                  #                   refseqfilename ='DATA/otus_join_nf_def/new_refseqs.fna',
                  version=1.9)

colnames(data.frame(OTUS@otu_table@.Data)) ###check what I have

####join sample data from mapfiles####
mapfile_in_use <- read.table('Y:/VPino_hpc/ubuntu_qiime/16SvpNS_output/mapfile_16SvpNS.txt',sep='\t',stringsAsFactors = F)
mapfile_in_use <- mapfile_in_use[order(mapfile_in_use[,1]),]#### need to order first

mapfile_in_use$Site <- c(rep(0:26,each = 6),rep('Blank',16),rep('Empty',10),rep('MC',4))
mapfile_in_use$Rep <- c(rep(1:3,length.out = 162),c(1,10,11:16,2:9),c(1,10,2:9),1:4)
mapfile_in_use$top <- 0 
mapfile_in_use$bottom <-5
mapfile_in_use$system <- c(rep(c(rep('Nat',3),rep('Crop',3)),27),rep('Blank',16),rep('Empty',10),rep('MC',4))
mapfile_in_use$index <-toupper(with(mapfile_in_use,paste(Site,system,top,bottom,'R1',sep='_'))) #0_NAT_0_5_R1
#####Import Chemical Data ####

Lab_data <-readRDS('RData/NSW_CHEMICAL_with_lab_values.RDS')

####JUST FOR NOW####
#copy the data of Site 12 nat 5 to 10 to Site 12 0 to 5 ... because that sample was not analysed#

Lab_data_corrected_tmp <-rbind(Lab_data$responses,Lab_data$responses[Lab_data$responses$index=='12_NAT_5_10_R1',]) 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$Depth <- '0-5' 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$index <- '12_NAT_0_5_R1' 
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$top <- '0'
Lab_data_corrected_tmp[nrow(Lab_data_corrected_tmp),]$bottom <- '5'

Sample_data <- merge(Lab_data_corrected_tmp,mapfile_in_use,by.x = 'index',by.y = 'index')

# write.csv(Sample_data,'DATA/Sample_data_export.csv')
####Now add that to the OTUS ####

#first order in the same order than the OTU table#

Sample_data<-Sample_data[match(colnames(OTUS_16S@otu_table@.Data),Sample_data$V1),]
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



#####Export biodiversity for NSW map####

bact_table<-as.matrix(t(OTUS_16S@otu_table))
class(bact_table)<-'matrix'
bact_table<-bact_table[order(colnames(OTUS_16S@otu_table)),]
biodiversity<-data.frame(Site=OTUS_16S@sam_data$Site,
                         pH_CaCl=as.numeric(OTUS_16S@sam_data$pH.Level..CaCl2.),
                         ECEC=as.numeric(OTUS_16S@sam_data$ECEC),
                         Tot_P=as.numeric(OTUS_16S@sam_data$Phosphorus.Colwell),
                         Total_N=as.numeric(OTUS_16S@sam_data$Total.Nitrogen),
                         System=OTUS_16S@sam_data$system.y,
                         x=OTUS_16S@sam_data$coords.x1,
                         y=OTUS_16S@sam_data$coords.x2,
                         index=OTUS_16S@sam_data$index,
                         index1=OTUS_16S@sam_data$V1)
biodiversity <- biodiversity[order(biodiversity$index1),]

test<-as.data.frame(do.call(rbind,by(bact_table,biodiversity$index,colMeans)))
test$index<-rownames(test)

require(plyr)
BIO_diver <- join(test,biodiversity)


require(vegan)
BIO_diver_DATA<-data.frame(BIO_diver[,187211:187220],
                           shannon_div=diversity(BIO_diver[,1:187210]),
                           simpson_div=diversity(BIO_diver[,1:187210],index = 'simpson'),
                           invsimpson_div=diversity(BIO_diver[,1:187210],index = 'invsimpson'))


BIO_diver_DATA<-BIO_diver_DATA[order(BIO_diver_DATA$Site),]
BIO_diver_DATA_samples <- BIO_diver_DATA[-c(1:27),]

require(naturalsort)
BIO_diver_DATA_samples <-BIO_diver_DATA_samples[naturalorder(BIO_diver_DATA_samples$Site),]
plot(BIO_diver_DATA_samples$x[],BIO_diver_DATA_samples$invsimpson_div)
plot(BIO_diver_DATA_samples$x[],BIO_diver_DATA_samples$simpson_div)
plot(BIO_diver_DATA_samples$x[],BIO_diver_DATA_samples$shannon_div)

# write.csv(BIO_diver_DATA,file='../../../../HPC/NSW_biodiversity_all.txt')
# ####Basic Analysis based on vignette example####
# 
# ###Prune by taxa
# Top_10_OTUs_names <- names(sort(taxa_sums(OTUS), TRUE)[1:10])
# Top_10_OTUs <- prune_taxa(Top_10_OTUs_names, OTUS)
# 
# #bar plots#
# plot_bar(Top_10_OTUs, "Site",fill='Class',facet_grid = ~System)
# 
# 
# Thermoleophilia_OTUs <- subset_taxa(OTUS,Class=='c__Thermoleophilia')
# 
# plot_bar(Thermoleophilia_OTUs, "Site",fill='Class',facet_grid = ~System)
# 
# Archaea_OTUs <- subset_taxa(OTUS,Kingdom=='k__Archaea')
# 
# plot_bar(Archaea_OTUs, "Site",fill='Species',facet_grid = ~System)
# 
# alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
# 
# (p <- plot_richness(Archaea_OTUs, "System", measures=alpha_meas))
# 
# 
# ## richness_estimates
# p + geom_boxplot(data=p$data, aes(x=System, y=value, color=System), alpha=0.1)
# 
# 
# variable <-Top_10_OTUs
# 
# PC_OTUS <- ordinate(variable, "PCoA", "bray")
# p <- plot_ordination(variable, PC_OTUS,type='species',color = "Family")
# p + geom_point(size = 6, alpha = 0.7)
# 
# ##
# plot_tree(Top_10_OTUs, color="Family", shape="System",label.tips="Family", size="Abundance")
# 
# ## 
# require(naturalsort)
# 
# test<-subset_samples(Top_10_OTUs,sample_names(Top_10_OTUs)%in%naturalsort(sample_names(Top_10_OTUs))[sort(c(seq(1,80,6),seq(1,80,6)+1,seq(1,80,6)+2))])
# 
# (p <- plot_heatmap(title = 'Ten most abundant OTUs distribution in Natural Systems',test, "NMDS", "bray", "Site",sample.order = naturalsort(sample_names(test))))
# 
# test<-subset_samples(Top_10_OTUs,sample_names(Top_10_OTUs)%in%naturalsort(sample_names(Top_10_OTUs))[sort(c(seq(4,80,6),seq(4,80,6)+1,seq(4,80,6)+2))])
# 
# (p <- plot_heatmap(test,title = 'Ten most abundant OTUs distribution in Crop Systems', "NMDS", "bray", "Site",sample.order = naturalsort(sample_names(test))))
# 
# 
# test<-subset_samples(OTUS,sample_names(OTUS)%in%naturalsort(sample_names(OTUS))[sort(c(seq(1,80,6),seq(1,80,6)+1,seq(1,80,6)+2))])
# 
# (p <- plot_heatmap(title = 'OTUs distribution in Natural Systems',test, "NMDS", "bray", "Site",sample.order = naturalsort(sample_names(test))))
# 
# test<-subset_samples(OTUS,sample_names(OTUS)%in%naturalsort(sample_names(OTUS))[sort(c(seq(4,80,6),seq(4,80,6)+1,seq(4,80,6)+2))])
# 
# (p <- plot_heatmap(title = 'OTUs distribution in Crop Systems',test, "NMDS", "bray", "Site",sample.order = naturalsort(sample_names(test))))
# 
# 
# 
# ##plot_sample_network #
# Network <- make_network(test, max.dist=0.9)
# plot_network(Network,test,color='Site',line_weight=0.4, label=NULL)
# 
# ------------------------------------------------
# library("reshape2")
# # Melt the species-data.frame, DF, to facet each CA axis separately
# mdf <- melt(p1$data[, c("CA1", "CA2", "Phylum", "Family", "Genus")], 
#             id=c("Phylum", "Family", "Genus") )
# # Select some special outliers for labelling
# LF <- subset(mdf, variable=="CA2" & value < -1.0)
# # build plot: boxplot summaries of each CA-axis, with labels
# p <- ggplot(mdf, aes(Phylum, value, color=Phylum)) + geom_boxplot() + 
#   facet_wrap(~variable, 2) + scale_colour_hue(guide = FALSE) +
#   theme_bw() + theme( axis.text.x = element_text(angle = -90, hjust = 0) )
# # Add the text label layer, and render ggplot graphic
# (p <- p + geom_text(aes(Phylum, value+0.1, color=Phylum, label=Family), 
#                     data=LF, vjust=0, size=2) )
# 
# 



######Diversity analysis#####
dat<- read.table("../../../../HPC/NSW_biodiversity_all.txt", header=T,sep=",") # soil data

print(ggplot(dat[dat$Site%in%0:26,],
             aes(x=Site*50,y=invsimpson_div,group=System,colour=System))+
        geom_point(size=3)+
        theme(strip.text=element_text(size=20,colour = 'white'),
              axis.text.x=element_text(size=20),
              axis.text.y=element_text(size=20),
              axis.title.x=element_text(size=20),
              axis.title.y=element_text(size=20),
              title=element_text(size=30),
              legend.text=element_text(size = 20))+
        ggtitle(label = 'Transect NS first 13 sites, Inverse Simpson diversity')+
        xlab(label = 'Distance in km (approx)')+
        geom_smooth(aes(colour=System),size=2))
ggsave(filename = '../../../../../../../Dropbox/Mario-Vane/Pedometrics_2015_Figures/Simpson_13sites.png',width = 16,height = 10)


