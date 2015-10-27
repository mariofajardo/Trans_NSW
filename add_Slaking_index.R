require(plyr)
require(reshape2)
STRUCTURE_DATA<-readRDS('../ASWAT/Files/Processed_DATA.rds')
CHEMICAL_DATA <-readRDS('RData/NSW_CHEMICAL_with_lab_values_13_8_2105.RDS')

View(CHEMICAL_DATA$responses)

# piecewise <- function(param){
#   (param[2]*(exp(-param[1]*(exp(-param[3]*log(t))))))
# }
# type<-'Gompertz log scale'
# 
# t_last<-7020
# t<-0:t_last

coefficients <-readRDS('../ASWAT/Files/coefficients_Gomp_logscale_logsample.RDS')

Struct<-ldply(names(STRUCTURE_DATA),function(x) {
  Struct_tmp<-data.frame(ID=x,SI_a=coefficients[[x]]$par[2],SI_b=coefficients[[x]]$par[1],SI_c=coefficients[[x]]$par[3])
})

CHEMICAL_DATA$responses$index_struct <-paste('Sample',CHEMICAL_DATA$responses$site,CHEMICAL_DATA$responses$system,CHEMICAL_DATA$responses$top,CHEMICAL_DATA$responses$bottom,sep='_') 

#Take this sample for paired statistical analyses since it wasn't analysed#
CHEMICAL_DATA$spectra <- CHEMICAL_DATA$spectra[!CHEMICAL_DATA$responses$index=='SAMPLE_12_CROP_0_5_R1',] 
CHEMICAL_DATA$responses <-CHEMICAL_DATA$responses[!CHEMICAL_DATA$responses$index=='SAMPLE_12_CROP_0_5_R1',]
#

Struct_chem <- merge(Struct,CHEMICAL_DATA$responses,by.x = 'ID',by.y = 'index_struct')
Struct_chem[,'Ca/Mg'] <-as.numeric(Struct_chem$Exc..Calcium)/as.numeric(Struct_chem$Exc..Magnesium)
Struct_chem[,'C/N'] <-as.numeric(Struct_chem$Total.Carbon)/as.numeric(Struct_chem$Total.Nitrogen)


colnames(Struct_chem)[-c(1:13,29:32)] <-c('Am_N','Nit_N','P','K','EC','pH_CaCl','pH_H2O','Ex_Al','Ex_Ca',
                                          'Ex_Mg','Ex_K','Ex_Na','N','TC','ECEC','Clay_perc','Ca/Mg','C/N')


identical(CHEMICAL_DATA$responses$index,Struct_chem$index)

Struct_chem_DATA<-list()
Struct_chem_DATA$spectra<-CHEMICAL_DATA$spectra
Struct_chem_DATA$responses<-Struct_chem

saveRDS(Struct_chem_DATA,'../usyd_spectral_lib/RData/Struct_chem_dataset.rds')

