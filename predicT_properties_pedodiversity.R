#####CEC####
pedodiversity_dataset <- readRDS('RData/pedodiversity_dataset.RDS') 
require(foreach)
require(doSNOW)
cl <-makeCluster(5)
registerDoSNOW(cl)
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_CEC_surf_cores.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_CEC_surf_cores.RData')
source('Katoomba_Properties_pedodiversity.R')

pedodiversity_dataset$predicted_CEC <- predictions$Mean
pedodiversity_dataset$predicted_CEC_sd <- predictions$Standard_deviation


#####pH####
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_pH_surf_cores.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_pH_surf_cores.RData')
source('Katoomba_Properties_pedodiversity.R')

pedodiversity_dataset$predicted_pH <- predictions$Mean
pedodiversity_dataset$predicted_pH_sd <- predictions$Standard_deviation

#####Clay####
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Clay_surf_transect.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Clay_surf_transect.RData')
source('Katoomba_Properties_pedodiversity.R')

pedodiversity_dataset$predicted_Clay <- predictions$Mean
pedodiversity_dataset$predicted_Clay_sd <- predictions$Standard_deviation

#####EC####
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]
models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_EC_surf_cores.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_EC_surf_cores.RData')
source('Katoomba_Properties_pedodiversity.R')

pedodiversity_dataset$predicted_EC <- predictions$Mean
pedodiversity_dataset$predicted_EC_sd <- predictions$Standard_deviation

#####TC####
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_TC_surf_cores.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_TC_surf_cores.RData')
source('Katoomba_Properties_pedodiversity.R')

pedodiversity_dataset$predicted_TC <- predictions$Mean
pedodiversity_dataset$predicted_TC_sd <- predictions$Standard_deviation


#####Slaking_a####
input <- pedodiversity_dataset[,grepl('X',colnames(pedodiversity_dataset))]

models <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_50_iter_models_Slaking_a_coef.RData')
validation <- readRDS('../usyd_spectral_lib/Katoomba_spectral_models/Katoomba_val_Slaking_a_coef.RData')
source('Katoomba_Properties_pedodiversity.R')
stopCluster(cl)

pedodiversity_dataset$predicted_SI_a <- predictions$Mean
pedodiversity_dataset$predicted_SI_a_sd <- predictions$Standard_deviation

saveRDS(pedodiversity_dataset,file = 'RData/pedodiversity_dataset_with_katoomba_preds.RDS')
