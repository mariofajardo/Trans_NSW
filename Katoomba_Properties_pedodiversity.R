require(Cubist)
require(plyr)
require(tripack)
require(spectroscopy)

if (nrow(input)<2) input<-rbind(input,input)
# input <- #needs to be a matrix which the server will provide
colnames(input)<-seq(350,2450)
#####process input####

if(names(validation)[2]=='slaking'){
  input <-trimSpec(input,wavlimits=c(500,2450),as.numeric(colnames(input)))
  input<-filter_sg(input,n=11,p=2,m=0)
  colnames(input)<-seq(500,2450)
} else {
  input <-trimSpec(input,wavlimits=c(500,2450),as.numeric(colnames(input)))
  # if(names(validation)[2]=='pH') input<-log(1/input)
  input<-filter_sg(input,n=11,p=2,m=0)
  colnames(input)<-seq(500,2450)
  input<-strip_spectra(input,datawavs = as.numeric(colnames(input)),
                       wavlimits=c(500,2450),
                       which=10)
  input <- snvBLC(input)
}
#####
# models <- #the server will provide them (Binary R object)
# validation <- #the server will provide it (Binary R object)

#####Check if its outside the first two PC convex hull ####
# pr_spectra<-prcomp(validation$spectra, center=T,scale=T) 
# # screeplot(pr_spectra)#visualize the PC
# # significance<-as.numeric(summary(pr_spectra)$importance[2,1:2]) #significance two first principal comp
# pr_scores <- pr_spectra$x 
# 
# ###do some convex.hull analysis###
# rand_tr <-tri.mesh(pr_scores[,1],pr_scores[,2])
# rand.ch <- convex.hull(rand_tr,plot.it=F)
# pr_poly <-cbind(x=c(rand.ch$y),y=c(rand.ch$y))
# # plot(pr_scores[,1],pr_scores[,2])
# # lines(c(rand.ch$x,rand.ch$x[1]),c(rand.ch$y,rand.ch$y[1]),col='blue')
# 
# scores_input <- predict(pr_spectra,input)

# logic_test<-in.convex.hull(rand_tr,scores_input[,1],scores_input[,2])
# points(scores_input[,1],scores_input[,2],col='red')

# samples_out <- which(logic_test==F)

# if(!any(logic_test)) stop('All of your samples are outside the prediction limits of the current models')


results <-foreach(model=models,.packages = 'Cubist',.combine = cbind)%dopar%{
  predict(model,input)
}

sd_predictions <- aaply(results,1,function(x) sd(x,na.rm = T))
mean_predictions <- rowMeans(results)

predictions <- data.frame(Mean=mean_predictions,Standard_deviation=sd_predictions)
if(nrow(predictions)==2) predictions <-predictions[1,] # workaround with single observations

rm(list = ls()[!ls() %in% c('predictions','samples_out','cl','pedodiversity_dataset','CORES','Chemical','usyd_DATA','sand','clay','soc','bd','Theta10','Theta1500','Log_Ks')])
gc()
#end#
