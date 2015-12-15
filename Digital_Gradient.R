#####CORES####
if(data_type!='pit'){
u_mat<-apply(fuzzy_data[[sample]]$Clustering$membership,2,rep,each=2)
nrowmat<-ncol(u_mat)

g<-list()
for (i in 1:ncol(u_mat)){
serie <- u_mat[,i]
u<-seq(1,length(serie),lag)
u1<-seq(lag,by=lag,along.with = u)

if(u1[length(u1)]>length(serie)) {
  u<-u[-length(u)]
  u1<-u1[-length(u1)]
  u1[length(u1)]<-length(serie)
}
u1[length(u1)]<-length(serie)
mean_lag_u<-mapply(function(x,y) mean(serie[x:y]),u,u1)
mean_lag_u1<-mean_lag_u[-1]
mean_lag_u<-mean_lag_u[-length(mean_lag_u)]
g[[i]]<-(mean_lag_u-mean_lag_u1)^2 #squared differences of all the segments pairs
}

dig_grad<-sqrt((colSums(do.call(rbind,g)))/2) 
dig_grad <- c(dig_grad,dig_grad[length(dig_grad)])
par(mfrow=c(3,1))


plot(seq(1,length(serie)),u_mat[,1],type='l',
     main=paste0('Classes membership sample ',sample),
     xlab='Depth',
     ylab='Membership',
     ylim = c(0,1))

lines(seq(1,length(serie)),u_mat[,2],col='blue')
if(nrowmat>=3){
lines(seq(1,length(serie)),u_mat[,3],col='red')
}
if(nrowmat>=4){
  lines(seq(1,length(serie)),u_mat[,4],col='green')
}
if(nrowmat>=5){
  lines(seq(1,length(serie)),u_mat[,5],col='purple')
}
if(nrowmat>=6){
  lines(seq(1,length(serie)),u_mat[,6],col='yellow')
}
if(nrowmat>=6){
  lines(seq(1,length(serie)),u_mat[,7],col='moccasin')
}
if(nrowmat>=6){
  lines(seq(1,length(serie)),u_mat[,8],col='orange2')
}
if(nrowmat>=6){
  lines(seq(1,length(serie)),u_mat[,9],col='olivedrab2')
}
if(nrowmat>=6){
  lines(seq(1,length(serie)),u_mat[,10],col='paleturquoise3')
}
if(nrowmat>=6){
  lines(seq(1,length(serie)),u_mat[,11],col='salmon4')
}
if(nrowmat>=6){
  lines(seq(1,length(serie)),u_mat[,12],col='tomato3')
}
if(nrowmat>=6){
  lines(seq(1,length(serie)),u_mat[,13],col='thistle2')
}


plot(seq(lag,length(serie),lag),dig_grad,type='l',main='Digital gradient',
     xlab='Depth',
     ylab='Digital gradient',
     xlim=c(0,length(serie)))

abline(h=.4,col='green')
abline(h=.3,col='yellow')
abline(h=.2,col='red')

boundaries<-which(dig_grad>threshold,dig_grad)*lag

#check contiguous boundaries#

if (length(boundaries)>1){
  index <- c()
  for (i in 1:(length(boundaries)-1)){
    if(boundaries[i+1]-boundaries[i]>lag) index<-c(index,T) else index <-c(index,F)
  }
  index <- c(index,T) #add the last boundary 
  if(boundaries[length(boundaries)] < length(serie)-lag){ ## check if the last boundary is too close to the final
  boundaries <- boundaries[index]
  } else {
  boundaries <- boundaries[index][-length(boundaries)] #if last boundary is too close to the final
  } 
  if(boundaries[1] < lag){ #if first boundary is too close from the beggining (Posible in Organic horizons)
  boundaries <- boundaries[-1]
  }
}
##

plot(seq(1,length(serie)),col='white',xlab = 'Depth',
     yaxt='n',
     ylab='',
     xlim=c(0,length(serie)))
for (z in boundaries) {
  abline(v = z)
}


upper_bound <- c(0,boundaries)
lower_bound <- c(boundaries,nrow(fuzzy_data[[sample]]$Clustering$membership)*2)

# create_hor <- function(up,down) {
#   Mode <- function(x) {
#     ux <- unique(x)
#     ux[which.max(tabulate(match(x, ux)))]
#   }
#   
#   hor$up<-rep(Mode(rep(fuzzy_data[[sample]]$Clustering$clustering,each=2)[up:down]),(down-up))
#   
# }
create_hor <- function(up,down) {
  
  hor$up<-rep(paste0('spd_',up,'_to_',down),(down-up))
}

hor<-list()
if(length(lower_bound)>1){
horizon<-unlist(mapply(create_hor,upper_bound,lower_bound))
}else{
  horizon <- rep('spd_0_to_100',(lower_bound-upper_bound))
}
horizon_class[[sample]] <-horizon 


# test_levels<-unique(no_out_details[[i]]$b.master_hor)
# limits<-sapply(test_levels,function(x){
#   sum(no_out_details[[i]]$b.master_hor==x)*5
# })
# 
# hor_lim <- c(limits[1],sum(limits[1:2])) 
# if (length(unique(no_out_details[[i]]$b.master_hor))==3) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3])) 
# if (length(unique(no_out_details[[i]]$b.master_hor))==4) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]))
# if (length(unique(no_out_details[[i]]$b.master_hor))==5) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5])) 
# if (length(unique(no_out_details[[i]]$b.master_hor))==6) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5]),sum(limits[1:6])) 
# if (length(unique(no_out_details[[i]]$b.master_hor))==7) hor_lim <- c(limits[1],sum(limits[1:2]),sum(limits[1:3]),sum(limits[1:4]),sum(limits[1:5]),sum(limits[1:6]),sum(limits[1:7])) 
# 
# plot(x=c(1:hor_lim[length(hor_lim)]),col='white',xlab = 'Depth',yaxt='n',ylab='',xlim=c(1,hor_lim[length(hor_lim)]),main = i)
# for (z in 1:length(hor_lim)) {
#   abline(v = hor_lim[z])
# }
# 
# text(hor_lim-limits/2,25,labels = unique(no_out_details[[i]]$b.master_hor))
# 


}



####PITS####
if(data_type=='pit'){
  u_mat<-apply(fuzzy_data[[sample]]$Clustering$membership,2,rep,each=5)
  nrowmat<-ncol(u_mat)
  
  g<-list()
  for (i in 1:ncol(u_mat)){
    serie <- u_mat[,i]
    u<-seq(1,length(serie),lag)
    u1<-seq(lag,by=lag,along.with = u)
    
    if(u1[length(u1)]>length(serie)) {
      u<-u[-length(u)]
      u1<-u1[-length(u1)]
      u1[length(u1)]<-length(serie)
    }
    u1[length(u1)]<-length(serie)
    mean_lag_u<-mapply(function(x,y) mean(serie[x:y]),u,u1)
    mean_lag_u1<-mean_lag_u[-1]
    mean_lag_u<-mean_lag_u[-length(mean_lag_u)]
    g[[i]]<-(mean_lag_u-mean_lag_u1)^2 #squared differences of all the segments pairs
  }
  
  dig_grad<-sqrt((colSums(do.call(rbind,g)))/2) 
  dig_grad <- c(dig_grad,dig_grad[length(dig_grad)])
  par(mfrow=c(3,1))
  
  
  plot(seq(1,length(serie)),u_mat[,1],type='l',
       main=paste0('Classes membership sample ',sample),
       xlab='Depth',
       ylab='Membership')
  
  lines(seq(1,length(serie)),u_mat[,2],col='blue')
  if(nrowmat>=3){
    lines(seq(1,length(serie)),u_mat[,3],col='red')
  }
  if(nrowmat>=4){
    lines(seq(1,length(serie)),u_mat[,4],col='green')
  }
  if(nrowmat>=5){
    lines(seq(1,length(serie)),u_mat[,5],col='purple')
  }
  if(nrowmat>=6){
    lines(seq(1,length(serie)),u_mat[,6],col='yellow')
  }
  
  
  plot(seq(lag,length(serie),lag),dig_grad,type='l',main='Digital gradient',
       xlab='Depth',
       ylab='Digital gradient',
       xlim=c(1,length(serie))) 
  
  abline(h=.6,col='green')
  abline(h=.5,col='yellow')
  abline(h=.4,col='red')
  
  boundaries<-which(dig_grad>threshold,dig_grad)*lag
  
  #check contiguous boundaries#
  
  if (length(boundaries)>1){
    index <- c()
  for (i in 1:(length(boundaries)-1)){
    if(boundaries[i+1]-boundaries[i]>lag ) index<-c(index,T) else index <-c(index,F)
  }
  
  boundaries <- boundaries[index]
  }
  ##
  
  plot(seq(1,length(serie)),col='white',xlab = 'Depth',
       yaxt='n',
       ylab='')
  for (z in boundaries) {
    abline(v = z)
  }
  upper_bound <- c(0,boundaries)
  lower_bound <- c(boundaries,nrow(fuzzy_data[[sample]]$Clustering$data)*5)
  
  create_hor <- function(up,down) {
    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    hor$up<-rep(Mode(rep(fuzzy_data[[sample]]$Clustering$clustering,each=5)[up:down]),(down-up))
  
  }
  hor<-list()
  horizon<-unlist(mapply(create_hor,upper_bound,lower_bound))
  horizon_class[[sample]] <-horizon 
}
