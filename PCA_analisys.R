Transect_NSW_surface_samples<-readRDS('RData/NSW_CHEMICAL_with_lab_values.RDS')

require(ggplot2)
require(reshape2)

check_properties<-melt(Transect_NSW_surface_samples$responses,id.vars = colnames(Transect_NSW_surface_samples$responses)[c(1:8,24:28)])
check_properties$value<-as.numeric(check_properties$value)

pdf('Plots/check_lab_results.pdf')
for (i in unique(check_properties$variable)){
  print(ggplot(check_properties[check_properties$variable==i,],aes(coords.x1,coords.x2,group=variable,colour=system))+
          geom_point(aes(size=value))+
          facet_grid(~system)+
          ggtitle(label = i)+
          coord_equal())
}
dev.off()


str(Transect_NSW_surface_samples$responses[,colnames(Transect_NSW_surface_samples$responses)[-c(1:8,24:28)]])


for (i in colnames(Transect_NSW_surface_samples$responses)[-c(1:8,24:28)]){
  Transect_NSW_surface_samples$responses[,i]<-as.numeric(Transect_NSW_surface_samples$responses[,i])
}

scaled_values<-scale(Transect_NSW_surface_samples$responses[,colnames(Transect_NSW_surface_samples$responses)[-c(1:8,24:28)]])
values<-Transect_NSW_surface_samples$responses[,colnames(Transect_NSW_surface_samples$responses)[-c(1:8,24:28)]]

components<-prcomp(values,scale. = T,center = T)

#CHECK
# scaled_values[1,] #first observation 
scale(values)[1,]

#17th observation coordinate in the third PC
sum(components$rotation[,3]*scale(values)[17,])

components$x[17,3]

check_PCA<-data.frame(scaled_values,components$x[,1:2])
check_PCA$PC1_rot<-components$rotation[,1]
check_PCA$PC2_rot<-components$rotation[,2]
check_PCA$coord_x<-as.numeric(Transect_NSW_surface_samples$responses$coords.x1)
check_PCA$coord_y<-as.numeric(Transect_NSW_surface_samples$responses$coords.x2)
check_PCA<-cbind(check_PCA,Transect_NSW_surface_samples$responses[,24:28])

print(ggplot(check_PCA,aes(coord_x,coord_y,colour=system))+
        geom_point(aes(size=PC1))+
        scale_size(guide = 'none')+
        facet_grid(~system)+
        ggtitle(label = paste0('first eigenvalue','  Var. exp = ',round(summary(components)$importance[2],2))))

check_PCA$Loadings<-rep(colnames(check_PCA)[1:15])
require(grid)
ggplot(check_PCA)+
  geom_point(aes(PC1,PC2))+
  facet_grid(~system)+
  ggtitle(label = paste0('First two eigenvectors for all measured properties ',' Var. exp = ',round(summary(components)$importance[6],2)))+
  geom_text(aes(PC1,PC2,label=site,vjust=2))+  
  geom_segment(aes(x=0,y=0,xend=PC1_rot*components$sdev[1]*sqrt(nrow(check_PCA)),
                   yend=PC2_rot*components$sdev[2]*sqrt(nrow(check_PCA))),
               arrow=arrow(length=unit(0.2,"cm")))+
  geom_text(aes(x=PC1_rot*components$sdev[1]*sqrt(nrow(check_PCA)),
                y=PC2_rot*components$sdev[2]*sqrt(nrow(check_PCA)),
                label=Loadings))





biplot(components,scale = F)
points(components$rotation[1,1],components$rotation[1,2],lwd = 4)



#####Pedometrics 2015####
# 
# for (i in unique(check_properties$variable)){
#   
#   print(ggplot(check_properties[check_properties$variable==i&check_properties$site%in%0:27,],aes(coords.x1,coords.x2,group=variable,colour=system))+
#           geom_point(aes(size=value))+
#           facet_grid(~system)+
#           scale_size(guide='none')+
#           ggtitle(label = i)+
#           coord_equal())
#   
# }
# 
# 
# library(sp)
# library(rgdal)
# 
# coord_map<-check_properties[,2:3]
# 
# coordinates(coord_map) <- ~ coords.x1 + coords.x2
# proj4string(coord_map) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# # writeOGR(Total_stocks, "Plots/Total_stocks.kml",layer='OC',driver="KML") 
# ######
# 
# 
# b <- 5
# plot_extent <- bbox(matrix(as.vector(bbox(coord_map) + c(-b,-b,b,b)),ncol=2,byrow=T)) #add border
# 
# 
# require(ggmap)
# require(grid)
# map_sat <- get_map(location = plot_extent,
#                    scale=4,
#                    maptype='satellite',
#                    messaging=T,
#                    source='google')
# map_polit <- get_map(location = plot_extent,
#                      scale=4,
#                      maptype='roadmap',
#                      messaging=T,
#                      source='google')
# 
# map_hybrid <- get_map(location = plot_extent,
#                       scale=4,
#                       maptype='hybrid',
#                       messaging=T,
#                       source='google')
# 
# for (i in unique(check_properties$variable)){
#  
# print(ggmap(map_hybrid, extent = 'device',maprange=F,legend='topright') +
#   geom_point(data=check_properties[check_properties$variable==i&check_properties$system=='Nat',],
#              aes(x=coords.x1, y=coords.x2,colour=value),size=10)+ 
#     scale_color_gradient(low = 'blue',high = 'red')+
#   scale_size(guide="none",range=c(5,20))+
#   theme(legend.key.size = unit(.8,'lines'),
#         legend.title = element_text(size = 20, face = 'bold'),
#         legend.text = element_text(size = 20),
#         axis.text.x = element_text(size=16),
#         axis.text.y = element_text(size=16),
#         axis.title.x = element_text(size=20),
#         axis.title.y = element_text(size=20,vjust=2),
#         title=element_text(size=35))+
#   ggtitle(paste0('Transect NSW :',i))+
#   geom_text(data=check_properties,
#             aes(x=coords.x1, y=coords.x2,label=site),size=4,color='white'))
# 
# ggsave(paste0('../../../../../../../Dropbox/Mario-Vane/Pedometrics_2015_Figures/map_property_',i,'.png'),width = 15,height = 10)
# }
# 
# 
# ggmap(map_hybrid, extent = 'device',maprange=F,legend='topright') +
#         geom_point(data=check_PCA[check_PCA$site%in%27:48,]
#                    , aes(x=coord_x, y=coord_y,colour=PC1,size=PC1),alpha=.8)+ 
#         scale_color_gradient(low = 'blue',high = 'red')+
#         scale_size(guide="none",range=c(5,20))+
#         theme(legend.key.size = unit(.8,'lines'),
#               legend.title = element_text(size = 20, face = 'bold'),
#               legend.text = element_text(size = 20),
#               axis.text.x = element_text(size=16),
#               axis.text.y = element_text(size=16),
#               axis.title.x = element_text(size=20),
#               axis.title.y = element_text(size=20,vjust=2),
#               title=element_text(size=35))+
#         ggtitle('First Principal component of soil properties transect E-W')+
#   geom_text(data=check_properties[check_properties$variable==i&check_properties$site%in%27:48,],
#             aes(x=coords.x1, y=coords.x2,label=site),size=4,color='white')
#         
# 
# ggsave('../../../../../../../Dropbox/Mario-Vane/Pedometrics_2015_Figures/PC1_E_W_.png',width = 15,height = 10)
# 
require(naturalsort)

check_properties$site<-as.numeric(check_properties$site)
for (i in unique(check_properties$variable)){
print(ggplot(check_properties[check_properties$site%in%0:26&check_properties$top==0&check_properties$variable==i,],
                    aes(x=site*50,y=value,group=system,colour=system))+
              geom_point(size=3)+
              theme(strip.text=element_text(size=20,colour = 'white'),
              axis.text.x=element_text(size=20),
              axis.text.y=element_text(size=20),
              axis.title.x=element_text(size=20),
              axis.title.y=element_text(size=20),
              title=element_text(size=30),
              legend.text=element_text(size = 20))+
        ggtitle(label = paste('Transect NS, Property: ',i))+
        xlab(label = 'Distance in km (approx)')+
        geom_smooth(aes(colour=system)))

# print(ggplot(check_properties[check_properties$site%in%27:48&check_properties$top==0&check_properties$variable==i,],
#              aes(x=(site-26)*50,y=value,group=variable,colour=system))+
#         geom_point(size=3)+
#         theme(strip.text=element_text(size=20,colour = 'white'),
#               axis.text.x=element_text(size=20),
#               axis.text.y=element_text(size=20),
#               axis.title.x=element_text(size=20),
#               axis.title.y=element_text(size=20),
#               title=element_text(size=30),
#               legend.text=element_text(size = 20))+
#         ggtitle(label = paste('Transect EW, Property: ',i))+
#         xlab(label = 'Distance in km (approx)')+
#         geom_smooth()+
#         scale_x_reverse())


ggsave(paste0('../../../../../../../Dropbox/Mario-Vane/Pedometrics_2015_Figures/property_NS ',i, '.png'),width = 15,height = 10)

}
      

