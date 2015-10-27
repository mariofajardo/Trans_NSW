data(iris)


components<-prcomp(as.data.frame(scale(iris[,1:4])))

summary(components)

components$sdev

plot(components$sdev,type='l')

screeplot(components,type='l',col=3)

components$rotation


sum(components$rotation[,1]^2)

#first observation

as.data.frame(scale(iris[,1:4]))[1,]

#17th observation coordinate in the third PC
sum(components$rotation[,3]*as.data.frame(scale(iris[,1:4]))[17,])

components$x[17,3]


biplot(components,scale = T)

points(components$rotation[1,1],components$rotation[1,2],lwd = 4)
points(components$rotation[2,1],components$rotation[2,2],lwd = 4)
points(components$rotation[3,1],components$rotation[3,2],lwd = 4)



components<-prcomp(iris[,1:4])

biplot(components,scale = F)
points(components$rotation[1,1],components$rotation[1,2],lwd = 4)
points(components$rotation[2,1],components$rotation[2,2],lwd = 4)
points(components$rotation[3,1],components$rotation[3,2],lwd = 4)

components<-prcomp(iris[,1:4],scale. = T)


biplot(components,scale = F)

points(components$rotation[1,1],components$rotation[1,2],lwd = 4)
points(components$rotation[2,1],components$rotation[2,2],lwd = 4)
points(components$rotation[3,1],components$rotation[3,2],lwd = 4)




#check the center of first variable i.e. mean
mean(scale(iris$Sepal.Length))

components$center

components$



