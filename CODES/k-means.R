library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(car)
library(rgl)  
library("plot3D")

data<-read.csv(choose.files(),row.names=NULL)
dataraw<-read.csv(choose.files())
row.names(data)<-data$X
data<-data[,-1]
data<-data[1:3]
data <- na.omit(data)
data <- scale(data)
head(data)
distance <- get_dist(data)

fviz_nbclust(data, kmeans, method = "wss")+
  geom_vline(xintercept =8, linetype = 2)+
  labs(subtitle = "Elbow method")+ geom_line(size=3)

gap_stat <- clusGap(data,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)

fviz_gap_stat(gap_stat)

set.seed(1)
km <- kmeans(data, centers = 8, nstart = 25)
fviz_cluster(km, data = data, geom = "point", stand=1, ellipse.type = "convex",
             ellipse.alpha = 0.2,
             
             ggtheme = theme_bw(),main = "K-means Cluster Analyze")


fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"),
          lab_size = 2, show_labels = FALSE)



k2 <- kmeans(data, centers = 5, nstart = 25)
str(k2)
fviz_cluster(k2, data = data)+theme_bw()


set.seed(123)

wss <- function(k) {
  kmeans(data, k, nstart = 10 )$tot.withinss
}

k.values <- 1:15
wss_values <- map_dbl(k.values, wss)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
fviz_nbclust(data, kmeans, method = "wss")

fviz_nbclust(data, kmeans, method = "wss", linecolor = "red", print.summary=1,
             verbose = 1)+ 
  geom_vline(xintercept =5, linetype = 2)+
  labs(subtitle = "Elbow method")+ geom_line(size=3)+theme_bw()


KM <- kmeans(data, centers = 8)

scatter3D(
  x = dataraw$pH.6.8,
  y = dataraw$ph.6.5, 
  z = dataraw$pH.7.1, 
  bg.col = c("white"), 
  ellipsoid.alpha = 0.2, 
  xlab = "Stiff/Soft_pH:6.8", 
  ylab = "Stiff/Soft_pH:6.5", 
  zlab = "Stiff/Soft_pH:7.1", 
  surface.col = colorRampPalette(c("blue", "yellow", "red"))(length(levels(factor(KM$cluster)))), 
  groups = factor(KM$cluster),
  grid = TRUE,
  surface = 0,
  ellipsoid = TRUE,
  phi = 0, bty = "g",
  pch = 20, cex = 3, ticktype = "detailed", main = "k-means cluster (centers=8)", 
  clab = c("logFC", "(Stiff/Soft)")
)





