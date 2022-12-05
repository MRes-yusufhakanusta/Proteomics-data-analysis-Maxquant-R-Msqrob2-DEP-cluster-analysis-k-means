

library(factoextra)
library(NbClust)
library(ggplot2)
df <- scale(USArrests)
head(df)

fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept =5, linetype = 2)+
  labs(subtitle = "Elbow method")


df<-read.csv(choose.files())
df2<-as.matrix(df)
row.names(df)<-df$X
colnames(df2)<-NULL
df<-df[2:4]
fviz_nbclust(df2, kmeans, method = "wss")+ geom_vline(xintercept =5, linetype = 2)+
  labs(subtitle = "Elbow method")+ geom_line(size=3)
fviz_nbclust(df2, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
