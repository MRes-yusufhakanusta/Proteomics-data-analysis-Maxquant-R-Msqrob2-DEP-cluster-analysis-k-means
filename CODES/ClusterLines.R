library(ggplot2)
library(DESeq2)
library(tidyverse)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(ggthemes)


data<-read.csv(choose.files(),row.names=NULL)
#dataraw<-read.csv(choose.files())
row.names(data)<-data$X
data<-data[,-1]
data <- na.omit(data)

for_clust<-data

max_itr <-  25
n_clust  <-  8# number of cluster 
set.seed(123) ## reproduce the cluster 
kmeans_out  <- kmeans(data,n_clust,iter.max = max_itr, nstart=1)

#kmeans_out  <- kmeans(data, centers = 3, nstart = 25)

data_with_cust_info <- data %>% 
  mutate(clust = paste("Cluster_", kmeans_out$cluster,sep = ""))

data_with_cust_info %>% 
  gather(key = "variable" , value = "value", -c(0,4)) %>%  ### 1 is the index of column 'geneName' and 7 is the index of column 'clust'
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num, lineType=row_num)) +   
  geom_point() + geom_line(alpha = 0.5 , aes(col = as.character(clust)))+
  theme_linedraw()+ 
  theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
  facet_wrap(~clust)+theme(text = element_text(size = 18))+
  xlab("") + ylab("Z-Score")+
  geom_line(stat = "summary", fun = "median", colour = "black", size = 1.5, 
           aes(group = 1))


write.csv(data_with_cust_info,"ClusterLines Results.csv")


data_with_cust_info %>% 
  gather(key = "variable" , value = "value", -c(0,4)) %>%  ### 1 is the index of column 'geneName' and 7 is the index of column 'clust'
  group_by(variable) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  variable , y = value , group = row_num, lineType=row_num)) +   
  geom_point() +  
  geom_line(alpha = 0.2 ,aes(col = as.character(clust))) + 
  theme_bw()+  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
  facet_wrap(~clust)+theme(text = element_text(size = 18))+
  xlab("") + ylab("Z-Score")

 # geom_line(stat = "summary", fun = "median", colour = "red", size = 1.5,aes(group = 1))
#####















