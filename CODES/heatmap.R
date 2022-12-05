library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(car)
library(rgl)  
library("plot3D")
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(ggthemes)
library("gplots")
library("pheatmap")
library(ComplexHeatmap)
library("gridtext")
library("RColorBrewer")
library("gplots")

data<-read.csv(file.choose())
data<-data[order(data$ph.6.5),]
row.names(data) <- data$X
data<-data[2:4]
data <- na.omit(data)
data_matrix <- data.matrix(data)
data_heatmap <- heatmap(data_matrix, Rowv=NA, Colv=NA, col = cm.colors(256), 
                       scale="column", margins=c(5,10))
heatmap(data_matrix)


library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_fun(seq(-3, 3))
Heatmap(data_matrix, name = "mat", col = col_fun)


f1 = colorRamp2(seq(min(data_matrix), max(data_matrix), length = 3), c("blue", "#EEEEEE", "red"))
Heatmap(data_matrix, name = "mat1", col = f1, column_title = "LAB color space")

col <- colorRampPalette(brewer.pal(10, "Paired"))(256)

Heatmap(data_matrix, name = "z-score",col = col, row_names_gp = gpar(fontsize = 10), 
        heatmap_width = unit(16, "cm"), 
        heatmap_height = unit(16, "cm"),
        border=TRUE,column_title = "Cluster 1")

Heatmap(data_matrix, name = "z-score", row_names_gp = gpar(fontsize = 5), 
        heatmap_width = unit(8, "cm"), 
        heatmap_height = unit(24, "cm"),
        border=TRUE,column_title = "Cluster 8")


Heatmap(data_matrix, 
        name = "logfoldchange", #title of legend
        column_title = "Conditions", 
        row_title = "Genes", 
        border=TRUE,
        show_column_dend = TRUE,
        row_dend_reorder = TRUE,
        column_dend_reorder = NULL,
        column_order=NULL,
        cluster_rows = TRUE) 



