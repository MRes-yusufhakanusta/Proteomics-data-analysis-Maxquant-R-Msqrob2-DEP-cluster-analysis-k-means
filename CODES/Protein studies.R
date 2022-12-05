
library("gplots")
library("pheatmap")

df<-read.csv(file.choose(), sep=",")
row.names(df)<-df$GeneName
df <- df[,2:6]
df <- data.matrix(df)
df <- df[order(df$Stiff.Soft),]

heatmap(df, Rowv=NA, Colv=NA, col = heat.colors(256), scale="row", margins=c(5,10))
heatmap(df, scale="row")
heatmap.2(df, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
pheatmap(df, cutree_rows = 4)

library(ComplexHeatmap)
library("gridtext")
set.seed(123)
Heatmap(df, 
        name = "logfoldchange", #title of legend
        column_title = "Conditions", 
        row_title = "Genes", 
        border=TRUE,
      show_column_dend = TRUE,
       row_dend_reorder = TRUE,
       column_dend_reorder = NULL,
        column_order=NULL,
        cluster_rows = TRUE, 
       row_km = 2,
        column_km = 4,
        row_names_gp = gpar(col = c(rep("red", 5), rep("blue", 5)))) 


 rowAnnotation(foo = anno_text(gt_render(sapply(LETTERS[1:10], strrep, 10), align_widths = TRUE), 
                  gp = gpar(box_col = "blue", box_lwd = 2), 
                  just = "right", 
                  location = unit(1, "npc")))



ha=HeatmapAnnotation(foo = letters[1:4],
                  annotation_label = gt_render("**Annotation** _one_",
                                               gp = gpar(box_col = "black")),
                  show_legend = FALSE)
Heatmap(df, top_annotation = ha)

##

library(ggplot2)
library(dplyr)
library(gganimate)

ggplot(df, aes(x=df$GeneName, y=df$Stiff.Soft.pH7.1, fill=df$Stiff.Soft.pH7.1))+
  geom_bar(stat="identity")+
  geom_text(aes(label=df$Stiff.Soft.pH7.1), vjust=-0.3, size=3.5)+
  theme_light()

ggplot(df, aes(x=Gene.names, y=logFC, fill=Condition)) + 
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+theme_light()+
  labs(title = "Upregulated and Downregulated Gene", subtitle = "ALL")+  xlab("Conditions") + 
  ylab("Log2 Fold Change")
  
  
  geom_text(aes(label=logFC), vjust=0, size=2.5)#value of all conditions





