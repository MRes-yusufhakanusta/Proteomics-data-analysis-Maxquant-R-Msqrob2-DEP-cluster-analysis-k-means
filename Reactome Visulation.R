library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library("RColorBrewer")

df = read.csv(file.choose(), header=TRUE)

col <- colorRampPalette(brewer.pal(10, "Paired"))(256)


S1<- ggplot(df, aes(x= ProteinFromGeneSet, y= HitGenes,
                    size=RatioOfProteinInPathway, color=P.value)) + geom_point(alpha = 1) +
  theme_linedraw()+theme(text = element_text(size = 15)) +
  ggtitle("Cluster 8") +
  xlab("Gene Count") + ylab("")

S1

S1+scale_y_continuous(
  
  # Features of the first axis
  name = "First Axis",
  
  # Add a second axis and specify its features
  sec.axis = sec_axis(~.*HitGenes, name="Second Axis")
)




S2<- ggplot(df, aes(x= reorder(Count,Description), y=Description,
                    size=BgRatio, color=qvalue, group=Group)) + geom_point(alpha = 1) +
  theme_linedraw()+theme(text = element_text(size = 15)) +
  facet_grid(.~ifelse(Significant==1, 'Upregulated', 'Downregulated'))+
  theme(panel.grid.major.y = element_line(linetype='dotted', color='blue'),
        panel.grid.major.x = element_blank())+guides(size = guide_legend(override.aes=list(shape=1)))+
  scale_color_continuous(low='blue', high='red')+ guides(size=FALSE)+
  ggtitle("Stiff vs Soft under pH:7.1") +
  xlab("Gene Count") + ylab("")

S2
